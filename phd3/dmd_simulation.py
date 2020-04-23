#!/usr/bin/env python3
"""
Author  ==>> Matthew R. Hennefarth
Date    ==>> April 16, 2020
"""

#Standard Library Imports
import logging
import os
import shutil
import json
import signal
import subprocess
import sys
import math
from timeit import default_timer as timer
import datetime
from subprocess import Popen, PIPE

#PHD3 Imports
import phd3.protein as protein
import phd3.utility.utilities as utilities
from phd3.setupjob import setupDMDjob
from phd3.utility.exceptions import ParameterError
from phd3.titrate import titrate_protein

logger=logging.getLogger(__name__)

__all__ = [
    'dmd_simulation'
]

class dmd_simulation:

    __slots__=["_submit_directory", "_scratch_directory", "_config", "_cores",
            "_time_to_run", "_timer_went_off",  "_start_time",
            "_parameter_file", "_raw_parameters", "_commands", "_src_files" ,"_titration"]

    def __init__(self, cores: int = 1, run_dir: str='./', time=-1, pro: protein.Protein=None, parameters: dict=None):

        logger.debug("Initializing variables")
        self._submit_directory = os.getcwd()
        self._scratch_directory = run_dir
        self._config = utilities.load_phd_config()
        self._cores = cores
        self._time_to_run = time
        self._timer_went_off = False
        self._start_time = 0

        # Want to make sure that we make the scratch directory!!
        try:
            logger.debug(f"Checking if scratch directory: {self._scratch_directory} is created")
            if not os.path.isdir(self._scratch_directory) and os.path.abspath(self._scratch_directory) != os.path.abspath(self._submit_directory):
                logger.info(f"Creating scratch directory: {self._scratch_directory}")
                os.mkdir(self._scratch_directory)

        except OSError:
            logger.warning("Could not make scratch directory : {self._scratch_directory}")
            logger.warning("I will run job in the current directory.")
            self._scratch_directory = './'

        utilities.setup_dmd_environ()

        if parameters is None:
            logger.debug("Checking for a dmdinput.json file")
            if os.path.isfile("dmdinput.json"):
                self._parameter_file = os.path.join(self._submit_directory, "dmdinput.json")

            else:
                logger.error("No parameters specified for the job!")
                raise FileNotFoundError("dmdinput.json")

            #Now we read in the parameters here to a dictionary
            try:
                logger.debug("Loading in parameters")
                with open(self._parameter_file, 'r') as inputfile:
                    self._raw_parameters = json.load(inputfile)

            except IOError:
                logger.exception("Could not open the parameter file correctly!")
                raise

        else:
            logger.debug("Using parameters passed")
            self._raw_parameters = parameters

        # Now we check to see if the parameters are indeed valid
        try:
            utilities.valid_dmd_parameters(self._raw_parameters)

        except ValueError:
            logger.exception("Missing a parameter definition!")
            raise

        except ParameterError:
            logger.exception("Invalid parameter specification")
            raise

        if self._raw_parameters["titr"]["titr on"]:
            self._titration = titrate_protein(self._raw_parameters["titr"])
            self._raw_parameters = self._titration.expand_commands(self._raw_parameters)            

        else:
            self._titration = None

        # TODO check for any exceptions raised from setupDMDjob
        if pro is None:
            if not os.path.isfile("initial.pdb"):
                logger.debug("initial.pdb not found, will try setting up from scratch")
                sj = setupDMDjob(parameters=self._raw_parameters)
                sj.full_setup()

        else:
            logger.debug("Will setup the protein for DMD")
            sj = setupDMDjob(parameters=self._raw_parameters, pro=pro)
            sj.full_setup()

        if os.path.isfile(self._raw_parameters["Echo File"]):
            with open(self._raw_parameters["Echo File"]) as echofile:
                lines = []
                for line in echofile:
                    lines.append(line)

                last_line = lines[-1].split()
                self._start_time = int(float(last_line[0]))
                logger.debug(f"Last recorded time: {self._start_time}")

        if self._raw_parameters["Remaining Commands"]:
            self._commands = self._raw_parameters["Remaining Commands"].copy()

            if self._raw_parameters["Commands"]:
                all_commands = self._raw_parameters["Commands"].copy()

                remove = len(self._commands)
                for i in range(remove):
                    all_commands.pop(list(self._commands.keys())[-1])

                time_elapsed = 0
                for step in all_commands:
                    if "Time" in all_commands[step].keys():
                        time_elapsed += all_commands[step]["Time"]

                    else:
                        time_elapsed += self._raw_parameters["Time"]

                logger.debug(f"Total time completed: {time_elapsed}")

                diff = self._start_time - time_elapsed
                logger.debug(f"We are off by: {diff}")

                if "Time" in self._commands[list(self._commands.keys())[0]].keys():
                    new_time = self._commands[list(self._commands.keys())[0]]["Time"] - diff

                else:
                    new_time = self._raw_parameters["Time"] - diff

                if new_time < 0:
                    logger.error("Somehow we moved onto a later step then what is reported in remaining calculations.")
                    raise ValueError("Invalid time")

                logger.debug(f"Setting new time for the first step to: {new_time}")
                self._commands[list(self._commands.keys())[0]]["Time"] = new_time

            else:
                logger.warning("Unknown how many steps prior to this one!")
                logger.warning("Will just start continue from where we left off then")

        elif self._raw_parameters["Commands"]:
            # There are no remaining commands, so continue like normal more or less
            self._commands = self._raw_parameters["Commands"].copy()

        else:
            self._commands = {"1": {}}
            logger.debug("Commands passed, using those")

        if os.path.abspath(self._scratch_directory) != os.path.abspath(self._submit_directory):
            logger.info(f"Copying files from {os.path.abspath(self._submit_directory)} to {os.path.abspath(self._scratch_directory)}")
            self._src_files = os.listdir(self._submit_directory)
            for file_name in self._src_files:
                full_file_name = os.path.join(self._submit_directory, file_name)
                dest_file_name = os.path.join(self._scratch_directory, file_name)
                if os.path.isfile(full_file_name):
                    shutil.copy(full_file_name, dest_file_name)

                elif os.path.isdir(full_file_name):
                    if os.path.isdir(dest_file_name):
                        shutil.rmtree(dest_file_name)

                    shutil.copytree(full_file_name, dest_file_name)

            os.chdir(os.path.abspath(self._scratch_directory))
            # Now we have to change the logger output so that we can properly save the output
            # TODO: have it use the formatter from an old handler instead!
            node_logger = logging.FileHandler("./tmpLog", 'a')
            formatter = logging.Formatter("%(asctime)-15s %(levelname)-8s %(message)s")
            node_logger.setFormatter(formatter)
            node_logger.setLevel(logging.INFO)
            # These are the old_handlers, we will save them when we go back
            # We don't remove them so that we can continue to get updates on the node
            old_handlers = logger.handlers[:]
            logger.addHandler(node_logger)

        # We can arm the timer
        if self._time_to_run != -1:
            logger.info("Starting the timer")
            signal.signal(signal.SIGALRM, self.calculation_alarm_handler)
            signal.alarm((self._time_to_run * 60 - 30) * 60)

        # We loop over the steps here and will pop elements off the beginning of the dictionary
        while len(self._commands.values()) != 0:
            logger.info("")
            if self._timer_went_off:
                logger.info("Timer went off, not continuing onto next command")
                break

            # Grab the next step dictionary to do
            steps = self._commands[list(self._commands.keys())[0]]
            logger.debug(f"On step: {steps}")
            updated_parameters = self._raw_parameters.copy()

            for changes in steps:
                logger.debug(f"Updating {changes}: changing {updated_parameters[changes]} to {steps[changes]}")
                updated_parameters[changes] = steps[changes]
            
            start = timer()
            if updated_parameters["titr"]["titr on"]:
                if self._titration is None:
                    logger.warning("Titration feature cannot be turned on in the middle of a run")
                    raise ValueError("Titration turned on")
                
                if os.path.isfile(updated_parameters["Movie File"]):
                    utilities.make_movie("initial.pdb", updated_parameters["Movie File"], "_tmpMovie.pdb")
                    #Append to movie
                    with open("_tmpMovie.pdb", 'r') as tmpMovie, open("movie.pdb", 'a') as movie:
                        for line in tmpMovie:
                            movie.write(line)

                    last_frame = utilities.last_frame("_tmpMovie.pdb")
                    
                    #Clean up our mess
                    logger.debug("Removing _tmpMovie.pdb file")
                    os.remove("_tmpMovie.pdb")
                    logger.debug(f"Removing {updated_parameters['Movie File']} file")
                    os.remove(updated_parameters["Movie File"])
                    
                else:
                    last_frame = utilities.load_pdb("initial.pdb")
                
                if os.path.isfile(updated_parameters['Restart File']):
                    logger.debug(f"Removing {updated_parameters['Restart File']} file")
                    os.remove(updated_parameters['Restart File'])

                #TODO check to see if any of the protonation states are invalids (ie, they affect statically held protonation
                #states defined by the user)
                updated_parameters["Custom protonation states"] = self._titration.evaluate_pkas(last_frame)

                sj = setupDMDjob(parameters=updated_parameters, pro=last_frame)
                
                #This will not do a quick dmd setup, so we should be able to expedite that slightly. Also no topparam file either
                #creates state, start and outConstr, inConstr
                sj.titrate_setup()

                #We don't want to use the restart file, so set last var to False
                self.run_dmd(updated_parameters, self._start_time, False)

            elif "Custom protonation states" in steps.keys():
                logger.warning("Why are you trying to change the protonation state in the middle of DMD?")

            elif "Frozen atoms" in steps.keys() or "Restrict Displacement" in steps.keys():
                logger.warning("Cannot freeze atoms or change displacement between atoms in the middle of a run.")
                logger.warning("This does not make any...ignoring these")

            else:
                # We can just run the job with no issues other than those raised from the above
                self.run_dmd(updated_parameters, self._start_time, True)

            end = timer()

            self.print_summary(updated_parameters['Time'], end-start)

            # Assuming we finished correctly, we pop off the last issue
            self._commands.pop(list(self._commands.keys())[0])
            # Update the new start time!
            self._start_time += updated_parameters["Time"]

        if self._commands:
            logger.info("Did not finish all of the commands, will save the remaining commands")

        else:
            logger.debug("Finished all commands...writing final dmdinput.json")

        logger.debug("Setting remaining commands to the rest of the commands")
       
        if self._titration is not None:
            logger.debug("Condensing any commands remaining from the titratable feature")
            self._raw_parameters = self._titration.condense_commands(self._raw_parameters)

        self._raw_parameters["Remaining Commands"] = self._commands
        with open("dmdinput.json", 'w') as dmdinput:
            logger.debug("Dumping to json")
            json.dump(self._raw_parameters, dmdinput, indent=4)


        if os.path.abspath(self._scratch_directory) != os.path.abspath(self._submit_directory):
            logger.info(
                f"Copying files from {os.path.abspath(self._scratch_directory)} to {os.path.abspath(self._submit_directory)}")
            self._src_files = os.listdir(self._scratch_directory)
            for file_name in self._src_files:
                full_file_name = os.path.join(self._scratch_directory, file_name)
                dest_file_name = os.path.join(self._submit_directory, file_name)
                if os.path.isfile(full_file_name):
                    shutil.copy(full_file_name, dest_file_name)

                # Want to remove and then copy over a directory and everything in it!
                elif os.path.isdir(full_file_name):
                    if os.path.isdir(dest_file_name):
                        shutil.rmtree(dest_file_name)

                    shutil.copytree(full_file_name, dest_file_name)

            os.chdir(os.path.abspath(self._submit_directory))
            # Now we swap back to the initial handlers
            # We want to get rid of all of our old handlers
            for hdlr in logger.handlers[:]:
                logger.removeHandler(hdlr)

            # This appends our temp file to our output file!
            # TODO, allow the ext://sys.stdout to be replaced by whatever the initial logger stream was...
            appendHandler = logging.StreamHandler(sys.stdout)
            formatter = logging.Formatter("%(message)s")
            appendHandler.setFormatter(formatter)
            logger.addHandler(appendHandler)
            try:
                with open('tmpLog') as logFile:
                    for line in logFile:
                        logger.info(line.rstrip())
                os.remove("tmpLog")

            except IOError:
                logger.critical("Error with appending node log file to initial log file")

            # Officially gets rid of everything
            for hdlr in logger.handlers[:]:
                logger.removeHandler(hdlr)

            # Now we go back to our initial handlers
            for hdlr in old_handlers:
                logger.addHandler(hdlr)


    def run_dmd(self, parameters, start_time: int, use_restart: bool):
        # Remake the start file with any changed parameters
        utilities.make_start_file(parameters, start_time)

        if use_restart:
            state_file = self._raw_parameters["Restart File"] if os.path.isfile(self._raw_parameters["Restart File"]) else "state"

        else:
            state_file = "state"

        logger.info(f"[Restart File]     ==>> {'False' if state_file == 'state' else 'True'}")

        #Now we execute the command to run the dmd
        try:
            with open("dmd.out", 'a') as dmd_out:
                logger.info(f"[Issuing command]  ==>> pdmd.linux -i dmd_start -s {state_file} -p param -c outConstr -m {self._cores} -fa")
                logger.info("...")
                with Popen(f"pdmd.linux -i dmd_start -s {state_file} -p param -c outConstr -m {self._cores} -fa",
                        stdout=PIPE, stderr=subprocess.STDOUT, universal_newlines=True, shell=True, env=os.environ) as shell:
                    while shell.poll() is None:
                        dmd_out.write(shell.stdout.readline().strip() + '\n')

        except OSError:
            logger.exception("Error calling pdmd.linux")
            raise

    @staticmethod
    def get_echo_data(echo_file):
        if not os.path.isfile(echo_file):
            logger.error(f"Echo file does not exist: {echo_file}")
            raise FileNotFoundError("Echo File")

        energies = []

        with open(echo_file, 'r') as echo:
            for line in echo:
                if line[0] == "#":
                    continue

                line = line.split()
                energies.append(line)

        return energies

    @staticmethod
    def get_average_potential_energy(echo_file):
        energies = dmd_simulation.get_echo_data(echo_file)

        energies = [float(line[4]) for line in energies]
        ave = sum(energies)/len(energies)

        stdev = 0
        for e in energies:
            stdev += (e - ave)**2

        stdev /= (len(energies)-1)
        stdev = math.sqrt(stdev)

        return [ave, stdev] 

    @staticmethod
    def get_average_kinetic_energy(echo_file):
        energies = dmd_simulation.get_echo_data(echo_file)

        energies = [float(line[5]) for line in energies]
        ave = sum(energies)/len(energies)

        stdev = 0
        for e in energies:
            stdev += (e - ave)**2

        stdev /= (len(energies)-1)
        stdev = math.sqrt(stdev)

        return [ave, stdev] 


    @staticmethod
    def get_average_temp_energy(echo_file):
        energies = dmd_simulation.get_echo_data(echo_file)

        energies = [float(line[1]) for line in energies]
        ave = sum(energies)/len(energies)

        stdev = 0
        for e in energies:
            stdev += (e - ave)**2

        stdev /= (len(energies)-1)
        stdev = math.sqrt(stdev)

        return [ave, stdev] 

    @staticmethod
    def get_average_pressure_energy(echo_file):
        energies = dmd_simulation.get_echo_data(echo_file)

        energies = [float(line[2]) for line in energies]
        ave = sum(energies)/len(energies)

        stdev = 0
        for e in energies:
            stdev += (e - ave)**2

        stdev /= (len(energies)-1)
        stdev = math.sqrt(stdev)

        return [ave, stdev] 

    def print_summary(self, sim_time, wall_time):
        pot_energy = self.get_average_potential_energy(self._raw_parameters['Echo File'])
        kinetic = self.get_average_kinetic_energy(self._raw_parameters['Echo File'])
        pressure = self.get_average_pressure_energy(self._raw_parameters['Echo File'])
        temperature = self.get_average_temp_energy(self._raw_parameters['Echo File'])
        
        logger.info(f"[Ave. Pot. Energy] ==>> {pot_energy[0]:.5f} ({pot_energy[1]:.5f}) kcal/mol")
        logger.info(f"[Ave. Kin. Energy] ==>> {kinetic[0]:.5f} ({kinetic[1]:.5f}) kcal/mol")
        logger.info(f"[Ave. Tot. Energy] ==>> {(pot_energy[0] + kinetic[0]) / 2.0:.5f} kcal/mol")
        logger.info(f"[Ave. Pressure   ] ==>> {pressure[0]:.5f} ({pressure[1]:.5f})")
        logger.info(f"[Ave. Temperature] ==>> {temperature[0]:.5f} ({temperature[1]:.5f})")
        logger.info(f"[Est. Phys. Time ] ==>> {sim_time*0.0000488882} ns")
        logger.info(f"Time elapsed during DMD simulation: {datetime.timedelta(seconds = int(wall_time))}")


    def calculation_alarm_handler(self, signum, frame):
        """
        Called if time is almost up in the dmd job! Will copy all files to a backup directory and write out the
        remaining commands to perform in the remaining_commands.json file.
        """
        logger.warning("Creating a backup directory!")
        self._timer_went_off = True

        logger.info("Placing the remaining commands into remaining_commands.json")
        with open("remaining_commands.json", 'w') as rc:
            json.dump(self._commands, rc)

        logger.debug("Checking if scratch directory is different from submit directory")
        if os.path.abspath(self._scratch_directory) != os.path.abspath(self._submit_directory):
            logger.warning("Creating dmd_backup in the submit directory")
            backup_dir = os.path.join(os.path.abspath(self._submit_directory), 'dmd_backup')
            if os.path.isdir(backup_dir):
                logger.warning("Removing backup directory already present in the submit directory")
                shutil.rmtree(backup_dir)

            logger.warning(f"Copying files from {os.path.abspath(self._scratch_directory)} to {backup_dir}")
            shutil.copytree(self._scratch_directory, backup_dir)

        logger.info("Turning off alarm")
        signal.alarm(0)