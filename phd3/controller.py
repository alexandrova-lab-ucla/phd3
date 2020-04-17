#!/usr/bin/env python3
"""
Author  ==>> Matthew R. Hennefarth
Date    ==>> April 16, 2020
"""

#Standard Library Imports
import logging
import os
import signal
import json
import shutil
import sys

#PHD Imports
import phd3.iteration as iteration
import phd3.bin.submitphd as submitphd

logger = logging.getLogger(__name__)

__all__ = [
        'controller'
    ]

class controller:

    def __init__(self, cores:int, time:int=-1, scratch:str="./"):
        self._cores = cores
        self._time = time
        self._scratch = scratch
        self._submit_directory = os.getcwd()
        self._iteration = 0
        self._parameters = {}
        self._curr_iterat = None

        self._stop = False

        #In case a signal is sent to stop!
        signal.signal(signal.SIGUSR1, self.alarm_handler)

        if not os.path.isfile("phdinput.json"):
            logger.error("PHDinput.json does not exist")
            raise FileNotFoundError("phdinput.jdon")
       
        #TODO verify phdinput.json is valid format
            
        logger.debug("Loading in input file")
        with open("phdinput.json", 'r') as f:
            self._parameters = json.load(f)

        logger.debug("Checking last iteration directory")
        dirs = [d for d in os.listdir() if os.path.isdir(d)]
        iterations = [int(d.split("_")[1]) for d in dirs if "Iteration_" in d]
        if not iterations:
            last_iteration_directory = -1

        else:
            last_iteration_directory = max(iterations)

        logger.debug("Checking last saved iteration")
        self.last_finished_iteration = -1
        if os.path.isfile("phd_energy"):
            energy_lines = []
            with open("phd_energy") as energy:
                for line in energy:
                    energy_lines.append(line)

            #Get rid of the first line
            if len(energy_lines) > 1:
                energy_lines = energy_lines[1:]
                self.last_finished_iteration = len(energy_lines) - 1

        #Then we assume we start the next after last_iteration_directory
        if last_iteration_directory != self.last_finished_iteration:
            #Then we need to start from the last_iteration_directory

            #We should throw an error if this is not the case
            assert(last_iteration_directory > self.last_finished_iteration)
            self._iteration = last_iteration_directory

        else:
            #They are the same, therefore start next iteration
            self._iteration = self.last_finished_iteration + 1

        logger.debug("Finding last to_next_iteration.pdb file")
        if self._iteration > 0 :
            if os.path.isfile(os.path.join(f"Iteration_{self._iteration-1}", "to_next_iteration.pdb")):
                self._parameters["last pdb"] = os.path.abspath(os.path.join(f"Iteration_{self._iteration-1}", "to_next_iteration.pdb"))

            else:
                logger.error("Last to_next_iteration.pdb not found!")
                raise FileNotFoundError("to_next_iteration.pdb")

        else:
            self._parameters["last pdb"] = os.path.abspath(self._parameters["pdb file"])

        #Now we move the necessary files to the new directory!!
        if os.path.abspath(self._scratch) != os.path.abspath(self._submit_directory):
            if not os.path.isdir(self._scratch):
                logger.debug("Making scratch directory")
                os.mkdir(self._scratch)
            
            if self._iteration == last_iteration_directory:
                if os.path.isdir(os.path.join(self._scratch, f"Iteration_{last_iteration_directory}")):
                    shutil.rmtree(os.path.join(self._scratch, f"Iteration_{last_iteration_directory}"))

                shutil.copytree(f"Iteration_{last_iteration_directory}", os.path.join(self._scratch, f"Iteration_{last_iteration_directory}"))
                
            #copy over the phdenergy file
            if os.path.isfile("phd_energy"):
                shutil.copy("phd_energy", self._scratch)

            logger.debug("Changing directory from {os.getcwd()} to {self._scratch}")            
            os.chdir(self._scratch)

        if self._time != -1:
            logger.info("Starting the timer")
            signal.signal(signal.SIGALRM, self.alarm_handler)
            signal.alarm((self._time* 60 - 55) * 60)

        while not self._stop and self._iteration <= self._parameters["Max Iterations"]:
            #This is the loop we stay in until we need to quit
            self._curr_iterat = iteration.iteration(f"Iteration_{self._iteration}", os.path.abspath("."), self._parameters, self._iteration, self._cores) 
            self._curr_iterat.continue_calculation()
            if self._stop:
                self._iteration -=1

            self._iteration += 1

        if self._time != -1:
            logger.info("Turning off timer")
            signal.alarm(0)

        if os.path.abspath(self._scratch) != os.path.abspath(self._submit_directory):
            dirs = [d for d in os.listdir() if "Iteration_" in d]
            for d in dirs:
                if os.path.isdir(os.path.join(self._submit_directory, d)):
                    shutil.rmtree(os.path.join(self._submit_directory, d))

                shutil.copytree(d, os.path.join(self._submit_directory, d))

            if os.path.isfile("phd_energy"):
                shutil.copy("phd_energy", self._submit_directory)
            
            logger.debug("Changing directory from {os.getcwd()} to {self._submit_directory}")            
            os.chdir(self._submit_directory)
        
        if self._stop and self._iteration <= self._parameters["Max Iterations"]:
            if self._parameters["Resubmit"]:
                submitphd.main(_cores = self._cores, _time=self._time)

        elif self._iteration > self._parameters["Max Iterations"]:
            logger.info("Finished with all QM/DMD cyles")

    def alarm_handler(self, signum, frame):
        #Alarm went off
        logger.info("Alarm went off!")
        self._stop = True

        #Propogate down to the current iteration
        self._curr_iterat.signal_alarm()

        if self._scratch != self._submit_directory:
            dirs = [d for d in os.listdir(self._scratch) if "Iteration_" in d]
            for d in dirs:
                if os.path.isdir(os.path.join(self._submit_directory, d)):
                    shutil.rmtree(os.path.join(self._submit_directory, d))

                shutil.copytree(os.path.join(self._scratch, d), os.path.join(self._submit_directory, d))

            if os.path.isfile(os.path.join(self._scratch, "phd_energy")):
                shutil.copy(os.path.join(self._scratch, "phd_energy"), self._submit_directory)
        
        signal.alarm(0)
