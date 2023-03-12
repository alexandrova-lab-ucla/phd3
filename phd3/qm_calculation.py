#!/usr/bin/env python3
"""
Author  ==>> Matthew R. Hennefarth
Date    ==>> April 16, 2020
"""

#Standard Library Impors
import os
import sys
import json
import shutil
import signal
import logging
from timeit import default_timer as timer
import datetime
from subprocess import Popen, PIPE, STDOUT
import time
import tarfile

#PHD3 Imports
from phd3.utility import utilities, exceptions, constants
from phd3.setupjob import setupTMjob 
from phd3.bin import submitturbomole

logger=logging.getLogger(__name__)

__all__ = [
    'TMcalculation'
]

class TMcalculation:
    def __init__(self, cores, run_dir: str='./', time=-1, coord=None, parameters: dict=None):
        logger.debug("Beginning Calculation")

        logger.debug("Initializing variables")
        self._submit_directory = os.path.abspath(os.getcwd())
        self._scratch_directory = os.path.abspath(run_dir)
        self._config = utilities.load_phd_config()
        self._turbodir = self._config["PATHS"]['TURBODIR']
        self._cores = cores
        self._time_to_run = time
        self._para_arch = "SMP"
        self._resub = False
        self._timer_went_off = False
        self._coord_file = ""

        #Start the general signal for sigUSR1
        #This will call function calculation_alarm_handler if SIGUSR1 is sent to the node
        signal.signal(signal.SIGUSR1, self.calculation_alarm_handler)

        # We need to make sure that we have some parameters for this job!
        if parameters is None:
            logger.debug("Checking for definput.json file")
            if os.path.isfile("definput.json"):
                self._parameter_file = os.path.join(self._submit_directory, "definput.json")

            else:
                logger.error("No parameters specified for the job!")
                raise FileNotFoundError("definput.json")

            # read in the parameters here to the dictionary!!
            try:
                # Read in the user supplied parameters from the input file
                logger.debug("Loading in parameters")
                with open(self._parameter_file, 'r') as inputfile:
                    self._raw_parameters = json.load(inputfile)

            except IOError:
                logger.exception("Could not open the parameter file correctly")
                raise

        else:
            logger.debug("Using parameters passed")
            self._raw_parameters = parameters

        # we want to make sure that turbomole environment is setup properly first now in case we have to run define or x2t!
        logger.debug("Setting up TURBOMOLE Environment")
        self.setup_turbomole_env_parallel(self._cores, self._para_arch, self._scratch_directory, self._turbodir)

        #Check to make sure that we have valid parameters before continuing!
        try:
            logger.debug("Ensuring proper parameters")
            utilities.valid_qm_parameters(self._raw_parameters)

        except exceptions.ParameterError:
            logger.error("Invalid parameters provided!")
            raise

        # Now we are in the submit directory, check for a coord file first!
        if coord is None:
            logger.debug("Looking for a coord file")
            if os.path.isfile("coord"):
                logger.debug(f"Found a coord file: {os.path.join(self._submit_directory, 'coord')}")
                self._coord_file = os.path.join(self._submit_directory, "coord")

            #Check to see if there is a .xyz to use instead
            else:
                logger.debug(f"Looking for an .xyz file in: {self._submit_directory}")
                files = [f for f in os.listdir(self._submit_directory) if os.path.isfile(os.path.join(self._submit_directory, f))]
                for f in files:
                    if os.path.splitext(f)[-1].lower() == ".xyz":
                        logger.info(f"Found an .xyz file to use: {f}, will convert to coord")
                        utilities.xyz_to_coord(f)
                        self._coord_file = os.path.join(self._submit_directory, "coord")
                        break

                #No such luck, can't do the job then!
                if self._coord_file == "":
                    logger.error("You don't seem to have provided a coord file, cannot proceed!")
                    raise FileNotFoundError("coord")

        else:
            self._coords = coord

        # At this point we will set up the turbomole job if it is not already setup (or more importantly, the control file!)
        if not os.path.isfile("control"):
            logger.debug("Control file does not exist. Attempting to run setupturbomole.py")
            stm = setupTMjob(self._raw_parameters)

        #Remove the stop file if it exists right now, otherwise we will have issues
        if os.path.isfile(os.path.join(self._submit_directory, "stop")):
            logger.warning("Stop file exists before job has run")
            logger.warning("Removing stop file and continuing")
            os.remove(os.path.join(self._submit_directory, "stop"))

        #We now want to transfer all of the files from this directory to the scratch directory if a diff. directory
        if self._scratch_directory != self._submit_directory:
            logger.info(f"Copying file from {self._submit_directory} to {self._scratch_directory}")
            self._scratch_directory = os.path.join(self._scratch_directory, os.path.basename(self._submit_directory))
            utilities.copy_directories(self._submit_directory, self._scratch_directory)

            os.chdir(self._scratch_directory)

        # Check for MOs in tar.gz format
        for mo_file in constants.MO_FILES:
            if os.path.isfile(f"{mo_file}.tar.gz"):
                logger.info(f"[Untaring] ==>> {mo_file}")
                with tarfile.open(f"{mo_file}.tar.gz", "r:gz") as tar:
                    def is_within_directory(directory, target):
                    	
                    	abs_directory = os.path.abspath(directory)
                    	abs_target = os.path.abspath(target)
                    
                    	prefix = os.path.commonprefix([abs_directory, abs_target])
                    	
                    	return prefix == abs_directory
                    
                    def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
                    
                    	for member in tar.getmembers():
                    		member_path = os.path.join(path, member.name)
                    		if not is_within_directory(path, member_path):
                    			raise Exception("Attempted Path Traversal in Tar File")
                    
                    	tar.extractall(path, members, numeric_owner=numeric_owner) 
                    	
                    
                    safe_extract(tar)

                os.remove(f"{mo_file}.tar.gz")

        switcher = {
            "numforce": self._numforce,
            "nf" : self._numforce,
            "forceopt": self._forceopt,
            "geo": self._geo,
            "singlepoint": self._singlepoint,
            "sp" : self._singlepoint,
            "trans": self._trans,
            "escf": self._escf,
            "woelfling" : self._woelfling,
            "egeo" : self._egeo,
            "eforceopt" : self._eforceopt,
            "enumforce" : self._enumforce,
            "enf" : self._enumforce
        }

        calc = switcher.get(self._raw_parameters["calculation"].lower(), None)
        try:
            logger.debug("Beginning the TURBOMOLE calulcation")

            # Arm the timer
            if self._time_to_run != -1:
                logger.info("Starting the timer")
                signal.signal(signal.SIGALRM, self.calculation_alarm_handler)
                signal.alarm(int((self._time_to_run * 60 - 45) * 60))
            
            start = timer()
            calc()
            end = timer()
            logger.info(f"Time elapsed: {datetime.timedelta(seconds = int(end -start))}")
            if self._time_to_run != -1:
                # Turn off the timer now!
                logger.info("Turning off timer")
                signal.alarm(0)

        except:
            logger.error(f"Failed executing the commands for the job type provided: {self._raw_parameters['calculation']}")
            raise

        #Now we copy the files back!
        if self._scratch_directory != self._submit_directory:
            utilities.copy_directories(self._scratch_directory, self._submit_directory)
            os.chdir(self._submit_directory)
            shutil.rmtree(self._scratch_directory)

        self.clean_up()
        logger.debug("Finished Calculation")

        if self._resub:
            logger.info("Resubmitting the job!")
            submitturbomole.main(_cores=self._cores, _time=self._time_to_run)

    def _woelfling(self):
        logger.debug("Woelfling transition state")
        self._run("frozen_woelfling-job > woelfling.out")
                      
    def _forceopt(self, ex=False):
        """ Executes the commands (jobex) to perform a geometry optimization. Will resubmit if not done.
        """
        logger.debug("Force Optimization!")
        self._raw_parameters['geo_iterations'] = 100000
                         
        if ex:
            self._egeo()
        
        else:                      
            self._geo()
        
        # This is where we do some light error checking
        if os.path.isfile("GEO_OPT_FAILED"):
            with open("GEO_OPT_FAILED") as error_file:
                for line in error_file:
                    if "stop file found" in line:
                        self._resub = True
                        break

                    elif "OPTIMIZATION DID NOT CONVERGE" in line:
                        self._resub = True
                        break

        # Now we determine the fate, do we resubmit? or not...
        if self._resub:
            # Check to see if too many iterations in general now
            try:
                i = 0
                with open("energy") as energy_file:
                    for line in energy_file:
                        i += 1

                if i > 2000:
                    self._resub = False

            except IOError:
                logger.error("No energy file, something is wrong!")
                self._resub = False

    def _eforceopt(self):
        self._forceopt(ex=True)
    
    def _egeo(self):
        logger.debug("Excited state geometry optimization job")
        logger.debug("Checking control file for proper itvc and exopt")
        
        control_lines = []
        with open("control", 'r') as control_file:
            for line in control_file:
                control_lines.append(line)
                         
        with open("control", "w+") as control_file:
            for line in control_lines:
                if "itrvec" in line:
                    logger.debug("Fixing itrvec in control file!")
                    control_file.write("    itrvec    0\n")
           
                else:
                    control_file.write(line)
                                       
        command = f"jobex -ex -c {str(self._raw_parameters['geo_iterations'])}"

        if self._raw_parameters["gcart"] is not None:
            command += f" -gcart {str(self._raw_parameters['gcart'])}"

        command += f" -np {self._cores}"

        self._run(command + " > jobex_ex.out")
                         
    def _geo(self):
        logger.debug("Geometry Optimization Job")
        logger.debug("Checking control file for proper itvc!")

        control_lines = []
        with open("control", 'r') as control_file:
            for line in control_file:
                control_lines.append(line)

        with open("control", "w+") as control_file:
            for line in control_lines:
                if "itrvec" in line:
                    logger.debug("Fixing itrvec in control file!")
                    control_file.write("    itrvec    0\n")

                else:
                    control_file.write(line)

        command = f"jobex -c {str(self._raw_parameters['geo_iterations'])}"
        if self._raw_parameters["rij"]:
            command += " -ri"

        if self._raw_parameters["gcart"] is not None:
            command += f" -gcart {str(self._raw_parameters['gcart'])}"

        command += f" -np {self._cores}"

        self._run(command + " > jobex.out")

    def _singlepoint(self):
        logger.debug("Single Point Job")
        command = "dscf" if not self._raw_parameters["rij"] else "ridft"
        self._run(command + f" -smpcpus {self._cores} > {command}.out")

    def _escf(self):
        logger.debug("Excited state job")
        if not os.path.isfile("energy"):
            logger.warning("Energy not found, you should perform a single point calculation first!")
            logger.warning("I will continue, maybe you bamboozled me...turbomole will decide our fates")

        command = f"escf -smpcpus {self._cores} > escf.out"
        self._run(command)

    def _enumforce(self):
        logger.debug("NumForce Job")
        if not os.path.isdir("numforce"):        
            self._run(f"dscf -smpcpus {self._cores} > dscf.out")
            self._run(f"escf -smpcpus {self._cores} > escf.out")
            self._run(f"egrad -smpcpus {self._cores} > egrad.out")

        elif os.path.isdir("numforce/KraftWerk"):
            logger.debug("Checking to see if Kraftwerk directory is cleaned up")
            kraftwerk_files = os.listdir("numforce/KraftWerk")

            deleteFiles = [f for f in kraftwerk_files if "ENVIRONMENT" in f or "lockhost." in f]
            lockFiles = [f for f in kraftwerk_files if "lock." in f]
            
            jobs = [f.split('.')[1] for f in lockFiles]
            
            deleteFiles.extend(lockFiles)
            # TODO add some try and except clauses around the remove and rmtree
            for e in deleteFiles:
                logger.debug(f"Removing file: {e}")
                os.remove(os.path.join("numforce/KraftWerk", e))

            for j in jobs:
                if os.path.isdir(os.path.join("numforce/KraftWerk", j)):
                    logger.debug(f"Removing directory: {j}")
                    shutil.rmtree(os.path.join("numforce/KraftWerk", j))

                if os.path.isfile(os.path.join("numforce/KraftWerk", j + ".log")):
                    logger.debug(f"Removing file: {j}")
                    os.remove(os.path.join("numforce/KraftWerk", j + ".log"))

                for k in kraftwerk_files:
                    if k.startswith(j + '.') and k.endswith('.err'):
                        if os.path.isfile(os.path.join("numforce/KraftWerk", k)):
                            logger.debug(f"Removing file:{k}")
                            os.remove(os.path.join("numforce/KraftWerk", k))

                        else:
                            logger.debug(f"File: {k}, must have been already removed")
                        
                        break

            logger.debug("Cleaned up the Kraftwerk directory!")

        ex_state = 1
                         
        with open("control") as control_file:
            for line in control_file:
                if line.startswith("$exopt"):
                    ex_state = int(line.split()[-1])
                    break
                         
        command = f"NumForce -ex {ex_state} -central"
        if self._raw_parameters["freeze_atoms"]:
            command += " -frznuclei"

        command += f" -mfile {constants.MFILE} > numforce.out"
        self._run(command)
        if not self._timer_went_off:
            logger.info("Trying to delete KraftWerk directory")
            try:
                shutil.rmtree("numforce/KraftWerk")

            except:
                logger.info("I guess there is no KraftWerk to delete")
                pass
                         
    def _numforce(self):
        """ Exceutes the commands (NumForce) to perform a numforce calculation """
        logger.debug("NumForce Job")
        if not os.path.isdir("numforce"):        
            if self._raw_parameters["rij"]:
                self._run(f"ridft -smpcpus {self._cores} > ridft.out")
                self._run(f"rdgrad -smpcpus {self._cores} > rdgrad.out")
            else:
                self._run(f"dscf -smpcpus {self._cores} > dscf.out")
                self._run(f"grad -smpcpus {self._cores} > grad.out")

        elif os.path.isdir("numforce/KraftWerk"):
            logger.debug("Checking to see if Kraftwerk directory is cleaned up")
            kraftwerk_files = os.listdir("numforce/KraftWerk")

            deleteFiles = [f for f in kraftwerk_files if "ENVIRONMENT" in f or "lockhost." in f]
            lockFiles = [f for f in kraftwerk_files if "lock." in f]
            
            jobs = [f.split('.')[1] for f in lockFiles]
            
            deleteFiles.extend(lockFiles)
            # TODO add some try and except clauses around the remove and rmtree
            for e in deleteFiles:
                logger.debug(f"Removing file: {e}")
                os.remove(os.path.join("numforce/KraftWerk", e))

            for j in jobs:
                if os.path.isdir(os.path.join("numforce/KraftWerk", j)):
                    logger.debug(f"Removing directory: {j}")
                    shutil.rmtree(os.path.join("numforce/KraftWerk", j))

                if os.path.isfile(os.path.join("numforce/KraftWerk", j + ".log")):
                    logger.debug(f"Removing file: {j}")
                    os.remove(os.path.join("numforce/KraftWerk", j + ".log"))

                for k in kraftwerk_files:
                    if k.startswith(j + '.') and k.endswith('.err'):
                        if os.path.isfile(os.path.join("numforce/KraftWerk", k)):
                            logger.debug(f"Removing file:{k}")
                            os.remove(os.path.join("numforce/KraftWerk", k))

                        else:
                            logger.debug(f"File: {k}, must have been already removed")
                        
                        break

            logger.debug("Cleaned up the Kraftwerk directory!")

        command = "NumForce -central"
        if self._raw_parameters["freeze_atoms"]:
            command += " -frznuclei"

        if self._raw_parameters["rij"]:
            command += " -ri"

        command += f" -mfile {constants.MFILE} > numforce.out"
        self._run(command)
        if not self._timer_went_off:
            logger.info("Trying to delete KraftWerk directory")
            try:
                shutil.rmtree("numforce/KraftWerk")

            except:
                logger.info("I guess there is no KraftWerk to delete")
                pass

    def _trans(self):
        """Checks for a Hessian file and then runs a jobex -trans job"""
        logger.debug("Transition State Search Job")

        logger.debug("Checking control file for proper itvc!")
        control_lines = []
        with open("control", 'r') as control_file:
            for line in control_file:
                control_lines.append(line)

        with open("control", "w+") as control_file:
            for line in control_lines:
                if "itrvec" in line:
                    logger.debug("Fixing itrvec in control file!")
                    control_file.write(f"   itrvec    {self._raw_parameters['stp']['itvc']}\n")
                else:
                    control_file.write(line)

        if not os.path.isfile("hessian"):
            logger.warning("Hessian not found, have to perform a numforce calculation first!")
            logger.warning("Continuing with the job to see if you bamboozled me")

        command = f"jobex -trans -c {str(self._raw_parameters['geo_iterations'])}"
        if self._raw_parameters["rij"]:
            command += " -ri"

        if self._raw_parameters["gcart"] is not None:
            command += f" -gcart {str(self._raw_parameters['gcart'])}"

        self._run(command + f" -np {self._cores} > trans.out")

    def _run(self, command):
        logger.info(f"[Issuing command] ==>> {command}")
        logger.info("...")
        with Popen(command, stdin=PIPE, stdout=PIPE, stderr=STDOUT, env=os.environ, shell=True, universal_newlines=True) as shell:
            while shell.poll() is None:
                logger.debug(shell.stdout.readline().strip())


    def calculation_alarm_handler(self, signum, frame):
        logger.warning("Creating stop file!")
        self._timer_went_off = True

        if self._scratch_directory != self._submit_directory:
            # This is just in case TURBOMOLE can't stop in time, we just copy everything over as a backup
            logger.warning("Creating a trun_backup in the submit directory")
            backup_dir = os.path.join(self._submit_directory, 'trun_backup')
            if os.path.isdir(backup_dir):
                logger.warning("Removing backup directory already present in submit directory")
                shutil.rmtree(backup_dir)
                         
            logger.warning(f"Copying files from {self._scratch_directory} to {backup_dir}")             
            try:
                shutil.copytree(self._scratch_directory, backup_dir)

            except:
                logger.warn("Roger, can't build here")
                time.sleep(5)
                if os.path.isdir(backup_dir):
                    shutil.rmtree(backup_dir)
                
                shutil.copytree(self._scratch_directory, backup_dir)
    

        with open("numforce/stop" if self._raw_parameters["calculation"].lower() == "numforce" else "stop" , 'w+') as stopfile:
            pass
        
        logger.info("Turning off alarm")
        signal.alarm(0)

    def clean_up(self):
        logger.debug("Cleaning up directory")
        files = [f for f in os.listdir('.') if os.path.isfile(os.path.join('.', f))]
        for f in files:
            if "slave" in f:
                try:
                    os.remove(f)

                except OSError:
                    logger.error(f"Could not remove file: {f}")

        if os.path.isdir('trun_backup'):
            logger.debug("Removing backup directory")
            shutil.rmtree('trun_backup')

    @staticmethod
    def get_energy(cycle=None, energyfile="energy"):
        if not os.path.isfile(energyfile):
            logger.error("No energy file")
            raise FileNotFoundError("energy")
        
        e_lines = []
        with open(energyfile, 'r') as e:
            for line in e:
                if "$" != line[0]:
                    e_lines.append(line)

        if cycle is None:
            step = e_lines[-1]

        else:
            try:
                step = e_lines[cycle-1]
            
            except IndexError:
                logger.exception("Invalid cycle number")
                raise

        step = step.split()
        return float(step[1])
    
    @staticmethod
    def create_MFILE(cores: int):
        logger.debug(f"Creating MFILE for {cores} cores")

        host_file_lines = ""

        #Check if we are on an SGE engine:
        if os.environ.get("PE_HOSTFILE") is not None:
            logger.debug("On SGE system, using $PE_HOSTFILE")
            with open(os.environ["PE_HOSTFILE"], 'r') as hostfile:
                for line in hostfile:
                    line = line.split()
                    try:
                        for i in range(int(line[1])):
                            host_file_lines += (line[0] + '\n')

                    except ValueError:
                        logger.error("Error in generating HOSTFILE")
                        raise

        #Check to see if we are on a SLURM system
        elif os.environ.get("SLURM_JOB_NODELIST") is not None:
            logger.debug("On SLURM system, generating MFILE")
            n_nodes = ""
            with Popen("scontrol show hostnames $SLURM_JOB_NODELIST", universal_newlines=True, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, bufsize=1) as shell:
                while shell.poll() is None:
                    n_nodes += shell.stdout.readline().strip()
            try:
                t_nodes = int(os.environ["SLURM_TASKS_PER_NODE"])

            except ValueError:
                logger.error("Error in getting tasks per node")
                raise

            n_nodes = n_nodes.split()

            for node in n_nodes:
                for i in range(t_nodes):
                    host_file_lines += (node + '\n')

        else:
            hostname = ""
            with Popen('hostname', universal_newlines=True, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, bufsize=1) as shell:
                while shell.poll() is None:
                    hostname += shell.stdout.readline().strip()
            logger.debug(f"Hostname: {hostname}")
            for i in range(cores):
                host_file_lines += (hostname + '\n')

        logger.debug("Writing out MFILE")
        try:
            with open(constants.MFILE, 'w+') as file:
                file.write(host_file_lines)

        except IOError:
            logger.error("IOError in writing MFILE")
            raise

    @staticmethod
    def setup_turbomole_env_parallel(cores: int, para_arch: str, scratch_directory: str, turbodir: str):
        TMcalculation.create_MFILE(cores)
        
        if "TURBOTMPDIR" in os.environ and "PARA_ARCH" in os.environ and "PARNODES" in os.environ:
            if os.environ["TURBOTMPDIR"] == scratch_directory and os.environ["PARA_ARCH"] == para_arch:
                if os.environ['PARNODES'] == str(cores):
                    #Already setup env. once...no need to again
                    os.environ["HOSTS_FILE"] = "./" + constants.MFILE
                    return

        logger.info(">>>> Setting up TURBOMOLE Environment >>>>")
        logger.info(f"[Cores]           ==>>   {str(cores)}")
        logger.info(f"[PARA_ARCH]       ==>>   {para_arch}")
        logger.info(f"[LOCAL_DIR]       ==>>   {scratch_directory}")
        logger.info(f"[TM_PAR_FORK]     ==>>   on")
        logger.info("")

        os.environ["TURBODIR"] = turbodir
        os.environ["PATH"] += os.pathsep + turbodir + '/scripts'
        os.environ["PATH"] += os.pathsep + turbodir

        os.environ["TURBOTMPDIR"] = scratch_directory

        os.environ["PARA_ARCH"] = para_arch
        os.environ["PARNODES"] = str(cores)
        os.environ["HOSTS_FILE"] = "./" + constants.MFILE
        # lets see if it still runs normally with this turned off
        os.environ["TM_PAR_FORK"] = "on"

        sysname = ""
        with Popen('sysname', universal_newlines=True, shell=True,
                   stdin=PIPE, stdout=PIPE, stderr=PIPE, bufsize=1) as shell:
            while shell.poll() is None:
                sysname += shell.stdout.readline().strip()
        os.environ["PATH"] += os.pathsep + turbodir + f'/bin/{sysname}'

    @staticmethod
    def scfiterfail():
        if not os.path.isfile("GEO_OPT_FAILED"):
            return False

        with open("GEO_OPT_FAILED", 'r') as errorfile:
            for line in errorfile:
                if "ERROR: your energy calculation did not converge" in line:
                    return True
        
        return False


