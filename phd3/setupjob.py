#!/usr/bin/env python3
"""
Author  ==>> Matthew R. Hennefarth
Date    ==>> April 16, 2020
"""

#Standard Library Imports
from subprocess import Popen, PIPE
import subprocess
import os
import json
import signal
import logging
import pkg_resources
import shutil
import random
import numpy as np

#PHD3 Imports
from phd3.utility import utilities, exceptions, constants
import phd3.protein.protein as protein
from phd3 import dmd_to_qm

__all__ = [
    'setupTMjob',
    'setupDMDjob',
    'setupPHDjob'
]

logger = logging.getLogger(__name__)

class setupTMjob:
    # states in which they want '\n' instead of '*\n' to exit
    blank_quit = ['scf', 'ff', 'internal', 'control', 'title', 'marij']

    _defineresponses = {  # lines that trigger change in state
        # stored as "phrase from define" : (lines till response, where we are in define)
        "IF YOU WANT TO READ DEFAULT-DATA FROM ANOTHER control-TYPE FILE,": (1, "control"),
        "INPUT TITLE OR": (1, "title"),
        "HIT >return< TO ACCEPT DEFAULT TITLE OR": (1, "title"),
        "TERMINATE MOLECULAR GEOMETRY SPECIFICATION": (5, "geometry"),
        "HIT >return< TO CONFIRM REMOVAL OF THESE INTERNAL COORDINATE": (2, "continue"),
        "IF YOU DO NOT WANT TO USE INTERNAL COORDINATES ENTER  no": (1, "cartesians"),
        "ENTER COMMAND OR HIT >return< TO GET BACK TO GEOMETRY MAIN MENU": (1, "internal"),
        "ENTER INTERNAL COORDINATE DEFINITION COMMAND": (12, "idef"),
        "BOND ANALYSIS MAIN MENU": (16, "iaut"),

        "GOBACK=& (TO GEOMETRY MENU !)": (1, "basis"),

        "THE COMMANDS  use  OR  eht  OR  *  OR q(uit) TERMINATE THIS MENU": (2, "occupation"),

        "DO YOU WANT THE DEFAULT PARAMETERS FOR THE EXTENDED HUECKEL CALCULATION ?": (2, "continue_y"),
        "ENTER THE MOLECULAR CHARGE": (1, "charge"),
        "ENTER THE ATOMIC CHARGE": (1, "charge"),
        "DO YOU WANT THE DEFAULT OCCUPATION ASSIGNMENT FOR ATOMS ?": (1, "continue_n"),

        "2 SUITED DEFINITIONS OF ATOMIC ORBITALS FOR cu": (6, "continue_y"),

        "DO YOU ACCEPT THIS OCCUPATION ?": (1, "open_shell"),
        "OCCUPATION NUMBER ASSIGNMENT MENU": (22, "open_shell"),
        "DO YOU REALLY WANT TO WRITE OUT NATURAL ORBITALS ?": (1, "continue_n"),

        "GO BACK TO OCCUPATION/ORBITAL ASSIGNMENT MENU": (1, "general"),
        "on:   TO SWITCH ON  DFT": (1, "dft_off"),
        "off:  TO SWITCH OFF DFT": (1, "dft_on"),
        "on: TO SWITCH ON  RI": (1, "rij_off"),
        "off: TO SWITCH OFF RI": (1, "rij_on"),
        "threshold for multipole neglect": (2, "marij"),
        "old :  to switch DFT-D2 correction on": (1, "dsp"),
        "change TRUST RADIUS": (2, "stp"),
        "0  for TS search": (0, "stp"),

        "SPIN ORBIT GENERALIZED SCF": (0, "scf"),
        "ENTER DESIRED ACCURACY OF SCF-ENERGY": (3, "scf_conv"),
        "ENTER NEW VALUE FOR MAXIMUM NUMBER OF SCF-ITERATIONS": (1, "scf_iter"),
        "electrostatic field definition menu": (13, "efield"),
        "*** specification of electrostatic field(s) ***": (4, "efield_specification"),
        "DO YOU STILL WANT TO CALCULATE THE PSEUDO-INVERSE B-MATRIX": (4, "continue_n"),
        "TO CONTINUE, ENTER <return>": (0, "continue"),
        "HIT >return< TO CONFIRM REMOVE OF THESE INTERNAL COORDINATE": (2, "continue"),
        "DO YOU WANT TO SWITCH OFF AUTOMATIC SHIFTING OF VIRTUAL ORBITALS ? DEFAULT=n": (0, "continue"),
        "CURRENTLY NO CLOSED SHELL SHIFT WILL BE APPLIED": (1, "scf_shift"),
        "ENTER START VALUE FOR DAMPING": (0, "scf_damp"),
        "ENTER INCREMENT FOR REDUCTION OF DAMPING": (0, "continue"),
        "ENTER MINIMAL DAMPING": (0, "continue"),
        "DO YOU WANT TO SWITCH IT ON AGAIN ?  DEFAULT=y": (0, "continue"),
        "ENTER VALUE FOR CLOSED SHELL SHIFT OR ENTER 0 TO SWITCH": (1, "scf_shift")
    }

    _timeout_err = "timeout error"
    _errorresponses = {  # lines that will trigger an error
        "M-1: Attention! Not enough linearly independent coord": (0, "ired error"),
        "DO YOU STILL WANT TO CALCULATE THE PSEUDO-INVERSE B-MATRIX": (5, "iaut error"),
        "ONLY AN INCOMPLETE SET OF    0 INTERNALS SPECIFIED": (4, "iaut error"),
        "nick    - REPEAT NICKNAME INPUT": (8, "basis set name error"),
        "CHOSEN OCCUPATION IMPOSSIBLE": (0, "multiplicity error"),
        "NO MORE DATA AVAILABLE FOR EHT !!": (3, "eht error"),
        "Please define more than one fragment!": (3, "frag error"),
        "NO ATOMS, NO MOLECULE, NOTHING !": (3, "no atoms error"),
        "WARNING! Lowest Eigenvalue": (0, "ired error"),
        "ired failed, retry with cartesian coords": (0, "ired error"),
        "REDCOR: STRANGE LOCAL GEOMETRY?": (0, "ired error"),
        "YOU DID NOT DEFINE A COMPLETE SET OF INTERNAL": (8, "iaut error"),
        "THE B*m*Bt-MATRIX IS SINGULAR !": (8, "iaut error"),
        "CARTESIAN COORDINATES AND VALUES OF INTERNAL COORDINATES DO  N O T   AGREE !": (3, "coord file has internals"),
        "!!! CHOSEN OCCUPATION IMPOSSIBLE !!!": (26, "Impossible Occupation"),
        _timeout_err: (1, _timeout_err)
    }
    _runFile = "tcommands.sh"  # name of file with the turbomole commands
    _define_input_file = "definput"  # input file name with user provided parameters

    def __init__(self, parameters: dict = None, directory="./", timeout: int = 15):
        """Initialized the object, overall main function for setting up a turbomole job

        Keyword arguments
        :param parameters: Can be passed a pre-loaded dictionary of parameters (default None)
        :param dir: Directory to setup the turbomole job in (default is current)
        """
        logger.debug("Initializing variables")
        # Initialize instance vars
        self.timeout = 15
        self._raw_parameters = {}  # user provided parameters
        self._initial_directory = ""  # Current working directory
        self._run_directory = ""
        self._MINN = ""  # stores info about minn functions being used
        self._shell = None  # shell that runs define commands
        self._define_logger = None  # Logger for define
        self._define_state = ""  # state in which define is in
        self._errstate = ''  # errstate of define, if any
        self._state_responses = {  # possible states we can be in within define
            "control": [],
            "title": [],
            "geometry": ["a coord"],
            "internal": [],
            "idef": [],
            "iaut": [],
            "cartesians": [],
            "basis": [],
            "charge": [],
            "occupation": [],
            "open_shell": [],
            "general": [],
            "dft_on": [], "dft_off": [],
            "rij_on": [], "rij_off": [],
            "marij": [],
            "dsp": [],
            "stp": [],
            "efield": [], "efield_specification": [],
            "scf": [], "scf_conv": [], "scf_iter": [], "scf_shift": [], "scf_damp": []
        }

        # Change directory to that with everything we need in it
        self._initial_directory = os.getcwd()
        logger.debug(f"Set initial directory to: {self._initial_directory}")
        self._run_directory = directory
        logger.debug(f"Set run directory to: {self._run_directory}")

        try:
            logger.debug(f"Changing to run directory: {self._run_directory}" )
            os.chdir(self._run_directory)

        except OSError:
            logger.exception(f"Failed moving to run directory: {self._run_directory}")
            raise

        if parameters is not None:
            try:
                logger.debug("Using parameter dictionary passed to setupTMjob")
                self._create_define_options(parameters)

            except exceptions.ParameterError:
                logger.error("Parameter values in dictionary passed to setupTMjob invalid")
                raise

        else:
            # Want to make sure that we have a parameters file
            logger.debug("Searching for definput.json file")
            if os.path.isfile(self._define_input_file + ".json"):
                self._define_input_file = self._define_input_file + ".json"
                try:
                    logger.debug("Found definput.json, loading in variables")
                    self._create_define_options()

                except exceptions.ParameterError:
                    logger.error(f"definput.json in directory: {self._run_directory} invalid parameters")
                    raise

                except IOError:
                    logger.error(f"Could not open definput.json in: {self._run_directory}")
                    raise

                except ValueError:
                    logger.error(f"Invalid json formatting for definput.json in: {self._run_directory}")
                    raise

            else:
                # Copy over a sample.json file
                logger.warn("Could not find a defjinput.json file")
                logger.info("Copying over a default definput.json file")
                try:
                    shutil.copy(pkg_resources.resource_filename('phd3.resources', 'definput.json'), './')

                except OSError:
                    logger.exception(f"Could not copy default definput.json from: {os.path.dirname(os.path.abspath(__file__))} to: {os.path.abspath(self._run_directory)}")
                    raise

                raise FileNotFoundError("definput.json")

        # Make sure that we have some coordinates/geometry to work with
        logger.debug("Checking for coord file")
        if not os.path.isfile("coord"):
            logger.debug(f"No coord file found in: {self._run_directory}")
            files = [f for f in os.listdir('./') if os.path.isfile(os.path.join('./', f))]
            logger.debug(f"The files in {self._run_directory} are: ")
            logger.debug(files)
            for f in files:
                if os.path.splitext(f)[-1].lower() == ".xyz":
                    logger.debug(f"Found an .xyz file to use: {f}, will convert to coord")
                    with Popen(f'x2t {f} > coord', shell=True, universal_newlines=True, env=os.environ) as shell:
                        while shell.poll() is None:
                            pass
                    break

            if not os.path.isfile("coord"):
                logger.error("There is no structure file (coord or .xyz). Cannot setup job!")
                raise FileNotFoundError("coord")

        # Check to see if a control file already exists, if so ABORT!
        if os.path.isfile("control"):
            logger.error("Control file already exists, proceeding will lead to errors")
            raise FileExistsError("control")

        # Begin to setup the job
        logger.debug("Beginning to setup up the job")
        try:
            self.predefine_process()

        except IOError:
            logger.error("IO Error encountered in predefine processing")
            raise

        self.execute_define()

        try:
            self.postdefine_process()

        except ValueError:
            logger.error("Value Error in post-define processing")
            raise

        except IOError:
            logger.error("IO Error in post-define processing")
            raise

        except IndexError:
            logger.error("Index Error in post-define processing")
            raise

        logger.debug("Finished setting up the job")
        
        logger.debug(f"Moving back to initial directory: {self._initial_directory}")
        try:
            os.chdir(self._initial_directory)
        except OSError:
            logger.exception(f"Could not return to initial directory: {self._initial_directory}")
            raise

        logger.info("[setup turbomole] ==>> SUCCESS")

    def _create_define_options(self, parameters: dict = None):
        """Loads in the parameters and pipes the appropriate commands
        into the corresponding dictionary

        :param parameters: parameters already loaded as a dictionary elsewhere (default None)
        :return: True if successful, false if it failed somewhere
        """
        if parameters is None:
            try:
                # Read in the user supplied parameters from the input file
                logger.debug(f"Trying to open file: {self._define_input_file}")
                with open(self._define_input_file, 'r') as inputfile:
                    self._raw_parameters = json.load(inputfile)

                logger.debug(f"Successfully loading in contents from: {self._define_input_file}")

            except IOError:
                logger.exception(f"Could not open file: {self._define_input_file} successfully")
                raise

            except ValueError:
                logger.exception("Error in file formatting")
                raise

        else:
            logger.debug(f"Using provided dictionary of parameters")
            self._raw_parameters = parameters
            # Saves a copy of the definput.json to the current directory for runturbomole to work properly
            try:
                logger.debug("Writing out parameters to definput.json")
                with open(f"{self._define_input_file}.json", 'w') as out:
                    json.dump(parameters, out, indent=4)

            except IOError:
                logger.warning("Could not write out parameters to definput.json")

        try:
            logger.debug("Checking if parameters provided are valid")
            utilities.valid_qm_parameters(self._raw_parameters)

        except exceptions.ParameterError:
            logger.error("Invalid parameter values provided")
            raise

        # Parse the parameters into the appropriate during/post define procedures
        # Geometry Specification Menu settings
        logger.debug("Parsing parameters to properly respond to define")
        logger.debug("Parsing geometry parameters")
        if self._raw_parameters["geometry"]["idef"]["idef_on"]:
            logger.debug("Turning on idef")
            self._state_responses["geometry"].append("idef")
            for bond in self._raw_parameters["geometry"]["idef"]["freeze_stretch"]:
                (atom1, atom2) = bond.split(',')
                logger.debug(f"Freezing bond between {atom1} {atom2}")
                self._state_responses["idef"].append(f"f stre {atom1} {atom2}")

            self._state_responses["internal"].append("")

        if self._raw_parameters["geometry"]["ired"]:
            logger.debug("Using ired")
            self._state_responses["geometry"].append("ired")

        elif self._raw_parameters["geometry"]["iaut"]["iaut_on"]:
            logger.debug("Using iaut")
            self._state_responses["geometry"].append("iaut")
            for bond in self._raw_parameters["geometry"]["iaut"]["bonds"]:
                (atom1, atom2) = bond.split(',')
                logger.debug(f"Adding bond: {atom1}-{atom2}")
                self._state_responses["iaut"].append(f"{atom1}-{atom2}")

            self._state_responses["internal"].extend(["imet", "irem d"])

        elif self._raw_parameters["geometry"]["cartesians"] and not self._raw_parameters["geometry"]["idef"]["idef_on"]:
            logger.debug("Turning on cartesian coordinates")
            self._state_responses["cartesians"].append("no")

        elif self._raw_parameters["geometry"]["cartesians"] and self._raw_parameters["geometry"]["idef"]["idef_on"]:
            logger.warning("You cannot do idef and cartesians, your idef won't work")
            logger.warning("Trying to do internal redundants instead")
            self._state_responses["geometry"].append("ired")

        else:
            logger.warning("You did not specify how you want to define the coordinates.")
            logger.warning("Attempting to use ired")
            self._state_responses["geometry"].append("ired")

        # Atom Attribute Definition Menu settings
        logger.debug("Parsing basis set")
        if self._raw_parameters["basis"]:
            for atom in self._raw_parameters["basis"].keys():
                if atom == "all":
                    logger.debug(f"Setting basis for all atoms to: {self._raw_parameters['basis'][atom]}")
                    self._state_responses["basis"].append(f"b {atom.lower()} {self._raw_parameters['basis'][atom]}")

                else:
                    logger.debug(f"Setting basis set for {atom.lower()} to: {self._raw_parameters['basis'][atom]}")
                    self._state_responses["basis"].append(
                        f"b \"{atom.lower()}\" {self._raw_parameters['basis'][atom]}")

        # Occupation menu settings
        logger.debug("Writing MOs to binary")
        self._state_responses["occupation"].append("atb")
        logger.debug("Switched on eht")
        self._state_responses["occupation"].append("eht")

        # Molecular Charge
        logger.debug(f"Molecular charge set to: {str(self._raw_parameters['charge'])}")
        self._state_responses["charge"].append(str(self._raw_parameters['charge']))

        # Open Shell parameters
        logger.debug("Parsing open shell parameters")
        if self._raw_parameters["open_shell"]["open_shell_on"]:
            logger.debug("Turning on open shell")
            logger.debug(f"Setting number of unpaired electrons to {str(self._raw_parameters['open_shell']['unpaired'])}")
            self._state_responses["open_shell"].extend(
                ["n", f"u {str(self._raw_parameters['open_shell']['unpaired'])}"])

        else:
            logger.debug("Using closed shell system")
            self._state_responses["open_shell"].append("y")

        # General Menu parameters
        logger.debug("Parsing general parameters")
        
        # DFT parameters
        logger.debug("Parsing DFT parameters")
        if self._raw_parameters["dft"]["dft_on"]:
            logger.debug("Turning DFT parameters")
            self._state_responses["general"].append("dft")
            self._state_responses["dft_off"].append("on")
            # Deal with the Minnesota Funcational
            if self._raw_parameters["dft"]["func"] in constants.MINN_FUNCS:
                self._MINN = self._raw_parameters["dft"]["func"]
                logger.debug(f"Using a Minnesota Functional: {self._MINN}")
                logger.debug("Using a dummy of tpss")
                self._raw_parameters["dft"]["func"] = "tpss"

            self._state_responses["dft_on"].extend(
                [f"func {self._raw_parameters['dft']['func']}", f"grid {self._raw_parameters['dft']['grid']}"])

        else:
            self._state_responses["dft_on"].append("off")

        # SCF parameters
        if self._raw_parameters["scf"]:
            self._state_responses["general"].append("scf")
            self._state_responses["scf"].append("iter")
            logger.debug(f"Setting scf iteration limit to: {str(self._raw_parameters['scf']['iter'])}")
            self._state_responses["scf_iter"].append(str(self._raw_parameters["scf"]["iter"]))
            self._state_responses["scf"].append("conv")
            logger.debug(f"Setting scf convergence to: f{str(self._raw_parameters['scf']['conv'])}")
            self._state_responses["scf_conv"].append(str(self._raw_parameters["scf"]["conv"]))
            if "orbital shift" in self._raw_parameters["scf"].keys():
                self._state_responses["scf"].append("shift")
                self._state_responses["scf_shift"].append(str(self._raw_parameters["scf"]["orbital shift"]))
            
            if "damp start" in self._raw_parameters["scf"].keys():
                self._state_responses["scf"].append("damp")
                self._state_responses["scf_damp"].append(str(self._raw_parameters["scf"]["damp start"]))
        
        # rij parameters
        if self._raw_parameters["rij"]:
            logger.debug("Switching on rij")
            self._state_responses["general"].append("ri")
            self._state_responses["rij_off"].append("on")

        # marij parameters
        if self._raw_parameters["marij"]:
            logger.debug("Switching on marij")
            self._state_responses["general"].append("marij")

        # dsp parameters
        if self._raw_parameters["dsp"]:
            logger.debug("Switching on dispersion correction")
            self._state_responses["general"].append("dsp")
            self._state_responses["dsp"].append("on")

        # stp parameters
        if self._raw_parameters["stp"]:
            self._state_responses["general"].append("stp")
            self._state_responses["stp"].append("itvc")
            logger.debug(f"Setting itvc to: {str(self._raw_parameters['stp']['itvc'])}")
            self._state_responses["stp"].append(str(self._raw_parameters['stp']['itvc']))
            logger.debug(f"Setting trust radius to: {str(self._raw_parameters['stp']['trad'])}")
            self._state_responses["stp"].append(f"trad {str(self._raw_parameters['stp']['trad'])}")

        if "efield" in self._raw_parameters.keys():
            if self._raw_parameters["efield"]["on"]:
                logger.debug("Electric Field is being applied to the system")
                self._state_responses["general"].append("e")
                self._state_responses["efield"].extend(["geofield", "man", "*"])

                assert(len(self._raw_parameters["efield"]["direction"])==3)
                for component in self._raw_parameters["efield"]["direction"]:
                    assert(type(component)==float)

                logger.debug(f"Direction: {self._raw_parameters['efield']['direction']}")

                assert(type(self._raw_parameters['efield']['magnitude']) == float)
                logger.debug(f"Magnitude: {self._raw_parameters['efield']['magnitude']}")

                self._state_responses["efield_specification"].append(f"{self._raw_parameters['efield']['direction'][0]} {self._raw_parameters['efield']['direction'][1]} {self._raw_parameters['efield']['direction'][2]} {self._raw_parameters['efield']['magnitude']}")


        logger.debug("Finished parsing parameters")

    ###################################################################################
    #                E X T E R N A L   D E F I N E   C O M M A N D S                  #
    ###################################################################################

    @staticmethod
    def freeze_coords(atoms):
        """within the coord file, it adds an 'f' at the end of each atom line so
        that TURBOMOLE knows to freeze the atom's cartesians

        :param atoms: List of atomnumbers [<atomnum>, <atomnum>,...]
        :return: True if successful, False if it failed
        """
        lines = []

        try:
            logger.debug("Opening coord file")
            with open("coord") as coord:
                for filelines in coord:
                    lines.append(filelines)

        except IOError:
            logger.exception("Could not open coord file")
            raise

        for i in range(len(lines)):
            if i in atoms:
                tmp = list(filter(None, lines[i].split()))
                if len(tmp) == 4:  # Checks to see if only x,y,z and atom type (ie not adding another f)
                    logger.debug(f"Freezing atom: {tmp}")
                    lines[i] = lines[i][:len(lines[i]) - 1] + " f\n"

                else:
                    logger.debug(f"Atom already frozen: {tmp}")

        try:
            logger.debug("Saving coord file")
            with open("coord", 'w') as tmp:
                for line in lines:
                    tmp.write(line)

        except IOError:
            logger.exception("Could not save coord file with frozen cartesians")
            raise

    @staticmethod
    def cosmo(epsilon, tmp):
        """Sticks in a line of the control file:
            $cosmo
                epsilon=<epsilon>

        :param epsilon: dielectric constant
        :param tmp: temporary control lines
        :return: True if successful, false otherwise
        """
        if epsilon < 1.0:
            logger.error("You specified a negative value for epsilon!!")
            raise ValueError

        for line in tmp:
            if "$end" in line:
                x = tmp.index(line)
                tmp.insert(x, f"$cosmo\n   epsilon={str(epsilon)}\n")
                return

        logger.error("Could not add cosmo to control file!")
        raise IndexError

    @staticmethod
    def weight(tmp):
        """Sticks in weight derivatives into the control file right after gridsize

        :param tmp: temporary control lines
        :return: True if successful, false otherwise
        """
        for lines in tmp:
            if "gridsize" in lines:
                x = tmp.index(lines)
                tmp.insert(x + 1, "   weight derivatives\n")

    ###################################################################################
    #                               B U I L D   J O B                                 #
    ###################################################################################

    def predefine_process(self):
        """ Performs any pre-define necessary processing

        :return: True if successful, false if an error occurred
        """
        logger.debug("Beginning pre-define processing")
        if self._raw_parameters["freeze_atoms"]:
            logger.debug("Freezing atom cartesians")
            atoms = self._raw_parameters["freeze_atoms"]
            try:
                self.freeze_coords(atoms)

            except IOError:
                logger.error("Exception encountered in freezing atom cartesians")
                raise

        logger.debug("Finished pre-define processing")

    def execute_define(self):
        """ Executes 'define' with supplied parameters

        :return: True if ran successfully, false if any error occurred, or if there is a state that is unfamiliar
        """
        logger.debug("Starting define")
       
        #Create the define logger
        if logger.isEnabledFor(logging.INFO):
            self._define_logger = logging.getLogger(__name__ + ".define")
            define_out = logging.FileHandler("define.out", 'w+')
            self._define_logger.addHandler(define_out)
            self._define_logger.setLevel(logging.INFO)
            self._define_logger.propagate = False

        else:
            self._define_logger = logger

        turbodir = utilities.load_phd_config()["PATHS"]['TURBODIR']
        utilities.setup_turbomole_env(turbodir)

        self._shell = Popen("define", shell=True, universal_newlines=True, stdin=PIPE, stdout=PIPE, stderr=PIPE,
                            bufsize=1, env = os.environ)

        logger.debug("Initializing timer")
        signal.signal(signal.SIGALRM, self.define_alarm_handler)

        while self._read_define_line():
            if self._define_state == "continue":
                self._write("")

            elif self._define_state == "continue_y":
                self._write("y")

            elif self._define_state == "continue_n":
                self._write("n")

            elif self._define_state in self._state_responses.keys():
                if len(self._state_responses[self._define_state]) != 0:
                    self._write(self._state_responses[self._define_state].pop(0))

                else:
                    self._write('*' if self._define_state not in self.blank_quit else '')

            else:
                self._define_logger.error(f"UNKNOWN DEFINE STATE: {self._define_state}")

        if self._errstate != "all done":
            self._define_logger.error(f"Error: {self._errstate}\n")
            logger.error(f"[Define]          ==>> ERROR TERMINATION ({self._errstate})")
            
            if self._errstate == "timeout error":
                if os.path.isfile("control"):
                    os.remove("control")
            
            raise exceptions.DefineError

        logger.info("[Define]          ==>> Ended Normally")
        logger.debug("Finished define")


    def postdefine_process(self):
        """ Adds any external options that are not available through define
        This includes cosmo, other functionals, weights, etc.

        :return: True if successful, false if an error occured
        """
        logger.debug("Beginning post-define processing")
        controllines = []
        try:
            logger.debug("Opening the control file")
            with open("control", 'r') as control:
                for line in control:
                    controllines.append(line)

        except IOError:
            logger.exception("Error reading in the control file for post-define editing")
            raise

        if self._raw_parameters["cosmo"] is not None:
            logger.debug("Adding cosmo")
            try:
                self.cosmo(self._raw_parameters["cosmo"], controllines)

            except ValueError:
                logger.exception("Invalid cosmo value")
                raise

            except IndexError:
                logger.exception("Could not add cosmo to control file")
                raise

        if self._MINN:
            logger.debug(f"Changing functional to: {self._MINN}")
            for line in controllines:
                if "functional" in line:
                    x = controllines.index(line)
                    controllines.pop(x)
                    controllines.insert(x, f"   functional {self._MINN}\n")

        if self._raw_parameters["weight"]:
            logger.debug("Adding weight derivatives")
            self.weight(controllines)

        if self._raw_parameters["denconv"]:
            logger.debug(f"Adding denconv to control: {self._raw_parameters['denconv']}")
            for line in controllines:
                if "$end" in line:
                    x = controllines.index(line)
                    controllines.insert(x, f"$denconv   {self._raw_parameters['denconv']}\n")
                    break

        try:
            logger.debug("Saving edited control file")
            with open("_control", 'w') as control:
                for line in controllines:
                    control.write(line)
            
            os.remove("control")
            os.rename("_control", "control")
        except IOError:
            logger.exception("Failed to save new control file in post-define editing")
            raise


    ###################################################################################
    #                      I O   H E L P E R   F U N C T I O N S                      #
    ###################################################################################

    def _write(self, output):
        """Responds to define based off of the state we are in. Saves responses to define.out'"""
        if self._shell.poll() is None:
            self._define_logger.info(">>>>" + output)
            self._shell.stdin.write(output + '\n')

    def _read_define_line(self):
        """Reads in the lines from define until it reachs a key line. Sets define_state to the necessary
        state in order to properly respond

        :return: True if a valid response line is reached. False if an error was encountered
        """
        while self._shell.poll() is None:
            line = self._read()
            self._define_logger.info(line)

            if line == "":
                continue
            elif "****  define : all done  ****" in line:
                self._errstate = "all done"
                return False

            elif "e n d   o f" in line:
                self._errstate = "all done"
                return False

            for possible_line in self._defineresponses.keys():
                if possible_line in line:

                    for wait in range(0, self._defineresponses[possible_line][0]):
                        self._define_logger.info(self._read())

                    self._define_state = self._defineresponses[possible_line][1]
                    return True

            for possible_error in self._errorresponses:
                if line == possible_error:
                    self._errstate = self._errorresponses[possible_error][1]
                    # abort define
                    self._write("qq")

                    while self._shell.poll() is None:
                        self._define_logger.info(self._read())

                    return False

        return False

    def _read(self):
        """ Reads in a line from define. Will return timeout_err if it takes too long to read

        :return: define line unless a timeout error has occurred, then returns timeout_err
        """
        # Check if timeout = 0
        if not self.timeout:
            return self._shell.stdout.readline().strip()

        # Set the alarm for timeout seconds
        signal.alarm(self.timeout)
        try:
            return self._shell.stdout.readline().strip()

        except exceptions.Alarm:
            return self._timeout_err

        finally:
            # turn off the alarm here
            signal.alarm(0)

    def __del__(self):
        """Makes sure the shell command is killed before destruction of object"""
        if not (self._shell is None) and self._shell.poll is None:
            self._shell.kill()

    def define_alarm_handler(self, signum, frame):
        """Raised if timeout in define"""
        raise exceptions.Alarm


class setupDMDjob:

    def __init__(self, parameters: dict=None, dir: str="./", pro: protein.Protein=None):
        logger.debug("Entered setupDMDjob")

        # Instance variables
        logger.debug("Initializing Variables")
        self._run_directory = dir
        logger.debug(f"Set run directory to: {self._run_directory}")
        self._initial_directory = os.getcwd()
        logger.debug(f"Set initial directory to: {self._initial_directory}")
        self._dmd_config = utilities.load_phd_config()
        os.environ["PATH"] += os.pathsep + self._dmd_config["PATHS"]["DMD_DIR"]
        self._protein = None

        try:
            logger.debug(f"Changing to run directory: {self._run_directory}" )
            os.chdir(self._run_directory)

        except OSError:
            logger.exception(f"Failed moving to run directory: {self._run_directory}")
            raise

        # Read in the parameters file, or use the passed parameters
        if parameters is None:
            logger.debug("Searching for a dmdinput.json file")
            if not os.path.isfile("dmdinput.json"):
                logger.error("Could not find dmdinput.json")
                logger.error("I will copy over a default dmdinput.json for you to edit")
                try:
                    shutil.copy(pkg_resources.resource_filename('phd3.resources', 'dmdinput.json'), './')

                except OSError:
                    logger.exception("Could not copy over default dmdinput.json!")
                    raise

                raise ValueError("dmdinput.json")

            try:
                with open('dmdinput.json') as inp:
                    self._raw_parameters = json.load(inp)

            except IOError:
                logger.exception("Error reading in dmdinput.json file!")
                raise

            except ValueError:
                logger.exception("dmdinput.json not formatted correctly")
                raise

        else:
            logger.debug("Using passed parameters")
            self._raw_parameters = parameters

        try:
            utilities.valid_dmd_parameters(self._raw_parameters)

        except ValueError:
            logger.exception("Missing a parameter definition!")
            raise ValueError("definition")

        except exceptions.ParameterError:
            logger.exception("Invalid parameter specification")
            raise

        if pro is None:
            logger.debug("Checking for a pdb")
            files = [f for f in os.listdir('./') if os.path.isfile(os.path.join('./', f))]
            logger.debug(f"The files in {self._run_directory} are :")
            logger.debug(files)
            for f in files:
                if os.path.splitext(f)[-1].lower() == ".pdb":
                    logger.debug(f"Found a pdb file to use: {f}")
                    try:
                        self._protein = utilities.load_pdb(f)

                    except IOError:
                        logger.debug("Ran into an issue with loading the pdb")
                        continue
                    break
            # We store the residue here and then we backtrack and figure out the correct atom to protonate later after
            # proper relabeling/renaming

        else:
            self._protein = pro

        if self._protein is None:
            raise ValueError("No Protein!")

        self._displacement = []
        if "Restrict Displacement" in self._raw_parameters.keys():
            for atom_pair in self._raw_parameters["Restrict Displacement"]:
                if type(atom_pair[0]) == list:
                    res1 = self._protein.get_atom(atom_pair[0])

                elif type(atom_pair[0]) == str:
                    pair_split = atom_pair[0].split(":")
                    res1 = self._protein.get_atom([pair_split[0], int(pair_split[1]), pair_split[2]])
            
                else:
                    logger.error("Invalid specification of displacement atom")
                    raise ValueError
                
                if type(atom_pair[1]) == list:
                    res2 = self._protein.get_atom(atom_pair[1])
                
                elif type(atom_pair[1]) == str:
                    pair_split = atom_pair[1].split(":")
                    res2 = self._protein.get_atom([pair_split[0], int(pair_split[1]), pair_split[2]])

                else:
                    logger.error("Invalid specification of displacement atom")
                    raise ValueError

                self._displacement.append([res1, res2, atom_pair[2]])

        self._static = []
        if "Frozen atoms" in self._raw_parameters.keys():
            for chain in self._raw_parameters["Frozen atoms"]["Chains"]:
                try:
                    protein_chain =self._protein.get_chain(chain)
                    for residue in protein_chain.residues:
                        self._static.extend(residue.atoms)

                except ValueError:
                    logger.exception("Could not find the chain!")
                    raise

            for residue in self._raw_parameters["Frozen atoms"]["Residues"]:
                try:
                    if type(residue) == list:
                        res = self._protein.get_residue(residue)

                    elif type(residue) == str:
                        residue_split = residue.split(":")
                        res = self._protein.get_residue([residue_split[0], int(residue_split[1])])
                    
                    else:
                        logger.error("Invalid specification of frozen residue")
                        raise ValueError

                    self._static.extend(res.atoms)

                except ValueError:
                    logger.exception("Could not find the residue!")
                    raise

            for atom in self._raw_parameters["Frozen atoms"]["Atoms"]:
                try:
                    if type(atom) == list:
                        self._static.append(self._protein.get_atom(atom))

                    elif type(atom) == str:
                        atom_split = atom.split(":")
                        self._static.append(self._protein.get_atom([atom_split[0], int(atom_split[1]), atom_split[2]]))

                    else:
                        logger.error("Invalid specification of frozen atom")
                        raise ValueError

                except ValueError:
                    logger.exception("Could not find the atom!")
                    raise

        # These holds all of the residues with weird protonation or deprotonation states
        self._protonate = []
        if "Custom protonation states" in self._raw_parameters.keys():
            for item in self._raw_parameters["Custom protonation states"]:
                if type(item[1]) == str:
                    tmp_item = item[0].split(":")
                    res_id = [tmp_item[0], int(tmp_item[1])]
                    self._protonate.append([self._protein.get_residue(res_id), item[1:]])
                
                elif type(item[1]) == int:
                    res_id = [item[0], item[1]]
                    self._protonate.append([self._protein.get_residue(res_id), item[2:]])
                
                else:
                    logger.error("Invalid specification of residue in protonation state")
                    raise ValueError


    def full_setup(self):

        logger.debug("Changing protein name to initial.pdb and writing out")
        self._protein.reformat_protein()
        self._protein.name = 'initial.pdb'
        self._protein.write_pdb()

        self.make_topparam()
        self.make_inConstr()
        utilities.make_state_file(self._raw_parameters, self._protein.name)
        self.short_dmd()
        utilities.make_start_file(self._raw_parameters)

        logger.info("[setup dmd]        ==>> SUCCESS")

    def titrate_setup(self):
        logger.debug("Skipping short dmd step")
        self._protein.reformat_protein()
        self._protein.name = 'initial.pdb'
        self._protein.write_pdb()

        self.make_inConstr()
        utilities.make_state_file(self._raw_parameters, self._protein.name)
        utilities.make_start_file(self._raw_parameters)
        logger.info("[titratable setup] ==>> SUCCESS")

    def short_dmd(self, keep_movie=False, time=1):
        try:
            with open("dmd_start_short", 'w') as dmdstart:
                dmdstart.write(f"THERMOSTAT     {self._raw_parameters['Thermostat']}\n")
                dmdstart.write(f"T_NEW          0.001\n")
                dmdstart.write(f"T_LIMIT        0.001\n")
                dmdstart.write(f"HEAT_X_C       1\n")
                dmdstart.write(f"RESTART_FILE   {self._raw_parameters['Restart File']}\n")
                dmdstart.write(f"RESTART_DT     1\n")
                dmdstart.write(f"ECHO_FILE      {self._raw_parameters['Echo File']}\n")
                dmdstart.write(f"ECHO_DT        1\n")
                dmdstart.write(f"MOVIE_FILE     {self._raw_parameters['Movie File']}\n")
                dmdstart.write(f"START_TIME     0\n")
                dmdstart.write(f"MOVIE_DT       1\n")
                dmdstart.write(f"MAX_TIME       {time}\n")

        except IOError:
            logger.exception("Error writing out dmd_start file")
            raise

        # Here we do a short run
        overlap = False
        try:
            with Popen(f"pdmd.linux -i dmd_start_short -s state -p param -c outConstr -m 1",
                    stdout=PIPE, stderr=subprocess.STDOUT, universal_newlines=True, shell=True, env=os.environ) as shell:
                while shell.poll() is None:
                    line = shell.stdout.readline().strip()
                    if "* Atoms overlap" in line:
                        overlap = True
                        logger.error(line)

                    else:
                        logger.debug(line)

        except OSError:
            logger.exception("Error calling pdmd.linux")
            raise

        if overlap:
            logger.error("Atoms are overlapping, checking for exact atoms!")
            for chain in self._protein.chains:
                for residues in chain.residues:
                    for i_atom in residues.atoms:
                        for n_res in chain.residues:
                            for n_atom in n_res.atoms:
                                if i_atom is n_atom:
                                    continue

                                if np.linalg.norm(i_atom.coords - n_atom.coords) < 0.5:
                                    logger.error(f"Check: {i_atom} and {n_atom} at residue {residues}!")


        if not os.path.isfile("movie"):
            logger.error("movie file was not made, dmd seems to be anrgy at your pdb")
            raise ValueError("initial.pdb")

        logger.debug("Finished the short DMD step successfully")

        utilities.make_movie("initial.pdb", "movie", "check.pdb")

        if not os.path.isfile("check.pdb"):
            logger.error("check.pdb not found, complex_M2P.linux did not run properly")
            raise ValueError("check.pdb")

        else:
            with open("check.pdb") as pdb:
                if len(pdb.readlines()) == 0:
                    logger.error("check.pdb is empty, complex_M2P.linux did not run properly")
                    raise ValueError("check.pdb")

        logger.debug("Was able to create a pdb from the short DMD run")
        logger.debug("Good to go, removing old file now")
        os.remove("check.pdb")
        if not keep_movie:
            os.remove(self._raw_parameters['Movie File'])
        
        os.remove(self._raw_parameters['Echo File'])
        os.remove(self._raw_parameters['Restart File'])
        
        os.remove("dmd_start_short")

    def make_topparam(self):
        try:
            logger.debug("Making topparam file")
            with open('topparam', 'w') as topparam_file:
                for residue in self._protein.sub_chain.residues:
                    try:
                        utilities.make_mol2(residue)

                    except OSError:
                        logger.error("Error in making mol2 file")
                        raise

                    topparam_file.write(f"MOL {residue.name} ./{residue.name}.mol2\n")

        except IOError:
            logger.exception("Error with writing topparam file!")
            raise

    def make_inConstr(self):
        try:
            with open('inConstr', 'w') as inConstr_file:
                if self._raw_parameters["Freeze Non-Residues"]:
                    logger.debug("Freeze Non-residues turned on, freezing residues")
                    for residue in self._protein.sub_chain.residues:
                        logger.debug(f"Freezing residue: {residue}")
                        inConstr_file.write(f"Static {residue.write_inConstr()}\n")

                for static_atom in self._static:
                    logger.debug(f"Freezing atom: {static_atom}")
                    inConstr_file.write(f"Static {static_atom.write_inConstr()}\n")

                for state in self._protonate:
                    logger.debug(f"Adding protonation state: {state[0]} and {state[1]}")
                    atom_id = ""
                    #TODO try and except for weird atoms or residues if it cannot find it
                    if len(state[1]) > 1:
                        #Then we had a number specify
                        logger.debug("Specified which atom specifically to use!")
                        
                        #For n-terminus
                        if state[1][1] == -1:
                            atom_id = "N"

                        #For c-terminus
                        elif state[1][1] == -2:
                            atom_id = "O"

                        elif state[1][0] == "protonate":
                            atom_id = constants.PROTONATED[state[0].name][state[1][1]-1]

                        elif state[1][0] == "deprotonate":
                            atom_id = constants.DEPROTONATED[state[0].name][state[1][1]-1]

                    else:
                        if state[1][0] == "protonate":
                            atom_id = constants.PROTONATED[state[0].name][0]

                        elif state[1][0] == "deprotonate":
                            atom_id = constants.DEPROTONATED[state[0].name][0]

                    if atom_id == "":
                        raise ValueError("Did not specify to protonate or deprotonate correctly")

                    try:
                        pro_atom = state[0].get_atom(atom_id[0])

                    except ValueError:
                        logger.exception("Could not find the correct atom to protonate or deprotonate in the residue")
                        logger.warning(f"{state}")
                        raise

                    if state[1][0] == "protonate":
                        inConstr_file.write(f"Protonate {pro_atom.write_inConstr()}\n")

                    else:
                        inConstr_file.write(f"Deprotonate {pro_atom.write_inConstr()}\n")

                if self._raw_parameters["Restrict Metal Ligands"]:
                    logger.debug("Restricting distance between atoms and metals!")
                    for metal in self._protein.metals:
                        logger.debug(f"Looking at metal: {metal}")
                        atoms_near_metal = self._protein.atoms_near_metal(metal, 3.1)
                        for atoms in atoms_near_metal:
                            logger.debug(f"Freezing atom: {atoms} since too close to a metal")
                            inConstr_file.write(f"Static {atoms.write_inConstr()}\n")

                            for bonded_atoms in atoms.bonds:
                                if bonded_atoms.element.lower() != "h" and bonded_atoms.element.lower() not in constants.METALS:
                                    logger.debug(f"Restricting motion of atom {bonded_atoms} and atom {atoms} by {0.05}")
                                    inConstr_file.write(
                                        f"AtomPairRel {bonded_atoms.write_inConstr()} {atoms.write_inConstr()} -{0.05} +{0.05}\n")

                for disp_atom in self._displacement:
                    logger.debug(
                        f"Restricting motion of atom: {disp_atom[0]} and atom {disp_atom[1]} by {disp_atom[2]}")
                    inConstr_file.write(
                        f"AtomPairRel {disp_atom[0].write_inConstr()} {disp_atom[1].write_inConstr()} -{disp_atom[2]} +{disp_atom[2]}\n")

        except IOError:
            logger.exception("Error opening inConstr file")
            raise

        logger.debug("Finished making the inConstr file!")

    def updated_parameters(self):
        #TODO make sure that this is correct
        new_parameters = self._raw_parameters.copy()

        # Update the custom protonation states
        for new_state, state in zip(new_parameters["Custom protonation states"], self._protonate):
            additional = False
            if len(new_state) == 4:
                additional = True

            new_state.clear()
            new_state.append(state[0].chain.name)
            new_state.append(state[0].number)
            new_state.append(state[1][0])
            if additional:
                new_state.append(state[1][1])

        # Update the frozen atoms
        new_parameters["Frozen atoms"]["Chains"].clear()
        new_parameters["Frozen atoms"]["Residues"].clear()
        new_parameters["Frozen atoms"]["Atoms"].clear()
        
        for a in self._static:
            new_parameters["Frozen atoms"]["Atoms"].append(a.label())

        # Update the displacement atoms
        new_parameters["Restrict Displacement"].clear()
        for state in self._displacement:
            new_parameters["Restrict Displacement"].append([state[0].label(), state[1].label(), state[2]])

        return new_parameters

    def update_from_movie(self):
        last_frame = utilities.load_movie("movie.pdb")[-1]
        linked_atoms = []
        
        for res in self._protein.sub_chain.residues:
            for atom in res.atoms:
                if atom.element.lower() not in constants.METALS:
                    linked_atoms.append([atom, last_frame.get_atom(atom.label())])
        
        #Now we reformat the last_frame so that it is consistent...
        last_frame.reformat_protein(relabel_protein=False)
        for atom_pairs in linked_atoms:
            atom_pairs[0].id = atom_pairs[1].id
            atom_pairs[0].number = atom_pairs[1].number

        for res in self._protein.sub_chain.residues:
            res.reorder_atoms()



class setupPHDjob:

    #call setupDMD first
    #get updated parameters, need a track residues/chains/atoms etc...
    #convert pdb to coord file (ie the chop)
    #go from coord to pdb again (should be no issue here)
    #save new phdinput file (expanded from what it was) with updated params
    def __init__(self):

        if not os.path.isfile("phdinput.json"):
            logger.error("No phdinput.json file")
            logger.error("Copying over a default file to fill in")
            try:
                shutil.copy(pkg_resources.resource_filename('phd3.resources', 'phdinput.json'), './')

            except OSError:
                logger.exception(f"Could not copy default phdinput.json from: {os.path.dirname(os.path.abspath(__file__))} to: {os.getcwd()}")
                raise

            raise FileNotFoundError("phdinput.json")

        with open("phdinput.json") as param_file:
            self._parameters = json.load(param_file)

        #Load in the protein
        protein = utilities.load_pdb(self._parameters["pdb file"])

        #Get the qm_chop info to track those residues and atoms as well
        track_residues = []
        track_multi = []
        for res in self._parameters["QM Chop"]["Residues"]:
            if "-" in res:

                res = res.split("-")
                res1 = res[0]
                res2 = res[1]

                res1_back = 'n'
                if res1[-1].isalpha():
                    res1_back = res1[-1]
                    res1 = res1[:-1]

                res1 = res1.split(":")
                res1 = protein.get_residue([res1[0], int(res1[1])])

                res2_back = 'n'
                if res2[-1].isalpha():
                    res2_back = res2[-1]
                    res2 = res2[:-1]
                
                res2 = res2.split(":")
                res2 = protein.get_residue([res2[0], int(res2[1])])

                track_multi.append([res1, res1_back, res2, res2_back])
                    
            else:
                res = res.split(":")
                track_residues.append(protein.get_residue([res[0], int(res[1])]))

        track_exclude_atoms = []
        if "Exclude Atoms" in self._parameters["QM Chop"].keys():
            for atom in self._parameters["QM Chop"]["Exclude Atoms"]:
                atom = atom.split(":")
                if len(atom) == 3:
                    track_exclude_atoms.append(protein.get_atom([atom[0], int(atom[1]), atom[2]]))
                elif len(atom) == 2:
                    track_exclude_atoms.extend(protein.get_residue([atom[0], int(atom[1])]).atoms)
                else:
                    logger.error(f"Invalid specification of Exclude atoms {atom}")

        track_substrate_chop = []
        if "Substrate Chop" in self._parameters["QM Chop"].keys():
            for chop in self._parameters["QM Chop"]["Substrate Chop"]:
                atom1 = chop.split("-")[0]
                atom1 = atom1.split(":")
                atom1 = protein.get_atom([atom1[0], int(atom1[1]), atom1[2]])

                atom2 = chop.split("-")[1]
                atom2 = atom2.split(":")
                atom2 = protein.get_atom([atom2[0], int(atom2[1]), atom2[2]])

                track_substrate_chop.append([atom1, atom2])

        track_exclude_side = []
        if "Exclude Side Chain" in self._parameters["QM Chop"].keys():
            for res in self._parameters["QM Chop"]["Exclude Side Chain"]:
                res = res.split(":")
                track_exclude_side.append(protein.get_residue([res[0], int(res[1])]))

        track_protonation = []
        if "Protonation" in self._parameters["QM Chop"].keys():
            for protonation_state in self._parameters["QM Chop"]["Protonation"]:
                res = protonation_state[0].split(":")
                res = protein.get_residue([res[0], int(res[1])])
                track_protonation.append([res] + protonation_state[1:])

        track_freeze = []
        if "Freeze Atoms" in self._parameters["QM Chop"].keys():
            for atom in self._parameters["QM Chop"]["Freeze Atoms"]:
                atom = atom.split(":")
                chain = atom[0]
                res_num = int(atom[1])
                atom_id = atom[2]
                track_freeze.append(protein.get_atom([chain, res_num, atom_id]))

        track_dummy = []
        if "Dummy H" in self._parameters["QM Chop"].keys():
            for atom in self._parameters["QM Chop"]["Dummy H"]:
                atom = atom.split(":")
                chain = atom[0]
                res_num = int(atom[1])
                atom_id = atom[2]
                track_dummy.append(protein.get_atom([chain, res_num, atom_id]))

        if os.path.isdir("dmd_setup"):
            shutil.rmtree("dmd_setup")

        os.mkdir("dmd_setup")
        os.chdir("dmd_setup")

        sdj = setupDMDjob(pro=protein, parameters = self._parameters["dmd params"])
        sdj.full_setup()

        self._parameters["dmd params"] = sdj.updated_parameters()
        logger.info(">>>> Running Short DMD >>>>")
        logger.info("...")
        sdj.short_dmd(keep_movie=True, time=50)
        utilities.make_movie("initial.pdb", self._parameters["dmd params"]["Movie File"], "movie.pdb")

        # Updates numbering of hydrogens on substrates
        sdj.update_from_movie()

        self._parameters["QM Chop"]["Residues"] = [res.label() for res in track_residues]
        self._parameters["QM Chop"]["Residues"] += [f"{res[0].label()}{res[1]}-{res[2].label()}{res[3]}" for res in track_multi]
        self._parameters["QM Chop"]["Exclude Atoms"] = [atom.label() for atom in track_exclude_atoms]
        self._parameters["QM Chop"]["Substrate Chop"] = [f"{atoms[0].label()}-{atoms[1].label()}" for atoms in track_substrate_chop]
        self._parameters["QM Chop"]["Exclude Side Chain"] = [res.label() for res in track_exclude_side]
        self._parameters["QM Chop"]["Protonation"] = [[res[0].label()] + res[1:] for res in track_protonation]
        self._parameters["QM Chop"]["Freeze Atoms"] = [atom.label() for atom in track_freeze]
        self._parameters["QM Chop"]["Dummy H"] = [atom.label() for atom in track_dummy]

        logger.info(">>>> Loading in Movie >>>>")
        logger.info("...")
        last_protein = utilities.load_movie("movie.pdb")[-1]

        os.chdir("..")
        if os.path.isdir("qm_setup"):
            shutil.rmtree("qm_setup")

        os.mkdir("qm_setup")
        os.chdir("qm_setup")

        logger.info(">>>> Converting to Coord >>>>")
        dmd_to_qm.protein_to_coord(last_protein, self._parameters["QM Chop"])
        stj = setupTMjob(parameters=self._parameters["qm params"])

        os.chdir("..")

        self._parameters["pdb file"] = "start.pdb"
        protein.write_pdb(name="start.pdb")

        #Write out the phdinput.json
        with open("phdinput.json", "w") as param_file:
            json.dump(self._parameters, param_file, indent=4)

        logger.info("Finished setting up job")
        logger.info("Check dmd_setup and qm_setup to ensure proper chop and parameter interpretation")
        logger.info("Also note that phd will use the start.pdb instead of the initial pdb provided")

