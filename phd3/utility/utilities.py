#!/usr/bin/env python3
"""
Author  ==>> Matthew R. Hennefarth
Date    ==>> April 16, 2020
"""

#Standard Library Imports
import os
import json
import pkg_resources
import logging
import shutil
import random
import copy
from logging.config import dictConfig
from subprocess import Popen, PIPE
import subprocess
from random import shuffle
import numpy as np
from datetime import datetime

#PHD3 Imports
from ..protein import atom, chain, residue, protein

#import phd3.protein.atom as atom
#import phd3.protein.chain as chain
#import phd3.protein.residue as residue
#import phd3.protein.protein as protein 
from . import constants
from .exceptions import ParameterError
#from phd3.utility import constants
#from phd3.utility.exceptions import ParameterError


__all__=[
    'load_pdb',
    'make_mol2',
    'make_start_file',
    'make_state_file',
    'make_movie',
    'setup_dmd_environ',
    'valid_qm_parameters',
    'load_movie',
    'setup_turbomole_env',
    'valid_dmd_parameters',
    'create_config',
    'load_phd_config',
    'copy_directories',
    'quote_me',
    'xyz_to_coord',
    'add_proton',
    'addH',
]

logger = logging.getLogger(__name__)

def setup_turbomole_env(turbodir: str):
    """
    Sets up the TURBOMOLE environment. Ensures that the paths are extended correctly in order to run any TURBOMOLE
    commands including jobex, define, NumForce. It also sets up the parallel environment for TURBOMOLE

    :param turbodir: Directory where TURBOMOLE is installed
    """
    logger.debug("Setting up the TURBOMOLE environment")

    if turbodir in os.environ["PATH"] and f"{turbodir}/scripts" in os.environ["PATH"]:
        #Should already be setup
        return
    
    os.environ["TURBODIR"] = turbodir
    os.environ["PATH"] += os.pathsep + turbodir + '/scripts'
    os.environ["PATH"] += os.pathsep + turbodir

    sysname = ""
    with Popen('sysname', universal_newlines=True, shell=True,
               stdin=PIPE, stdout=PIPE, stderr=PIPE, bufsize=1) as shell:
        while shell.poll() is None:
            sysname += shell.stdout.readline().strip()
    
    os.environ["PATH"] += os.pathsep + turbodir + f'/bin/{sysname}'
    logger.debug("Finished setting up the TURBOMOLE environment")

def valid_qm_parameters(parameters: dict):
    """
    Ensures that the dictionary provided has the proper parameters for _create_define_options to run within the class
    SetupTurbomole. Can be used by other scripts to ensure that they have the proper parameters (but not necessarily
    that define will run properly) without do too much extra work

    :param parameters: dictionary that contains all of the parameters needed to setup a turbomole job. Should be modeled
    after the definput.json file
    :return: True if everything goes well. Will raise exceptions of type ParameterError if there is an issue/error
    encountered
    """
    logger.debug("Checking if parameters are valid")
    logger.debug("Checking geometry parameters")
    if parameters["geometry"]["idef"]["idef_on"]:
        logger.debug("Checking idef")
        for bond in parameters["geometry"]["idef"]["freeze_stretch"]:
            try:
                atom1, atom2 = bond.split(',')
                assert (atom1.isdigit() and int(atom1) >= 0)
                assert (atom2.isdigit() and int(atom2) >= 0)
            except ValueError:
                raise ParameterError(f"incorrect formatting of idef bonds : {bond}")

        if "freeze_dihedral" in parameters["geometry"]["idef"].keys():
            for dihedral in parameters["geometry"]["idef"]["freeze_dihedral"]:
                try:
                    atoms = dihedral.split(',')
                    for a in atoms:
                        assert(a.isdigit() and int(a) >= 0)

                except ValueError:
                    raise ParameterError(f"incorrect formatting of idef dihedral : {dihedral}")

    if parameters["geometry"]["iaut"]['iaut_on']:
        logger.debug("Checking iaut")
        for bond in parameters["geometry"]["iaut"]["bonds"]:
            try:
                atom1, atom2 = bond.split(',')
                assert (atom1.isdigit())
                assert (atom2.isdigit())
            except ValueError:
                raise ParameterError(f"incorrect formatting of iaut bonds : {bond}")

    if parameters["basis"]:
        logger.debug("Checking basis")
        for atom in parameters["basis"].keys():
            if not parameters["basis"][atom] in constants.AVAILABLE_BASIS:
                raise ParameterError(f"invalid basis set for: {atom} | {parameters['basis'][atom]}")

    try:
        logger.debug("Checking charge")
        assert (type(parameters["charge"]) == int)

    except:
        logger.error("Invalid charge type provided!")
        raise ParameterError("invalid charge provided")

    if parameters["open_shell"]["open_shell_on"]:
        logger.debug("Checking unpaired electrons")
        try:
            assert (type(parameters["open_shell"]["unpaired"]) == int)

        except:
            logger.error("Invalid value for unpaired electrons")
            raise ParameterError("invalid number of unpaired electrons")

    if parameters["dft"]["dft_on"]:
        logger.debug("Checking dft")
        if parameters["dft"]["func"] not in constants.AVAILABLE_FUNCS and parameters["dft"]["func"] not in constants.MINN_FUNCS:
            logger.error("Functional provided not available")
            raise ParameterError(f"invalid functional provided: {parameters['dft']['func']}")
        if parameters["dft"]["grid"] not in constants.AVAILABLE_GRIDS:
            logger.error("Grid provided is not available")
            raise ParameterError(f"invalid grid provided: {parameters['dft']['grid']}")

    if parameters["scf"]:
        if "damp start" in parameters["scf"]:
            if type(parameters["scf"]["damp start"]) != int and type(parameters["scf"]["damp start"]) != float:
                logger.error("Invalid damp start value, must be int or float")
                raise ParameterError(f"Invalid scf damp value provided: {parameters['scf']['damp start']}")
            
            if parameters["scf"]["damp start"] < 0 :
                logger.error("Cannot have a negative damp start value, must be >=0")
                raise ParametersError(f"Invalid scf damp value provided: {parameters['scf']['damp start']}")

        if "orbital shift" in parameters["scf"]:
            if type(parameters["scf"]["orbital shift"]) != int and type(parameters["scf"]["orbital shift"]) != float:
                logger.error("Invalid orbital shift value, must be int or float")
                raise ParameterError(f"Invalid orbital shift value: {parameters['scf']['orbital shift']}")

            if parameters["scf"]["orbital shift"]< 0:
                logger.error(f"Invlid orbital shift value, must be >= 0")
                raise ParameterError(f"Invalid orbital shift value: {parameters['scf']['orbital shift']}")


    if parameters["stp"]:
        try:
            logger.debug("Checking itvc")
            assert (type(parameters["stp"]["itvc"]) == int)

        except:
            logger.error("Invalid itvc value")
            raise ParameterError("invalid input for itvc")

        try:
            logger.debug("Checking trust radius")
            assert (type(parameters["stp"]["trad"]) == float)

        except:
            logger.error("Invalid trust radius value")
            raise ParameterError("invalid input for trad")

    if parameters["cosmo"]:
        logger.debug("Checking cosmo")
        try:
            assert (type(parameters["cosmo"]) == float or type(parameters["cosmo"]) == int)

        except:
            logger.error("Invalid cosmo value")
            raise ParameterError("invalid input for cosmo")

    if parameters["freeze_atoms"]:
        try:
            logger.debug("Checking frozen atoms")
            assert (type(atom) == int for atom in parameters["freeze_atoms"])

        except:
            logger.error("Invalid list of frozen cartesian atoms")
            raise ParameterError("invalid list of frozen cartesian atoms")

    if parameters["gcart"]:
        logger.debug("Checking gcart")
        try:
            assert (type(parameters["gcart"]) == int)

        except:
            logger.error("Invalid gcart value")
            raise ParameterError("invalid value for gcart")

    try:
        assert (type(parameters["weight"]) == bool)

    except:
        logger.error("Invalid value for weight")
        raise ParameterError("invalid value for weight")

    if parameters["calculation"] == "":
        logger.warning("No calculation type set!")

    else:
        if parameters["calculation"] not in constants.AVAILABLE_CALCULATIONS:
            logger.error("Invalid value for calculation")
            raise ParameterError("invalid calculation type")

        if parameters["calculation"] == "trans" and parameters["stp"]["itvc"] == 0:
            logger.error("Cannot perform trans calculation with itvc = 0!")
            raise ParameterError("calculation and itvc do not match")

        if parameters["calculation"] == "geo" or parameters["calculation"] == "forceopt":
            if parameters["stp"]["itvc"] != 0:
                logger.error("Cannot perform geo calculation with itvc not 0!")
                raise ParameterError("calculation and itvc do not match")

def create_config():
    # TODO: have a function that verifies the accuracy of the config file!
    logger.debug(f"creating .phd3 in the home directory: {os.path.expanduser('~')}")
    home = os.path.expanduser("~")
    try:
        os.mkdir(os.path.join(home, '.phd3'))

    except FileExistsError:
        logger.debug(".phd3 directory already exists, continuing")

    logger.debug("Placing default phd_config.json in the .turbopy directory")
    shutil.copy(pkg_resources.resource_filename('phd3.resources', 'phd_config.json'), os.path.join(home, '.phd3'))
    shutil.copy(pkg_resources.resource_filename('phd3.resources', 'logger_config.json'), os.path.join(home, '.phd3'))
    logger.info(f"Please ensure that {os.path.join(home, '.phd3/phd_config.json')} has the correct values.")


def load_phd_config():
    """Loads in the PHD Config File

    :return: Dictionary of the config file
    """
    home = os.path.expanduser("~")
    
    # Check in the .config directory first, and then in the home directory
    path_to_config = os.path.join(home, ".config/phd3/phd_config.json")
    
    if not os.path.isdir(os.path.join(home, ".config/phd3")) and not os.path.isfile(path_to_config):
        path_to_config = os.path.join(home, '.phd3/phd_config.json')

        if not os.path.isdir(os.path.join(home, '.phd3')):
            create_config()
            raise ValueError(".phd3 missing")

        if not os.path.isfile(path_to_config):
            logger.error("Config file not in .phd3 directory")
            logger.warning("Copying default logger over now")
            shutil.copy(pkg_resources.resource_filename('phd3.resources', 'phd_config.json'), os.path.join(home, '.phd3'))
            raise ValueError("phd_config.json file missing")

    try:
        logger.debug("Loading in phd_config")
        with open(path_to_config, 'r') as turbopy_config_file:
            config = json.load(turbopy_config_file)

    except IOError:
        logger.critical("Couldn't find config file!")
        raise

    return config

phd_config = load_phd_config()

def load_logger_config():
    """Loads in the Logger Config File"""
    home = os.path.expanduser("~")
    path_to_config = os.path.join(home, ".config/phd3/logger_config.json")
    
    if not os.path.isdir(os.path.join(home, ".config/phd3")) and not os.path.isfile(path_to_config):
        path_to_config = os.path.join(home, '.phd3/logger_config.json')
        if not os.path.isdir(os.path.join(home, '.phd3')):
            create_config()
            raise ValueError(".phd3 missing")

        if not os.path.isfile(path_to_config):
            print("Logger file not in .phd3 directory")
            print("Copying default logger over now")
            shutil.copy(pkg_resources.resource_filename('phd3.resources', 'logger_config.json'), os.path.join(home, '.phd3'))

    try:
        logger.debug("Loading in logger_config")
        with open(path_to_config, 'r') as config_file:
            logging_config = json.load(config_file)

    except IOError:
        print("CRITICAL: Error in loading in config file!!")
        raise
    
    try:
        dictConfig(logging_config)
    
    except ValueError:
        print("CRITICAL: Error in the logging_config.json file")
        raise IOError

    logger.debug("Loaded in logger parameters")
    logger.debug("Logger has started!")

def quote_me():
    return constants.QUOTES[random.randint(1,10000) % len(constants.QUOTES)]

def setup_dmd_environ():
    """ setups of an os.environ so that DMD can be run """
    logger.debug("Setting up DMD Environment")
    os.environ["PATH"] =  phd_config["PATHS"]["DMD_DIR"] + os.pathsep + os.environ["PATH"]
    return os.environ

def valid_dmd_parameters(parameters: dict):
    """Checks to see if the dmd parameters passed are valid at all"""
    logger.debug("Checking if parameters are valid")
    if "Thermostat" not in parameters.keys():
        raise ParameterError("No Thermostat")

    elif parameters["Thermostat"] not in constants.THERMOSTATS:
        raise ValueError("Invalid thermostat")

    if "Initial Temperature" not in parameters.keys():
        raise ParameterError("No Initial Temperature")

    else:
        try:
            assert(type(parameters["Initial Temperature"]) is float and parameters["Initial Temperature"] > 0)

        except ValueError:
            raise ParameterError("Incorrect initial temperature provided")

    if "Final Temperature" not in parameters.keys():
        raise ParameterError("No Final Temperature")

    else:
        try:
            assert (type(parameters["Final Temperature"]) is float and parameters["Final Temperature"] > 0)

        except ValueError:
            raise ParameterError("Incorrect Final temperature provided")

    if "HEAT_X_C" not in parameters.keys():
        raise ParameterError("No HEAT_X_C")

    else:
        try:
            assert (type(parameters["HEAT_X_C"]) is float and parameters["HEAT_X_C"] > 0)

        except ValueError:
            raise ParameterError("Incorrect HEAT_X_C provided")

    if "Echo File" not in parameters.keys():
        raise ParameterError("You did not specify what to label the echo file")

    if "Movie File" not in parameters.keys():
        raise ParameterError("You did not specify what to call the movie file")

    if 'Restart File' not in parameters.keys():
        raise ParameterError("You did not specify what to call the restart file")

    if "dt" not in parameters.keys():
        raise ParameterError("You did not specify the time in which save data")

    else:
        try:
            assert(type(parameters["dt"]) is int and parameters["dt"] > 0)

        except ValueError:
            raise ParameterError(f"Invalid value provided for dt: {parameters['dt']}")

    if "Time" not in parameters.keys():
        raise ParameterError("You did not specify the default length of time to be for the DMD simulation")

    else:
        try:
            assert (type(parameters["Time"]) is int and parameters["Time"] > 0)

        except ValueError:
            raise ParameterError(f"Invalid value provided for Time: {parameters['Time']}")

    if "titr" in parameters.keys():
        try:
            assert(type(parameters["titr"]["titr on"]) == bool)

        except ValueError:
            raise ParameterError("You did not correctly type the titr on parameters, must be a bool")

    if "Freeze Non-Residues" not in parameters.keys():
        raise ValueError("Missing Freeze Non-Residues parameters")

    else:
        try:
            assert(type(parameters["Freeze Non-Residues"]) == bool)

        except ValueError:
            raise ParameterError("Freeze Non-Residues MUST be a bool")

    if "Restrict Metal Ligands" not in parameters.keys():
        pass

    else:
        try:
            assert (type(parameters["Restrict Metal Ligands"]) == bool)

        except ValueError:
            raise ParameterError("Restrict Metal Ligands MUST be a bool")

    if "Custom protonation states" not in parameters.keys():
        raise ValueError("Missing Custom protonation state parameters")

    else:
        for state in parameters["Custom protonation states"]:
            try:
                assert(type(state[0]) == str)
                if ":" in state[0]:
                    assert(len(state) == 2 or len(state) == 3)
                    assert(type(state[1]) == str)
                    if len(state) == 3:
                        assert(type(state[2]) == int)

                else:
                    assert(len(state) == 3 or len(state) == 4)
                    assert(type(state[1]) == int)
                    assert(type(state[2]) == str)
                    if len(state) == 4:
                        assert(type(state[3]) == int)

            except ValueError:
                raise ParameterError(f"Invalid specification of protonation states: {state}")

    if "Frozen atoms" in parameters.keys():
        if "Chains" not in parameters["Frozen atoms"].keys():
            raise ValueError("Missing Chains in frozen atom definition")

        else:
            for chain in parameters["Frozen atoms"]["Chains"]:
                try:
                    assert(type(chain) == str)
                    assert(chain.isalpha())

                except ValueError:
                    raise ParameterError(f"Incorrect value provided for frozen chain: {chain}")

        if "Residues" not in parameters["Frozen atoms"].keys():
            raise ValueError("Missing Residue in frozen atom definition")

        else:
            for residue in parameters["Frozen atoms"]["Residues"]:
                try:
                    if type(residue) == list:
                        assert(type(residue[0]) == str and residue[0].isalpha())
                        assert(type(residue[1]) == int and residue[1] > 0)

                    elif type(residue) == str:
                        tmp = residue.split(":")
                        assert(len(tmp) == 2)
                        assert(type(tmp[0]) == str and tmp[0].isalpha())
                        assert(int(tmp[1]) > 0)

                    else:
                        raise ValueError

                except ValueError:
                    raise ParameterError(f"Incorrect value provided for frozen residue: {residue}")

        if "Atoms" not in parameters["Frozen atoms"].keys():
            raise ValueError("Missing Atom in frozen atom definition")

        else:
            for atom in parameters["Frozen atoms"]["Atoms"]:
                try:
                    if type(atom) == list:
                        assert(type(atom[0]) is str and atom[0].isalpha())
                        assert(type(atom[1]) is int and atom[1] > 0)
                        assert(type(atom[2]) is str)

                    elif type(atom) == str:
                        tmp = atom.split(":")
                        assert(len(tmp) == 3)
                        assert(tmp[0].isalpha())
                        assert(int(tmp[1]) > 0)

                except ValueError:
                    raise ParameterError(f"Incorrect value provided for frozen atom: {atom}")

    if "Restrict Displacement" in parameters.keys():
        for id in parameters["Restrict Displacement"]:
            try:
                assert(len( id) == 3)
                assert(type(id[2]) is float and id[2] > 0)

                if type(id[0]) == list:
                    assert(len(id[0]) == 3)
                    assert(type(id[0][0]) is str and id[0][0].isalpha())
                    assert(type(id[0][1]) is int and id[0][1] > 0)
                    assert(type(id[0][2]) is str)

                elif type(id[0]) == str:
                    tmp = id[0].split(":")
                    assert(len(tmp) == 3)
                    assert(tmp[0].isalpha())
                    assert(int(tmp[1]) > 0)
                
                else:
                    raise ValueError
                
                if type(id[1]) is list:
                    assert(len(id[1]) == 3)
                    assert (type(id[1][0]) is str and id[1][0].isalpha())
                    assert (type(id[1][1]) is int and id[1][1] > 0)
                    assert (type(id[1][2]) is str)

                elif type(id[1]) is str:
                    tmp = id[1].split(":")
                    assert(len(tmp) == 3)
                    assert(tmp[0].isalpha())
                    assert(int(tmp[1]) > 0)

                assert(id[0] != id[1])

            except ValueError:
                raise ParameterError(f"Invalid definition of Restrict Displacement: {id}")

    if "Commands" not in parameters.keys():
        raise ValueError("Missing commands definition")

    else:
        for command in parameters["Commands"]:
            try:
                assert(type(parameters["Commands"][command]) is dict)

            except ValueError:
                raise ParameterError(f"Command provided is not a dictionary: {command}")


def load_pdb(file: str):
    logger.debug(f"Finding file: {file}")

    chains = []

    resNum = 0
    chainLet = ""
    lineNumber = 1
    try:
        with open(file, 'r') as pdb:
            for line in pdb:
                if "ATOM" == line[0:4] or "HETATM" == line[0:6]:
                    tmpAtom = atom.Atom(line)
                    if chainLet != line[21:22]:
                        chainLet = line[21:22]
                        chains.append(chain.Chain(chainLet))
                        resNum = 0

                    if resNum != int(line[22:26]):
                        resNum = int(line[22:26])
                        chains[-1].add_residue(residue.Residue(line))

                    chains[-1].residues[-1].add_atom(tmpAtom)
                lineNumber += 1
    except IOError:
        logger.exception(f"Error opening {file}")
        raise
    except ValueError:
        logger.exception("Something bad happened!")
        raise

    logger.debug("Successfully loaded in the file!")
    return protein.Protein(file, chains)


def make_mol2(res: residue, reformat: bool=True):
    # TODO: try and except wrap
    # TODO: add logger stuff to this
    successful = False
    with open(f"{res.name}.pdb", 'w') as mol2_file:
        for atom in res.atoms:
            mol2_file.write(atom.pdb_line())

        mol2_file.write('TER\nENDMDL')
    # Now we execute the babel command here -- try openbabel 3 first
    with Popen(f"obabel -i pdb {res.name}.pdb -o mol2 -O {res.name}.mol2", stdin=PIPE, stdout=PIPE, stderr=PIPE,
               universal_newlines=True, shell=True, bufsize=1, env=os.environ) as shell:
        while shell.poll() is None:
            logger.debug(shell.stdout.readline().strip())
            output = shell.stderr.readline().strip()
            logger.debug(output)
            if "1 molecule converted" in output:
                successful = True

    if not successful:  # try with openbabel 2 syntax
        with Popen(f"babel {res.name}.pdb {res.name}.mol2", stdin=PIPE, stdout=PIPE, stderr=PIPE,
                   universal_newlines=True, shell=True, bufsize=1, env=os.environ) as shell:
            while shell.poll() is None:
                logger.debug(shell.stdout.readline().strip())
                output = shell.stderr.readline().strip()
                logger.debug(output)
                if "1 molecule converted" in output:
                    successful = True

    if not successful:
        logger.error("Could not create {residue.name} mol2 file!")
        raise OSError("mol2_file")

    if reformat:
        # Now we fix the mol2 file
        no_atoms = []
        hatoms = []
        polarhatoms = []
        atom_section = False
        bond_section = False

        mol_file_lines = []

        with open(f"{res.name}.mol2") as mol_file:
            for line in mol_file:
                if "ATOM" in line:
                    atom_section = True
                    bond_section = False

                elif "BOND" in line:
                    atom_section = False
                    bond_section = True

                elif atom_section:
                    if line[47:49] == 'N.' or line[47:49] == 'O.':
                        no_atoms.append(int(line[:7]))

                    elif line[47:49] == 'H ':
                        hatoms.append(int(line[:7]))

                elif bond_section:
                    if int(line[6:12]) in no_atoms and int(line[12:18]) in hatoms:
                        polarhatoms.append(int(line[12:18]))

                    elif int(line[12:18]) in no_atoms and int(line[6:12]) in hatoms:
                        polarhatoms.append(int(line[6:12]))

                mol_file_lines.append(line)

        atom_section = False
        atomline = 0
        with open(f"{res.name}.mol2", 'w+') as mol_file:
            for line in mol_file_lines:
                if "ATOM" in line:
                    atom_section = True

                elif "BOND" in line:
                    atom_section = False

                if atom_section:
                    if atomline in hatoms and atomline not in polarhatoms:
                        line = line[:47] + 'Eh' + line[49:]

                    atomline += 1

                mol_file.write(line)

    os.remove(f"{res.name}.pdb")
    logger.info(f"Successfuly made: {res.name} mol2")


def make_start_file(parameters: dict, start_time: int =0):
    logger.debug("Making the Start File")
    try:
        with open("dmd_start", 'w') as dmdstart:
            dmdstart.write(f"THERMOSTAT     {parameters['Thermostat']}\n")
            dmdstart.write(f"T_NEW          {parameters['Initial Temperature']}\n")
            dmdstart.write(f"T_LIMIT        {parameters['Final Temperature']}\n")
            dmdstart.write(f"HEAT_X_C       {parameters['HEAT_X_C']}\n")
            dmdstart.write(f"RESTART_FILE   {parameters['Restart File']}\n")
            dmdstart.write(f"RESTART_DT     {parameters['dt']}\n")
            dmdstart.write(f"ECHO_FILE      {parameters['Echo File']}\n")
            dmdstart.write(f"ECHO_DT        {parameters['dt']}\n")
            dmdstart.write(f"MOVIE_FILE     {parameters['Movie File']}\n")
            dmdstart.write(f"START_TIME     {start_time}\n")
            dmdstart.write(f"MOVIE_DT       {parameters['dt']}\n")
            dmdstart.write(f"MAX_TIME       {parameters['Time'] + start_time}\n")

    except IOError:
        logger.exception("Error writing out dmd_start file")
        raise

    logger.debug("made the start file!")

def copy_restart_velocities(state_file, restart_file):
    restart_velocities = []

    with open(restart_file, "r+") as res_f:
        read_velocity = False
        for line in res_f.readlines():
            split_line = line.split()

            if read_velocity == True:
                if split_line[0] == 'REACTIONED':
                    read_velocity = False
                    break

                restart_velocities += [[int(split_line[9]), split_line[5][:-2], split_line[6][:-2], split_line[7][:-2]]]

            else:
                if split_line[0] == '#Format:':
                    if split_line[1] == "AtomIndex":
                        read_velocity = True

    new_state_content = ''
    with open(state_file, "r+") as sta_f:
        change_velocity = False
        atom_line_offset = 0
        for idx, line in enumerate(sta_f):
            split_line = line.split()
            if change_velocity == True:
                if int(split_line[9]) == restart_velocities[idx - atom_line_offset][0]:
                    new_state_content += '  ' + split_line[0] + '  ' + split_line[1] + '  ' + split_line[2] + '  ' + split_line[3] + '  ' + split_line[4] + '  ' + restart_velocities[idx - atom_line_offset][1] + '  ' + restart_velocities[idx - atom_line_offset][2] + '  ' + restart_velocities[idx - atom_line_offset][3] + '  ' + split_line[8] + ' ' + split_line[9]
                    new_state_content += '\n'
                elif int(split_line[9]) == (restart_velocities[idx - atom_line_offset][0] + 1): # Any given residue should only vary by one atom, and residue numbers should only be offset by one number
                    atom_line_offset -= 1
                    new_state_content += '  ' + split_line[0] + '  ' + split_line[1] + '  ' + split_line[2] + '  ' + split_line[3] + '  ' + split_line[4] + '  ' + restart_velocities[idx - atom_line_offset][1] + ' ' + restart_velocities[idx - atom_line_offset][2] + '  ' + restart_velocities[idx - atom_line_offset][3] + '  ' + split_line[8] + ' ' + split_line[9]
                    new_state_content += '\n'
                elif (int(split_line[9]) + 1) == restart_velocities[idx - atom_line_offset][0]:
                    new_state_content += '  ' + split_line[0] + '  ' + split_line[1] + '  ' + split_line[2] + '  ' + split_line[3] + '  ' + split_line[4] + '  ' + '0.0000000000000' + '  ' + '0.0000000000000' + '  ' + '0.0000000000000' + '  ' + split_line[8] + ' ' + split_line[9]
                    new_state_content += '\n'
                    atom_line_offset += 1
                else: # Except if the last residue has a protonation state change immediately before the substrate
                    if int(split_line[9]) < restart_velocities[idx - atom_line_offset][0]: # Any given residue should only vary by one atom, and residue numbers should only be offset by one number
                        atom_line_offset -= 1
                        new_state_content += '  ' + split_line[0] + '  ' + split_line[1] + '  ' + split_line[2] + '  ' + split_line[3] + '  ' + split_line[4] + '  ' + restart_velocities[idx - atom_line_offset][1] + '  ' + restart_velocities[idx - atom_line_offset][2] + '  ' + restart_velocities[idx - atom_line_offset][3] + '  ' + split_line[8] + ' ' + split_line[9]
                        new_state_content += '\n'
                    elif int(split_line[9]) > restart_velocities[idx - atom_line_offset][0]:
                        new_state_content += '  ' + split_line[0] + '  ' + split_line[1] + '  ' + split_line[2] + '  ' + split_line[3] + '  ' + split_line[4] + '  ' + '0.0000000000000' + '  ' + '0.0000000000000' + '  ' + '0.0000000000000' + '  ' + split_line[8] + ' ' + split_line[9]
                        new_state_content += '\n'
                        atom_line_offset += 1
            else:
                new_state_content += line
                if split_line[0] == 'ATOMS':
                    change_velocity = True
                    atom_line_offset = idx + 1

    with open(state_file, "w+") as sta_of:
        sta_of.write(new_state_content)


def make_state_file(parameters: dict, pdbName):
    if os.path.exists("restart_velocity_cycle"):
        if os.path.exists("_last_restart"):
            #shutil.copy("_last_restart", "_last_old_restart")
            os.remove("_last_restart")

        shutil.copy("restart_velocity_cycle", "_last_restart")
        os.remove("restart_velocity_cycle")


    logger.debug("Calling complex.linux")
    try:
        # There is an issue here with complex.linux not actually running
        # It is coming from the mol2 of the substrate...interesting...
        # There is a Segmentation fault (core dumped) error that occurs, asking Jack if he knows what the issue is
        # Jack thinks it is the segfault mike wrote in his HACK ALERT section
        # I have emailed the Dohkyan group regarding it...its only for certain pdbs...
        with Popen(
                f"{os.path.join(phd_config['PATHS']['DMD_DIR'], 'complex.linux')} -P {phd_config['PATHS']['parameters']} -I {pdbName} -T topparam -D 200 -p param -s state -C inConstr -c outConstr 2> complex.linux.err",
                stdout=PIPE, stderr=subprocess.STDOUT, universal_newlines=True, shell=True, env=os.environ) as shell:
            while shell.poll() is None:
                logger.debug(shell.stdout.readline().strip())

        attempts = 100
        if not os.path.isfile("state"):
            logger.warning("Could not make state file first time around, could be a segfault error")
            logger.warning("Going to reorder the bond list in the mol2 files!")
            attempts = 1

        while not os.path.isfile("state") and attempts < 100:
            logger.warning(f"complex fix attempt {attempts}")
            mol2_files = []
            with open("topparam") as topparm:
                for line in topparm:
                    mol2_files.append(line.split()[2])

            for mol2 in mol2_files:
                with open(mol2, 'r') as mf:
                    bonds = []
                    save = []
                    bond_section = False
                    for line in mf:
                        if not bond_section and "BOND" in line:
                            save.append(line)
                            bond_section = True
                            continue

                        if bond_section:
                            bonds.append(line.split())

                        else:
                            save.append(line)

                    # now we reorder based upon the attempt!
                    shuffle(bonds)

                with open(mol2, 'w+') as mf:
                    for line in save:
                        mf.write(line)

                    for bond in bonds:
                        mf.write(f"{bond[0]}\t{bond[1]}\t{bond[2]}\t{bond[3]}\n")

            with Popen(
                    f"{os.path.join(phd_config['PATHS']['DMD_DIR'], 'complex.linux')} -P {phd_config['PATHS']['parameters']} -I {pdbName} -T topparam -D 200 -p param -s state -C inConstr -c outConstr 2> complex.linux.err",
                    stdout=PIPE, stderr=subprocess.STDOUT, universal_newlines=True, shell=True,
                    env=os.environ) as shell:
                while shell.poll() is None:
                    logger.debug(shell.stdout.readline().strip())

            attempts += 1

        if not os.path.isfile("state"):
            logger.critical("could not create state file, something is very wrong!")
            raise ValueError("state file")

    except OSError:
        logger.exception("Error calling complex-1.linux!")
        raise

    logger.debug("Made the state file!")

    if os.path.exists("state"):
        '''
        if os.path.exists("_last_new_state"):
            new_count = 1
            new_count_name = "_last_new_state_" + str(new_count)
            while os.path.exists(new_count_name) == True:
                new_count += 1
                new_count_name = "_last_new_state_" + str(new_count)
            shutil.copy("state", new_count_name)
        else:
            shutil.copy("state", "_last_new_state")
        '''

        # Copy velocities from restart into the new state file
        if os.path.exists("_last_restart"):
            shutil.copy("state", "_last_state")
            copy_restart_velocities("state", "_last_restart")


def make_movie(initial_pdb, movie_file, output_pdb, protonate=[]):
    """

    :param initial_pdb: name of the initial pdb for the dmd run
    :param movie_file: name of the movie file created from dmd
    :param output_pdb: name of the output pdb that is generated from the movie file
    :return:
    """
    try:
        logger.debug("Creating movie file")
        with Popen(
                f"{os.path.join(phd_config['PATHS']['DMD_DIR'], 'complex_M2P.linux')} {phd_config['PATHS']['parameters']} {initial_pdb} topparam {movie_file} {output_pdb} inConstr > m2p.err",
                stdout=PIPE, stderr=subprocess.STDOUT, bufsize=1, universal_newlines=True, shell=True,
                env=os.environ) as shell:
            while shell.poll() is None:
                logger.debug(shell.stdout.readline().strip())
    except OSError:
        logger.exception("Error calling complex_M2P.linux")
        raise

    if not os.path.isfile(output_pdb):
        raise FileNotFoundError(output_pdb)

    if os.path.isfile("m2p.err"):
        os.remove("m2p.err")

    if protonate:
        proteins = load_movie(output_pdb)
        
        if protonate[0] == "all":
            os.remove(output_pdb)
            for pro in proteins:
                pro = addH(pro)
                pro.write_pdb(name=output_pdb, append=True)
        
        else:
            for p in protonate:
                try:
                    index = int(p)
                    pro = proteins[index]
                
                except (ValueError, IndexError):
                    logger.error(f"Skipping {p}")
                    continue

                pro = addH(pro, remove_h=False)
                pro.write_pdb(name=f"{output_pdb.split('.')[0]}_{index}.pdb")


def load_movie(movie_file:str):
    if not os.path.isfile(movie_file):
        logger.error(f"File does not exist: {movie_file}")
        raise ValueError(movie_file)

    #ENDMDL is what seperates the proteins in the movie file
    proteins = []
    try:
        with open(movie_file, 'r') as mf:
            chains = []
            resNum = 0
            chainLet = ""

            protein_number = 0

            for line in mf:
                try:
                    if "ATOM" == line[0:4] or "HETATM" == line[0:6]:
                        tmpAtom = atom.Atom(line)
                        if chainLet != line[21:22]:
                            chainLet = line[21:22]
                            chains.append(chain.Chain(chainLet))
                            resNum = 0

                        if resNum != int(line[22:26]):
                            resNum = int(line[22:26])
                            chains[-1].add_residue(residue.Residue(line))
               
                        chains[-1].residues[-1].add_atom(tmpAtom)

                    elif "ENDMDL" in line:
                        if chains:
                            proteins.append(protein.Protein(f"{movie_file.split('.')[0]}_{protein_number:0>4d}", chains.copy()))   

                        else:
                            logger.warn("Empty chain while loading in movie")
                            protein_number -= 1

                        chains = []
                        resNum = 0
                        chainLet = ""
                        protein_number += 1

                except ValueError:
                    logger.exception(f"Error reading in model {protein_number} in {movie_file}")
                    raise

    except IOError:
        logger.exception(f"Error opening {file}")
        raise

    logger.debug("Successfully loaded in the file!")
    if os.path.isfile("protonation_states.json"):
        logger.debug("Found protonation_states file")
        logger.debug("Loading in protonation states for each structure")
        with open("protonation_states.json") as protonation_states:
            states = json.load(protonation_states)

        for struct, state in zip(proteins, states):
            struct.protonation_states = state

    return proteins

def last_frame(movie_file):
    return load_movie(movie_file)[-1]

def print_header():
    main_logger = logging.getLogger("phd3")

    main_logger.info(datetime.now())
    main_logger.info("")
    main_logger.info("==============================================================================")
    main_logger.info("")
    main_logger.info("                      _______  __   __  ______    _______ ")
    main_logger.info("                     |       ||  | |  ||      |  |       |")
    main_logger.info("                     |    _  ||  |_|  ||  _    | |___    |")
    main_logger.info("                     |   |_| ||       || | |   |  ___|   |")
    main_logger.info("                     |    ___||       || |_|   | |___    |")
    main_logger.info("                     |   |    |   _   ||       |  ___|   |")
    main_logger.info("                     |___|    |__| |__||______|  |_______|")
    main_logger.info("")
    main_logger.info("")
    main_logger.info("--------------------[Protein Hybrid Discrete Dynamics/DFT]--------------------")
    main_logger.info("")
    main_logger.info("[Version]            ==>>    1.1.0")
    main_logger.info("")
    main_logger.info("[Idea and Director]  ==>>    Anastassia N. Alexandrova ")
    main_logger.info("[Idea and Director]  ==>>    Manuel Sparta")
    main_logger.info("[Program Developer]  ==>>    Matthew R. Hennefarth")
    main_logger.info("[Titrate Developer]  ==>>    David J. Reilley")
    main_logger.info("")
    main_logger.info("[Research Group]     ==>>    Alexandrova Research Group")
    main_logger.info("                     ==>>    University of California, Los Angeles")
    main_logger.info("")
    main_logger.info("---------------------------------[References]---------------------------------")
    main_logger.info("")
    main_logger.info(">>>> Quantum Mechanics/Discrete Molecular Dynamics (QM/DMD) >>>>")
    main_logger.info("==>> M. Sparta, D. Shirvanyants, F. Ding,") 
    main_logger.info("     N. V. Dokholyan, A. N. Alexandrova")
    main_logger.info("     Hybrid Dynamics Simulation Engine for Metalloproteins")
    main_logger.info("     Biophys J. 103: 4 (2012)")
    main_logger.info("==>> N. M. Gallup, A. N. Alexandrova")
    main_logger.info("     Use of QM/DMD as a Multiscale Approach to Modeling Metalloenzymes")
    main_logger.info("     Methods Enzymol. 577 (2016)")
    main_logger.info("")
    main_logger.info(">>>> Titratable DMD >>>>")
    main_logger.info("==>> D. J. Reilley, A. N. Alexandrova")
    main_logger.info("     Manuscript in Preparation (2020).")
    main_logger.info("")
    main_logger.info(">>>> piDMD >>>>")
    main_logger.info("==>> D. Shirvanyants, F. Ding, D. Tsao,")
    main_logger.info("     S. Ramachandran, N. V. Dokholyan")
    main_logger.info("     DMD: an Efficient and Versatile Simulation Method for Fine")
    main_logger.info("     Protein Characterization")
    main_logger.info("     J. Phys. Chem. 116: 29 (2012)")
    main_logger.info("")
    main_logger.info(">>>> TURBOMOLE >>>>")
    main_logger.info("==>> R. Ahlrichs, M. Baer, M. Haeser, H. Horn,")
    main_logger.info("     C. Koelmel")
    main_logger.info("     Electronic structure calculations on workstation")
    main_logger.info("     computers: the program system TURBOMOLE")
    main_logger.info("     Chem. Phys. Lett. 162: 165 (1989)")
    main_logger.info("")
    main_logger.info("==============================================================================")
    main_logger.info("")

def copy_directories(src, dst):
    if not os.path.isdir(dst):
        logger.info(f"Making directory {dst}")
        os.mkdir(dst)

    logger.info(f"Copying files from {os.path.abspath(src)} to {os.path.abspath(dst)}")
    src_files = os.listdir(src)
    for file_name in src_files:
        skip = False
        for ignore in constants.IGNORE_FILES:
            if ignore in file_name:
                skip = True
                break

        if skip or file_name.startswith("."):
            continue
        
        src_file_name = os.path.join(src, file_name)
        dst_file_name = os.path.join(dst, file_name)
        try:
            if os.path.isfile(src_file_name):
                shutil.copy(src_file_name, dst_file_name)

            elif os.path.isdir(src_file_name):
                if os.path.isdir(dst_file_name):
                    shutil.rmtree(dst_file_name)

                shutil.copytree(src_file_name, dst_file_name)

        except:
            logger.warn(f"{file_name} suddently vanished in this air...")

def xyz_to_coord(xyz_file):
    with Popen(f"x2t {xyz_file} > coord", shell=True, universal_newlines=True, stdin=PIPE,
            stdout=PIPE, stderr=PIPE, bufsize=1, env=os.environ) as shell:
        while shell.poll() is None:
            logger.info(shell.stdout.readline().strip())
            logger.info(shell.stderr.readline().strip())


def add_proton(atm, ID="HY"):
    vectors = []
    for a in atm.bonds:
        vectors.append(a.coords - atm.coords)
        vectors[-1] = vectors[-1]/np.linalg.norm(vectors[-1])

    #If this is a carboynl, it adds the hydrogen linearly...not ideal, but should be good enough...
    direction = -1.1* sum(vectors)

    new_proton = atom.Atom(element = "H", coords = atm.coords + direction, id=ID, number=20)

    #Update the bond lists etc...
    atm.add_bond(new_proton)
    atm.residue.add_atom(new_proton)

    return new_proton


def addH(protein):
    phd_config = load_phd_config()
    chimera_path = phd_config["PATHS"]["chimera"]

    protein.reformat_protein()
    protein.write_pdb("_temp.pdb", exclude_sub_chain=True, hydrogens=False)  # hydrogens=False is Jack's edit

    with open("chimeraaddh.com", "w") as comfile:
        comfile.write("open _temp.pdb\n")
        comfile.write("addh\n")
        comfile.write("write format pdb 0 addh.pdb\n")
        comfile.write("stop\n")

    try:
        with Popen(f"{os.path.join(chimera_path, 'chimera')} --nogui chimeraaddh.com",
                stdout=PIPE, stderr=subprocess.STDOUT, universal_newlines=True, shell=True, env=os.environ) as shell:
            while shell.poll() is None:
                logger.debug(shell.stdout.readline().strip())

    except OSError:
        logger.exception("Error calling chimera")
        raise

    if not os.path.isfile("addh.pdb"):
        logger.exception("Could not call chimera, check your path")
        raise OSError("chimera")

    logger.debug("Removing chimeraddh.com")
    os.remove("chimeraaddh.com")
    logger.debug("Removing _temp.pdb")
    os.remove("_temp.pdb")

    new_protein = load_pdb("addh.pdb")
    new_protein.chains.append(protein.sub_chain)
    new_protein.reformat_protein()

    logger.debug("Removing addh.pdb")
    os.remove("addh.pdb")

    #Remove all epsilon hydrogens on the histidines
    for chain in new_protein.chains:
        for res in chain.residues:
            if res.name.upper() == "HIS":
                epsilon_nitrogen = res.get_atom("NE2")
                deleted_hydrogen = False
                for atom in epsilon_nitrogen.bonds:
                    if atom.element.lower() == "h":
                        #DELETE THIS ATOM, no good two-timer
                        epsilon_nitrogen.bonds.remove(atom)
                        res.atoms.remove(atom)
                        del atom
                        deleted_hydrogen = True

                if not deleted_hydrogen:
                    logger.warn(f"Could not find his ({res}) 2HNE atom, searching again!")
                    for a in res.atoms:
                        if a.id.lower() == "2hne" or a.id.lower() == "he2":
                            #bastard got through
                            logger.warn("Found the hydrogen!")
                            for b in a.bonds:
                                b.bonds.remove(a)

                            res.atoms.remove(a)
                            del a
                            break

                #Now we add a proton to the ND1 atom
                delta_nitrogen = res.get_atom("ND1")
                skip = False
                for a in res.atoms:
                    if a.id.lower() == "hd1":
                        skip = True

                for a in delta_nitrogen.bonds:
                    if a.element.lower() == "h":
                        skip = True

                if not skip:
                    logger.debug(f"Adding delta nitrogen to {res}")
                    add_proton(delta_nitrogen, ID="HD1")

    #For any added protons (fix the labeling and numbering)
    new_protein.reformat_protein()

    return new_protein


