#!/usr/bin/env python3

"""Calculated the Free Energy correction from a Numforce calculation"""

import math
import os
import logging
import subprocess

import phd3.utility.constants as constant
from phd3.utility import utilities

__all__ = [
    'gcorrtrans',
    'gcorrvib',
    'free_energy_correction'
]


logger = logging.getLogger(__name__)

def gcorrtrans(atoms = None, temp:float = 298.15):
    """
    calculated the free energy correction of a molecules Free energy from only the translational energy. Typically for
    single atom systems.

    :param atoms: List of the atomic elements in the system
    :param temp: Temperature at which the correction should be calculated at
    :return: Correction to the Gibbs Free Energy, from translation only, in Hartrees
    """

    if atoms is None:
        if not os.path.isfile("coord"):
            logger.error("No coord file present, cannot calculate translational contribution!")
            raise ValueError("no coord file")

        logger.debug("Opening the coord file")
        atoms = []
        coord_section = False
        try:
            with open("coord") as coordfile:
                for line in coordfile:
                    if "$coord" in line:
                        coord_section = True
                        continue

                    elif "$" in line:
                        break

                    if coord_section:
                        line = line.split()
                        atoms.append(line[3])

        except IOError:
            logger.exception("Error opening the coord file")
            raise

    total_mass = 0
    logger.debug("Calculating total mass of the system")
    for atom in atoms:
        try:
            total_mass += constant.ATOM_MASS[atom]

        except:
            logger.exception("Atom is missing its mass definition in turbopy, sad day")
            raise

    return -3.166811e-6 * temp * math.log(constant.Kb * temp / 101325 * math.pow(2 * math.pi * total_mass * 1.6605e-27 * constant.Kb * temp/math.pow(6.626e-34, 2), 1.5))


def gcorrvib(temp:float = 298.15):
    """
    Calculates the vibrational correction to the free energy using the spectrum data in the vibspectrum output file
    :param temp: Temperature at which to do the correction, default is room temp
    :return: returns the Gibbs Free Energy correction in Hartrees
    """
    logger.debug("Calculating the vibrational correction to the Free Energy Manually")
    spectrumStarted = False
    freq = []

    if not os.path.isfile("vibspectrum"):
        logger.error("vibspectrum does not exist")
        raise ValueError("vibspectrum")

    try:
        with open('vibspectrum') as inputfile:
            for line in inputfile:
                if line.split()[0] == "1":
                    spectrumStarted = True
                if spectrumStarted:
                    if len(line.split()) == 6:
                        frequency = float(line.split()[2])
                        if frequency < 0.1:
                            logger.info(f"Imaginary frequency detected: {frequency}")

                        elif frequency < 10000:
                            freq.append(frequency)

    except IOError:
        logger.exception("Could not open 'vibspectrum'")
        raise

    logger.debug("Converting frequencies to a unit in Hartrees")
    for i, w, in enumerate(freq):
        freq[i] = w/219474.63 #To hartrees

    gcorr = 0.0

    kbT = 3.166811e-6 * temp
    for e in freq:
        gcorr += e / 2 + kbT * math.log1p(-math.exp(-e/kbT))

    return gcorr

def free_energy_correction(temp:float = 298.15):
    vibscale = 1
    gcorr = 0

    frozen_coords = False
    if not os.path.isfile("coord"):
        logger.exception("No coord file!")
        raise ValueError("coord")

    try:
        with open('coord') as coordFile:
            numatoms = -1
            for line in coordFile:
                if line[0] == '$' and "$coord" not in line:
                    logger.debug("Finished in the coordinate section")
                    break
                if len(line.split()) > 4:
                    logger.info("Frozen Coordinates detected; calculating vibrational free energy correction")
                    frozen_coords = True
                    gcorr = gcorrvib(temp)
                    break
                numatoms += 1

    except:
        logger.exception("Error reading the coord file!")
        raise

    if not frozen_coords:
        if numatoms == 1:
            logger.info("Single atom detected; calculating translational free energy correction")
            gcorr = gcorrtrans()

        else:
            with open("vibspectrum") as vibfile:
                for line in vibfile:
                    if not line.startswith(('$', '#')) and float(line[20:].split()[0]) < 0:
                        logger.warning(f"Imaginary frequency detected: {line[20:].split()[0]}")

            logger.debug("Setting up environment to call freeh")
            config = utilities.load_phd_config()
            utilities.setup_turbomole_env(config["PATHS"]["TURBODIR"])

            try:
                logger.debug("Calling freeh")
                freehrun = subprocess.run(['freeh'], input=f'\n{vibscale}\n\nq\n', stdout=subprocess.PIPE, universal_newlines=True)

            except:
                logger.exception("Error in calling 'freeh' program!")
                raise ValueError

            chempotline = False
            linenum=0
            for line in freehrun.stdout.split('\n'):
                if chempotline:
                    linenum+= 1

                elif 'chem.pot.' in line.split():
                    chempotline = True

                if linenum == 3:
                    try:
                        logger.debug("Trying to grab the gcorr value from freeh output")
                        gcorr = float(line.split()[5])/2625.50

                    except ValueError:
                        logger.exception("Error in converting str to float in freeh output")
                        raise

                    break
    if not gcorr:
        raise ValueError("Could not get a free energy correction")

    return gcorr
