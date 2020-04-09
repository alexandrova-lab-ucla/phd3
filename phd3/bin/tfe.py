#!/usr/bin/env python3

import sys
import logging
import argparse

from phd3 import free_energy
from phd3.utility import utilities


def main():

    try:
        utilities.load_logger_config()

    except ValueError:
        print("CRITICAL: Created .phd3 in the root.")
        sys.exit(1)

    logger = logging.getLogger(__name__)

    # This is the logger for the free_energy script. We now need to remove all of the handlers and append the typical sys.std to it as this is a command line script
    internal_logger = logging.getLogger("phd3.free_energy")
   
    for hdlr in internal_logger.handlers[:]:
        internal_logger.removeHandler(hdlr)
    
    appendHandler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter("\t%(message)s")
    appendHandler.setFormatter(formatter)
    internal_logger.addHandler(appendHandler)

    internal_logger.setLevel(logging.INFO) 

    # Now we have switched out the logger for the sys.stoud so that everything at info or above is output to the console
    logger.debug("Logger is setup!")

    logger.debug("Creating argparser")
    parser = argparse.ArgumentParser(description="Calculates the free energy correction following a Vibrational Frequency calculation")
    parser.add_argument("-T", dest="temp", type=float, required=False, default=298.15,
                        help="The temperature at which to calculate the free energy correction")

    args = parser.parse_args()

    if args.temp < 0:
        logger.error("Cannot have a negative temperature!")
        sys.exit(1)

    try:
        gcorr = free_energy.free_energy_correction(args.temp)

    except:
        #TODO Make this except clause better...not as general
        logger.exception("Error while calculating free energy")
        sys.exit(1)

    print(f"\tFree energy correction (Hartrees): {gcorr}")

    logger.debug("Finished with Calculation")


if __name__ == "__main__":
    main()
