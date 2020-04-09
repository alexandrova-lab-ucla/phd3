#!/usr/bin/env python3

"""Command line script to run a turbomole calculation"""

import argparse
import logging


from phd3.utility import utilities, exceptions
from phd3 import TMcalculation

def main():
    # Load the logger
    utilities.load_logger_config()
    logger = logging.getLogger(__name__)

    # Sets up the argument parser for the necessary vars
    logger.debug("Parsing arguments")
    parser = argparse.ArgumentParser(description="Runs a TURBOMOLE job")
    parser.add_argument('-n', nargs=1, dest="cores", type=int, required=True,
                        help='number of cores available to TURBOMOLE')
    parser.add_argument('-t', nargs=1, dest="time", type=int, default=[-1], required=False,
                        help='time available to run job (hours)')
    parser.add_argument('-s', nargs=1, dest="scratch_directory", default='./', type=str, required=False,
                        help='scratch directory to run the job in')

    args = parser.parse_args()

    logger.debug("Proper arguments passed")

    try:
        job = TMcalculation(cores=args.cores[0], time=args.time[0], run_dir=args.scratch_directory[0])

    except FileNotFoundError as e:
        if str(e) == "definput.json":
            logger.error("There were no parameters for this job!")

        elif str(e) == "coord":
            logger.error("You do not have a structure file!")

        else:
            logger.exception("Unknown FileNotFoundError encountered!")
            raise

    except exceptions.ParameterError:
        logger.error("Parameters passed where invalid! Check to ensure that the value in definput are valid!")

    except:
        logger.exception("Unknown exception encountered while trying to run job")
        raise

if __name__ == "__main__":
    main()
