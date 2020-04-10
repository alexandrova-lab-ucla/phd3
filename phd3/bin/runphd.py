#!/usr/bin/env python3

import logging
import argparse

import phd3.utility.utilities as utilities
from phd3 import controller

logger = logging.getLogger(__name__)

def main():
    # Load the logger
    utilities.load_logger_config()
    logger = logging.getLogger(__name__)

    # Sets up the argument parser for the necessary vars
    logger.debug("Parsing arguments")
    parser = argparse.ArgumentParser(description="Runs a QM/DMD job")
    parser.add_argument('-n', nargs=1, dest="cores", type=int, required=True,
                        help='number of cores available')
    parser.add_argument('-t', nargs=1, dest="time", type=int, default=[-1], required=False,
                        help='time available to run job (hours)')
    parser.add_argument('-s', nargs=1, dest="scratch_directory", default='./', type=str, required=False,
                        help='scratch directory to run the job in')

    args = parser.parse_args()

    logger.debug("Proper arguments passed")

    try:
        utilities.print_header()
        job = controller(cores=args.cores[0], time=args.time[0], scratch=args.scratch_directory[0])

    except:
        logger.exception("Unknown exception encountered while trying to run job")
        raise

if __name__ == "__main__":
    main()
