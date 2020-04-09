#!/usr/bin/env python3

import logging
import sys
import os


from phd3.utility import utilities
import phd3.dmd_simulation


def main():
    try:
        utilities.load_logger_config()

    except ValueError:
        print("CRITICAL: Created .dmdpy in the root")
        sys.exit(1)

    logger = logging.getLogger(__name__)

    logger.debug("Logger initialized")


    logger.debug("Parsing arguments")
    parser = argparse.ArgumentParser(description="Runs a DMD simulation")
    parser.add_argument('-n', nargs=1, dest="cores", type=int, required=True,
                        help='number of cores available')
    parser.add_argument('-t', nargs=1, dest="time", type=int, default=[-1], required=False,
                        help='time available to run job (hours)')
    parser.add_argument('-s', nargs=1, dest="scratch_directory", default='./', type=str, required=False,
                        help='scratch directory to run the job in')

    args = parser.parse_args()

    logger.debug("Proper arguments passed")

    # TODO have it receive argument for the number of cores, time to run, and a scratch directory

    try:
        logger.debug("Attempting to begin the calculation")
        c = phd3.dmd_simulation(cores=args.cores[0], time=args.time[0], run_dir=args.scratch_directory[0])

    except:
        logger.exception("Check the error")
        logger.error("Error on the calculation")
        sys.exit(1)


if __name__ == "__main__":
    main()
