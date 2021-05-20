#!/usr/bin/env python3

import argparse
import logging
import os
import sys


from phd3.utility import utilities


def main():
    try:
        utilities.load_logger_config()

    except ValueError:
        print("CRITICAL: Created .phd3 in the root")
        sys.exit(1)

    logger = logging.getLogger(__name__)

    logger.debug("Parsing Arguments")

    parser = argparse.ArgumentParser(description="Converts a piDMD movie file to a pdb")

    parser.add_argument("pdbfile", type=str, nargs=1, help="initial pdb file (typically initial.pdb)")
    parser.add_argument("moviefile", type=str, nargs=1, help="Movie file from piDMD")
    parser.add_argument("-o", dest="outputFile", default=["movie.pdb"], type=str, nargs=1, required=False, help="Output for movie file")
#    parser.add_argument("--p", dest="protonate", nargs='*', default=[], required=False, help="Protonates the movie file")
    
    args = parser.parse_args()

    if not os.path.isfile(args.pdbfile[0]):
        logger.error(f"pdb file provided does not exist: {args.pdbfile[0]}")
        sys.exit(1)

    if not os.path.isfile(args.moviefile[0]):
        logger.error(f"movie file provided does not exist: {args.moviefile[0]}")
        sys.exit(1)

    logger.debug("Passing parameters...")
    try:
        logger.debug("Passing args to the function in utilities")
        utilities.make_movie(args.pdbfile[0], args.moviefile[0], args.outputFile[0])
        #utilities.make_movie(args.pdbfile[0], args.moviefile[0], args.outputFile[0], protonate=args.protonate)

    except OSError:
        logger.error("Error creating movie file")
        sys.exit(1)

    #TODO check for the output file actually...
    logger.info("Successfully converted movie file!")
    logger.debug("Finished")


if __name__ == "__main__":
    main()
