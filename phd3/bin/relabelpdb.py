#!/usr/bin/env python3

import logging
import argparse
import sys
import os


from phd3.utility import utilities


def main():
    # Load in the logger
    try:
        utilities.load_logger_config()

    except ValueError:
        print("CRITICAL: Created .phd3 in the root")
        sys.exit(1)

    logger = logging.getLogger(__name__)

    logger.debug("Parsing Arguments")

    parser = argparse.ArgumentParser(description="Relabels a pdb to a new atom labeling scheme")
    parser.add_argument("pdbfile", type=str, nargs=1, help="PDB file to edit")
    parser.add_argument("scheme", type=str, nargs=1, help="Scheme to convert to")
    parser.add_argument("-o", dest="outputFile", default="", type=str, nargs=1, required=False, help="PDB file to write to")

    args = parser.parse_args()

    if not os.path.isfile(args.pdbfile[0]):
        logger.error(f"Could not find pdb: {args.pdbfile[0]}")
        sys.exit(1)

    logger.debug("Proper Arguments Passed")

    try:
        protein = utilities.load_pdb(args.pdbfile[0])

    except IOError:
        logger.error("Error in loading in the PDB file provided")
        sys.exit(1)

    logger.debug("Relabeling")
    try:
        protein.relabel(args.scheme[0])

    except ValueError:
        logger.error("Error in relabeling the protein")
        sys.exit(1)

    if args.outputFile != "":
        logger.debug(f"Changing protein name to {args.outputFile[0]}")
        protein.name = args.outputFile[0]

    try:
        protein.write_pdb()

    except IOError:
        logger.error("Error writing the protein to pdb file")
        sys.exit(1)

    logger.debug("Finished")


if __name__ == "__main__":
    main()
