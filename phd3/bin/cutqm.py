#!/usr/bin/env python3

import logging
import os
import sys
import json
import argparse

from phd3.utility import utilities
import phd3.dmd_to_qm as dmd_to_qm

def main():
    try:
        utilities.load_logger_config()

    except ValueError:
        print("CRITICAL: Created .phd3 in the root")
        sys.exit(1)

    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(description="Converts pdb -> coord, or reinstalls coord in pdb")

    parser.add_argument("--reinstall", dest="install", required=False, action="store_true", help="Reinstall coord into pdb")

    args = parser.parse_args()

    if not os.path.isfile("phdinput.json"):
        logger.error("No phdinput.json file")
        sys.exit(1)

    with open("phdinput.json", 'r') as optionFile:
        parameters = json.load(optionFile)

    protein = utilities.load_pdb(parameters["pdb file"])
    
    if not args.install:
        logger.info("[Cutting]...")
        dmd_to_qm.protein_to_coord(protein, parameters["QM Chop"])

    else:
        logger.info("Reinstalling...")
        protein = dmd_to_qm.coord_to_protein(protein, parameters["QM Chop"])
        protein.write_pdb("reinstalled.pdb")


    logger.info("[Finished]")

if __name__ == "__main__":
    main()
