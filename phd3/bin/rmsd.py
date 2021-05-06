#!/usr/bin/env python3

import logging
import os
import sys
import json
import argparse

from phd3.utility import utilities

def main():
    try:
        utilities.load_logger_config()

    except ValueError:
        print("CRITICAL: Created .phd3 in the root")
        sys.exit(1)

    logger = logging.getLogger(__name__)

    proteins = utilities.load_movie("movie.pdb")
    rmsd = []
    for protein in proteins:
        rmsd.append(proteins[0].aa_rmsd(protein))

    print(rmsd)


    logger.info("[Finished]")

if __name__ == "__main__":
    main()
