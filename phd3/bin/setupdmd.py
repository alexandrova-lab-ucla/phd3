#!/usr/bin/env python3

import logging
import sys


import phd3.utility.utilities as utilities
from phd3.setupjob import setupDMDjob
from phd3.utility.exceptions import ParameterError


def main():

    try:
        utilities.load_logger_config()

    except ValueError:
        print("CRITICAL: Created .phd3 in the root")
        sys.exit(1)

    logger = logging.getLogger(__name__)

    try:
        sdj = setupDMDjob()
        sdj.full_setup()

    except OSError:
        logger.exception("OSError encountered, likely an issue with moving files around")
        logger.critical("If the issue persists, contact the developers")
        sys.exit(1)

    except ValueError as e:
        if "dmdinput.json" in str(e):
            logger.error("Error with the dmdinput.json file")
            sys.exit(1)

        elif "definition" in str(e):
            logger.error("Please provide the correct parameter definition!")
            sys.exit(1)

        elif "No Protein" in str(e):
            logger.error("No protein or pdb was provided")
            sys.exit(1)

        else:
            logger.error("Unknown ValueError encountered")
            raise

    except ParameterError:
        logger.error("Please make sure your parameters are correct")
        sys.exit(1)

    except:
        logger.exception("Unknown exception encountered, quitting")
        sys.exit(1)


if __name__ == "__main__":
    main()

