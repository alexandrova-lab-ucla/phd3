#!/usr/bin/env python3

"""Sets up a turbomole calculation by running define and editing the control file"""

import logging
import sys
import argparse

from phd3 import setupTMjob
from phd3.utility import utilities, exceptions

sys.settrace

def main():
    # Begin the logger (__main__)
    try:
        utilities.load_logger_config()

    except ValueError:
        print("CRITICAL: Created .phd3 in the root.")
        sys.exit(1)

    logger = logging.getLogger(__name__)

    logger.debug("Parsing Arguments")
    parser = argparse.ArgumentParser(description="Setups a turbomole job")
    
    parser.add_argument("-t",dest="timeout", action='store_false', default=True, required=False, help="Turns off timeout")
    
    args = parser.parse_args()
    
    timeout_time = 0
    if args.timeout:
        logger.debug("Timeout is turned on")
        timeout_time = 15
    
    # Start to do the dew
    try:
        tm = setupTMjob(timeout = timeout_time)

    except FileExistsError as e:
        if str(e) == 'control':
            logger.info("Either delete control file, or perform define by hand!")

        else:
            logger.exception("Unknown FileExistsError encountered!")

    except FileNotFoundError as e:
        if str(e) == "definput.json":
            logger.error("Please resubmit once you are done editing the parameters")

        elif str(e) == "coord":
            logger.info("Please provide a structure file for this job!")

        else:
            logger.exception("Unknown FileNotFoundError encountered!")

    except exceptions.DefineError:
        logger.info("Check define.out for the exact issue!")

    except exceptions.ParameterError as e:
        logger.info(e)

    except:
        logger.exception("Exception encountered in running setupjob.py")

if __name__ == "__main__":
    main()
