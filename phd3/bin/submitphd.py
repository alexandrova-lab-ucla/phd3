#!/usr/bin/env python3

"""Submits a TURBOMOLE job to a queuing system using a template submit script"""
# Give it some workybois and some time

import argparse
import datetime
import os
import sys
import logging
from subprocess import Popen, PIPE


import jinja2


from phd3.utility import utilities, constants

def main(_cores: int=None, _time: int=None, _nodes: int=1, _sub: bool=True):
    # Load in the logger
    utilities.load_logger_config()
    logger = logging.getLogger(__name__)

    # Load in the configuration for this user
    logger.debug("Calling load_phd_config")
    config = utilities.load_phd_config()

    # Parse the configuration!
    logger.debug("Assigning variables for users queueing")
    max_nodes = config['QUEUING']['max_nodes']
    node_types = config['QUEUING']['node_types']
    max_time = config['QUEUING']['max_time']
    high_priority_nodes = config['QUEUING']['high_priority_nodes']

    # Creates a list of (num_cores, max_time)
    avail_nodes = list(zip(node_types, max_time))

    # Assign the vars passed
    logger.debug("Checking for any passed args")
    cores = _cores
    time = _time
    nodes = _nodes
    sub = _sub

    # If no args passed, check for command line args
    if cores is None or time is None:
        # Arg Parser Here
        logger.debug("Parsing args")
        parser = argparse.ArgumentParser(description="Submits a QM/DMD job to a queuing system")
        parser.add_argument("-N", dest="nodes", type=int, required=False, default=1,
                            help=f"The number of nodes to run the QM/DMD job on.")
        parser.add_argument("-n", dest="cores", type=int, required=True,
                            help=f"The number of cores per node to run the QM/DMD job with.")
        parser.add_argument("-t", dest="time", type=int, required=True,
                            help="The amount of time to submit the job for (hours).")
        parser.add_argument("-sub", dest="sub", action='store_true',
                            help="Submit the job to the queueing system.")
        args = parser.parse_args()

        if args.nodes > max_nodes or args.cores > max(node_types) or args.time > max(max_time):
            logger.error("Invalid args parsed based off of what user has specified in config file")
            parser.print_help()
            sys.exit(1)

        cores = args.cores
        time = args.time
        nodes = args.nodes
        sub = args.sub

    # Finds the appropriate node to use
    logger.debug("Checking to ensure user has enough resources")
    big_enough_nodes = list(filter(lambda i: i[0] >= cores, avail_nodes))
    enough_time_nodes = list(filter(lambda i: i[1] >= time, big_enough_nodes))

    possible_nodes = sorted(enough_time_nodes, key=lambda x: x[0])
    if not possible_nodes:
        logger.error("User does not have access to resources")
        raise ValueError

    # Grab the minimum processing node
    node_type = min(possible_nodes)[0]
    logger.debug(f"Using node type: {node_type}")

    high_priority = True if node_type in high_priority_nodes else False

    # Create the jinja2 dictionary
    job_params = {
        "date": datetime.datetime.now(),
        "submit_dir": os.path.realpath(os.getcwd()),
        "node_type": node_type,
        "cores": cores,
        "nodes": nodes,
        "high_priority": ",highp" if high_priority else "",
        "user": os.environ["USER"],
        "job_name": os.path.basename(os.getcwd()),
        "run_script": os.path.join(os.path.dirname(__file__), "runphd.py"),
        "time": time
    }

    # Make the submit script for the queuing system
    home = os.path.expanduser("~")
    path_to_template = os.path.join(home, ".config/phd3")
    
    if not os.path.isdir(path_to_template):
        path_to_template = os.path.join(home, ".phd3")

    logger.debug(
        f"Finding file in: {path_to_template} ")
    templateLoader = jinja2.FileSystemLoader(searchpath=path_to_template)
    templateEnv = jinja2.Environment(loader=templateLoader)
    TEMPLATE_FILE = "submit.j2"
    try:
        template = templateEnv.get_template(TEMPLATE_FILE)

    except jinja2.exceptions.TemplateNotFound:
        logger.error(f"Template file not found in {path_to_template}")
        raise

    logger.debug("Found template file!")
    
    # Now we dump text to a file
    outputText = template.render(job_params)
    try:
        logger.debug(f"Writing out submit file: {constants.SUBMIT_FILE_NAME}")
        with open(constants.SUBMIT_FILE_NAME, "w+") as submit_file:
            submit_file.write(outputText)

    except IOError:
        logger.exception("Could not write out submit file!")
        sys.exit(1)

    # Submit the Job!!
    if sub:
        logger.info(">>>> Submitting job to queue >>>>")
        with Popen(f"{config['QUEUING']['submit']} {constants.SUBMIT_FILE_NAME}", shell=True, universal_newlines=True,
                   stdin=PIPE, stdout=PIPE, stderr=PIPE, bufsize=1, env=os.environ) as shell:
            while shell.poll() is None:
                logger.info(shell.stdout.readline().strip())
                logger.info(shell.stderr.readline().strip())


if __name__ == "__main__":
    main()
