#!/usr/bin/env python3
"""
Author  ==>> Matthew R. Hennefarth
Date    ==>> April 16, 2020
"""

#Standard Library Imports
import logging
import os
import signal
import json
import shutil
import sys
from timeit import default_timer as timer

#PHD Imports
from . import iteration
from .bin import submitphd

logger = logging.getLogger(__name__)

__all__ = [
        'controller'
    ]

class controller:

    def __init__(self, cores:int, time:int=-1, scratch:str="./"):
        self._cores = cores
        self._time = time
        self._scratch = scratch
        self._submit_directory = os.getcwd()
        self._iteration = 0
        self._parameters = {}
        self._curr_iterat = None
        self._start = timer()
        self._stop = False


        if not os.path.isfile("phdinput.json"):
            logger.error("PHDinput.json does not exist")
            raise FileNotFoundError("phdinput.jdon")
       
        #TODO verify phdinput.json is valid format
            
        logger.debug("Loading in input file")
        with open("phdinput.json", 'r') as f:
            self._parameters = json.load(f)

        logger.debug("Checking last iteration directory")
        dirs = [d for d in os.listdir() if os.path.isdir(d)]
        iterations = [int(d.split("_")[1]) for d in dirs if "Iteration_" in d]
        if not iterations:
            last_iteration_directory = -1

        else:
            last_iteration_directory = max(iterations)

        logger.debug("Checking last saved iteration")
        self.last_finished_iteration = -1
        if os.path.isfile("phd_energy"):
            energy_lines = []
            with open("phd_energy") as energy:
                for line in energy:
                    energy_lines.append(line)

            #Get rid of the first line
            if len(energy_lines) > 1:
                energy_lines = energy_lines[1:]
                self.last_finished_iteration = len(energy_lines) - 1

        #Then we assume we start the next after last_iteration_directory
        if last_iteration_directory != self.last_finished_iteration:
            #Then we need to start from the last_iteration_directory

            #We should throw an error if this is not the case
            assert(last_iteration_directory > self.last_finished_iteration)
            self._iteration = last_iteration_directory

        else:
            #They are the same, therefore start next iteration
            self._iteration = self.last_finished_iteration + 1

        if "Default MOs" in self._parameters.keys():
            self._parameters["Default MOs"] = os.path.abspath(self._parameters["Default MOs"])

        logger.debug("Finding last to_next_iteration.pdb file")
        if self._iteration > 0 :
            if os.path.isfile(os.path.join(f"Iteration_{self._iteration-1}", "to_next_iteration.pdb")):
                self._parameters["last pdb"] = os.path.abspath(os.path.join(f"Iteration_{self._iteration-1}", "to_next_iteration.pdb"))
                self._parameters["Default MOs"] = os.path.abspath(f"Iteration_{self._iteration-1}/Optimization")
            
            else:
                logger.error("Last to_next_iteration.pdb not found!")
                raise FileNotFoundError("to_next_iteration.pdb")

        else:
            self._parameters["last pdb"] = os.path.abspath(self._parameters["pdb file"])

        while not self._stop and self._iteration <= self._parameters["Max Iterations"]:
            #This is the loop we stay in until we need to quit
            self._curr_iterat = iteration.iteration(self, f"Iteration_{self._iteration}", os.path.abspath("."), self._parameters, self._iteration, self._cores, scratch=self._scratch) 
            self._curr_iterat.continue_calculation()
            
            #Checks to see if the iteration is telling us to stop
            if self._curr_iterat.stop:
                self._stop = self._curr_iterat.stop

            if self._stop:
                self._iteration -=1

            self._iteration += 1

        if self._stop and self._iteration <= self._parameters["Max Iterations"]:
            if self._parameters["Resubmit"]:
                submitphd.main(_cores = self._cores, _time=self._time)

        elif self._iteration > self._parameters["Max Iterations"]:
            logger.info("Finished with all QM/DMD cyles")

    def time_left(self):
        if self._time == -1:
            return -1

        now = timer()
        time_elapsed = int(now-self._start)/3600.0
        if self._time < time_elapsed:
            return 0


        return self._time - time_elapsed 

