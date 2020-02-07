#!/usr/bin/env python3

# This will be created on each iteration

import logging
import os
import shutil
import numpy as np
import scipy.cluster

import dmdpy

logger = logging.getLogger(__name__)

class iteration:

    def __init__(self, directory:str, root_dir:str, parameters:dict):
        self.directory = directory
        self.root_directory = root_dir
        self.parameters = parameters

        if not os.path.isdir(root_dir):
            raise ValueError("root directory does not exist")

        if not os.path.isdir(os.path.abspath(self.directory)):
            try:
                logger.info("Making new iteration directory")
                os.mkdir(self.directory)
          
            except FileExistsError:
                logger.exception("Could not create directory!")
                raise

            self.dmd_structures = None

        else:
            logger.info("Directory exists, loading in the appropriate files")
            
    def continue_calculation(self):
        """
        This will look to see what the next step is in the process of performing the calculation, most likely in the
        initialization of the files it will check to see what the last step was and put that into a string
        """

        logger.debug(f"Changing to iteration directory: {self.directory}")
        os.chdir(self.directory)
    
    def performDMD(self):
        #Check to see if there is any current dmd directory, otherwise delete and restart it!
        if os.path.isdir(os.path.join(os.path.abspath(self.directory), "dmd")):
            logger.info("Removing old DMD directory and will restart with the last protein structure")
            shutil.rmtree(os.path.join(self.directory, "dmd")
  
        logger.debug("Moving the last step pdb to the new dmd directory")
        shutil.copy(self.parameters["last pdb"], os.path.abspath(self.directory))
        
        #TODO update the cores and time with the appropriate value from the controller?????
        dmd_simulation = dmdpy.calculation(cores = 8, time = 12, parameters = self.parameters["dmd params"])
        
        #convert movie to a movie.pdb file
        dmdpy.utility.utilities.make_movie("initial.pdb", self.parameters["dmd params"]["Movie File"], "movie.pdb")
        
        #TODO Need to implement in dmdpy, have it return a list of proteins!!!!!!
        self.dmd_structures = dmdpy.utility.utulities.load_movie("movie.pdb")

        #TODO update the next step so that it was "cluster"


    def cluster(self):
        if self.dmd_structures is None or len(self.dmd_structures) == 0:
            logger.error("No DMD structures generated!")
            raise ValueError("No DMD structures")

        # Create a zero  matrix of the appropriate size, will fill in with the RMSD between structures
        distance_matrix = np.array([np.array([0 for i in self.dmd_structures]) for j in self.dmd_structures])

        # This is a symmetrix matrix
        for row, first_structure in enumerate(self.dmd_structures):
            for col in range(row, len(self.dmd_structure)):
                distance_matrix[row][col] = first_structure.rmsd(self.dmd_structure[col])
                distance_matrix[col][row] = distance_matrix[row][col] 

        #TODO have it actually cluster and get some amount of structures
        possible_structures = self.dmd_structures[:3]
        # Will use the above method as just a quick way of doing this

        

    def qm_singlepoints(self):
        pass

    def qm_optimization(self):
        pass

    def finish_iteration(self):
        pass



    #Dictionary that references the appropriate class method
    steps = {
        "performDMD" : performDMD,
        "cluster" : cluster,
        "qm_singlepoints" : qm_singlepoints,
        "qm_optimization" : qm_optimization,
        "finish_iteration" : finish_iteration
    }




