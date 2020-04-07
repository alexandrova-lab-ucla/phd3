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

    def __init__(self, directory:str, root_dir:str, parameters:dict, cores=1):
        self.directory = directory
        self.root_directory = root_dir
        self.parameters = parameters
        self.cores = cores
        self.stop = False #This checks to see if we have to stop!
        self.final_dmd_average_energy = [0, 0]
        
        # Assign this to the correct next function step
        self.next_step = self.performDMD

        if not os.path.isdir(root_dir):
            logger.error("Root directory for iteration does not exist")
            raise ValueError("root directory does not exist")

        if not os.path.isdir(os.path.abspath(self.directory)):
            try:
                logger.info("Making new iteration directory")
                os.mkdir(self.directory)
          
            except FileExistsError:
                logger.exception("Could not create directory!")
                raise

            self.dmd_structures = None
            self.sp_PDB_structures = None
            self.pdb_winner = None

        else:
            logger.info("Directory exists, loading in the appropriate files")
            self.load_from_directory()
            
    def continue_calculation(self):
        """
        This will look to see what the next step is in the process of performing the calculation, most likely in the
        initialization of the files it will check to see what the last step was and put that into a string
        
        Maybe return the next step, or step to restart on to the main controller?
        """

        logger.debug(f"Changing to iteration directory: {self.directory}")
        os.chdir(self.directory)
   
        #TODO some checks in the directory to update self.next_step

        #Now we do the next step
        while self.next_step is not None and not self.stop:
            self.next_step()

        if self.stop:
            logger.debug("Timer went off, we are ending this iteration early")
            #TODO add some saving things here

        #Now we go back to original directory
        logger.debug(f"Changing to iteration directory: {self.root_directory}")
        os.chdir(self.root_directory)

    def performDMD(self):
        if os.path.isdir("dmd"):
            # Want to continue (check to see if we have first converged, if that is an option as well)
            logger.info("Removing old DMD directory and will restart with the last protein structure")
            logger.info("In the future, we can check this directory and restart internally")
            shutil.rmtree("dmd")

            #TODO get rid of this, and load them in from the directory???
            dmd_average_energies = []
            self.dmd_steps = 0
  
        else:
            dmd_average_energies = []
            self.dmd_steps = 0

        logger.debug("Making dmd directory")
        os.mkdir("dmd")
        logger.debug("Moving to dmd directory")
        os.chdir("dmd")

        logger.debug("Moving the last step pdb to the new dmd directory")
        shutil.copy(self.parameters["last pdb"], ".")
      
        #This is where we add any annelaing or equilibration stuff here
        #TODO add annealing/equilibration

        if self.parameters["DMD CONVERGE"]:
            while not self.stop and self.dmd_steps <= self.parameters["MAX DMD STEPS"]:
                #TODO try and except the following
                os.mkdir(f"dmdstep_{self.dmd_steps}")
                if os.path.isfile(self.parameters["dmd params"]["Restart File"]):
                    logger.debug("Copying restart file")
                    shutil.copy(self.parameters["dmd params"]["Restart File"], f"dmdstep_{self.dmd_steps}/")
                
                shutil.copy("initial.pdb", f"dmdstep_{self.dmd_steps}/")
                os.chdir(f"dmdstep_{self.dmd_steps}")

                calc = dmdpy.calculation(cores = self.cores, parameters = self.parameters["dmd params"])

                curr_energy = calc.get_average_energy()

                dmd_simulations.append(curr_energy)
                #if converged then break
                if len(dmd_simulations) < 2:
                    prev_energy = curr_energy
                    continue

               
                if abs(curr_energy[0]-prev_enery)[0] < prev_energy[1] and abs(curr_energy[0]-prev_energy[0]) < curr_energy[1]:
                    logger.debug("We have converged")
                    if prev_energy[0] < curr_energy[0]:
                        logger.debug("Using previous energy")
                        dmd_simulations.pop()
                        self.final_dmd_average_energy = prev_energy
                        os.chdir(f"../dmdstep_{self.dmd_steps - 1}")

                    else:
                        logger.debug("Using current energy")
                        self.final_dmd_average_energy = curr_energy
            
                    #TODO try and except the following
                    shutil.copy(self.parameters["dmd params"]["Restart File"], "../")
                    shutil.copy(self.parameters["dmd params"]["Echo File"], "../")
                    shutil.copy(self.parameters["dmd params"]["Movie File"], "../")
                    os.chdir("..")
                    break

                logger.debug("Not converged yet")
                
                logger.debug("Copying files up")
                #TODO try and except the following
                shutil.copy(self.parameters["dmd params"]["Restart File"], "../")
                shutil.copy(self.parameters["dmd params"]["Echo File"], "../")
                shutil.copy(self.parameters["dmd params"]["Movie File"], "../")
                os.chdir("..")
                self.dmd_steps += 1

        else:
            calc = dmdpy.calculation(cores = self.cores, parameters = self.parameters["dmd params"])
            self.final_dmd_average_energy = calc.get_average_energy()
            print(self.final_dmd_average_energy)            

        if self.stop:
            logger.info("Stop signal received while in DMD loop")
            #TODO save this information?
            #create a dmd_energy file that has all dmd_average energies
            #could also save this to the phd_control file...
            logger.info("Returning")
            return

        if self.dmd_steps > self.parameters["MAX DMD STEPS"]:
            logger.info("Max DMD steps taken, continuing on to QM calculation")
            logger.info("Finding lowest energy DMD")
            min_step = min(dmd_simulations, key = lambda step: step[0])
            self.final_dmd_average_energy = min_step
            min_step_index = dmd_simulations.index(min_step) + 1
            logger.info(f"Lowest DMD simulation: {min_step_index}")
            logger.info(f"Energy: {min_step}")

            logger.debug(f"Moving to dmdstep_{min_step_index}")
            os.chdir(f"dmdstep_{min_step_index}")
            
            logger.debug("Copying files up")
            shutil.copy(self.parameters["dmd params"]["Restart File"], "../")
            shutil.copy(self.parameters["dmd params"]["Echo File"], "../")
            shutil.copy(self.parameters["dmd params"]["Movie File"], "../")

        elif self.parameters["DMD CONVERGE"]:
            logger.info(f"DMD converged in: {self.dmd_steps} steps")

        #convert movie to a movie.pdb file
        dmdpy.utility.utilities.make_movie("initial.pdb", self.parameters["dmd params"]["Movie File"], "movie.pdb")

        self.dmd_structures = dmdpy.utility.utilities.load_movie("movie.pdb")

        self.dmd_structures[0].write_pdb()

        os.remove("movie.pdb")
        os.chdir("../")

        self.next_step = self.cluster

    def cluster(self):
        if self.dmd_structures is None or len(self.dmd_structures) == 0:
            logger.error("No DMD structures generated!")
            raise ValueError("No DMD structures")

        # Create a zero  matrix of the appropriate size, will fill in with the RMSD between structures
        distance_matrix = np.array([np.array([0.0 for i in self.dmd_structures]) for j in self.dmd_structures])

        # This is a symmetrix matrix
        for row, first_structure in enumerate(self.dmd_structures):
            for col in range(row, len(self.dmd_structures)):
                distance_matrix[row][col] = first_structure.aa_rmsd(self.dmd_structures[col])
                distance_matrix[col][row] = distance_matrix[row][col] 


        print(distance_matrix)

        #TODO have it actually cluster and get some amount of structures
        self.sp_PDB_structures = self.dmd_structures[:3]

        self.next_step = self.qm_singlepoints

    def qm_singlepoints(self):
        if self.sp_PDB_structures is None or len(self.sp_PDB_structures) == 0:
            logger.error("No PDB structures for single point calculations")
            raise ValueError("No single point structures")

        #TODO Ensure that scoring values are saved in some dictionary?
        #TODO set self.pdb_winner to the correct winning pdb structure

        self.next_step = self.qm_optimization

    def qm_optimization(self):
        if self.pdb_winner is None:
            logger.error("No PDB structure for qm optimization!")
            raise ValueError("No qm optimization pdb structure")


        self.next_step = self.finish_iteration

    def finish_iteration(self):
        #Do some finishing stuff here...
        #clean up this directory
        #write out the data (energies)
        self.next_step = None

    def load_from_directory(self):
        pass


