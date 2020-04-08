#!/usr/bin/env python3

# This will be created on each iteration

import logging
import os
import shutil
import numpy as np
import scipy.cluster
import sys
from timeit import default_timer as timer
import datetime

import dmdpy
import turbopy

import pdb_to_coord

logger = logging.getLogger(__name__)

H_TO_KCAL = 627.509

class iteration:

    def __init__(self, directory:str, root_dir:str, parameters:dict, iteration:int, cores=1):
        self.directory = directory
        self.iter_number = iteration
        self.root_directory = root_dir
        self.parameters = parameters
        self.cores = cores
        self.stop = False #This checks to see if we have to stop!
        self.final_dmd_average_energy = [0, 0]
        self.dmd_sp_energies = []
        self.pdb_winner = None
        self.qm_sp_energies = []
        self.scoring_energies = []
        self.qm_final_energy = 0.0
        self.to_next_iteration = None

        # Assign this to the correct next function step
        self.next_step = self.performDMD

        if not os.path.isdir(root_dir):
            logger.error("Root directory for iteration does not exist")
            raise ValueError("root directory does not exist")
        
        logger.info("")
        logger.info("+---------------------------------------------------------+")
        logger.info(f"|                   ITERATION      {self.iter_number:0>2d}                     |")
        logger.info("+---------------------------------------------------------+")
        logger.info("")
        logger.info("")

        if not os.path.isdir(os.path.abspath(self.directory)):
            try:
                logger.debug("Making new iteration directory")
                os.mkdir(self.directory)
          
            except FileExistsError:
                logger.exception("Could not create directory!")
                raise

            self.dmd_structures = None
            self.sp_PDB_structures = None
            self.pdb_winner = None

        else:
            logger.info("Directory exists, loading in previous runs")
            logger.info("")
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
            start = timer()
            self.next_step()
            end = timer()
            logger.info(f"Time elapsed: {datetime.timedelta(seconds = int(end -start))}")
            logger.info("")
            logger.info("")

        if self.stop:
            logger.debug("Timer went off, we are ending this iteration early")
            #TODO add some saving things here

        #Now we go back to original directory
        logger.debug(f"Changing to iteration directory: {self.root_directory}")
        os.chdir(self.root_directory)

    def performDMD(self):
        logger.info("====================[Beginning DMD]====================")
        logger.info("")

        #This tracks to see if we are done with the dmd simulations!
        finished = False
        dmd_average_energies = []
        self.dmd_steps = 1
        if os.path.isdir("dmd"):
            # Want to continue (check to see if we have first converged, if that is an option as well)
            logger.info("dmd Directory found")
            logger.info("Loading in previous DMD simulation data")
            logger.info("")

            logger.debug("Moving to dmd directory")
            os.chdir("dmd")
            if self.parameters["DMD CONVERGE"]:
                dirs = [d for d in os.listdir() if "dmdstep_" in d]
                self.dmd_steps = len(dirs)+1
                #Load in the previous energies
                
                logger.info("      Attempting to converge DMD simulations")
                logger.info("[RUN #]---------[Ave Pot. Energy]-------------")

                for i, d in enumerate(dirs):
                    ave, stdev = dmdpy.calculation.get_average_energy(os.path.join(d, self.parameters["dmd params"]["Echo File"]))
                    dmd_average_energies.append([ave, stdev])
                    logger.info(f"[Run {i+1}]\t==>> {ave:.5f} ({stdev:.5f})\t\t(finished)")
            
                if len(dmd_average_energies) >= 2:
                    if abs(dmd_average_energies[-1][0]-dmd_average_energies[-2][0]) < dmd_average_energies[-2][1] and abs(dmd_average_energies[-1][0]-dmd_average_energies[-2][0]) < dmd_average_energies[-1][1]:
                        logger.info("Convergence previously achieved")
                        logger.info(f"Absolute delta in energy: {abs(dmd_average_energies[-1][0]-dmd_average_energies[-2][0]):.5f}")
                        finished = True
                        self.dmd_steps -= 1

            else:
                if os.path.isfile(self.parameters["dmd params"]["Echo File"]) and os.path.isfile(self.parameters["dmd params"]["Movie File"]):
                    logger.info("[RUN #]---------[Ave Pot. Energy]-------------")
                    ave, stdev = dmdpy.calculation.get_average_energy(self.parameters["dmd params"]["Echo File"])
                    logger.info(f"[Run 0]\t==>> {ave:.5f} ({stdev:.5f})\t\t(finished)")
                    finished = True

  
        else:
            logger.debug("Making dmd directory")
            os.mkdir("dmd")
            logger.debug("Moving to dmd directory")
            os.chdir("dmd")
            logger.debug("Moving the last step pdb to the new dmd directory")
            shutil.copy(self.parameters["last pdb"], "./dmdStart.pdb")
      
        if not finished:
            #TODO set the dmdpy logger to dmd.out in this directory!
            dmd_logger = logging.getLogger("dmdpy")
            dmd_logger.setLevel(logging.CRITICAL)

            #This is where we add any annelaing or equilibration stuff here
            #TODO add annealing/equilibration
            if self.parameters["DMD CONVERGE"]:
                
                while not self.stop and self.dmd_steps <= self.parameters["MAX DMD STEPS"]:
                    #TODO try and except the following

                    if os.path.isdir(f"dmdstep_{self.dmd_steps-1}"):
                        logger.debug(f"Copying {self.dmd_steps-1} files to new dmd start")

                        shutil.copytree(f"dmdstep_{self.dmd_steps-1}", f"dmdstep_{self.dmd_steps}")
                        logger.debug("Removing echo and movie file from new directory")
                        os.remove(os.path.join(f"dmdstep_{self.dmd_steps}", self.parameters["dmd params"]["Echo File"]))
                        os.remove(os.path.join(f"dmdstep_{self.dmd_steps}", self.parameters["dmd params"]["Movie File"]))

                    else:
                        logger.info("      Attempting to converge DMD simulations")
                        logger.info("[RUN #]---------[Ave Pot. Energy]-------------")
                        os.mkdir(f"dmdstep_{self.dmd_steps}")
                        shutil.copy("dmdStart.pdb", f"dmdstep_{self.dmd_steps}/")

                    logger.debug(f"Changing direcetory from {os.getcwd()} to dmdstep_{self.dmd_steps}")
                    os.chdir(f"dmdstep_{self.dmd_steps}")

                    logger.debug("Begining dmdpy")
                    start = timer()
                    calc = dmdpy.calculation(cores = self.cores, parameters = self.parameters["dmd params"])
                    
                    #Harvest the time that the dmd run took
                    end = timer()

                    curr_energy = dmdpy.calculation.get_average_energy(self.parameters["dmd params"]["Echo File"])
                    logger.info(f"[Run {self.dmd_steps}]\t==>> {curr_energy[0]:.5f} ({curr_energy[1]:.5f})")
                    logger.info(f"Time elapsed during DMD simulation: {datetime.timedelta(seconds = int(end -start))}")
                    
                    dmd_average_energies.append(curr_energy)
                    #if converged then break
                    if len(dmd_average_energies) >= 2:
                        if abs(dmd_average_energies[-1][0]-dmd_average_energies[-2][0]) < dmd_average_energies[-2][1] and abs(dmd_average_energies[-1][0]-dmd_average_energies[-2][0]) < dmd_average_energies[-1][1]:
                            logger.info("")
                            logger.info("Converence achieved!")
                            logger.info(f"Absolute delta in energy: {abs(curr_energy[0]-prev_energy[0]):.5f}")
                            if prev_energy[0] < curr_energy[0]:
                                logger.debug("Using previous energy")
                                dmd_average_energies.pop()
                                self.final_dmd_average_energy = dmd_average_energies[-2] 
                                os.chdir(f"../dmdstep_{self.dmd_steps - 1}")

                            else:
                                logger.debug("Using current energy")
                                self.final_dmd_average_energy = curr_energy
                
                            #TODO try and except the following
                            shutil.copy("topparam", "../")
                            shutil.copy("inConstr", "../")
                            shutil.copy(self.parameters["dmd params"]["Restart File"], "../")
                            shutil.copy(self.parameters["dmd params"]["Echo File"], "../")
                            shutil.copy(self.parameters["dmd params"]["Movie File"], "../")
                            os.chdir("..")
                            break

                        else:
                            logger.debug("Not converged yet")

                    logger.debug("Copying files up")
                    #TODO try and except the following
                    shutil.copy("topparam", "../")
                    shutil.copy("inConstr", "../")
                    shutil.copy("initial.pdb", "../")
                    shutil.copy(self.parameters["dmd params"]["Restart File"], "../")
                    shutil.copy(self.parameters["dmd params"]["Echo File"], "../")
                    shutil.copy(self.parameters["dmd params"]["Movie File"], "../")
                    for f in [a for a in os.listdir() if a.endswith(".mol2")]:
                        shutil.copy(f, "..")
                
                    os.chdir("..")
                    self.dmd_steps += 1

            else:
                logger.info("[RUN #]-------[Ave Pot. Energy]-------------")
                start = timer()
                calc = dmdpy.calculation(cores = self.cores, parameters = self.parameters["dmd params"])
                end = timer()
                self.final_dmd_average_energy = dmdpy.calculation.get_average_energy()
                logger.info(f"[Run 0]\t==>> {self.final_dmd_average_energy[0]:.5f} ({self.final_dmd_average_energy[1]:.5f})")
                logger.info(f"Time elapsed during DMD simulation: {datetime.timedelta(seconds = int(end -start))}")

        if self.stop:
            logger.info("Stop signal received while performing DMD simulation")
        
        else:
            if self.dmd_steps > self.parameters["MAX DMD STEPS"]:
                logger.info("")
                logger.info("Max DMD steps taken")
                logger.info("Finding lowest energy DMD simulation")
                min_step = min(dmd_average_energies, key = lambda step: step[0])
                self.final_dmd_average_energy = min_step
                min_step_index = dmd_average_energies.index(min_step) + 1
                logger.info(f"Lowest DMD simulation: {min_step_index}")
                logger.info(f"[Run {min_step_index}]\t==>> {min_step[0]:.5f} ({min_step[1]:.5f})")

                logger.debug(f"Moving to dmdstep_{min_step_index}")
                os.chdir(f"dmdstep_{min_step_index}")
                
                logger.debug("Copying files up")
                shutil.copy(self.parameters["dmd params"]["Restart File"], "../")
                shutil.copy(self.parameters["dmd params"]["Echo File"], "../")
                shutil.copy(self.parameters["dmd params"]["Movie File"], "../")
                shutil.copy("topparam", "../")
                shutil.copy("initial.pdb", "../")
                shutil.copy("inConstr", "../")
                for f in [a for a in os.listdir() if a.endswith(".mol2")]:
                    shutil.copy(f, "..")

                os.chdir("../")

            elif self.parameters["DMD CONVERGE"]:
                logger.info("")
                logger.info(f"DMD converged in: {self.dmd_steps} steps")

            #convert movie to a movie.pdb file
            logger.debug("Making movie.pdb")
            dmdpy.utility.utilities.make_movie("initial.pdb", self.parameters["dmd params"]["Movie File"], "movie.pdb")

            logger.debug("Loading in movie.pdb")
            self.dmd_structures = dmdpy.utility.utilities.load_movie("movie.pdb")

            os.remove("movie.pdb")
            os.chdir("../")

            self.next_step = self.cluster
        
        logger.info("")
        logger.info("====================[Finished  DMD]====================")

    def cluster(self):
        logger.info("================[Beginning  Clustering]================")
        logger.info("")
        if self.dmd_structures is None or len(self.dmd_structures) == 0:
            logger.error("No DMD structures generated!")
            raise ValueError("No DMD structures")

        # Create a zero  matrix of the appropriate size, will fill in with the RMSD between structures
        distance_matrix = np.array([np.array([0.0 for i in self.dmd_structures]) for j in self.dmd_structures])

        # This is a symmetrix matrix
        for row, first_structure in enumerate(self.dmd_structures):
            for col in range(row, len(self.dmd_structures)):
                continue
                distance_matrix[row][col] = first_structure.aa_rmsd(self.dmd_structures[col])
                distance_matrix[col][row] = distance_matrix[row][col] 


        #TODO have it actually cluster and get some amount of structures
        self.sp_PDB_structures = self.dmd_structures[:3]
        #Saving the structures in case we have to restart

        for struct in self.sp_PDB_structures:
            struct.write_pdb(f"{struct.name}.pdb")

        rows = []
        for struct in self.sp_PDB_structures:
            n = int(struct.name.split("_")[1])
            rows.append(n)

        echo_lines = []
        with open("dmd/echo", 'r') as echo:
            for line in echo:
                if line[0] != "#":
                    echo_lines.append(line)

        for r in rows:
            line = echo_lines[r]
            line = line.split()
            self.dmd_sp_energies.append(float(line[4]))

        self.next_step = self.qm_singlepoints
        logger.info("")
        logger.info("=================[Finished Clustering]=================")

    def qm_singlepoints(self):
        logger.info("===================[Beginning QM SP]===================")
        logger.info("")

        #This is where we fix the logger???
        turbo_logger = logging.getLogger("turbopy")
        turbo_logger.setLevel(logging.CRITICAL)

        logger.info("[movie_####]---------[QM Energy (Hart)]------------[DMD Energy (kcal)]------------")
        if self.sp_PDB_structures is None or len(self.sp_PDB_structures) == 0:
            logger.error("No PDB structures for single point calculations")
            raise ValueError("No single point structures")
    
        #Make all of the directories
        for index, struct in enumerate(self.sp_PDB_structures):
            if not os.path.isdir(f"sp_{struct.name}"):
                os.mkdir(f"sp_{struct.name}")
                
            else:
                logger.debug(f"Directory for sp: {struct.name} exists")

            logger.debug(f"Changing directory from {os.getcwd()} to sp_{struct.name}")
            os.chdir(f"sp_{struct.name}")

            #TODO we should really make this adaptable for orca, lets say...
            if os.path.isfile("energy"):
                with open("energy", 'r') as energyFile:
                    for i,l in enumerate(energyFile):
                        pass

                if i >= 2:
                    logger.info(f"Singlepoint completed for {struct.name}")
                    self.qm_sp_energies.append(turbopy.calculation.get_energy(cycle=1))
                    logger.info(f"[{struct.name}] ==>>\t{self.qm_sp_energies[-1]:.5f} \t\t\t{self.dmd_sp_energies[index]:.5f}")
                    os.chdir("..")
                    continue
                
            elif not os.path.isfile("coord"):
                logger.debug("Creating coord file")

                pdb_to_coord.protein_to_coord(struct, self.parameters["QM Chop"])

            qm_params = self.parameters["qm params"].copy()
            qm_params["calculation"] = 'sp'
            start = timer()
            sp = turbopy.calculation(self.cores, parameters=qm_params)
            end = timer()

            if os.path.isfile("energy"):
                with open("energy", 'r') as energyFile:
                    for i,l in enumerate(energyFile):
                        pass
                
                if i < 2:
                    logger.info(f"Singlepoint failed for {struct.name}")
                    logger.info("Trying again")
                    sp = turbopy.calculation(self.cores, parameters=qm_params)
                    end = timer()

            else:
                logger.error("NO ENERGY FILE FOUND")
                logger.error("THERE COULD BE AN ERROR WITH TURBOMOLE")
                raise OSError("Turbomole")

            try:
                self.qm_sp_energies.append(turbopy.calculation.get_energy(cycle=1))
            
            except IndexError:
                logger.error("Singlepoint could not converge in {qm_params['scf']['iter']*2} scf cycles")
                self.qm_sp_energies.append(0.0)

            logger.info(f"[{struct.name}] ==>>\t{self.qm_sp_energies[-1]:.5f} \t\t\t{self.dmd_sp_energies[index]:.5f}")
            logger.info(f"Time elapsed during QM SP calculation: {datetime.timedelta(seconds = int(end -start))}")
            logger.debug(f"Changing directory from {os.getcwd()} to {os.path.abspath('..')}")
            os.chdir("../")

            if self.stop:
                logger.debug("Stop received while doing singlepoint calculations")
                return

        #We are done with all sp calculations now
        min_dmd = min(self.dmd_sp_energies)
        min_qm = min(self.qm_sp_energies)

        #We are keeping the order the same here...so everything should be good
        adjusted_dmd = [d-min_dmd for d in self.dmd_sp_energies]
        adjusted_qm = [(q-min_qm) * H_TO_KCAL for q in self.qm_sp_energies]
   
        full_scores = list(zip(adjusted_dmd, adjusted_qm))
        #TODO change the 0.5 to some value that is stored in .phd
        adjusted_scores = [s[0]*0.5 + s[1]*0.5 for s in full_scores]
        winning_score = min(adjusted_scores)
        winning_index = adjusted_scores.index(winning_score)

        self.scoring_energies = list(zip(self.dmd_sp_energies, self.qm_sp_energies, adjusted_scores))

        logger.info("")
        logger.info(f"[movie_####]------[DMD Inst. Pot (kcal)]------[QM Energy (Hart)]------[Score]------")
        for e, struct in zip(self.scoring_energies, self.sp_PDB_structures):
            out_string = f"[{struct.name}] ==>>\t{e[0]:.5f}      {e[1]:.5f}      {e[2]:.5f}"
            if e[2] == winning_score:
                out_string += "      --w---w--"

            logger.info(out_string)

        logger.info("")
        logger.info(f" >>>> [WINNER] >>>> {self.sp_PDB_structures[winning_index].name}")

        self.pdb_winner = self.sp_PDB_structures[winning_index]
        self.pdb_winner.write_pdb("dmdWinner.pdb")

        logger.info("")
        logger.info("===================[Finished  QM SP]===================")
        self.next_step = self.qm_optimization

    def qm_optimization(self):
        if self.pdb_winner is None:
            logger.error("No PDB structure for qm optimization!")
            raise ValueError("No qm optimization pdb structure")

        finished = False
        qm_params = self.parameters["qm params"].copy()
        qm_params["calculation"] = 'geo'

        total_cycles = 0

        if os.path.isdir("Optimization"):
            logger.debug("Optimization directory exists")

        else:
            os.mkdir("Optimization")
        
        logger.debug("Changing directory to 'Optimization'")
        os.chdir("Optimization")

        if os.path.isfile("energy"):
            with open('energy', 'r') as energyFile:
                for i,l in enumerate(energyFile):
                    pass
                
                i -= 1
                total_cycles = i

            if total_cycles >= qm_params["geo_iterations"]:
                logger.info(f"Previous run completed max number of cycles: {total_cycles}")
                finished = True

            else:
                logger.info(f"Restarting geometry optimization from cycle: {total_cycles}")
                qm_params['geo_iterations'] -= total_cycles

        if not finished:
            if not os.path.isfile("coord"):
                logger.debug("Creating coord file")
                pdb_to_coord.protein_to_coord(self.pdb_winner, self.parameters["QM Chop"])
  
            logger.info("      Attempting to converge active site geometry")
            start = timer()
            geo = turbopy.calculation(self.cores, parameters=qm_params)
            end = timer()
            if os.path.isfile("GEO_OPT_FAILED"):
                logger.info(f"Failed to converge in {qm_params['geo_iterations']} geometry cycles ({qm_params['geo_iterations'] + total_cycles} total cycles)")

            elif os.path.isfile("GEO_OPT_CONVERGED"):
                logger.info("Active site optimization converged!")

            else:
                logger.error("TURBOMOLE OPTIMIZATION FAILED TO RUN PROPERLY!")
                logger.error("THERE COULD BE AN ERROR WITH TURBOMOLE")
                raise OSError("Turbomole")

            self.qm_final_energy = turbopy.calculation.get_energy()
            logger.info("")
            logger.info(" >>>> [Final QM Energy] >>>> {self.qm_final_energy}")
            logger.info(f"Time elapsed during QM Opt calculation: {datetime.timedelta(seconds = int(end -start))}")

        if self.stop:
            logger.debug("Stop received while doing optimization calculation")
        
        #Now we reinstall the coords into the protein!
        self.to_next_iteration = pdb_to_coord.coord_to_pdb(self.pdb_winner)

        #Now we reinstall into the protein
        logger.debug(f"Changing directory to {os.path.abspath('..')}")
        os.chdir("..")

        self.to_next_iteration.write_pdb("to_next_iteration.pdb")
        self.next_step = self.finish_iteration

    def finish_iteration(self):
        if self.to_next_iteration is None:
            logger.error("No final PDB structure!")
            raise ValueError("No final PDB structure")

        #Do some finishing stuff here...
        #clean up this directory
        #write out the data (energies)
        self.next_step = None

    def load_from_directory(self):
        pass


