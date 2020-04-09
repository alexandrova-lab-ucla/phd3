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

from phd3.dmd_simulation import dmd_simulation
import phd3.qm_calculation as qm_calculation
import phd3.utility.utilities as utilities

import phd3.pdb_to_coord as pdb_to_coord

logger = logging.getLogger(__name__)

__all__ = [
        'iteration'
        ]

class iteration:

    def __init__(self, directory:str, root_dir:str, parameters:dict, iteration:int, cores=1):
        self.directory = os.path.abspath(directory)
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

        self.dmd_winner_energy = 0.0
        # Assign this to the correct next function step
        self.next_step = self.performDMD

        if not os.path.isdir(root_dir):
            logger.error("Root directory for iteration does not exist")
            raise ValueError("root directory does not exist")
        
        logger.info("")
        logger.info("+---------------------------------------------------------+")
        logger.info(f"|                     ITERATION      {self.iter_number:0>2d}                   |")
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
            
    def continue_calculation(self):
        logger.debug(f"Changing to iteration directory: {self.directory}")
        os.chdir(self.directory)
   
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

        #Now we go back to original directory
        logger.debug(f"Changing to iteration directory: {self.root_directory}")
        os.chdir(self.root_directory)

    def performDMD(self):
        logger.info("====================[Beginning DMD]====================")
        logger.info("")

        #Setup the dmd logger!!!!
        dmd_loggers = [logging.getLogger("phd3.setupjob"), logging.getLogger("phd3.dmd_simulation")]
        dmd_out = logging.FileHandler("dmdpy.out", 'a')
        
        for l in dmd_loggers:
            l.addHandler(dmd_out)
            l.setLevel(logging.INFO)
            l.propagate = False
        
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
                    ave, stdev = dmd_simulation.get_average_energy(os.path.join(d, self.parameters["dmd params"]["Echo File"]))
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
                    ave, stdev = dmd_simulation.get_average_energy(self.parameters["dmd params"]["Echo File"])
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

            #This is where we add any annelaing or equilibration stuff here
            if self.parameters["Equilibrate"]["Equilibrate On"] and not dmd_average_energies:
                if os.path.isdir("equilibrate"):
                    logger.info("Already equilibrated")

                else:
                    logger.info("    Equilibrating structure")
                    logger.info(f"        Starting Temp ==>> {self.parameters['Equilibrate']['Initial Temperature']}")
                    logger.info(f"        Ending Temp   ==>> {self.parameters['dmd params']['Initial Temperature']}")
                    logger.info(f"        Time          ==>> {self.parameters['Equilibrate']['Time']}")
                    logger.info("")
                    os.mkdir("equilibrate")
                    logger.debug("Changing directory to 'equilibrate'")
                    shutil.copy("dmdStart.pdb", "equilibrate")
                    os.chdir("equilibrate")
                    temp_params = self.parameters["dmd params"].copy()

                    #update params
                    temp_params["Initial Temperature"] = self.parameters['Equilibrate']['Initial Temperature']
                    temp_params["Final Temperature"] = self.parameters['dmd params']['Initial Temperature']
                    temp_params["Time"] = self.parameters['Equilibrate']['Time']

                    logger.debug("Begining dmd")
                    start = timer()
                    calc = dmd_simulation(cores = self.cores, parameters = temp_params)
                    end = timer()
                    logger.info(f"Time elapsed during equilibration: {datetime.timedelta(seconds = int(end -start))}")
                    logger.info("")

                    logger.debug("Changing directory to {os.getcwd()}")
                    os.chdir("..")


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
                        if os.path.isdir("equilibrate") and os.path.isdir(f"equilibrate/{self.parameters['dmd params']['Restart File']}"):
                            logger.debug("Copying equilibrate files to new dmd start")
                            shutil.copytree("equilibrate", f"dmdstep_{self.dmd_steps}")
                            logger.debug("Removing echo and movide file from new directory")
                            os.remove(os.path.join(f"dmdstep_{self.dmd_steps}", self.parameters["dmd params"]["Echo File"]))
                            os.remove(os.path.join(f"dmdstep_{self.dmd_steps}", self.parameters["dmd params"]["Movie File"]))

                        else:
                            os.mkdir(f"dmdstep_{self.dmd_steps}")
                            shutil.copy("dmdStart.pdb", f"dmdstep_{self.dmd_steps}/")
                        
                        logger.info("      Attempting to converge DMD simulations")
                        logger.info("[RUN #]---------[Ave Pot. Energy]-------------")

                    logger.debug(f"Changing direcetory from {os.getcwd()} to dmdstep_{self.dmd_steps}")
                    os.chdir(f"dmdstep_{self.dmd_steps}")

                    logger.debug("Begining dmd")
                    start = timer()
                    calc = dmd_simulation(cores = self.cores, parameters = self.parameters["dmd params"])
                    
                    #Harvest the time that the dmd run took
                    end = timer()

                    curr_energy = dmd_simulation.get_average_energy(self.parameters["dmd params"]["Echo File"])
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
                            for f in [d for d in os.listdir() if os.path.isfile(d)]:
                                shutil.copy(f, "../")
                 
                            os.chdir("..")
                            break

                        else:
                            logger.debug("Not converged yet")

                    os.chdir("..")
                    self.dmd_steps += 1

            else:
                if os.path.isdir("equilibrate"):
                    logger.debug("Copying equilibrate files to new dmd start")
                    for f in [d for d in os.listdir('equilibrate') if os.path.isfile(d)]:
                        shutul.copy(os.path.join("equilibrate", f), "./")

                    logger.debug("Removing echo and movie file from directory")
                    os.remove(self.parameters["dmd params"]["Echo File"])
                    os.remove(self.parameters["dmd params"]["Movie File"])
                
                logger.info("[RUN #]-------[Ave Pot. Energy]-------------")
                start = timer()
                calc = dmd_simulation(cores = self.cores, parameters = self.parameters["dmd params"])
                end = timer()
                self.final_dmd_average_energy = dmd_simulation.get_average_energy()
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
                for f in [d for d in os.listdir() if os.path.isfile(d)]:
                    shutil.copy(f, "../")

                os.chdir("../")

            elif self.parameters["DMD CONVERGE"]:
                logger.info("")
                logger.info(f"DMD converged in: {self.dmd_steps} steps")

            #convert movie to a movie.pdb file
            logger.debug("Making movie.pdb")
            utilities.make_movie("initial.pdb", self.parameters["dmd params"]["Movie File"], "movie.pdb")

            logger.debug("Loading in movie.pdb")
            self.dmd_structures = utilities.load_movie("movie.pdb")

            os.remove("movie.pdb")
            os.chdir("../")

            self.next_step = self.cluster

        for l in dmd_loggers:
            del l

        del dmd_loggers

        logger.info("")
        logger.info("====================[Finished  DMD]====================")

    def cluster(self):
        logger.info("================[Beginning  Clustering]================")
        logger.info("")
        if self.dmd_structures is None or len(self.dmd_structures) == 0:
            logger.error("No DMD structures generated!")
            raise ValueError("No DMD structures")

        #TODO check to see if there are any movie_####.pdb files...otherwise we do not need to recluster

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

            #Redirect the output to turbopy.out so that it doesn't clutter the phd output
            tm_loggers = [logging.getLogger("phd3.setupjob"), logging.getLogger("phd3.qm_calculation")]
            tm_out = logging.FileHandler("turbopy.out", 'w+')
            for l in tm_loggers:
                l.addHandler(tm_out)
                l.setLevel(logging.INFO)
                l.propagate = False

            #TODO we should really make this adaptable for orca, lets say...
            if os.path.isfile("energy"):
                with open("energy", 'r') as energyFile:
                    for i,l in enumerate(energyFile):
                        pass

                if i >= 2:
                    self.qm_sp_energies.append(qm_calculation.TMcalculation.get_energy(cycle=1))
                    logger.info(f"[{struct.name}] ==>>\t{self.qm_sp_energies[-1]:.5f} \t\t\t{self.dmd_sp_energies[index]:.5f}          (finished)")
                    os.chdir("../")
                    continue
                
            elif not os.path.isfile("coord"):
                logger.debug("Creating coord file")
                sys.exit(0)
                pdb_to_coord.protein_to_coord(struct, self.parameters["QM Chop"])

            qm_params = self.parameters["qm params"].copy()
            qm_params["calculation"] = 'sp'
            start = timer()
            sp = qm_calculation.TMcalculation(self.cores, parameters=qm_params)
            end = timer()

            if os.path.isfile("energy"):
                with open("energy", 'r') as energyFile:
                    for i,l in enumerate(energyFile):
                        pass
                
                if i < 2:
                    logger.info(f"Singlepoint failed for {struct.name}")
                    logger.info("Trying again")
                    sp = qm_calculation.TMcalculation(self.cores, parameters=qm_params)
                    end = timer()

            else:
                logger.error("NO ENERGY FILE FOUND")
                logger.error("THERE COULD BE AN ERROR WITH TURBOMOLE")
                raise OSError("Turbomole")

            try:
                self.qm_sp_energies.append(qm_calcualtion.TMcalculation.get_energy(cycle=1))
            
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

            #clean up the loggers
            for l in tm_loggers:
                del l

            del tm_loggers

        logger.info("")
        logger.info("                        >>>> Scoring Structures <<<<")

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
        logger.info(f"[movie_####]------[DMD Inst. Pot (kcal)]------[QM Energy (Hart)]-------------[Score]------")
        for e, struct in zip(self.scoring_energies, self.sp_PDB_structures):
            out_string = f"[{struct.name}] ==>>\t{e[0]:.5f}                {e[1]:.5f}               {e[2]:.5f}"
            if e[2] == winning_score:
                out_string += "      --w---w--"
                self.dmd_winner_energy = e[0]

            logger.info(out_string)

        logger.info("")
        logger.info(f" >>>> [WINNER] >>>> {self.sp_PDB_structures[winning_index].name}")

        self.pdb_winner = self.sp_PDB_structures[winning_index]
        self.pdb_winner.write_pdb("dmdWinner.pdb")

        logger.info("")
        logger.info("===================[Finished  QM SP]===================")
        self.next_step = self.qm_optimization

    def qm_optimization(self):
        logger.info("===================[Beginning QM OP]===================")
        logger.info("")

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
        
        #Setup the tm loggers
        tm_loggers = [logging.getLogger("phd3.setupjob"), logging.getLogger("phd3.qm_calculation")]
        tm_out = logging.FileHandler("turbopy.out", 'w+')
        for l in tm_logger:
            l.addHandler(tm_out)
            l.setLevel(logging.INFO)
            l.propagate = False

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
            geo = qm_calculation.TMcalculation(self.cores, parameters=qm_params)
            end = timer()
            if os.path.isfile("GEO_OPT_FAILED"):
                logger.info(f"Failed to converge in {qm_params['geo_iterations']} geometry cycles ({qm_params['geo_iterations'] + total_cycles} total cycles)")

            elif os.path.isfile("GEO_OPT_CONVERGED"):
                logger.info("Active site optimization converged!")

            else:
                logger.error("TURBOMOLE OPTIMIZATION FAILED TO RUN PROPERLY!")
                logger.error("THERE COULD BE AN ERROR WITH TURBOMOLE")
                raise OSError("Turbomole")

            self.qm_final_energy = qm_calculation.TMcalculation.get_energy()
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
      
        #Cleanup the loggers
        for l in tm_loggers:
            del l

        del tm_loggers

        logger.info("")
        logger.info("===================[Finished  QM OP]===================")


    def finish_iteration(self):
        logger.debug("Saving energies to phd_energy")
        if self.to_next_iteration is None:
            logger.error("No final PDB structure!")
            raise ValueError("No final PDB structure")

        if os.path.isfile(os.path.join(self.root_directory, "phd_energy")):
            #Just append
            outfile_lines = []
            with open(os.path.join(self.root_directory, "phd_energy"), 'r') as outfile:
                for line in outfile:
                    outfile_lines.append(line)

            first_line = outfile_lines.pop(0)
            prev_iter = [int(i.split()[0]) for i in outfile_lines][-1]

            while prev_iter + 1 != self.iter_number and outfile_lines:
                outfile_lines.pop()
                prev_iter = [int(i.split()[0]) for i in outfile_lines][-1]
            
            outfile_lines.append(f"{self.iter_number:0>2d}            {self.qm_final_energy:.5f}            {self.dmd_winner_energy:.5f}\n")

            with open(os.path.join(self.root_directory, "phd_energy"), 'w') as outfile:
                outfile.write(first_line)
                for line in outfile_lines:
                    outfile.write(line)

        else:
            with open(os.path.join(self.root_directory, "phd_energy")) as outfile:
                outfile.write("[Iteration]        [QM Energy (Hart)]        [DMD Energy (kcal)]\n")
                outfile.write(f"{self.iter_number:0>2d}            {self.qm_final_energy:.5f}            {self.dmd_winner_energy:.5f}\n")
        

        #TODO clean up the iteration directory

        self.next_step = None


    def signal_alarm(self):
        #This is called if the controller ran out of time
        self._stop = True

        #Create a stop file for the geometry optimization...otherwise, just have to wait for it to stop itself
        if self.next_step == self.qm_optimization:
            with open(os.path.join(self.directory, "Optimization/stop"), 'w') as stopfile:
                pass




