#!/usr/bin/env python3
"""
Author  ==>> Matthew R. Hennefarth
Date    ==>> April 16, 2020
"""

#Standard Library Imports
import logging
import os
import shutil
import numpy as np
import sys
import time
from timeit import default_timer as timer
import datetime
import itertools
import json
from multiprocessing import sharedctypes, Process

#3rd party libraries
import hdbscan
import scipy.cluster

#phd libraries
from phd3.dmd_simulation import dmd_simulation
import phd3.qm_calculation as qm_calculation
import phd3.utility.utilities as utilities
from phd3.setupjob import setupTMjob
import phd3.utility.constants as constants
import phd3.dmd_to_qm as dmd_to_qm
import phd3.protein as protein
import phd3.utility.exceptions as exceptions
from .bin import submitturbomole

logger = logging.getLogger(__name__)

__all__ = [
        'iteration'
        ]

class iteration:

    """Upon creati"""
    def __init__(self, controller, directory:str, root_dir:str, parameters:dict, iteration:int, cores=1, scratch="./"):
        
        #Iteration directory, iteration number, root directory of job, number of cores to use, and parameters saved here
        self.directory = os.path.abspath(directory)
        self.iter_number = iteration
        self.root_directory = root_dir
        self.parameters = parameters
        self.cores = cores
        self.scratch = scratch
        self.controller = controller
        
        #Set to True if a timer went off, or we need to stop the iteration
        self.stop = False
        
        #Final (ave pot. energy, stdev) from the DMD simulations
        self.final_dmd_average_energy = [0, 0]
        
        #Winning protein structure after scoring
        self.pdb_winner = None

        #QM sp energies,. same order as self.dmd_structures
        self.qm_sp_energies = []
        self.scoring_energies = []
        self.qm_final_energy = 0.0

        #Final, updated protein after QM optimization
        self.to_next_iteration = None
        
        #Stores the dmd structures from the final movie file
        #Then cleared at end of clustering to free up memory
        self.dmd_structures = None
        #Holds the pdb structures that will move on SP analysis
        self.sp_PDB_structures = None

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
            
            self.stop = self.timer_went_off()

        if self.stop:
            logger.debug("Timer went off, we are ending this iteration early")
        
        if self.next_step is None:
            logger.info("")
            logger.info("+---------------------------------------------------------+")
            logger.info(f"|                 ITERATION      {self.iter_number:0>2d} COMPLETE              |")
            logger.info("+---------------------------------------------------------+")
            logger.info("")

        #Now we go back to original directory
        logger.debug(f"Changing to iteration directory: {self.root_directory}")
        os.chdir(self.root_directory)

    def timer_went_off(self):
        if self.controller.time_left() == -1:
            return False

        elif self.controller.time_left() < 0.75:
            return True
        
        return self.stop

    def performDMD(self):
        logger.info("====================[Beginning DMD]====================")
        logger.info("")

        if self.iter_number == -1:
            logger.info("On iteration 0, we skip the DMD portion")
            logger.info("and go straight to QM Optimization!")
            self.pdb_winner = [utilities.load_pdb(self.parameters["last pdb"]), 0.0]
            self.pdb_winner[0].name = "movie_0000"
            self.next_step = self.qm_optimization
            logger.info("")
            logger.info("====================[Finished  DMD]====================")
            return

        #This tracks to see if we are done with the dmd simulations!
        finished = False
        converged = False
        dmd_average_energies = []
        self.dmd_steps = 1
        if os.path.isdir("dmd"):
            # Want to continue (check to see if we have first converged, if that is an option as well)
            logger.info("dmd directory found")
            logger.info(">>>> Previous DMD simulations >>>>")

            logger.debug("Moving to dmd directory")
            os.chdir("dmd")
            if self.parameters["DMD CONVERGE"]:
                dirs = [d for d in os.listdir() if "dmdstep_" in d]

                #Make sure that the directories are sorted
                dirs.sort(key = lambda i: int(i.split("_")[-1]))
                rm_dirs = False
                for i, d in enumerate(dirs):
                    if rm_dirs:
                        shutil.rmtree(d)
                        continue

                    elif not os.path.isfile(os.path.join(d, self.parameters['dmd params']['Echo File'])):
                        logger.info(f"Did not complete run {i+1}, removing directory and continuing")
                        shutil.rmtree(d)
                        #Want to remove any directories (if they exist) after this one! and start from here
                        rm_dirs = True
                        continue

                    echo_lines = dmd_simulation.get_echo_data(os.path.join(d, self.parameters['dmd params']['Echo File']))

                    last_time = int(float(echo_lines[-1][0]))
                    if last_time < self.parameters['dmd params']['Time']:
                        logger.info(f"Did not complete run {i+1}, removing directory and continuing")
                        rm_dirs = True
                        shutil.rmtree(d)

                    else:
                        ave, stdev = dmd_simulation.get_average_potential_energy(os.path.join(d, self.parameters["dmd params"]["Echo File"]))
                        dmd_average_energies.append([ave, stdev])
                        logger.info(f"[Run {i+1:0>2d}] ==>> {ave:0<12.5f} ({stdev:0<9.5f})(finished)") 
            
                if len(dmd_average_energies) >= 2:
                    #Check to see if we have previously converged
                    if self.dmd_converged(dmd_average_energies):
                        logger.info("Convergence previously achieved")
                        logger.info(f"Absolute delta in energy: {abs(dmd_average_energies[-1][0]-dmd_average_energies[-2][0]):.5f}")
                        converged = True
                        finished = True
                        self.final_dmd_average_energy = dmd_average_energies[-1] if dmd_average_energies[-1][0] < dmd_average_energies[-2][0] else dmd_average_energies[-2]
               
                #Check to see how many directories we have actually gone through now that we have removed them...
                dirs = [d for d in os.listdir() if "dmdstep_" in d]
                self.dmd_steps = len(dirs)+1

            else:
                #We are not converging....check the echo file
                if os.path.isfile(self.parameters["dmd params"]["Echo File"]) and os.path.isfile(self.parameters["dmd params"]["Movie File"]):
                    ave, stdev = dmd_simulation.get_average_potential_energy(self.parameters["dmd params"]["Echo File"])
                    logger.info(f"[Run 01] ==>> {ave:0<12.5f} ({stdev:0<9.5f}) (finished)")
                    finished = True
  
        else:
            logger.debug("Making dmd directory")
            os.mkdir("dmd")
            logger.debug("Moving to dmd directory")
            os.chdir("dmd")
            #This is set by the controller.py class, no need to worry too much about this...
            logger.debug("Moving the last step pdb to the new dmd directory")
            shutil.copy(self.parameters["last pdb"], "./dmdStart.pdb")
      
        if not finished:

            #This is where we add any annelaing or equilibration stuff here
            if self.parameters["Equilibrate"]["Equilibrate On"] and not dmd_average_energies:
                if os.path.isdir("equilibrate") and os.path.isfile(f"equilibrate/{self.parameters['dmd params']['Echo File']}"):
                    logger.info("Already Equilibrated")

                else:
                    logger.info(">>>> Equilibrating structure >>>>")
                    logger.info(f"[Starting Temp.]   ==>> {self.parameters['Equilibrate']['Initial Temperature']}")
                    logger.info(f"[Ending Temp.]     ==>> {self.parameters['dmd params']['Initial Temperature']}")
                    logger.info(f"[Time]             ==>> {self.parameters['Equilibrate']['Time']}")
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

                    calc = dmd_simulation(cores = self.cores, parameters = temp_params, run_dir=self.scratch)
                    self.stop = self.timer_went_off()
                    logger.debug("Changing directory to {os.getcwd()}")
                    os.chdir("..")

            #If we are converging DMD, then we enter here into a nice loop
            if self.parameters["DMD CONVERGE"]:
                logger.info("")
                logger.info("Attempting to converge DMD simulations")
                
                while not self.stop and self.dmd_steps <= self.parameters["MAX DMD STEPS"]:
                    #TODO try and except the following
                    #This allows us to use the previous step restart file if it exists
                    if os.path.isdir(f"dmdstep_{self.dmd_steps-1}"):
                        logger.debug(f"Copying {self.dmd_steps-1} files to new dmd start")

                        shutil.copytree(f"dmdstep_{self.dmd_steps-1}", f"dmdstep_{self.dmd_steps}")
                        logger.debug("Removing echo and movie file from new directory")
                        os.remove(os.path.join(f"dmdstep_{self.dmd_steps}", self.parameters["dmd params"]["Echo File"]))
                        os.remove(os.path.join(f"dmdstep_{self.dmd_steps}", self.parameters["dmd params"]["Movie File"]))

                    elif os.path.isdir("equilibrate") and os.path.isfile(f"equilibrate/{self.parameters['dmd params']['Restart File']}"):
                        logger.debug("Copying equilibrate files to new dmd start")
                        shutil.copytree("equilibrate", f"dmdstep_{self.dmd_steps}")
                        logger.debug("Removing echo and movide file from new directory")
                        os.remove(os.path.join(f"dmdstep_{self.dmd_steps}", self.parameters["dmd params"]["Echo File"]))
                        os.remove(os.path.join(f"dmdstep_{self.dmd_steps}", self.parameters["dmd params"]["Movie File"]))

                    else:
                        os.mkdir(f"dmdstep_{self.dmd_steps}")
                        shutil.copy("dmdStart.pdb", f"dmdstep_{self.dmd_steps}/")
                        

                    logger.debug(f"Changing direcetory from {os.getcwd()} to dmdstep_{self.dmd_steps}")
                    os.chdir(f"dmdstep_{self.dmd_steps}")

                    logger.info("")
                    logger.info(f">>>> Run {self.dmd_steps:0>2d} >>>>")
                    calc = dmd_simulation(cores = self.cores, parameters = self.parameters["dmd params"], run_dir=self.scratch)
                    self.stop = self.timer_went_off()
                    curr_energy = dmd_simulation.get_average_potential_energy(self.parameters["dmd params"]["Echo File"])
                    dmd_average_energies.append(curr_energy)

                    #if converged then break
                    if len(dmd_average_energies) >= 2:
                        logger.info(f"[Delta Last Sim. ] ==>> {abs(curr_energy[0] - dmd_average_energies[-2][0]):.5f}")
                        if self.dmd_converged(dmd_average_energies):
                            logger.info("")
                            logger.info("Converence achieved!")
                            converged = True
                            logger.info(f"Absolute delta in energy: {abs(curr_energy[0]-dmd_average_energies[-2][0]):.5f}")
                            if dmd_average_energies[-2][0] < curr_energy[0]:
                                logger.debug("Using previous energy")
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

            #Not converging the dmd simulations
            else:
                if os.path.isdir("equilibrate"):
                    logger.debug("Copying equilibrate files to new dmd start")
                    for f in [d for d in os.listdir('equilibrate') if os.path.isfile(d)]:
                        shutul.copy(os.path.join("equilibrate", f), "./")

                    logger.debug("Removing echo and movie file from directory")
                    os.remove(self.parameters["dmd params"]["Echo File"])
                    os.remove(self.parameters["dmd params"]["Movie File"])
                
                logger.info("")
                logger.info(">>>> Run 01 >>>>")
                calc = dmd_simulation(cores = self.cores, parameters = self.parameters["dmd params"], run_dir=self.scratch)
                self.final_dmd_average_energy = dmd_simulation.get_average_potential_energy()
                self.stop = self.timer_went_off()


        if self.stop:
            logger.info("Stop signal received while performing DMD simulation")
        
        else:
            if converged:
                logger.info("")
                logger.info(f"DMD converged in: {self.dmd_steps} steps")

            elif self.dmd_steps > self.parameters["MAX DMD STEPS"]:
                logger.info("")
                logger.info("Max DMD steps taken")
                logger.info("Finding lowest Ave. Pot. Energy DMD simulation")
                min_step = min(dmd_average_energies, key = lambda step: step[0])
                self.final_dmd_average_energy = min_step
                min_step_index = dmd_average_energies.index(min_step) + 1
                logger.info(f"[Lowest DMD simulation] ==>> {min_step_index}")
                logger.info(f"[Ave. Potential Energy] ==>> {min_step[0]:.5f} ({min_step[1]:.5f})")
                logger.debug(f"Moving to dmdstep_{min_step_index}")
                os.chdir(f"dmdstep_{min_step_index}")
                
                logger.debug("Copying files up")
                for f in [d for d in os.listdir() if os.path.isfile(d)]:
                    shutil.copy(f, "../")

                os.chdir("../")

            #If we are supposed to converge, but for some reason left the cycle previously...major error most likely
            #if this actually happens
            elif self.parameters["DMD CONVERGE"]:
                logger.warn("Unknown reason why we stopped the DMD cycle")

            #Print out the summary of all the dmd cycles...could add some additional information here, but
            #not really necessary
            logger.info("")
            logger.info(">>>> DMD Summary >>>>")
            logger.info("[RUN ##]---------[Ave Pot. Energy]---------[Est. Phys Time (ns)]---------")
            for index, energy in enumerate(dmd_average_energies):
                out = f"[Run {index+1:0>2d}]      {energy[0]:0<12.5f} ({energy[1]:.5f})\t\t{self.parameters['dmd params']['Time']*0.0000488882}"
                if energy == self.final_dmd_average_energy:
                    out += "    --w---w--"
                
                logger.info(out)

            logger.info("")
            logger.info("Loading in trajectory")
            #convert movie to a movie.pdb file
            logger.debug("Making movie.pdb")
            logger.info("...")
            utilities.make_movie("initial.pdb", self.parameters["dmd params"]["Movie File"], "movie.pdb")
            logger.debug("Loading in movie.pdb")
            #returns list of protein structures with movie_#### as the name of the protein
            self.dmd_structures = utilities.load_movie("movie.pdb")

            dmd_structures_energies = [float(line[4]) for line in dmd_simulation.get_echo_data(self.parameters['dmd params']['Echo File'])]
            #Now we are done
            #A list of (protein, dmd energy), this way we can always get the dmd energy...could make that a member of the protein class actually...
            self.dmd_structures = list(zip(self.dmd_structures, dmd_structures_energies))

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

        #Need these here because when we print out summary, need these vars. to be at least initialized so that it
        #doesn't throw an error
        min_structures = []
        centroid_structures = []

        if [d for d in os.listdir() if "movie_" in d]:
            logger.info("Loading in previous clustering results")
            self.sp_PDB_structures = []
            for f in [d for d in os.listdir() if d.startswith("movie_") and os.path.isfile(d)]:
                f = f.split("_")
                f = f[1].split(".")
                num = int(f[0])
                self.sp_PDB_structures.append(self.dmd_structures[num])

        else:
            # Create an uninitialized  matrix of the appropriate size, will fill in with the RMSD between structures
            logger.info("Computing RMSD between all pairs of structures")
            logger.info(f"[Cores]             ==>> {self.cores}")
            logger.info("...")
            start = timer()
           
            #If I want to speed this up, I have to look at numba or numpy optimizations in the aa_rmsd func.
            if self.cores > 1:
                #There is a pretty large speed up from using multicore chunking of the jobs
                #This works fairly well!
                indexes = [(i, j) for i, j in itertools.product(range(0, len(self.dmd_structures)), range(0, len(self.dmd_structures)))]
                index = np.array([a for a in indexes if a[0] < a[1]])
                
                #Now we have some chunks to pass to the function
                chunked_index = np.array_split(index, self.cores)
                
                #This is used for no data racing, etc..
                result = np.ctypeslib.as_ctypes(np.zeros((len(self.dmd_structures), len(self.dmd_structures) )))
                shared_array = sharedctypes.RawArray(result._type_, result)
                def computeRMSD(vals):
                    tmp = np.ctypeslib.as_array(shared_array)
                    for i, j in vals:
                        tmp[i, j] = self.dmd_structures[i][0].aa_rmsd(self.dmd_structures[j][0])
                        tmp[j, i] = tmp[i, j]

                #start all the procs
                procs = [Process(target=computeRMSD, args=(vals,)) for vals in chunked_index]
                [p.start() for p in procs]
                [p.join() for p in procs]

                result = np.ctypeslib.as_array(shared_array)
                distance_matrix = result 
                del procs                

            else:
                distance_matrix = np.zeros([len(self.dmd_structures), len(self.dmd_structures)], dtype=float)
                #We just use the slow method of this...
                for row, first_structure in enumerate(self.dmd_structures):
                    for col, second_structure in enumerate(self.dmd_structures[row+1:]):
                        distance_matrix[row][col] = first_structure[0].aa_rmsd(second_structure[0])
                        distance_matrix[col][row] = distance_matrix[row][col] 

            end = timer()
            logger.info(f"Time elapsed during RMSD calculation: {datetime.timedelta(seconds = int(end -start))}")
            
            cluster_method = "HDBSCAN"
            
            logger.info("")
            logger.info(">>>> Cluster Parameters >>>>")
            logger.info(f"[Frames Collected]  ==>> {len(self.dmd_structures)}")
            logger.info(f"[Cluster Method]    ==>> {cluster_method}")
            logger.info(f"[Min. Cluster Size] ==>> {50}")
            logger.info(f"[Min. Samples]      ==>> {25}")
            logger.info("...")

            #These clustering methods are incredibly fast in comparison to the agglomerative clustering
            #Really no need to time these as they are practically instantaneous
            start = timer()
            output = hdbscan.HDBSCAN(metric='precomputed', core_dist_n_jobs=self.cores, min_cluster_size=50, min_samples=25)
            output.fit(distance_matrix)
            end = timer()
            
            num_clusters = len([x for x in set(output.labels_) if x >= 0])
            logger.info(f"[Num. of Clusters]  ==>> {num_clusters}")
            logger.info(f"Time elapsed during HDBSCAN clustering: {datetime.timedelta(seconds = int(end -start))}")
            
            if num_clusters > self.parameters["Max Clusters"] or num_clusters <= 0:
                logger.info("")
                logger.info("Too many clusters")
                logger.info("Trying to cluster with single linkage hierachy")
                logger.info("")
                logger.info(">>>> Cluster Parameters >>>>")
                logger.info(f"[Max Num. Clusters] ==>> {self.parameters['Max Clusters']}")
                logger.info("[Cluster Method]    ==>> Single Linkage Heirarchy")
                logger.info("[Cluster Criterion] ==>> Max Cluster")
                logger.info("...") 
                start = timer()
                y = [distance_matrix[j][i] for j in range(len(distance_matrix)) for i in range(j+1, len(distance_matrix))]
                Z = scipy.cluster.hierarchy.linkage(y)
                cluster_labels = scipy.cluster.hierarchy.fcluster(Z, t=self.parameters['Max Clusters'], criterion='maxclust')
                end = timer()
                
                num_clusters = len(set(cluster_labels))
                logger.info(f"[Num. of Clusters]  ==>> {num_clusters}")
                logger.info(f"Time elapsed during Single Linkage clustering: {datetime.timedelta(seconds = int(end -start))}")


            else:
                cluster_labels = output.labels_ 

            logger.info(f"[Minimum Energy]    ==>> {'true' if self.parameters['Cluster Energy'] else 'false'}")
            logger.info(f"[Cluster Centroid]  ==>> {'true' if self.parameters['Cluster Centroid'] else 'false'}")

            #Reformat so that we can work the various clusters
            clustered_structures = {}
            for index in set(cluster_labels):
                #These are considered noise by the HDBSCAN
                #So we don't necessarily care for them
                if index == -1:
                    continue

                clustered_structures[index] = []
                for struct in zip(self.dmd_structures, cluster_labels):
                    if struct[1] == index:
                        clustered_structures[index].append(struct[0])

            #Now we can access the clusters by the keys in this dictionary
            # This is a dictionary, with values of a list of tuples (protein, dmd energy)
            old_dmd_structures = self.dmd_structures.copy()
            self.dmd_structures = clustered_structures
            
            #Get the minimum of each cluster
            self.sp_PDB_structures = []
            min_structures = []
            if self.parameters["Cluster Energy"]:
                for index in self.dmd_structures.keys():
                    min_structures.append(min(self.dmd_structures[index], key = lambda i: i[1]))
            
            self.sp_PDB_structures.extend(min_structures)

            #Get the centroid of each cluster
            centroid_structures = []
            if self.parameters["Cluster Centroid"]:
                for index in self.dmd_structures.keys():

                    #Want to find the structure that minimizes RMSD values with the rest of the cluster
                    sub_distance_matrix = np.zeros((len(self.dmd_structures[index]), len(self.dmd_structures[index])))

                    #Extract the sub distance matrix for this cluster only
                    #Too slow to reconstruct the distance matrix from scratch
                    for row, first_structure in enumerate(self.dmd_structures[index]):
                        first_structure_index = old_dmd_structures.index(first_structure)
                        for col, second_structure in enumerate(self.dmd_structures[index][row+1:]):
                            second_structure_index = old_dmd_structures.index(second_structure)
                            sub_distance_matrix[row][col] = distance_matrix[first_structure_index][second_structure_index]
                            sub_distance_matrix[col][row] = distance_matrix[row][col]

                    #Find the one with the lowest average row sum...
                    rows_sums = [i.sum() for i in distance_matrix]
                    with_rows = list(zip(self.dmd_structures[index], rows_sums))
                    winner = min(with_rows, key = lambda i: i[1])
                    centroid_structures.append(winner[0])

            self.sp_PDB_structures.extend(centroid_structures)
            
            if not self.sp_PDB_structures:
                logger.info("No clustering method specified, using minimum energy from each cluster")
                for index in self.dmd_structures.keys():
                    min_structures.append(min(self.dmd_structures[index], key = lambda i: i[1]))
       
                self.sp_PDB_structures.extend(min_structures)

            #As we leave the final part, free up the memory...no longer needed
            del distance_matrix
            self.dmd_structures.clear()
            old_dmd_structures.clear()

        #Print out choices, and save the files
        logger.info("")
        logger.info(">>>> Final Structures >>>>")
        logger.info("[structure]-------[DMD Energy (kcal)]------------")
        #Saving the structures in case we have to restart
        for struct in self.sp_PDB_structures:
            out = f"[{struct[0].name}]            {struct[1]:.5f}"
            if struct in min_structures:
                out += "    (min energy)"

            if struct in centroid_structures:
                out += "    (centroid)"

            logger.info(out)
            struct[0].write_pdb(f"{struct[0].name}.pdb")

        self.next_step = self.qm_singlepoints
        logger.info("")
        logger.info("=================[Finished Clustering]=================")

    def qm_singlepoints(self):
        logger.info("===================[Beginning QM SP]===================")

        if self.sp_PDB_structures is None or len(self.sp_PDB_structures) == 0:
            logger.error("No PDB structures for single point calculations")
            raise ValueError("No single point structures")
    
        #Make all of the directories
        for index, [struct, dmd_energy] in enumerate(self.sp_PDB_structures):
            logger.info("")
            logger.info(f">>>> {struct.name} >>>>")
            
            if not os.path.isdir(f"sp_{struct.name}"):
                os.mkdir(f"sp_{struct.name}")
                
            else:
                logger.debug(f"Directory for sp: {struct.name} exists")

            logger.debug(f"Changing directory from {os.getcwd()} to sp_{struct.name}")
            os.chdir(f"sp_{struct.name}")

            #TODO we should really make this adaptable for orca, lets say...
            #Make a base setupQM object, have orca and TM inherit (overlap on interface funcs)
            #Same for the executation (running the job)
            if os.path.isfile("energy"):
                i = 0
                with open("energy", 'r') as energyFile:
                    for i,l in enumerate(energyFile):
                        pass
                
                #recall that i goes 0 -> 1 -> 2, so if i = 2, there are 3 lines => an energy!
                if i >= 2:
                    self.qm_sp_energies.append(qm_calculation.TMcalculation.get_energy(cycle=1))
                    logger.info(f"[QM Energy]       ==>> {self.qm_sp_energies[-1]:.5f} Hart (finished)")
                    os.chdir("../")
                    continue
                
            elif not os.path.isfile("coord"):
                logger.debug("Creating coord file")
                dmd_to_qm.protein_to_coord(struct, self.parameters["QM Chop"])

            if not os.path.isfile("control"):
                logger.debug("Setting up TM job")

                #we set timeout to 0 so that we don't have a weird define error from timout in the middle of a job
                try:
                    sj = setupTMjob(parameters = self.parameters["qm params"], timeout=0)

                except exceptions.DefineError:
                    #Sometimes define ends abnormally and just deleting the control file and redoing works
                    if os.path.isfile("control"):
                        os.remove("control")

                    sj = setupTMjob(parameters = self.parameters['qm params'], timeout=0)
                
                dirs = [d for d in os.listdir("../") if os.path.isdir(os.path.join("../", d)) and "sp_movie_" in d]
                dirs.sort(reverse = True, key=lambda i: int(i.split("_")[-1]))
                
                found_mos = False 
                for sp_directory in dirs:
                    if os.path.abspath(os.path.join("..", sp_directory)) == os.path.abspath(os.getcwd()):
                        continue

                    sp_directory = os.path.join("../", sp_directory)
                    found_mos = False
                    for mo_file in constants.MO_FILES:
                        if os.path.isfile(os.path.join(sp_directory, mo_file)):
                            logger.debug(f"Copying over {mo_file} file")
                            found_mos = True
                            shutil.copy(os.path.join(sp_directory, mo_file), f"./{mo_file}")

                    if found_mos:
                        break
                
                if not found_mos:
                    #TODO allow user to specify a default mos, alpha, beta to use!
                    pass

            if os.path.isdir("trun_backup"):
                for mo_file in constants.MO_FILES:
                    if os.path.isfile(os.path.join("trun_backup", mo_file)):
                        logger.debug(f"Copying up {mo_file} file from trun_backup")
                        shutil.copy(os.path.join("trun_backup", mo_file), f"./{mo_file}")

            qm_params = self.parameters["qm params"].copy()
            qm_params["calculation"] = 'sp'
            start = timer()
            sp = qm_calculation.TMcalculation(self.cores, parameters=qm_params, run_dir=self.scratch, time=self.controller.time_left())
            end = timer()
            self.stop = self.timer_went_off()

            if os.path.isfile("energy"):
                i = 0
                with open("energy", 'r') as energyFile:
                    for i,l in enumerate(energyFile):
                        pass
                
                if i < 2 and not self.stop:
                    logger.info(f"Singlepoint failed for {struct.name}")
                    logger.info("Trying again")
                    sp = qm_calculation.TMcalculation(self.cores, parameters=qm_params, run_dir=self.scratch, time=self.controller.time_left())
                    end = timer()
                    self.stop=self.timer_went_off()

            else:
                logger.error("NO ENERGY FILE FOUND")
                logger.error("THERE COULD BE AN ERROR WITH TURBOMOLE")
                raise OSError("Turbomole")

            try:
                self.qm_sp_energies.append(qm_calculation.TMcalculation.get_energy(cycle=1))
            
            except IndexError:
                logger.error("Singlepoint could not converge in {qm_params['scf']['iter']*2} scf cycles")
                logger.error("Structure may be very weird, continuing to next structure")
                self.qm_sp_energies.append(0.0)

            logger.info(f"[QM Energy] ==>> {self.qm_sp_energies[-1]:.5f} Hart")
            logger.info(f"Time elapsed during QM SP calculation: {datetime.timedelta(seconds = int(end -start))}")
            logger.debug(f"Changing directory from {os.getcwd()} to {os.path.abspath('..')}")
            os.chdir("../")

            if self.stop:
                logger.debug("Stop received while doing singlepoint calculations")
                return

        logger.info("")
        logger.info(">>>> Scoring Structures >>>>")

        #We are done with all sp calculations now, onto scoring the structures
        min_dmd = min([d[1] for d in self.sp_PDB_structures])
        min_qm = min(self.qm_sp_energies)

        #We are keeping the order the same here...so everything should be good
        adjusted_dmd = [d-min_dmd for d in [i[1] for i in self.sp_PDB_structures]]
        adjusted_qm = [(q-min_qm) * constants.H_TO_KCAL for q in self.qm_sp_energies]
   
        full_scores = list(zip(adjusted_dmd, adjusted_qm))
        #TODO change the 0.5 to some value that is stored in .phd
        adjusted_scores = [s[0]*0.5 + s[1]*0.5 for s in full_scores]
        winning_score = min(adjusted_scores)
        winning_index = adjusted_scores.index(winning_score)

        self.scoring_energies = list(zip([d[0] for d in self.sp_PDB_structures], [d[1] for d in self.sp_PDB_structures], self.qm_sp_energies, adjusted_scores))

        logger.info(f"[movie_####]------[DMD Inst. Pot (kcal)]------[QM Energy (Hart)]-----------[Score]------")
        for struct in self.scoring_energies:
            out_string = f"[{struct[0].name}]         {struct[1]:0<12}                {struct[2]:0<14.6f}            {struct[3]:.3f}"
            if struct[3] == winning_score:
                out_string += "      --w---w--"

            logger.info(out_string)

        logger.info("")
        logger.info(f"[WINNER]     ==>> {self.sp_PDB_structures[winning_index][0].name}")
        logger.info(f"[QM  Energy] ==>> {self.qm_sp_energies[winning_index]:.5f} Hartrees")
        logger.info(f"[DMD Energy] ==>> {self.sp_PDB_structures[winning_index][1]:.5f} kcal/mol")

        self.pdb_winner = self.sp_PDB_structures[winning_index]
        self.pdb_winner[0].reformat_protein()
        self.pdb_winner[0].write_pdb("dmdWinner.pdb")

        logger.info("")
        logger.info("===================[Finished  QM SP]===================")
        self.next_step = self.qm_optimization

    def qm_optimization(self):
        logger.info("==================[Beginning  QM OPT]==================")
        logger.info("")

        if self.pdb_winner is None:
            logger.error("No PDB structure for qm optimization!")
            raise ValueError("No qm optimization pdb structure")

        finished = False
        qm_params = self.parameters["qm params"].copy()
        qm_params["calculation"] = 'geo'

        #Keeps track of total QM optimization steps taken
        total_cycles = 0

        if os.path.isdir("Optimization"):
            logger.debug("Optimization directory exists")

        else:
            os.mkdir("Optimization")
        
        logger.debug("Changing directory to Optimization")
        os.chdir("Optimization")
       
        if os.path.isdir("trun_backup"):
            logger.debug("Overriding directory with trun_backup")
            for f in [f for f in os.listdir("./trun_backup")]:
                f = os.path.join("./trun_backup", f)
                if os.path.isfile(f):
                    if "coord" in f or "energy" in f or "gradient" in f:
                        shutil.copy(f, "./")

            shutil.rmtree("trun_backup")

        if os.path.isfile("energy"):
            with open('energy', 'r') as energyFile:
                i = 0
                for i,l in enumerate(energyFile):
                    pass
                
                i -= 1
                total_cycles = i if i > 0 else 0

            if total_cycles >= qm_params["geo_iterations"]:
                logger.info(f"Previous run completed max number of cycles: {total_cycles}")
                finished = True

            else:
                logger.info(f"Restarting geometry optimization from cycle: {total_cycles}")
                qm_params['geo_iterations'] -= total_cycles

        if not finished:
            logger.info("      Attempting to converge active site geometry")
            if not os.path.isfile("coord"):
                logger.debug("Creating coord file")
                dmd_to_qm.protein_to_coord(self.pdb_winner[0], self.parameters["QM Chop"])
 
            if not os.path.isfile("control"):
                logger.debug("Setting up TM job")
                
                try:
                    sj = setupTMjob(parameters = qm_params, timeout=0)

                except exceptions.DefineError:
                    #Sometimes define throws an error even if everything worked, just retry it once to see if that fixes the issue
                    if os.path.isfile("control"):
                        os.remove("control")

                    sj = setupTMjob(parameters = qm_params, timeout=0)

                possible_sp_dir = "../sp_" + self.pdb_winner[0].name
                found_mos = False
                if os.path.isdir(possible_sp_dir):
                    for mo_file in constants.MO_FILES:
                        if os.path.isfile(os.path.join(possible_sp_dir, mo_file)):
                            logger.debug(f"Copying over {mo_file} file")
                            found_mos = True
                            shutil.copy(os.path.join(possible_sp_dir, mo_file), f"./{mo_file}")

                else:
                    #Check for any of these directories actually...
                    dirs = [d for d in os.listdir("../") if os.path.isdir(os.path.join("../", d)) and "sp_movie_" in d]
                    dirs.sort(reverse = True, key=lambda i: int(i.split("_")[-1]))
                    
                    for sp_directory in dirs:
                        sp_directory = os.path.join("../", sp_directory)
                        for mo_file in constants.MO_FILES:
                            if os.path.isfile(os.path.join(possible_sp_directory, mo_file)):
                                logger.debug(f"Copying over {mo_file} file")
                                found_mos = True
                                shutil.copy(os.path.join(sp_directory, mo_file), f"./{mo_file}")

                        if found_mos:
                            break
                
                if not found_mos:
                    #TODO allow user to specify a default mos, alpha, beta to use!
                    pass

            start = timer()
            geo = qm_calculation.TMcalculation(self.cores, parameters=qm_params, time=self.controller.time_left(), run_dir=self.scratch)
            end = timer()
            logger.info(f"Time elapsed during QM Opt calculation: {datetime.timedelta(seconds = int(end -start))}")
            #Get information on if the timer went off
            self.stop = geo._timer_went_off
            
            if geo.scfiterfail():
                #We need to quit and save everything
                self.stop = True
                logger.info(f"SCF iterations exceeded")
                self.parameters["Resubmit"] = False
                self.next_step = None
                return

            elif os.path.isfile("GEO_OPT_CONVERGED"):
                logger.info("Active site optimization converged")

            elif os.path.isfile("energy") and os.path.isfile("GEO_OPT_FAILED"):
                with open("energy", 'r') as energyfile:
                    i = -1
                    for l in energyfile:
                        i += 1
                    
                    total_cycles = i - 1
                    if total_cycles == -2:
                        self.stop = True
                        logger.warn("Energy file is empty")
                        self.next_step = None
                        self.parameters["Resubit"] = False
                        return

                logger.info(f"Failed to converge in a total of {total_cycles} geometry cycles")

            else:
                logger.error("TURBOMOLE OPTIMIZATION FAILED TO RUN PROPERLY!")
                logger.error("THERE COULD BE AN ERROR WITH TURBOMOLE")
                raise OSError("Turbomole")

        if self.stop:
            logger.debug("Stop received while doing optimization calculation")
        
        else:
            self.qm_final_energy = qm_calculation.TMcalculation.get_energy()
            logger.info("")
            logger.info(f"[Final QM Energy] ==>> {self.qm_final_energy}")
        
        #Now we reinstall the coords into the protein!
        self.to_next_iteration = dmd_to_qm.coord_to_protein(self.pdb_winner[0])


        #Now we reinstall into the protein
        logger.debug(f"Changing directory to {os.path.abspath('..')}")
        os.chdir("..")

        if not self.stop:
            self.to_next_iteration.write_pdb("to_next_iteration.pdb")
        
        self.next_step = self.finish_iteration

        # Check if we do a forceopt submit:
        if qm_params['geo_iterations'] <= total_cycles and self.parameters["qm params"]["calculation"] == "forceopt" and not os.path.isfile("Optimization/GEO_OPT_CONVERGED"):
            # Change calculation type to forceopt
            os.chdir("Optimization")
            with open("definput.json", "w") as definput:
                json.dump(self.parameters["qm params"], definput)

            logger.info(">>>> Submitting QM Optimization, forceopt >>>>")
            submitturbomole.main(_cores=self.cores, _time=self.controller._time, _nodes=1, _sub=True)
            os.chdir("..")

        logger.info("")
        logger.info("===================[Finished QM OPT]===================")


    def finish_iteration(self):
        logger.info("=====================[Saving Data]=====================")
        logger.info("")
        
        logger.debug("Saving energies to phd_energy")
        if self.to_next_iteration is None:
            logger.error("No final PDB structure!")
            raise ValueError("No final PDB structure")

        outfile_lines = []
        if os.path.isfile(os.path.join(self.root_directory, "phd_energy")):
            #Just append
            with open(os.path.join(self.root_directory, "phd_energy"), 'r') as outfile:
                for line in outfile:
                    outfile_lines.append(line)

            first_line = outfile_lines.pop(0)
            if outfile_lines:
                prev_iter = [int(i.split()[0]) for i in outfile_lines][-1]

                while prev_iter + 1 != self.iter_number and outfile_lines:
                    prev_iter = [int(i.split()[0]) for i in outfile_lines][-1]
                    outfile_lines.pop(0)
            
        else:
            first_line = ("[Iteration]        [QM Energy (Hart)]        [DMD Energy (kcal)]\n")

        outfile_lines.append(f"{self.iter_number:0>2d}                 {self.qm_final_energy:.5f}                 {self.pdb_winner[1]:.5f}\n")
        with open(os.path.join(self.root_directory, "phd_energy"), 'w') as outfile:
            outfile.write(first_line)
            for line in outfile_lines:
                outfile.write(line)

        #TODO clean up the iteration directory

        logger.info(">>>> Cleaning up iteration directory >>>>")
        for d in os.listdir("."):
            if "sp_movie" in d:
                logger.info(f"[Removing] ==>> {d}")
                shutil.rmtree(d)

        for d in os.listdir("dmd"):
            if "equilibrate" in d:
                logger.info(f"[Removing] ==>> {d}")
                shutil.rmtree(f"dmd/{d}")

            elif d.startswith("dmdstep_"):
                logger.info(f"[Removing] ==>> {d}")
                shutil.rmtree(f"dmd/{d}")


        logger.info("")
        logger.info("=====================[Data  Saved]=====================")
        self.next_step = None


    def signal_alarm(self):
        #This is called if the controller ran out of time
        self.stop = True

        #Create a stop file for the geometry optimization...otherwise, just have to wait for it to stop itself
        if self.next_step == self.qm_optimization:
            logger.info("MAKING STOP FILE")
            with open(os.path.join(self.directory, "Optimization/stop"), 'w') as stopfile:
                pass

    
    @staticmethod
    def dmd_converged( dmd_average_energies):
        if abs(dmd_average_energies[-1][0]-dmd_average_energies[-2][0]) < dmd_average_energies[-2][1]:
            if abs(dmd_average_energies[-1][0]-dmd_average_energies[-2][0]) < dmd_average_energies[-1][1]:
                return True

        return False


