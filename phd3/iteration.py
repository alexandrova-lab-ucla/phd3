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
from sklearn import cluster

from phd3.dmd_simulation import dmd_simulation
import phd3.qm_calculation as qm_calculation
import phd3.utility.utilities as utilities
from phd3.setupjob import setupTMjob
import phd3.utility.constants as constants
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

        if self.iter_number == -1:
            logger.info("On iteration 0, we skip the DMD portion")
            logger.info("and go straight to QM Optimization!")
            self.pdb_winner = [utilities.load_pdb(self.parameters["last pdb"]), 0.0]
            self.pdb_winner[0].name = "movie_0000"
            self.next_step = self.qm_optimization
            logger.info("")
            logger.info("====================[Finished  DMD]====================")
            return

        def print_energy(energy, kinetic, pressure, temperature, time):
            logger.info(f"[Ave. Pot. Energy] ==>> {energy[0]:.5f} ({energy[1]:.5f}) kcal/mol")
            logger.info(f"[Ave. Kin. Energy] ==>> {kinetic[0]:.5f} ({kinetic[1]:.5f}) kcal/mol")
            logger.info(f"[Ave. Tot. Energy] ==>> {(energy[0] + kinetic[0]) / 2.0:.5f} kcal/mol")
            logger.info(f"[Ave.  Pressure  ] ==>> {pressure[0]:.5f} ({pressure[1]:.5f})")
            logger.info(f"[Ave. Temperature] ==>> {temperature[0]:.5f} ({temperature[1]:.5f})")
            logger.info(f"[Est. Phys. Time ] ==>> {self.parameters['dmd params']['Time']*0.0000488882} ns")
            logger.info(f"Time elapsed during DMD simulation: {datetime.timedelta(seconds = time)}")

        #This tracks to see if we are done with the dmd simulations!
        finished = False
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
                self.dmd_steps = len(dirs)+1
                #Make sure that the directories are sorted
                dirs.sort(key = lambda i: int(i.split("_")[-1]))
                for i, d in enumerate(dirs):
                    echo_lines = []
                    if not os.path.isfile(os.path.join(d, self.parameters['dmd params']['Echo File'])):
                        logger.info(f"Did not complete run {i+1}, removing directory and continuing")
                        shutil.rmtree(d)
                        continue

                    with open(os.path.join(d, self.parameters['dmd params']['Echo File'])) as echo_file:
                        for l in echo_file:
                         echo_lines.append(l)

                    last_time = int(float(echo_lines[-1].split()[0]))
                    if last_time < self.parameters['dmd params']['Time']:
                        logger.info(f"Did not complete run {i+1}, removing directory and continuing")
                        shutil.rmtree(d)

                    else:
                        ave, stdev = dmd_simulation.get_average_potential_energy(os.path.join(d, self.parameters["dmd params"]["Echo File"]))
                        dmd_average_energies.append([ave, stdev])
                        logger.info(f"[Run {i+1:0>2d}] ==>> (finished)") 
            
                if len(dmd_average_energies) >= 2:
                    if abs(dmd_average_energies[-1][0]-dmd_average_energies[-2][0]) < dmd_average_energies[-2][1] and abs(dmd_average_energies[-1][0]-dmd_average_energies[-2][0]) < dmd_average_energies[-1][1]:
                        logger.info("Convergence previously achieved")
                        logger.info(f"Absolute delta in energy: {abs(dmd_average_energies[-1][0]-dmd_average_energies[-2][0]):.5f}")
                        finished = True
                        self.dmd_steps -= 1

            else:
                if os.path.isfile(self.parameters["dmd params"]["Echo File"]) and os.path.isfile(self.parameters["dmd params"]["Movie File"]):
                    ave, stdev = dmd_simulation.get_average_potential_energy(self.parameters["dmd params"]["Echo File"])
                    logger.info(f"[Run 01] ==>> (finished)")
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
                    logger.info(">>>> Equilibrating structure >>>>")
                    logger.info(f"[Starting Temp.]  ==>> {self.parameters['Equilibrate']['Initial Temperature']}")
                    logger.info(f"[Ending Temp.]    ==>> {self.parameters['dmd params']['Initial Temperature']}")
                    logger.info(f"[Time]            ==>> {self.parameters['Equilibrate']['Time']}")
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

                    logger.debug("Changing directory to {os.getcwd()}")
                    os.chdir("..")


            if self.parameters["DMD CONVERGE"]:
                logger.info("")
                logger.info(">>>> Attempting to converge DMD simulations >>>>")
                
                while not self.stop and self.dmd_steps <= self.parameters["MAX DMD STEPS"]:
                    #TODO try and except the following
                    if os.path.isdir(f"dmdstep_{self.dmd_steps-1}"):
                        logger.debug(f"Copying {self.dmd_steps-1} files to new dmd start")

                        shutil.copytree(f"dmdstep_{self.dmd_steps-1}", f"dmdstep_{self.dmd_steps}")
                        logger.debug("Removing echo and movie file from new directory")
                        os.remove(os.path.join(f"dmdstep_{self.dmd_steps}", self.parameters["dmd params"]["Echo File"]))
                        os.remove(os.path.join(f"dmdstep_{self.dmd_steps}", self.parameters["dmd params"]["Movie File"]))

                    elif os.path.isdir("equilibrate") and os.path.isdir(f"equilibrate/{self.parameters['dmd params']['Restart File']}"):
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
                    start = timer()
                    calc = dmd_simulation(cores = self.cores, parameters = self.parameters["dmd params"])
                    #Harvest the time that the dmd run took
                    end = timer()

                    curr_energy = dmd_simulation.get_average_potential_energy(self.parameters["dmd params"]["Echo File"])
                    kinetic_energy = dmd_simulation.get_average_kinetic_energy(self.parameters["dmd params"]["Echo File"])
                    press = dmd_simulation.get_average_pressure_energy(self.parameters["dmd params"]["Echo File"]) 
                    temp = dmd_simulation.get_average_temp_energy(self.parameters["dmd params"]["Echo File"])
                    #Prints nicely out what just happened
                    print_energy(curr_energy, kinetic_energy, press, temp, int(end-start))
                    
                    dmd_average_energies.append(curr_energy)
                    #if converged then break
                    if len(dmd_average_energies) >= 2:
                        if abs(dmd_average_energies[-1][0]-dmd_average_energies[-2][0]) < dmd_average_energies[-2][1] and abs(dmd_average_energies[-1][0]-dmd_average_energies[-2][0]) < dmd_average_energies[-1][1]:
                            logger.info("")
                            logger.info("Converence achieved!")
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
                
                start = timer()
                calc = dmd_simulation(cores = self.cores, parameters = self.parameters["dmd params"])
                end = timer()
                self.final_dmd_average_energy = dmd_simulation.get_average_potential_energy()
                kinetic_energy = dmd_simulation.get_average_kinetic_energy(self.parameters["dmd params"]["Echo File"])
                press = dmd_simulation.get_average_pressure_energy(self.parameters["dmd params"]["Echo File"]) 
                temp = dmd_simulation.get_average_temp_energy(self.parameters["dmd params"]["Echo File"])
                #Prints nicely out what just happened
                print_energy(curr_energy, kinetic_energy, press, temp, int(end-start))

        if self.stop:
            logger.info("Stop signal received while performing DMD simulation")
        
        else:
            if self.dmd_steps > self.parameters["MAX DMD STEPS"]:
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

            elif self.parameters["DMD CONVERGE"]:
                logger.info("")
                logger.info(f"DMD converged in: {self.dmd_steps} steps")

            logger.info("")
            logger.info(">>>> DMD Summary >>>>")
            logger.info("[RUN ##]---------[Ave Pot. Energy]---------[Est. Phys Time (ns)]---------")
            for index, energy in enumerate(dmd_average_energies):
                out = f"[Run {index+1:0>2d}]      {energy[0]:.5f} ({energy[1]:.5f})\t\t{self.parameters['dmd params']['Time']*0.0000488882}"
                if energy == self.final_dmd_average_energy:
                    out += "    --w---w--"
                
                logger.info(out)

            #convert movie to a movie.pdb file
            logger.debug("Making movie.pdb")
            utilities.make_movie("initial.pdb", self.parameters["dmd params"]["Movie File"], "movie.pdb")

            logger.debug("Loading in movie.pdb")
            self.dmd_structures = utilities.load_movie("movie.pdb")

            dmd_structures_energies = []
            if not os.path.isfile(self.parameters['dmd params']['Echo File']):
                logger.error("No echo file! Cannot properly cluster")
                raise FileNotFoundError("Echo File")

            with open(self.parameters['dmd params']['Echo File']) as echo_file:
                for line in echo_file:
                    if line[0] != '#':
                        line = line.split()
                        dmd_structures_energies.append(float(line[4]))

            #Now we are done
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

        min_structures = []
        centroid_structures = []

        if [d for d in os.listdir() if "movie_" in d]:
            logger.info("Loading in previous clustering results")
            logger.info("")
            self.sp_PDB_structures = []
            for f in [d for d in os.listdir() if d.startswith("movie_") and os.path.isfile(d)]:
                f = f.split("_")
                f = f[1].split(".")
                num = int(f[0])
                self.sp_PDB_structures.append(self.dmd_structures[num])

        else:
            # Create a zero  matrix of the appropriate size, will fill in with the RMSD between structures
            distance_matrix = np.array([np.array([0.0 for i in self.dmd_structures]) for j in self.dmd_structures])

            logger.info("Computing RMSD between all pairs of structures")
            logger.info("")
            # This is a symmetrix matrix
            for row, first_structure in enumerate(self.dmd_structures):
                for col in range(row, len(self.dmd_structures)):
                    distance_matrix[row][col] = first_structure[0].aa_rmsd(self.dmd_structures[col][0])
                    distance_matrix[col][row] = distance_matrix[row][col] 
            
            #Sets up
            assert(self.parameters["Number of Clusters"] > 2)
            agg = cluster.AgglomerativeClustering(n_clusters=self.parameters["Number of Clusters"], affinity='precomputed', linkage="average" )

            logger.info(">>>> Cluster Parameters >>>>")
            logger.info(f"[Frames Collected]  ==>> {len(self.dmd_structures)}")
            logger.info(f"[Num. of Clusters]  ==>> {self.parameters['Number of Clusters']}")
            logger.info("[Cluster Method]    ==>> Agglomerative Clustering")
            logger.info(f"[Minimum Energy]    ==>> {'true' if self.parameters['Cluster Energy'] else 'false'}")
            logger.info(f"[Cluster Centroid]  ==>> {'true' if self.parameters['Cluster Centroid'] else 'false'}")
            logger.info("")
            
            #Actually cluster
            output = agg.fit(distance_matrix)

            #Reformat so that we can work the various clusters
            clustered_structures = {}
            for index in range(self.parameters["Number of Clusters"]):
                clustered_structures[index] = []
                for struct in zip(self.dmd_structures, output.labels_):
                    if struct[1] == index:
                        clustered_structures[index].append(struct[0])

            #Now we can access the clusters by the keys in this dictionary
            # This is a dictionary, with values of a list of tuples (protein, dmd energy)
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
                    distance_matrix = np.array([np.array([0.0 for i in self.dmd_structures[index]]) for j in self.dmd_structures[index]])

                   #Create another distance matrix...then find the sum of each row and take the minimum 
                    for row, first_structure in enumerate(self.dmd_structures[index]):
                        for col in range(row, len(self.dmd_structures[index])):
                            distance_matrix[row][col] = first_structure[0].aa_rmsd(self.dmd_structures[index][col][0])
                            distance_matrix[col][row] = distance_matrix[row][col] 

                    rows_sums = [sum(i) for i in distance_matrix]
                    with_rows = list(zip(self.dmd_structures[index], rows_sums))
                    winner = min(with_rows, key = lambda i: i[1])
                    centroid_structures.append(winner[0])

            self.sp_PDB_structures.extend(centroid_structures)
            
            if not self.sp_PDB_structures:
                logger.info("No clustering method specified, using minimum energy from each cluster")
                for index in self.dmd_structures.keys():
                    min_structures.append(min(self.dmd_structures[index], key = lambda i: i[1]))
       
                self.sp_PDB_structures.extend(min_structures)

        #Print out choices, and save the files
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
        logger.info("")

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
            if os.path.isfile("energy"):
                with open("energy", 'r') as energyFile:
                    for i,l in enumerate(energyFile):
                        pass

                if i >= 2:
                    self.qm_sp_energies.append(qm_calculation.TMcalculation.get_energy(cycle=1))
                    logger.info(f"[QM Energy (Hart)] ==>> {self.qm_sp_energies[-1]:.5f}")
                    os.chdir("../")
                    continue
                
            elif not os.path.isfile("coord"):
                logger.debug("Creating coord file")
                pdb_to_coord.protein_to_coord(struct, self.parameters["QM Chop"])

            if not os.path.isfile("control"):
                logger.debug("Setting up TM job")
                sj = setupTMjob(parameters = self.parameters["qm params"], timeout=0)

                dirs = [d for d in os.listdir("../") if os.path.isdir(os.path.join("../", d)) and "sp_movie_" in d]
                dirs.sort(reverse = True, key=lambda i: int(i.split("_")[-1]))
                
                for sp_directory in dirs:
                    if os.path.abspath(os.path.join("..", sp_directory)) == os.path.abspath(os.getcwd()):
                        continue

                    sp_directory = os.path.join("../", sp_directory)
                    found_mos = False
                    if os.path.isfile(os.path.join(sp_directory, "mos")):
                        logger.debug("Copying over mos files")
                        found_mos = True
                        shutil.copy(os.path.join(sp_directory, "mos"), "./mos")

                    if os.path.isfile(os.path.join(sp_directory, "alpha")):
                        logger.debug("Copying over alpha files")
                        found_mos = True
                        shutil.copy(os.path.join(sp_directory, "alpha"), "./alpha")
                    
                    if os.path.isfile(os.path.join(sp_directory, "beta")):
                        logger.debug("Copying over beta files")
                        found_mos = True
                        shutil.copy(os.path.join(sp_directory, "beta"), "./beta")

                    if found_mos:
                        break
                
                if not found_mos:
                    #TODO allow user to specify a default mos, alpha, beta to use!
                    pass

            qm_params = self.parameters["qm params"].copy()
            qm_params["calculation"] = 'sp'
            start = timer()
            #sp = qm_calculation.TMcalculation(self.cores, parameters=qm_params)
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

            logger.info(f"[QM Energy (Hart)] ==>> {self.qm_sp_energies[-1]:.5f}")
            logger.info(f"Time elapsed during QM SP calculation: {datetime.timedelta(seconds = int(end -start))}")
            logger.debug(f"Changing directory from {os.getcwd()} to {os.path.abspath('..')}")
            os.chdir("../")

            if self.stop:
                logger.debug("Stop received while doing singlepoint calculations")
                return

        logger.info("")
        logger.info(">>>> Scoring Structures >>>>")

        #We are done with all sp calculations now
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

        logger.info(f"[movie_####]------[DMD Inst. Pot (kcal)]------[QM Energy (Hart)]--------[Score]------")
        for struct in self.scoring_energies:
            out_string = f"[{struct[0].name}]         {struct[1]:0<12}                {struct[2]:0<12}            {struct[3]:.3f}"
            if struct[3] == winning_score:
                out_string += "      --w---w--"

            logger.info(out_string)

        logger.info("")
        logger.info(f"[WINNER]     ==>> {self.sp_PDB_structures[winning_index][0].name}")
        logger.info(f"[QM  Energy] ==>> {self.qm_sp_energies[winning_index]:.5f} Hartrees")
        logger.info(f"[DMD Energy] ==>> {self.sp_PDB_structures[winning_index][1]:.5f} kcal/mol")

        self.pdb_winner = self.sp_PDB_structures[winning_index]
        self.pdb_winner[0].write_pdb("dmdWinner.pdb")

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
                pdb_to_coord.protein_to_coord(self.pdb_winner[0], self.parameters["QM Chop"])
 
            if not os.path.isfile("control"):
                logger.debug("Setting up TM job")
                sj = setupTMjob(parameters = qm_params, timeout=0)

                possible_sp_dir = "../sp_" + self.pdb_winner[0].name
                found_mos = False
                if os.path.isdir(possible_sp_dir):
                    if os.path.isfile(os.path.join(possible_sp_directory), "mos"):
                        logger.debug("Copying over mos files")
                        found_mos = True
                        shutil.copy(os.path.join(possible_sp_directory), "mos", "./mos")

                    if os.path.isfile(os.path.join(possible_sp_directory), "alpha"):
                        logger.debug("Copying over alpha files")
                        found_mos = True
                        shutil.copy(os.path.join(possible_sp_directory), "alpha", "./alpha")
                    
                    if os.path.isfile(os.path.join(possible_sp_directory), "beta"):
                        logger.debug("Copying over beta files")
                        found_mos = True
                        shutil.copy(os.path.join(possible_sp_directory), "beta", "./beta")

                else:
                    #Check for any of these directories actually...
                    dirs = [d for d in os.listdir("../") if os.path.isdir(os.path.join("../", d)) and "sp_movie_" in d]
                    dirs.sort(reverse = True, key=lambda i: int(i.split("_")[-1]))
                    
                    for sp_directory in dirs:
                        sp_directory = os.path.join("../", sp_directory)
                        if os.path.isfile(os.path.join(sp_directory), "mos"):
                            logger.debug("Copying over mos files")
                            found_mos = True
                            shutil.copy(os.path.join(sp_directory), "mos", "./mos")

                        if os.path.isfile(os.path.join(sp_directory), "alpha"):
                            logger.debug("Copying over alpha files")
                            found_mos = True
                            shutil.copy(os.path.join(sp_directory), "alpha", "./alpha")
                        
                        if os.path.isfile(os.path.join(sp_directory), "beta"):
                            logger.debug("Copying over beta files")
                            found_mos = True
                            shutil.copy(os.path.join(sp_directory), "beta", "./beta")

                        if found_mos:
                            break
                
                if not found_mos:
                    #TODO allow user to specify a default mos, alpha, beta to use!
                    pass

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
            logger.info("[Final QM Energy] ==>> {self.qm_final_energy}")
            logger.info(f"Time elapsed during QM Opt calculation: {datetime.timedelta(seconds = int(end -start))}")

        if self.stop:
            logger.debug("Stop received while doing optimization calculation")
        
        #Now we reinstall the coords into the protein!
        self.to_next_iteration = pdb_to_coord.coord_to_pdb(self.pdb_winner[0])

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
            
            outfile_lines.append(f"{self.iter_number:0>2d}            {self.qm_final_energy:.5f}            {self.pdb_winner[1]:.5f}\n")

            with open(os.path.join(self.root_directory, "phd_energy"), 'w') as outfile:
                outfile.write(first_line)
                for line in outfile_lines:
                    outfile.write(line)

        else:
            with open(os.path.join(self.root_directory, "phd_energy")) as outfile:
                outfile.write("[Iteration]        [QM Energy (Hart)]        [DMD Energy (kcal)]\n")
                outfile.write(f"{self.iter_number:0>2d}            {self.qm_final_energy:.5f}            {self.pdb_winner[1]:.5f}\n")
        

        #TODO clean up the iteration directory

        self.next_step = None


    def signal_alarm(self):
        #This is called if the controller ran out of time
        self._stop = True

        #Create a stop file for the geometry optimization...otherwise, just have to wait for it to stop itself
        if self.next_step == self.qm_optimization:
            with open(os.path.join(self.directory, "Optimization/stop"), 'w') as stopfile:
                pass




