#!/usr/bin/env python3
"""
Author  ==>> David J. Reilley
Author  ==>> Matthew R. Hennefarth
Date    ==>> April 16, 2020
"""

#Utility Paths
import phd3.utility.utilities as utilities
config = utilities.load_phd_config()

#Standard Library Imports
import logging
import os
import shutil
import pkg_resources

#3rd Party Libraries
import propka
from propka import run

#Titrate/PHD3
from . import montecarlo
from ..utility import constants, exceptions

logger = logging.getLogger(__name__)

#Turn off the propka logging
propka_logger = logging.getLogger("propka")
propka_logger.setLevel(logging.CRITICAL)
propka_logger.propagate = False


__all__ = [
        'titrate_protein'
        ]

PROTON_PARTNER_CUTOFF =3.5

class titrate_protein:

    __slots__ = ['_updated_protonation', '_pH', '_buried_cutoff', '_partner_dist', "_step", "_history", "_default_protonation_states"]

    @staticmethod
    def expand_commands(parameters):
        logger.debug("Expanding parameters")

        if parameters["Commands"]:
            expanded_command = {}
            for command in parameters["Commands"].keys():
                updated_params = parameters.copy()
                for key in parameters["Commands"][command].keys():
                    updated_params[key] = parameters["Commands"][command][key]

                steps = int(updated_params["Time"]/updated_params["titr"]["dt"])
                for s in range(steps):
                    expanded_command[f"{command}:{s}"] = parameters["Commands"][command].copy()
                    expanded_command[f"{command}:{s}"]["Time"] = updated_params["titr"]["dt"]

            parameters["Commands"] = expanded_command

        else:
            steps = int(parameters["Time"] / parameters["titr"]["dt"])
            for s in range(steps):
                parameters["Commands"][f":{s}"] = {"Time" : parameters["titr"]["dt"]}
        
        return parameters
    

    @staticmethod
    def condense_commands(parameters):
        logger.debug("Condesing commands")

        command_label = ""
        time = 0
        condensed_commands = {}
        for command in parameters["Remaining Commands"].keys():
            command_split = [c for c in command.split(":") if c != ""]
            if len(command_split) == 1:
                #we have only one command
                if command_label == "":
                    command_label = "A"
                    condensed_commands[command_label] = parameters["Remaining Commands"][command].copy()

                time += parameters["Remaining Commands"][command]["Time"]
            

            elif command_label == command_split[0]:
                time += parameters["Remaining Commands"][command]["Time"]

            else:
                #new command
                if command_label in condensed_commands:
                    condensed_commands[command_label]["Time"] = time

                time = parameters["Remaining Commands"][command]['Time']
                command_label = command_split[0]
                condensed_commands[command_label] = parameters["Remaining Commands"][command].copy()

        condensed_commands[command_label]["Time"] = time

        parameters["Commands"].clear()
        parameters["Commands"] = condensed_commands
        parameters["Remaining Commands"].clear()
        return parameters

    def __init__(self, parameters, initial_protonation = []):

        #Parameters for titration
        self._pH = parameters["pH"]
        self._buried_cutoff = parameters["Buried Cutoff"]
        self._partner_dist = parameters["Partner Distance"]
        
        if "Fixed States" in parameters.keys():
            self._default_protonation_states = parameters["Fixed States"]
        
        else:
            self._default_protonation_states = []

        self._history = []

        # Sets initial protonation states
        self._updated_protonation = initial_protonation.copy()        
        # Make save directory
        if os.path.isdir("save"):
            pkas = [f for f in os.listdir("save") if ".pka" in f]
            f = [int(f.split(".")[0]) for f in pkas]
            self._step = max(f)
            if os.path.isfile("inConstr"):
                shutil.copy("inConstr", f"save/{self._step}.inConstr")
                self._step += 1

        else:
            self._step = 0

    def history(self):
        return self._history

    def evaluate_pkas(self, protein):
        #First we transform protein to Standard
        protein.relabel(format="Standard")
        
        protein.write_pdb("_propka_inp.pdb")

        #Then we call propka from the import

        if not hasattr(run, "single"):
            logger.error("Propka not properly installed, or wrong version")
            raise exceptions.Propka_Error

        try:
            my_molecule = run.single("_propka_inp.pdb")

        except:
            logger.error("Error running propka")
            raise exceptions.Propka_Error

        if os.path.isfile("_propka_inp.propka_input"):
            logger.debug("Removing file: _propka_inp.propka_input")
            os.remove("_propka_inp.propka_input")


        logger.info("[propka]           ==>> SUCCESS")

        #And run MSMS to generate the SAS if called for
        if self._buried_cutoff == "sas":
            try:
                pdb_to_xyzrn_cmd = config["PATHS"]["MSMS_DIR"] + 'pdb_to_xyzrn _propka_inp.pdb > _msms_inp.xyzrn'
                msms_cmd = config["PATHS"]["MSMS_DIR"] + 'msms.x86_64Linux2.2.6.1 -if _msms_inp.xyzrn -af _msms_out'
                os.system(pdb_to_xyzrn_cmd)
                os.system(msms_cmd)
            except:
                logger.error("Error running MSMS")
                raise exceptions.MSMS_Error
        
        #Now we move onto davids actual script for evaluation of the protons and what not


        #SAVE THE DATA
        if not os.path.isdir("save"):
            os.mkdir("save")

        if os.path.isfile("_propka_inp.pka"):
            shutil.copy("_propka_inp.pka", f"save/{self._step}.pka")
        
        if os.path.isfile("inConstr") and self._step > 0:
            shutil.copy("inConstr", f"save/{self._step-1}.inConstr")
        
        if os.path.isfile("_msms_out.area"):
            shutil.copy("_msms_out.area", f"save/{self._step}.area")
            
        self._step += 1

        with open("_propka_inp.pdb", 'r') as in_pdb:
            #Get all of the titratable residues as a list
            titratable_residues = montecarlo.process_pdb(in_pdb.readlines())
        
        ## REMOVE THE RESIDUES HERE THAT ARE STATIC

        remove_residues = []
        logger.debug("Removing static protonation state residues...")
        for i, titratable_residue in enumerate(titratable_residues):
            residue = protein.get_residue([titratable_residue.chain, int(titratable_residue.res_num)])
            for default_state in self._default_protonation_states:
                default_state_residue = protein.get_residue(default_state)
                if residue == default_state_residue:
                    remove_residues.append(i)

        remove_residues.reverse()
        [titratable_residues.pop(i) for i in remove_residues]


        #Define connections between residues
        montecarlo.define_connections(titratable_residues, PROTON_PARTNER_CUTOFF)
        
        if titratable_residues[0].chain:
            chains = True

        else:
            chains = False

        #propka output, titratable residues, and if we have multiple chains...
        calc_pKa_data = montecarlo.calc_pKa_total_pdb("_propka_inp.pka", titratable_residues, chains)
        for res in titratable_residues:
            res.assign_pKa(calc_pKa_data)

        solv_data = montecarlo.find_solv_shell("_propka_inp.pka", chains)

        if self._buried_cutoff == "sas":
            msms_data = montecarlo.store_sas_area("_msms_out.area", chains)
        
        titr_stack = [] # Construct the stack form of all_titr_res for use in find_solv_shell
        for res in titratable_residues:
            titr_stack += [res]

        #Need to make a copy so that we don't accidently screw up out list
        #Should check this...i think we want a copy of a list, but it pointing to the same res in all_titr_res
        all_networks = montecarlo.define_aa_networks(titr_stack)
        if self._buried_cutoff == "sas":
            all_networks = montecarlo.find_network_solvent_access(all_networks, msms_data, self._buried_cutoff, self._partner_dist)
        else:
            all_networks = montecarlo.find_network_solvent_access(all_networks, solv_data, self._buried_cutoff, self._partner_dist)
        
        #Now we do monte carlo
        montecarlo.MC_prot_change(all_networks, self._pH)
        for residue in titratable_residues:
            residue.update_prots()
        
        protonation_changes = []
       
        remove = []
        switch_his = []

        for residue in titratable_residues:
            if residue.change[0] != "None":
                change = []
                if residue.change[0] == "Add":
                    protonate = True

                elif residue.change[0] == "Remove":
                    protonate = False

                else:
                    logger.warn("Unknown residue change command")
                    continue
               
                change.extend([residue.chain, int(residue.res_num)])
                
                # We check to see if the residue is static in regard to the protonation state...
                # change_residue = protein.get_residue(change)
                # ignore = False
                # for default_state in self._default_protonation_states:
                    # default_state_residue = protein.get_residue(default_state)
                    # if change_residue == default_state_residue:
                        # logger.warn(f"Cannot change protonation state of static residue {default_state_residue}")
                        # ignore = True

                # if ignore:
                    # continue

                change.append("protonate" if protonate else "deprotonate" )
                
                #Use -1 and -2 to deprotonate the C and N Terminus
                if residue.ter_name == 'N+' and not protonate:
                    #check if we are deprotonating the n-terminus
                    change.append(-1)
        
                elif residue.ter_name == 'C-' and protonate:
                    #change.append(-2)
                    logger.warn("Tried to protonate the c-terminus")
                    logger.warn("This has been turned off permanently due to DMD issues")
                    continue

                else:
                    if protonate:
                        if residue.amino_acid.upper() not in constants.PROTONATED_STANDARD.keys():
                            logger.debug("Cannot protonate residue {residue.amino_acid}")
                            remove.append(change[:2])
                            continue
                        
                        elif residue.amino_acid.upper() == "HIS" and residue.change_heteroatom[0] == "ND1":
                            switch_his.append(change[:2])
                            continue

                        for index, pairs in enumerate(constants.PROTONATED_STANDARD[residue.amino_acid.upper()]):
                            if residue.change_heteroatom[0] == pairs[0]:
                                change.append(index+1)
                                break

                    else:
                        if residue.amino_acid.upper() not in constants.DEPROTONATED_STANDARD.keys():
                            logger.debug("Cannot deprotonate residue {residue.amino_acid}")
                            remove.append(change[:2])
                            continue
                        
                        elif residue.amino_acid.upper() == "HIS" and residue.change_heteroatom[0] == "NE2":
                            switch_his.append(change[:2])
                            continue
                        
                        for index, pairs in enumerate(constants.DEPROTONATED_STANDARD[residue.amino_acid.upper()]):
                            if residue.change_heteroatom[0] == pairs[0]:
                                change.append(index+1)
                                break

                    if len(change) != 4:
                        logger.warn(f"Cannot change protonation state of atom {residue.change_heteroatom[0]} in res {residue.amino_acid} {residue.res_num}")
                        continue

                protonation_changes.append(change)

        for change in protonation_changes:
            residue = protein.get_residue(change[:2])
            for i,current in enumerate(self._updated_protonation):
                current_residue = protein.get_residue(current[:2])
                if residue == current_residue:
                    self._updated_protonation[i] = change
                    #current = change
                    if residue.name.upper() == "HIS":
                        if change[2] == "protonate":
                            pass
                        
                        else:
                            self._updated_protonation.append([change[0], change[1], "protonate"])
                    
                    break
                
            else:
                if residue.name.upper() == "HIS":
                    if change[2] == "protonate":
                        pass
                    
                    else:
                        self._updated_protonation.append([change[0], change[1], "protonate"])
                
                self._updated_protonation.append(change)
        

        if switch_his + remove:
            for switch in switch_his + remove:
                switch_res = protein.get_residue(switch[:2])
                remove_index = []
                for i in range(len(self._updated_protonation)):
                    residue = protein.get_residue(self._updated_protonation[i][:2])
                    if residue == switch_res:
                        remove_index.append(i)

                remove_index.reverse()
                [self._updated_protonation.pop(i) for i in remove_index]

        # Assign protonations to self._updated_protonation
        # List -> Tuple -> Set -> List to get rid of duplicates
        self._updated_protonation = [list(item) for item in set(tuple(row) for row in self._updated_protonation)]
        for state in self._updated_protonation:
            residue = protein.get_residue(state[:2])
            if residue.name.upper() == "HIS" and state[2] == "deprotonate":
                for other_state in self._updated_protonation:
                    other_residue = protein.get_residue(other_state[:2])
                    if residue == other_residue and other_state[2] == "protonate":
                        break

                else:
                    print("BRUHHHH!!!!!!!")
                    raise exceptions.Propka_Error

        self._history.append(self._updated_protonation.copy())
        return self._updated_protonation


    def get_new_protonation_states(self):
        return self._updated_protonation if self._updated_protonation is not None else []

