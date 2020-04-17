#!/usr/bin/env python3
"""
Author  ==>> David J. Reilley
Author  ==>> Matthew R. Hennefarth
Date    ==>> April 16, 2020
"""

#Standard Library Imports
import logging
import os
import pkg_resources

#Rd Party Libraries
import propka.molecular_container

#Titrate/PHD3
import phd3.titrate.montecarlo as montecarlo
import phd3.utility.constants as constants

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

    __slots__ = ['_updated_protonation', '_pH', '_buried_cutoff', '_partner_dist']

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
        for command in parameters["Commands"].keys():
            if command_label == command.split(":")[0]:
                time += parameters["Commands"][command]["Time"]

            else:
                #new command
                if command_label in condensed_commands:
                    condensed_commands[command_label]["Time"] = time

                time = parameters["Commands"][command]['Time']
                command_label = command.split(":")[0]
                condensed_commands[command_label] = parameters["Commands"][command].copy()

        parameters["Commands"] = condensed_commands
        return parameters

    def __init__(self, parameters):

        #Parameters for titration
        self._pH = parameters["pH"]
        self._buried_cutoff = parameters["Buried Cutoff"]
        self._partner_dist = parameters["Partner Distance"]

        self._updated_protonation = None

    def evaluate_pkas(self, protein):
        #First we transform protein to Standard
        protein.relabel(format="Standard")
        
        protein.write_pdb("_propka_inp.pdb")

        #Holds the default options for propka to use as an imported module
        class option:
            def __init__(self):
                self.keep_protons = False
                self.protonate_all = False
                self.parameters = pkg_resources.resource_filename("propka", "propka.cfg")
                self.chains = []
                self.titrate_only = None
                self.display_coupled_residues = False

        #Then we call propka from the import
        my_molecule = propka.molecular_container.Molecular_container("_propka_inp.pdb", option())
        my_molecule.calculate_pka()
        my_molecule.write_pka()
        if os.path.isfile("_propka_inp.propka_input"):
            logger.debug("Removing file: _propka_inp.propka_input")
            os.remove("_propka_inp.propka_input")

        logger.info("[propka]           ==>> SUCCESS")
        #Now we move onto davids actual script for evaluation of the protons and what not
        
        with open("_propka_inp.pdb", 'r') as in_pdb:
            #Get all of the titratable residues as a list
            titratable_residues = montecarlo.process_pdb(in_pdb.readlines())
        
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
        
        titr_stack = [] # Construct the stack form of all_titr_res for use in find_solv_shell
        for res in titratable_residues:
            titr_stack += [res]

        #Need to make a copy so that we don't accidently screw up out list
        #Should check this...i think we want a copy of a list, but it pointing to the same res in all_titr_res
        all_networks = montecarlo.define_aa_networks(titr_stack)
        all_networks = montecarlo.find_network_solvent_access(all_networks, solv_data, self._buried_cutoff, self._partner_dist)
        
        #Now we do monte carlo
        montecarlo.MC_prot_change(all_networks, self._pH)
        for residue in titratable_residues:
            residue.update_prots()


        #Any cleanup if necessary
        #Also, may want to print out new protein with the protonation states...but I don't know why we would really care fo that
        if self._updated_protonation is not None:
            self._updated_protonation.clear()
        
        else:
            self._updated_protonation = []
        
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
                change.append("protonate" if protonate else "deprotonate" )
                
                #Use -1 and -2 to deprotonate the C and N Terminus
                if residue.ter_name == 'N+' and not protonate:
                    #check if we are deprotonating the n-terminus
                    change.append(-1)
        
                elif residue.ter_name == 'C-' and protonate:
                    change.append(-2)

                else:
                    if protonate:
                        if residue.amino_acid.upper() not in constants.PROTONATED_STANDARD.keys():
                            logger.debug("Cannot protonate residue {residue.amino_acid}")
                            continue
                        
                        for index, pairs in enumerate(constants.PROTONATED_STANDARD[residue.amino_acid.upper()]):
                            if residue.change_heteroatom[0] == pairs[0]:
                                change.append(index+1)
                                break

                    else:
                        if residue.amino_acid.upper() not in constants.DEPROTONATED_STANDARD.keys():
                            logger.debug("Cannot deprotonate residue {residue.amino_acid}")
                            continue
                        
                        for index, pairs in enumerate(constants.DEPROTONATED_STANDARD[residue.amino_acid.upper()]):
                            if residue.change_heteroatom[0] == pairs[0]:
                                change.append(index+1)
                                break

                    if len(change) != 4:
                        logger.warn(f"Cannot change protonation state of atom {residue.change_heteroatom[0]} in res {residue.amino_acid}")

                self._updated_protonation.append(change)

        #Assign protonations to self._updated_protonation
        return self._updated_protonation


    def get_new_protonation_states(self):
        return self._updated_protonation if self._updated_protonation is not None else []



