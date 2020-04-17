#!/usr/bin/env python3
"""
Author  ==>> David J. Reilley
Date    ==>> April 16, 2020
"""


#Standard Library Imports
import logging
import os
import numpy as np
import random
from itertools import combinations

from . import titrate_data

__all__ = [
        'process_pdb',
        'find_network_solvent_access',
        'define_connections',
        'calc_pKa_total_pdb',
        'find_solv_shell',
        'define_aa_networks',
        'MC_prot_change'
        ]

logger = logging.getLogger(__name__)

class titr_res:
    def __init__(self, amino_acid, res_num, chain, atoms, ter=''):
        self.name = amino_acid + res_num + chain
        self.amino_acid = amino_acid
        self.res_num = res_num
        self.chain = chain
        if ter:
            ter_res = ter
        
        else:
            ter_res = amino_acid
        
        self.ter_name = ter_res
        self.full_name = ter_res + res_num + chain # Like self.name but with 'C' or 'N' for C and N-terminal residues
        self.atoms = atoms # List of lists with atom type followed by cartesian coordinates for every atom in the residue
        self.prots = [] # List of lists with titratable protons, stored in same format as self.atoms
        self.heteroatoms = [] # List of lists with heteroatoms bound to the titratable protons, stored in the same format as self.atoms
        self.prot_state = [] # List with two entries, with the first a string + or - for protonated or deprotonated, and the second a number corresponding to its form
        self.partners = [] # List of titr_res instances, containing nearby amino acids that this amino acid can interact with
        self.pKa = None
        self.new_prot_state = ['-', 0] # New protonation state decided by script, same formatting as self.prot_state
        self.change = ['None'] # List with two entries, the first saying whether there is a change to protonation state and what it is, the second as a list specifying the added or removed proton(s) if there is a change
        self.change_heteroatom = [] # List: If there is an added proton, name of the heteroatom it is bound to (1st entry) and a series of reference atoms for hydrogen placement (following entries)
        self.covalent_link = False
        self.change_prob = 0.0
        self.change_roll = 0.0

    def define_prot_state(self):
        titr_atoms = []
        titr_het_atoms = []

        for atom in self.atoms: # Find protons and heteroatoms that define the titration state of the residue
            atom_name = self.ter_name + ':' + atom[0]
            if atom_name in titrate_data.titr_form_prots:
                self.prots.append(atom)
                titr_atoms.append(atom[0])

            elif atom_name in titrate_data.titr_form_heteroatoms:
                self.heteroatoms.append(atom)
                titr_het_atoms.append(atom[0])

        titr_atoms_ID = self.ter_name # Sort present titratable protons and heteroatoms to define protonation state
        titr_atoms.sort()
        titr_het_atoms.sort()
        for atom in titr_atoms:
            titr_atoms_ID += ':' + atom

        for atom in titr_het_atoms:
            titr_atoms_ID += ':' + atom

        #print(self.full_name)
        prot_state = titrate_data.prots2titr_form[titr_atoms_ID]
        prot_state = prot_state.split(':')
        self.prot_state = [prot_state[1], prot_state[2]]

    def update_ter(self, ter_res):
        self.ter = ter_res

    def assign_pKa(self, pKa_data): # Looks up pKa from dictionary pKa_data and saves this to class instance
        if self.full_name in pKa_data:
            self.pKa = pKa_data[self.full_name]
        
        # If for some reason PropKa3.1 didn't calculate the pKa for this residue,
        #set the pKa so that the protonation state can't change and notify the user
        #print('pKa not found for ' + self.full_name + ', skipping for now but the formatting for this residue should be checked')
        elif self.prot_state[0] == '-':
            self.pKa = 1.0
        
        else:
            self.pKa = 20.0

        if self.covalent_link: # Do not allow residues involved in covalent interactions with other residues change protonation state
            if self.prot_state[0] == '-':
                self.pKa = 1.0

            else:
                self.pKa = 20.0

    def create_new_prot(self, new_atom_name):
        ref_coords = np.zeros(3, dtype=float)
        heteroatom_coords = np.zeros(3, dtype=float)
        reference_atoms = self.change_heteroatom[1:]

        for atom in self.atoms: # Generate vector pointing from average of reference atoms to the heteroatom and project that from the heteroatom at a distance of 1 angstrom for the proton position
            if atom[0] == self.change_heteroatom[0]:
                heteroatom_coords += atom[1]
            
            elif atom[0] in reference_atoms:
                ref_coords += atom[1] / len(reference_atoms)

        prot_vec = heteroatom_coords - ref_coords
        new_prot_pos = prot_vec/np.linalg.norm(prot_vec) + heteroatom_coords

        self.atoms.append([new_atom_name, new_prot_pos])

    def update_prots(self):
        if self.change[0] == 'Remove':
            for index, atom in enumerate(self.atoms):
                if atom[0] == self.change[1][0]:
                    self.atoms.pop(index)
                    return

        elif self.change[0] == 'Add':
            for new_atom_name in self.change[1]:
                self.create_new_prot(new_atom_name)


def assess_pdb_line_format(pdb_in_lines): # Gives parsing scripts the correct positions for things like atom coordinates as pdb formatting can vary
    pdb_line_atom_num = 1
    pdb_line_atom_type = 2
    pdb_line_res_type = 3
    pdb_line_res_num = 0
    pdb_line_x_coord = 0
    pdb_line_element = 0
    pdb_line_chain = None
    for line in pdb_in_lines:
        if len(line) > 20:
            if line.startswith('ATOM') or line.startswith('HETATM'): # Assess based on the first line about an atom in a titratable residue in the pdb
                split_line = line.split()
                if split_line[pdb_line_res_type] in titrate_data.titr_amino_acids_3let:
                    if split_line[4] == 'A' or split_line[4] == 'B': # if the 5th column holds the chain rather than the residue number
                        pdb_line_chain = 4
                        pdb_line_res_num = 5
                        pdb_line_x_coord = 6
                        pdb_line_element = 11
                    else:
                        pdb_line_res_num = 4
                        pdb_line_x_coord = 5
                        pdb_line_element = 10
                    current_res = split_line[pdb_line_res_type] + split_line[pdb_line_res_num]
                    break

    return pdb_line_atom_num, pdb_line_atom_type, pdb_line_res_type, pdb_line_res_num, pdb_line_x_coord, pdb_line_element, pdb_line_chain

def process_pdb(pdb_lines):
    all_titr_res = [] # List of titratable residues
    current_res = 0 # Place holders for iterating
    current_res_allatoms = []
    current_res_name = ''
    hit_ter = True
    ter_res = [] # List of terminal residue names
    ter_type = [] # And whether they are C or N terminal

    # Assess pdb line format
    pdb_line_atom_num, pdb_line_atom_type, pdb_line_res_type, pdb_line_res_num, pdb_line_x_coord, pdb_line_element, pdb_line_chain = assess_pdb_line_format(pdb_lines)

    # Find all titratable amino acids in pdb and initialize class instances
    for i in range(len(pdb_lines)):
        line = pdb_lines[i]
        if i != 0: # Define the previous line to identify N-terminal residues
                prev_line = pdb_lines[i - 1]

        if len(line) > 20:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                split_line = line.split()
                if split_line[pdb_line_res_type] in titrate_data.titr_amino_acids_3let:
                    if split_line[pdb_line_res_num] != current_res:
                        if current_res != 0:
                            current_titr_res = titr_res(current_res_name, current_res, current_chain, current_res_allatoms)
                            current_titr_res.define_prot_state()
                            all_titr_res.append(current_titr_res)

                        current_res = split_line[pdb_line_res_num] # Reset variables for the next titratable residue
                        current_res_name = split_line[pdb_line_res_type]
                        if pdb_line_chain is not None:
                            current_chain = split_line[pdb_line_chain]
                        
                        else:
                            current_chain = ''
                        
                        current_res_allatoms = []

                    coord = np.array([float(split_line[pdb_line_x_coord]), float(split_line[pdb_line_x_coord + 1]), float(split_line[pdb_line_x_coord + 2])])
                    current_res_allatoms.append([split_line[pdb_line_atom_type], coord])

                if hit_ter:
                    if split_line[pdb_line_res_type] in titrate_data.amino_acids_3let:
                        if pdb_line_chain is not None:
                            ter_res.append(split_line[pdb_line_res_type] + split_line[pdb_line_res_num] + split_line[pdb_line_chain])
                        
                        else:
                            ter_res.append(split_line[pdb_line_res_type] + split_line[pdb_line_res_num])
                        
                        ter_type.append('N+')
                        hit_ter = False

        if len(line) >= 4:
            if line.startswith('TER'):
                split_line = prev_line.split()
                if split_line[pdb_line_res_type] in titrate_data.amino_acids_3let:
                    if pdb_line_chain is not None:
                        ter_res.append(split_line[pdb_line_res_type] + split_line[pdb_line_res_num] + split_line[pdb_line_chain])
                    
                    else:
                        ter_res.append(split_line[pdb_line_res_type] + split_line[pdb_line_res_num])
                    
                    ter_type.append('C-')
                    hit_ter = True

    # Enter the last titratable residue in
    current_titr_res = titr_res(current_res_name, current_res, current_chain, current_res_allatoms)
    current_titr_res.define_prot_state()
    all_titr_res.append(current_titr_res)

    current_res = 0

    # Loop back through the pdb to save the C and N terminal residues
    for line in pdb_lines:
        if len(line) > 20:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                split_line = line.split()
                if pdb_line_chain is not None:
                    current_res_full_name = split_line[pdb_line_res_type] + split_line[pdb_line_res_num] + split_line[pdb_line_chain]
                
                else:
                    current_res_full_name = split_line[pdb_line_res_type] + split_line[pdb_line_res_num]

                if current_res_full_name in ter_res:
                    if split_line[pdb_line_res_num] != current_res:
                        if current_res != 0:
                            current_titr_res = titr_res(current_res_name, current_res, current_chain, current_res_allatoms, ter=ter_type[0])
                            current_titr_res.define_prot_state()
                            all_titr_res += [current_titr_res] # Store terminal residues that are titratable redundantly
                            ter_res.pop(0) # Clear the 0th entry because they are found sequentially
                            ter_type.pop(0)

                        current_res = split_line[pdb_line_res_num] # Reset variables for the next titratable residue
                        current_res_name = split_line[pdb_line_res_type]
                        if pdb_line_chain is not None:
                            current_chain = split_line[pdb_line_chain]
                        
                        else:
                            current_chain = ''
                        
                        current_res_allatoms = []
                    
                    coord = np.array([float(split_line[pdb_line_x_coord]), float(split_line[pdb_line_x_coord + 1]), float(split_line[pdb_line_x_coord + 2])]) 
                    current_res_allatoms.append([split_line[pdb_line_atom_type], coord])

    # Enter the last terminal residue in
    current_titr_res = titr_res(current_res_name, current_res, current_chain, current_res_allatoms, ter=ter_type[0])
    current_titr_res.define_prot_state()
    all_titr_res.append(current_titr_res)

    return all_titr_res 

def define_connections(all_titr_res, int_cutoff):
    all_heteroatoms = [[res, heteroatom] for res in all_titr_res for heteroatom in res.heteroatoms]

    # Sort list of heteroatoms by x position
    all_heteroatoms.sort(key=lambda x: x[1][1][0])

    for index, atom in enumerate(all_heteroatoms[:-1]):
        next_atom = all_heteroatoms[index+1]

        j = index+1
        while next_atom[1][1][0] - atom[1][1][0] <= int_cutoff:
            if atom[0] is not next_atom[0]:
                if next_atom[0] not in atom[0].partners:
                    if np.linalg.norm(next_atom[1][1] - atom[1][1]) <= int_cutoff:
                        next_atom[0].partners.append(atom[0])
                        atom[0].partners.append(next_atom[0])
                        
                        if atom[0].amino_acid + ':' + atom[1][0] in titrate_data.cov_link_heteroatoms:
                            if next_atom[0].amino_acid + ':' + next_atom[1][0] in titrate_data.cov_link_heteroatoms:
                                atom[0].covalent_link = True
                                next_atom[0].covalent_link = True

            if next_atom is all_heteroatoms[-1]:
                break
            
            else:
                j += 1
                next_atom = all_heteroatoms[j]

def calc_pKa_total_pdb(propka_file, all_titr_res, chains):
    if not os.path.isfile(propka_file):
        logger.error(f"File does not exist: {propka_file}")
        raise FileNotFoundError(f"{propka_file}")
    
    calc_pKa_data = {}

    with open(propka_file, 'r') as propka_in:
        start = False
        for line in propka_in:
            line = line.split()
            if not start:
                if len(line) > 2:
                    if line[0] == 'Group' and line[1] == 'pKa':
                        start = True

            else:
                if len(line) > 4:
                    if line[0] == "Free" and line[1] == 'energy':
                        break
                
                    else:
                        if chains:
                            calc_pKa_data[line[0] + line[1] + line[2]] = float(line[3])

                        else:
                            calc_pKa_data[line[0] + line[1]] = float(line[3])
    
    return calc_pKa_data

def find_solv_shell(propka_file, chains):
    if not os.path.isfile(propka_file):
        logger.error(f"File does not exist: {propka_file}")
        raise FileNotFoundError(f"{propka_file}")
    
    # Dictionary that stores pKa by residue name consistent with the established titr_res class
    solv_data = {}
    with open(propka_file, 'r') as propka_file:
        start = False
        for line in propka_file:
            line = line.split()
            
            if not start:
                if len(line) > 9:
                    if line[0] == 'RESIDUE' and line[1] == 'pKa':
                        start = True

            else:
                if len(line) > 15 and line[5] == '%':
                    #store the buried % as a decimal
                    if chains:
                        solv_data[line[0] + line[1] + line[2]] = (float(line[4]) / 100)

                    else:
                        solv_data[line[0] + line[1]] = (float(line[4]) / 100)

                elif len(line) > 8:
                    if line[0] == "Coupled" and line[1] == 'residues':
                        break
    
    return solv_data

def define_aa_networks(all_titr_res):
    all_networks = []
    current_network = []
    already_in_network = []
    already_in_current_network = []

    while all_titr_res: # Empty stack of all residues to trace all networks
        if all_titr_res[0] not in already_in_network: # Don't consider amino acids already part of a network
            current_network.append(all_titr_res[0])
            current_partners = all_titr_res[0].partners
            already_in_current_network.append(all_titr_res[0])
        
        all_titr_res.pop(0)

        while current_partners: # If there are partners to the current residue, it is part of a large network
            if current_partners[0] not in already_in_current_network: # Takes care of loops, when multiple residues interact with eachother it add duplicates to the current_partner list
                current_network.append(current_partners[0])
                current_partners += current_partners[0].partners # Trace out any chain of partners
                already_in_network.append(current_partners[0])
                already_in_current_network.append(current_partners[0])

            current_partners.pop(0) # This partner has now been checked

        if current_network:
            all_networks.append(current_network)

        current_network = []
        already_in_current_network = []

    return all_networks

def find_network_solvent_access(old_networks, run_solv_data, solv_cutoff, solv_prob):
    new_networks = []

    for network in old_networks:
        for atom in network:
            if atom.full_name in run_solv_data:
                if run_solv_data[atom.full_name] < solv_cutoff:
                    new_networks.append(["Solv", network])
                    break

        else:
            new_networks.append(['NoSolv', network])

    return new_networks

def MC_prot_change(networks, solution_pH):
    for network in networks:
        if network[0] == 'Solv':
            for residue in network[1]:
                MC_decide_solv(residue, solution_pH)
        
        else:
            MC_decide_nosolv(network[1])


def MC_decide_solv(residue, resevoir_pH):
    prob_add = (10.0**(residue.pKa - resevoir_pH)) / (1.0 + 10.0**(residue.pKa - resevoir_pH))
    MC_prot_state_roll = float(random.randint(0, 1000000) / 1000000.0)

    if MC_prot_state_roll <= prob_add:
        if residue.prot_state[0] == '-': # If this is a change to protonation state, take proper action
            possible_prot_states = titrate_data.old_titr_form2new_titr_form[residue.ter_name + ':' + 'Add' + ':' + str(residue.prot_state[1])]
            MC_prot_form_roll = random.randint(1, len(possible_prot_states))
            residue.new_prot_state = possible_prot_states[MC_prot_form_roll - 1][0]
            residue.change = ['Add', possible_prot_states[MC_prot_form_roll - 1][1]]
            residue.change_heteroatom = titrate_data.hydrogen2boundheteroatom[residue.ter_name + ':' + possible_prot_states[MC_prot_form_roll - 1][1][0] + ':' + str(residue.prot_state[1])] # Selects the first hydrogen if two are added (which only affects N-terminus)
    else:
        if residue.prot_state[0] == '+': # If this is a change to protonation state, take proper action
            possible_prot_states = titrate_data.old_titr_form2new_titr_form[residue.ter_name + ':' + 'Remove' + ':' + str(residue.prot_state[1])]
            MC_prot_form_roll = random.randint(1, len(possible_prot_states))
            residue.new_prot_state = possible_prot_states[MC_prot_form_roll - 1][0]
            residue.change = ['Remove', possible_prot_states[MC_prot_form_roll - 1][1]]
            residue.change_heteroatom = titrate_data.hydrogen2boundheteroatom[residue.ter_name + ':' + possible_prot_states[MC_prot_form_roll - 1][1][0] + ':' + str(possible_prot_states[MC_prot_form_roll - 1][0][1])] # Selects the first hydrogen if two are added (which only affects N-terminus)


    residue.change_prob = prob_add
    residue.change_roll = MC_prot_state_roll

def MC_decide_nosolv(network):
    available_prots = 0

    for residue in network: # Tally the number of mobile protons in the system
        if residue.prot_state[0] == '+':
            available_prots += 1

    prot_combos = combinations(network, available_prots) # Enumerate all possible total protonation states of the system
    probabilities = []
    partition_func = 0.0
    for combo in prot_combos: # Construct the partition function for the system
        temp_prob = 0.0
        for residue in combo:
            temp_prob += residue.pKa

        temp_prob = 10.0**(temp_prob)
        probabilities.append([temp_prob, combo])
        partition_func += temp_prob

    MC_prot_state_roll = float(random.randint(0, 1000000) / 1000000.0) # Roll the dice and decide the state
    current_prob_total = 0.0
    unchanged_protonated_residues = []
    for state in probabilities:
        previous_prob_total = current_prob_total
        current_prob_total += state[0] / partition_func
        if MC_prot_state_roll < current_prob_total and MC_prot_state_roll > previous_prob_total: # If the MC roll falls within the current state probability range, choose this state
            for residue in state[1]:
                if residue.prot_state[0] == '-': # If this involves a protonation state change to this residue, update it
                    possible_prot_states = titrate_data.old_titr_form2new_titr_form[residue.ter_name + ':' + 'Add' + ':' + str(residue.prot_state[1])]
                    MC_prot_form_roll = random.randint(1, len(possible_prot_states))
                    residue.new_prot_state = possible_prot_states[MC_prot_form_roll - 1][0]
                    residue.change = ['Add', possible_prot_states[MC_prot_form_roll - 1][1]]
                    residue.change_heteroatom = titrate_data.hydrogen2boundheteroatom[residue.ter_name + ':' + possible_prot_states[MC_prot_form_roll - 1][1][0] + ':' + str(residue.prot_state[1])] # Selects the first hydrogen if two are added (which only affects N-terminus)
                
                if residue.prot_state[0] == '+': # If this is an unchanged residue, record that
                    unchanged_protonated_residues.append(residue)
                
                residue.change_prob = current_prob_total
                residue.change_roll = MC_prot_state_roll

    for residue in network: # Go back around and update the residues losing protons
        if residue.prot_state[0] == '+' and residue not in unchanged_protonated_residues: # Remove protons from newly deprotonated residues
            possible_prot_states = titrate_data.old_titr_form2new_titr_form[residue.ter_name + ':' + 'Remove' + ':' + str(residue.prot_state[1])]
            MC_prot_form_roll = random.randint(1, len(possible_prot_states))
            residue.new_prot_state = possible_prot_states[MC_prot_form_roll - 1][0]
            residue.change = ['Remove', possible_prot_states[MC_prot_form_roll - 1][1]]
            residue.change_heteroatom = titrate_data.hydrogen2boundheteroatom[residue.ter_name + ':' + possible_prot_states[MC_prot_form_roll - 1][1][0] + ':' + str(possible_prot_states[MC_prot_form_roll - 1][0][1])] # Selects the first hydrogen if two are added (which only affects N-terminus)

