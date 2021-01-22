#!/usr/bin/env python3
"""
Author  ==>> Matthew R. Hennefarth
Date    ==>> April 16, 2020
"""

#Standard Library Imports
import logging
import subprocess
import os
import numpy as np
from subprocess import Popen, PIPE
import sys
sys.setrecursionlimit(10000)

#PHD3 Imports
import phd3.utility.utilities as utilities 
import phd3.utility.constants as constants
import phd3.protein as pro


logger = logging.getLogger(__name__)

_all__ = [
        'protein_to_coord',
        'coord_to_protein'
        ]

#Now we recursively delete the atoms
def remove_bonds_from_list(atom, remove_list, residue=None):
    proceed = False
   
    for a in atom.bonds:
        if residue is None:
            if a not in remove_list:
                proceed = True
                break
    
        else:
            if a not in remove_list and a in residue.atoms:
                proceed = True
                break

    if proceed:
        if residue is None:
            for a in atom.bonds:
                remove_list.append(a)
                remove_bonds_from_list(a, remove_list)

        else:
            while atom.bonds:
                a = atom.bonds.pop()
                if a.residue is residue:
                    if atom in a.bonds:
                        a.bonds.remove(atom)
                
                    remove_list.append(a)
                    remove_bonds_from_list(a, remove_list, residue)


def add_to_cut_list(atom_to_keep, atom_to_replace, cut_list, remove_list):
    try:
        atom_to_keep.bonds.remove(atom_to_replace)

    except ValueError:
        logger.warn(f"{atom_to_replace} not in {atom_to_keep} bonds")

    try:
        atom_to_replace.bonds.remove(atom_to_keep)

    except ValueError:
        logger.warn(f"{atom_to_keep} not in {atom_to_replace} bonds")

    # prevents the recursion from deleting our "frozen" atom
    for b in atom_to_replace.bonds:
        try:
            b.bonds.remove(atom_to_replace)
        
        except ValueError:
            logger.warn(f"{atom_replace} not in {b} bonds")
    
    remove_bonds_from_list(atom_to_replace, remove_list, atom_to_replace.residue)
    atom_to_replace.bonds = []

    cut_list.append([atom_to_keep, atom_to_replace])


def protein_to_coord(initial_protein, chop_params):

    logger.debug("Protonating the protein")

    #Easiest way to create a copy of the protein I have right now:
    #TODO actually implement copy functions in the protein library
    initial_protein.write_pdb("_to_copy.pdb")
    protein = utilities.load_pdb("_to_copy.pdb")
    os.remove("_to_copy.pdb")

    protein = utilities.addH(protein)

    #We fill with all possible atoms, and then remove what is not in the QM
    atoms = []

    #These will be chopped normally
    normal_residues = []

    #These are at the end of a list of continuous residues
    linked_residues_end = []

    logger.debug("Including all substrate/metals")
    for res in protein.sub_chain.residues:
        atoms.extend(res.atoms)

    if "Residues" in chop_params.keys():
        for res in chop_params["Residues"]:
            if "-" in res:
                #We are dealing with a range of residues
                try:
                    res1, res2 = res.split("-")
                
                except:
                    logger.error("Invalid specification of multiple residues")
                    raise
                
                res1 = res1.split(":")
                res2 = res2.split(":")
                
                #Specify the cut type
                cutleft = 'n'
                cutright = 'n'

                #Cannot pop the last thing...oh well
                if res1[-1][-1].isalpha():
                    cutleft = res1[-1][-1].lower()
                    res1[-1] = res1[-1][:-1]

                if res2[-1][-1].isalpha():
                    cutright = res2[-1][-1]
                    res2[-1] = res2[-1][:-1]

                chain1 = res1[0]
                res_num_1 = int(res1[1])
                
                res_num_2 = int(res2[1])
                
                if res_num_1 > res_num_2:
                    t = res_num_2
                    res_num_2 = res_num_1
                    res_num_1 = t
                    del t

                res_num = res_num_1
                while res_num <= res_num_2:
                    atoms.extend(protein.get_residue([chain1, res_num]).atoms)
                    res_num += 1

                #can later include a specification of what to include here
                linked_residues_end.append([protein.get_residue([chain1, res_num_1]), protein.get_residue([chain1, res_num_2]), cutleft.lower(), cutright.lower()])


            else:
                #Normal residue specification
                normal_residues.append(protein.get_residue(res))
                atoms.extend(protein.get_residue(res).atoms)

    #Strictly removed
    remove_atoms = []
    #[a1 -> freeze, a2->change to hydrogen and freeze]
    chop_atoms = []

    # Gets rid of any HDUM atoms (like for O2 and HO- ligands...)
    for a in atoms:
        if a.id.upper() == "HDUM":
            remove_atoms.append(a)

    #For something like O2, can just place the hydrogens in the "Exclude Atoms" for the QM region...
    if "Exclude Atoms" in chop_params.keys():
        for exclude in chop_params["Exclude Atoms"]:
            exclude = exclude.split(":")
            assert(len(exclude) == 3)
            remove_atoms.append(protein.get_atom([exclude[0], int(exclude[1]), exclude[2]]))
                
    if "Substrate Chop" in chop_params.keys():
        atoms_to_remove = []
        #We make all of the chops first, then we will recursively add the atoms bonded to remove_atoms to the removeList
        for chop in chop_params["Substrate Chop"]:
            chop = chop.split("-")
            
            # Actually retrives atom objects
            keepatom = protein.get_atom(chop[0])
            removeatom = protein.get_atom(chop[1])

            try:
                keepatom.bonds.remove(removeatom)
            
            except:
                pass

            try:
                removeatom.bonds.remove(keepatom)

            except:
                pass

            #This way we don't premptively remove the removeatom (which will be replaced by a hydrogen)
            for b in removeatom.bonds:
                b.bonds.remove(removeatom)
                
            atoms_to_remove.append(removeatom)
            chop_atoms.append([keepatom, removeatom])

            # Can now quickly specify a different type of cut for residues
            # Though, this won't work for including sequential residues (linked)
            if keepatom.residue.name in constants.AMINO_ACID_RESIDUES:
                # First we cut n and alpha carbons
                for a in keepatom.residue.get_atom("N").bonds:
                    if a.id == "C":
                        keepatom.residue.get_atom("N").bonds.remove(a)

                for a in keepatom.residue.get_atom("C").bonds:
                    if a.id == "N":
                        keepatom.residue.get_atom("C").bonds.remove(a)

                # Then we remove from normal residue list
                if keepatom.residue in normal_residues:
                    normal_residues.remove(keepatom.residue)

                if "Exclude Side Chain" in chop_params.keys():
                    res_id = f"{chop[0][0]}:{chop[0][1]}"
                    if res_id in chop_params["Exclude Side Chain"]:
                        chop_params["Exclude Side Chain"].remove(res_id)

        for a in atoms_to_remove:
            if a.residue.name in constants.AMINO_ACID_RESIDUES:
                remove_bonds_from_list(a, remove_atoms, a.residue)
            
            else:
                remove_bonds_from_list(a, remove_atoms, a.residue)
            
            a.bonds = []

    for res in normal_residues:
        if res.name in constants.AMINO_ACID_RESIDUES:
            #Get rid of all backbone atoms...
            if res.name == "GLY":
                logger.warn("Tried to include glycine as a normal residue in the QM region")
                remove_atoms.extend(res.atoms)
                
            else:
                if res.name == "PRO":
                    logger.warn("Trying to include proline as a normal residue in QM region")
                    logger.warn("Issues/exceptions/infinite recursions can occur")
                        
                #We chop the residue from the rest of protein
                for a in res.get_atom("N").bonds:
                    if a.id == "C":
                        res.get_atom("N").bonds.remove(a)

                for a in res.get_atom("C").bonds:
                    if a.id == "N":
                        res.get_atom("C").bonds.remove(a)
            
                add_to_cut_list(res.get_atom("CB"), res.get_atom("CA"), chop_atoms, remove_atoms)
        
        else:
            logger.error("Cannot specify non-amino acid residues in the 'Residue' section of the chop")
            logger.error("Skipping residue: {str(res)}")

    for linked_residues in linked_residues_end:
        nterm = linked_residues[0]
       
        n_atom = nterm.get_atom("N")
        if linked_residues[2] == 'c':
            for a in n_atom.bonds:
                if a.id == "C":
                    atoms.append(a)
                    # I don't know why it gets added to this list, but there must be a good reason...
                    if a in remove_atoms:
                        remove_atoms.remove(a)
                    chop_atoms.append([n_atom, a])
                    break

            #If we don't find a "C", then we have an n-terminus...and therefore kinda pointless to cut...
            else:
                logger.warn(f"Could not cut {linked_residues[0]} properly, could be N-terminus")

        elif linked_residues[2] == 'a' or linked_residues[2] == 'n':
            for a in n_atom.bonds:
                if a.id == 'C':
                    n_atom.bonds.remove(a)
                    try:
                        a.bonds.remove(n_atom)

                    except ValueError:
                        logger.warn(f"{n_atom} not in bond list for {a}")

            if linked_residues[2] == 'a':
                add_to_cut_list(nterm.get_atom("C"), nterm.get_atom("CA"), chop_atoms, remove_atoms)
                if "Exclude Side Chain" in chop_params.keys():
                    if nterm.label() in chop_params["Exclude Side Chain"]:
                        chop_params["Exclude Side Chain"].remove(nterm.label())

            else:
                add_to_cut_list(nterm.get_atom("CA"), nterm.get_atom("N"), chop_atoms, remove_atoms)

        else:
            logger.error("Invalid cut specification for residue {str(nterm)}")
            raise ValueError("Cut specification for {str(nterm)}")

        cterm = linked_residues[1]
        c_atom = cterm.get_atom("C")

        if linked_residues[3] == 'n':
            #chop between c and n to make aldehyde
            for a in cterm.get_atom("C").bonds:
                if a.id == "N":
                    atoms.append(a)
                    try:
                        chop_atoms.append([cterm.get_atom("C"), a])
                    
                    except ValueError:
                        logger.warn(f"{c_atom} not in bond list for {a}")

                    break

            else:
                logger.warn(f"Could not cut {linked_residues[1]} properly, could be a C-terminus")

        elif linked_residues[3] == 'c' or linked_residues[3] == 'a':
            # Cut residues from rest of protein
            for a in c_atom.bonds:
                if a.id == "N":
                    c_atom.bonds.remove(a)
                    try:
                        a.bonds.remove(c_atom)

                    except ValueError:
                        logger.warn(f"{c_atom} not in bond list for {a}")

            if linked_residues[3] == 'c':
                add_to_cut_list(cterm.get_atom("CA"), cterm.get_atom("C"), chop_atoms, remove_atoms)

            else:
                add_to_cut_list(cterm.get_atom("N"), cterm.get_atom("CA"), chop_atoms, remove_atoms)
                if "Exclude Side Chain" in chop_params.keys():
                    if cterm.label() in chop_params["Exclude Side Chain"]:
                        chop_params["Exclude Side Chain"].remove(cterm.label())

        else:
            logger.error("Invalid cut specification for residue {str(cterm)}")
            raise ValueError(f"Cut specification for {str(cterm)}")

    if "Exclude Side Chain" in chop_params.keys():
        for sidechain in chop_params["Exclude Side Chain"]:
            res = protein.get_residue(sidechain)

            if res.name == "GLY":
                logger.warn("Cannot exclude the sidechain of glycine!")
                continue
            
            add_to_cut_list(res.get_atom("CA"), res.get_atom("CB"), chop_atoms, remove_atoms)

    if "Protonation" in chop_params.keys():
        for protonation_state in chop_params["Protonation"]:
            state = protonation_state[1:]

            res = protein.get_residue(protonation_state[0])

            if len(state) > 1:
                if state[0] == "protonate":
                    heavy_atom, hydrogen_atom = constants.PROTONATED[res.name][state[1]]
                
                elif state[0] == "deprotonate":
                    heavy_atom, hydrogen_atom = constants.DEPROTONATED[res.name][state[1]]

                else:
                    logger.error("Invalid protonate/deprotonate specified")
                    raise ValueError("protonate/deprotonate")

            else:
                if state[0] == "protonate":
                    heavy_atom, hydrogen_atom = constants.PROTONATED[res.name][0]
                
                elif state[0] == "deprotonate":
                    heavy_atom, hydrogen_atom = constants.DEPROTONATED[res.name][0]

                else:
                    logger.error("Invalid protonate/deprotonate specified")
                    raise ValueError("protonate/deprotonate")
            
            if state[0] == "deprotonate":
                for atom in res.atoms:
                    if atom.id == hydrogen_atom:
                        remove_atoms.append(atom)
                        break 
                
            if state[0] == "protonate":
                #Need to actually compute some angles and place in an atom...
                atom_to_protonate = res.get_atom(heavy_atom)
               
                proton = utilities.add_proton(atom_to_protonate)
                atoms.append(proton)

    if "Freeze Atoms" in chop_params.keys():
        for a in chop_params["Freeze Atoms"]:
            protein.get_atom(a).freeze = True

    if "Dummy H" in chop_params.keys():
        for a in chop_params["Dummy H"]:
            remove_atoms.append(protein.get_atom(a))

    #Make sure that the chops are not connected to each other
    for chop in chop_atoms:
        #The first one just gets frozen
        #The second one becomes hydrogen, and frozen
        chop[0].freeze = True
        chop[1].freeze = True
        chop[1].element = 'h'
        chop[1].id = "HX"

        #Now we compute the distance between these two
        dis = chop[1].coords - chop[0].coords
        dis = dis / np.linalg.norm(dis)
        dis = dis * 1.1
        chop[1].coords = chop[0].coords + dis

    for a in remove_atoms:
        if a in atoms:
            atoms.remove(a)
    
    with open("coord", 'w') as final:
        final.write("$coord\n")
        for a in atoms:
            final.write(a.coord_line())
        final.write("$end\n")

    #This holds the information to go backwards after all of this!
    with open("label", 'w') as lb:
        for a in atoms:
            lb.write(a.label() + '\n')


def coord_to_protein(initial_protein, chop_params):
    #Will update coords from the coord file in the protein passed
    #Needs to have a label and a coord file present in directory
    if not os.path.isfile("coord"):
        logger.error("No coord file found")
        raise FileNotFoundError("coord")

    if not os.path.isfile("label"):
        logger.error("No label file found")
        raise FileNotFoundError("label")

    protein = utilities.addH(initial_protein)

    coord_lines = []
    labels = []
    with open("coord", 'r') as coord:
        for line in coord:
            if line.startswith("$coord"):
                continue
        
            elif line[0] == '$':
                break

            coord_lines.append(line.rstrip())

    with open("label", 'r') as lblfile:
        for line in lblfile:
            labels.append(line.rstrip())

    assert(len(coord_lines) == len(labels))

    coord_lines = [i.split() for i in coord_lines]
    coord_lines = [[np.array( [float(i[0]), float(i[1]), float(i[2])]), i[3]] for i in coord_lines if len(i) >= 4]

    #coord_lines is (coords, element)
    for atom, label in zip(coord_lines, labels):
        atomID = label.split(":")[2]

        #HX are dummy from chops, HY are protons for protonation
        if atomID == "HX" or atomID == "HY":
            continue

        try:
            if protein.get_atom(label).element != "Zn" and protein.get_atom(label).element.lower() == atom[1].lower():
                protein.get_atom(label).coords = atom[0] / constants.A_TO_BOHR
           
            elif protein.get_atom(label).id.lower() == atom[1].lower():
                protein.get_atom(label).coords = atom[0] / constants.A_TO_BOHR

            else:
                logger.error("Label file is out of order from coord file!")
                logger.error(label)
                raise ValueError("Label file")

        except ValueError:
            logger.error("Cannot find atoms in the protein!")
            raise
   
                          
    if "Dummy H" in chop_params.keys():
        for atom in chop_params["Dummy H"]:
            atom = protein.get_atom(atom)
            if atom.bonds:
                b = atom.bonds[0]
                d = b.coords - atom.coords
                distance = np.linalg.norm(d)
                d = d/distance
                if distance > 1.5:
                    atom.coords = b.coords - (1.1*d)

                

    #Now we are tech. done!
    return protein

