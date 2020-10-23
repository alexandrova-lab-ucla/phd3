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
        'addH',
        'protein_to_coord',
        'coord_to_protein'
        ]

def add_proton(atom, ID="HY"):
    vectors = []
    for a in atom.bonds:
        vectors.append(a.coords - atom.coords)
        vectors[-1] = vectors[-1]/np.linalg.norm(vectors[-1])

    #If this is a carboynl, it adds the hydrogen linearly...not ideal, but should be good enough...
    direction = -1.1* sum(vectors)

    new_proton = pro.atom.Atom(element = "H", coords = atom.coords + direction, id=ID, number=20)
    
    #Update the bond lists etc...
    atom.add_bond(new_proton)
    atom.residue.add_atom(new_proton)
    
    return new_proton


def addH(protein):
    phd_config = utilities.load_phd_config()
    chimera_path = phd_config["PATHS"]["chimera"]

    protein.reformat_protein()
    protein.remove_h()
    protein.write_pdb("_temp.pdb", exclude_sub_chain=True)

    with open("chimeraaddh.com", "w") as comfile:
        comfile.write("open _temp.pdb\n")
        comfile.write("addh\n")
        comfile.write("write format pdb 0 addh.pdb\n")
        comfile.write("stop\n")

    try:
        with Popen(f"{os.path.join(chimera_path, 'chimera')} --nogui chimeraaddh.com",
                stdout=PIPE, stderr=subprocess.STDOUT, universal_newlines=True, shell=True, env=os.environ) as shell:
            while shell.poll() is None:
                logger.debug(shell.stdout.readline().strip())

    except OSError:
        logger.exception("Error calling chimera")
        raise

    if not os.path.isfile("addh.pdb"):
        logger.exception("Could not call chimera, check your path")
        raise OSError("chimera")

    logger.debug("Removing chimeraddh.com")
    os.remove("chimeraaddh.com")
    logger.debug("Removing _temp.pdb")
    os.remove("_temp.pdb")

    new_protein = utilities.load_pdb("addh.pdb")
    new_protein.chains.append(protein.sub_chain)

    new_protein.reformat_protein()
    logger.debug("Removing addh.pdb")
    os.remove("addh.pdb")

    #Remove all epsilon hydrogens on the histidines
    for chain in new_protein.chains:
        for res in chain.residues:
            if res.name.upper() == "HIS":
                epsilon_nitrogen = res.get_atom("NE2")
                deleted_hydrogen = False
                for atom in epsilon_nitrogen.bonds:
                    if atom.element.lower() == "h":
                        #DELETE THIS ATOM, no good two-timer
                        epsilon_nitrogen.bonds.remove(atom)
                        res.atoms.remove(atom)
                        del atom
                        deleted_hydrogen = True
                
                if not deleted_hydrogen:
                    logger.warn(f"Could not find his ({res}) 2HNE atom, searching again!")
                    for a in res.atoms:
                        if a.id.lower() == "2hne" or a.id.lower() == "he2":
                            #bastard got through
                            logger.warn("Found the hydrogen!")
                            for b in a.bonds:
                                b.bonds.remove(a)
                            
                            res.atoms.remove(a)
                            del a
                            break
                
                #Now we add a proton to the ND1 atom
                delta_nitrogen = res.get_atom("ND1")
                skip = False
                for a in res.atoms:
                    if a.id.lower() == "hd1":
                        skip = True

                for a in delta_nitrogen.bonds:
                    if a.element.lower() == "h":
                        skip = True

                if not skip:
                    logger.debug(f"Adding delta nitrogen to {res}")
                    add_proton(delta_nitrogen, ID="HD1")
               
    #For any added protons (fix the labeling and numbering)
    new_protein.reformat_protein()

    return new_protein

def protein_to_coord(initial_protein, chop_params):

    logger.debug("Protonating the protein")

    #Easiest way to create a copy of the protein I have right now:
    #TODO actually implement copy functions in the protein library
    initial_protein.write_pdb("_to_copy.pdb")
    protein = utilities.load_pdb("_to_copy.pdb")
    os.remove("_to_copy.pdb")

    protein = addH(protein)

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
                res = res.split(":")
                assert(len(res) == 2)

                chain = res[0]
                try:
                    res_num = int(res[1])
            
                except ValueError:
                    logger.error("Invalid specification of residue")
                    raise
                
                normal_residues.append(protein.get_residue([chain, res_num]))
                atoms.extend(protein.get_residue([chain, res_num]).atoms)

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
                
    #Now we recursively delete the atoms
    def remove_bonds_from_list(atom, residue=None):
        proceed = False
       
        for a in atom.bonds:
            if residue is None:
                if a not in remove_atoms:
                    proceed = True
                    break
        
            else:
                if a not in remove_atoms and a in residue.atoms:
                    proceed = True
                    break

        if proceed:
            if residue is None:
                for a in atom.bonds:
                    remove_atoms.append(a)
                    remove_bonds_from_list(a)

            else:
                while atom.bonds:
                    a = atom.bonds.pop()
                    if a.residue is residue:
                        if atom in a.bonds:
                            a.bonds.remove(atom)
                    
                        remove_atoms.append(a)
                        remove_bonds_from_list(a, residue)

    if "Substrate Chop" in chop_params.keys():
        atoms_to_remove = []
        #We make all of the chops first, then we will recursively add the atoms bonded to remove_atoms to the removeList
        for chop in chop_params["Substrate Chop"]:
            chop = chop.split("-")
            keepatom = chop[0].split(":")
            removeatom = chop[1].split(":")
            
            assert(len(keepatom) == 3)
            assert(len(removeatom) == 3)

            keepatom = protein.get_atom([keepatom[0], int(keepatom[1]), keepatom[2]])
            removeatom = protein.get_atom([removeatom[0], int(removeatom[1]), removeatom[2]])

            keepatom.bonds.remove(removeatom)
            removeatom.bonds.remove(keepatom)

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

                if "Exclude Sidechain" in chop_params.keys():
                    res_id = f"{chop[0][0]}:{chop[0][1]}"
                    if res_id in chop_params["Exclude Sidechain"]:
                        chop_params["Exclude Sidechain"].remove(res_id)

        for a in atoms_to_remove:
            if a.residue.name in constants.AMINO_ACID_RESIDUES:
                remove_bonds_from_list(a, a.residue)
            
            else:
                remove_bonds_from_list(a, a.residue)
            
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
            
                if res.get_atom("CB") in res.get_atom("CA").bonds:
                    res.get_atom("CA").bonds.remove(res.get_atom("CB"))
                    res.get_atom("CB").bonds.remove(res.get_atom("CA"))

                else:
                    logger.warn(f"CB not bonded to CA in residue {res}")

                #We do this so we don't remove the CA atom
                for b in res.get_atom("CA").bonds:
                    b.bonds.remove(res.get_atom("CA"))

                remove_bonds_from_list(res.get_atom("CA"), res)
                res.get_atom("CA").bonds = []

                chop_atoms.append([res.get_atom("CB"), res.get_atom("CA")])
        
        else:
            logger.error("Cannot specify non-amino acid residues in the 'Residue' section of the chop")
            logger.error("Skipping residue: {str(res)}")

    for linked_residues in linked_residues_end:
        nterm = linked_residues[0]
        
        if linked_residues[2] == 'c':
            for a in nterm.get_atom("N").bonds:
                if a.id == "C":
                    atoms.append(a)
                    chop_atoms.append([nterm.get_atom("N"), a])

            #If we don't find a "C", then we have an n-terminus...and therefore kinda pointless to cut...

        elif linked_residues[2] == 'n':
            #cut residue from other residues first
            for a in nterm.get_atom("N").bonds:
                if a.id == "C":
                    nterm.get_atom("N").bonds.remove(a)
                    a.bonds.remove(nterm.get_atom("N"))

            if nterm.get_atom("N") in nterm.get_atom("CA").bonds:
                nterm.get_atom("CA").bonds.remove(nterm.get_atom("N"))
                nterm.get_atom("N").bonds.remove(nterm.get_atom("CA"))
        
            else:
                logger.warn(f"N not in CA bonds in residue {nterm}")

            #ensures that we don't remove the N atom...
            for b in nterm.get_atom("N").bonds:
                b.bonds.remove(nterm.get_atom("N"))

            #Cut this stuff out
            remove_bonds_from_list(nterm.get_atom("N"))
            nterm.get_atom("N").bonds = []
            chop_atoms.append([nterm.get_atom("CA"), nterm.get_atom("N")])


        else:
            logger.error("Invalid cut specification for residue {str(nterm)}")
            raise ValueError("Cut specification for {str(nterm)}")

        cterm = linked_residues[1]

        if linked_residues[3] == 'c':
            #chop between CA and C to make CRH2
            for a in cterm.get_atom("C").bonds:
                if a.id == "N":
                    cterm.get_atom("C").bonds.remove(a)
                    a.bonds.remove(cterm.get_atom("C"))

            cterm.get_atom("C").bonds.remove(cterm.get_atom("CA"))
            cterm.get_atom("CA").bonds.remove(cterm.get_atom("C"))

            for b in cterm.get_atom("C").bonds:
                b.bonds.remove(cterm.get_atom("C"))

            remove_bonds_from_list(cterm.get_atom("C"))
            cterm.get_atom("C").bonds = []

            chop_atoms.append([cterm.get_atom("CA"), cterm.get_atom("C")])

        elif linked_residues[3] == 'n':
            #chop between c and n to make aldehyde
            for a in cterm.get_atom("C").bonds:
                if a.id == "N":
                    atoms.append(a)
                    chop_atoms.append([cterm.get_atom("C"), a])

        else:
            logger.error("Invalid cut specification for residue {str(cterm)}")
            raise ValueError(f"Cut specification for {str(cterm)}")

    if "Exclude Sidechain" in chop_params.keys():
        for sidechain in chop_params["Exclude Sidechain"]:
            sidechain = sidechain.split(":")
            assert(len(sidechain) == 2)
            chain = sidechain[0]
            resNum = int(sidechain[1])

            res = protein.get_residue([chain, resNum])

            if res.name == "GLY":
                logger.warn("Cannot exclude the sidechain of glycine!")
                continue

            res.get_atom("CA").bonds.remove(res.get_atom("CB"))
            res.get_atom("CB").bonds.remove(res.get_atom("CA"))

            for a in res.get_atom("CB").bonds:
                a.bonds.remove(res.get_atom("CB"))

            remove_bonds_from_list(res.get_atom("CB"))
            chop_atoms.append([res.get_atom("CA"), res.get_atom("CB")])

    if "Protonation" in chop_params.keys():
        for protonation_state in chop_params["Protonation"]:
            resID = protonation_state[0]
            state = protonation_state[1:]
  
            chain = resID.split(":")[0]
            resNum = int(resID.split(":")[1])

            res = protein.get_residue([chain, resNum])

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
               
                proton = add_proton(atom_to_protonate)
                atoms.append(proton)

    if "Freeze Atoms" in chop_params.keys():
        for a in chop_params["Freeze Atoms"]:
            a = a.split(":")
            chain = a[0]
            res_num = int(a[1])
            atom_id = a[2]
            protein.get_atom([chain, res_num, atom_id]).freeze = True

    if "Dummy H" in chop_params.keys():
        for a in chop_params["Dummy H"]:
            a = a.split(":")
            chain = a[0]
            res_num = int(a[1])
            atom_id = a[2]
            remove_atoms.append(protein.get_atom([chain, res_num, atom_id]))

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

    protein = addH(initial_protein)

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
        label = label.split(":")
        chain = label[0]
        resNum = int(label[1])
        atomID = label[2]

        #HX are dummy from chops, HY are protons for protonation
        if atomID == "HX" or atomID == "HY":
            continue

        try:
            if protein.get_atom([chain, resNum, atomID]).element != "Zn" and protein.get_atom([chain, resNum, atomID]).element.lower() == atom[1].lower():
                protein.get_atom([chain, resNum, atomID]).coords = atom[0] / constants.A_TO_BOHR
           
            elif protein.get_atom([chain, resNum, atomID]).id.lower() == atom[1].lower():
                protein.get_atom([chain, resNum, atomID]).coords = atom[0] / constants.A_TO_BOHR

            else:
                logger.error("Label file is out of order from coord file!")
                logger.error(f"Chain: {chain}, Res: {resNum}, atom: {atomID}")
                raise ValueError("Label file")

        except ValueError:
            logger.error("Cannot find atoms in the protein!")
            raise
   
                          
    if "Dummy H" in chop_params.keys():
        for atom in chop_params["Dummy H"]:
            atom = atom.split(":")
            atom = protein.get_atom([atom[0], int(atom[1]), atom[2]])
            if atom.bonds:
                b = atom.bonds[0]
                d = b.coords - atom.coords
                distance = np.linalg.norm(d)
                d = d/distance
                if distance > 1.5:
                    atom.coords = b.coords - (1.1*d)

                

    #Now we are tech. done!
    return protein

