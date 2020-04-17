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

#PHD3 Imports
import phd3.utility.utilities as utilities 
import phd3.utility.constants as constants
import phd3.protein

logger = logging.getLogger(__name__)

_all__ = [
        'addH',
        'protein_to_coord',
        'coord_to_protein'
        ]

def addH(protein):
    protein.write_pdb("_temp.pdb")

    with open("chimeraaddh.com", "w") as comfile:
        comfile.write("open _temp.pdb\n")
        comfile.write("addh\n")
        comfile.write("write format pdb 0 addh.pdb\n")
        comfile.write("stop\n")

    try:
        with Popen("chimera --nogui chimeraaddh.com",
                stdout=PIPE, stderr=subprocess.STDOUT, universal_newlines=True, shell=True, env=os.environ) as shell:
            while shell.poll() is None:
                logger.debug(shell.stdout.readline().strip())

    except OSError:
        logger.exception("Error calling chimera")
        raise
   
    logger.debug("Removing chimeraddh.com")
    os.remove("chimeraaddh.com")
    logger.debug("Removing _temp.pdb")
    os.remove("_temp.pdb")

    pro = utilities.load_pdb("addh.pdb")

    pro.reformat_protein()
    logger.debug("Removing addh.pdb")
    os.remove("addh.pdb")

    return pro

def protein_to_coord(protein, chop_params):

    logger.debug("Protonating the protein")
    
    #TODO can instead check substrates and exclude, except if we need to add a hydrogen for waters...
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

                if res1[-1].isalpha():
                    cutleft = res1.pop().lower()

                if res2[-1].isalpha():
                    cutright = res2.pop().lower()

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
                linked_residues_end.append([protein.get_residue([chain1, res_num_1]), protein.get_residue([chain1, res_num_2]), cutleft, cutright])


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

    if "Exclude Atoms" in chop_params.keys():
        for exclude in chop_params["Exclude Atoms"]:
            exclude = exclude.split(":")
            assert(len(exclude) == 3)
            remove_atoms.append(protein.get_atom([exclude[0], int(exclude[1], exclude[2])]))

    #Now we recursively delete the atoms
    def remove_bonds_from_list(atom):
        proceed = False

        for a in atom.bonds:
            if a not in remove_atoms:
                proceed = True
                break

        if proceed:
            for a in atom.bonds:
                remove_atoms.append(a)
                remove_bonds_from_list(a)


    for a in remove_atoms:
        atoms.remove(a)
    
    if "Substrate Chop" in chop_params.keys():
        
        remove_atoms = []
        #We make all of the chops first, then we will recursively add the atoms bonded to remove_atoms to the removeList
        for chop in chop_params["Substrate Chop"]:
            chop = chop.split("-")
            keepatom = chop[0].split(":")
            removeatom = chop[1].split(":")
            
            assert(len(keepatom) == 3)
            assert(len(removeatom) == 3)

            keepatom = protein.get_atom(keepatom[0], int(keepatom[1]), keepatom[2])
            removeatom = protein.get_atom(removeatom[0], int(removeatom[1]), removeatom[2])

            keepatom.bonds.remove(removeatom)
            removeatom.bonds.remove(keepatom)

            #This way we don't premptively remove the removeatom (which will be replaced by a hydrogen)
            for b in removeatom.bonds:
                b.bonds.remove(removeatom)
                
            remove_atoms.append(removeatom)
            chop_atoms.append(keepatom, removeatom)

        for a in remove_atoms:
            remove_bonds_from_list(a)
            a.bonds = []

        del remove_atoms

    for res in normal_residues:
        if res.name in constants.AMINO_ACID_RESIDUES:
            #Get rid of all backbone atoms...
            if res.name == "GLY":
                logger.warn("Tried to include glycine as a normal residue in the QM region")
                remove_atoms.extend(res.atoms)
            
            else:
                #We chop the residue from the rest of protein
                for a in res.get_atom("N").bonds:
                    if a.id == "C":
                        res.get_atom("N").bonds.remove(a)

                for a in res.get_atom("C").bonds:
                    if a.id == "N":
                        res.get_atom("C").bonds.remove(a)

                res.get_atom("CA").bonds.remove(res.get_atom("CB"))
                res.get_atom("CB").bonds.remove(res.get_atom("CA"))

                #We do this so we don't remove the CA atom
                for b in res.get_atom("CA").bonds:
                    b.bonds.remove(res.get_atom("CA"))

                remove_bonds_from_list(res.get_atom("CA"))
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
                    chop_atoms.append(nterm.get_atom("N"), a)

            #If we don't find a "C", then we have an n-terminus...and therefore kinda pointless to cut...

        elif linked_residues[2] == 'n':
            #cut residue from other residues first
            for a in nterm.get_atom("N").bonds:
                if a.id == "C":
                    nterm.get_atom("N").bonds.remove(a)
                    a.bonds.remove(nterm.get_atom("N"))

            nterm.get_atom("CA").bonds.remove(nterm.get_atom("N"))
            nterm.get_atom("N").bonds.remove(nterm.get_atom("CA"))

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

            remove_bronds_from_list(cterm.get_atom("C"))
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
            raise ValueError("Cut specification for {str(cterm)}")

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

            remove_bronds_from_list(res.get_atom("CB"))
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
                remove.append(res.get_atom(hydrogen_atom))

            if state[0] == "protonate":
                #Need to actually compute some angles and place in an atom...
                atom_to_protonate = res.get_atom(heavy_atom)
                
                vectors = []
                for  a in atom_to_protonate.bonds:
                    vectors.append(a.coords - atom_to_protonate.coords)
                    vectors[-1] = vectors[-1]/np.linalg.norm(vectors[-1])

                #If this is a carboynl, it adds the hydrogen linearly...not ideal, but should be good enough...
                direction = -1.1* sum(vectors)

                #Label them with HX so that when we load back in, it does so normally (ie ignores these hydrogens because they shouldn't be in the protein!)
                atoms.append(protein.atom.Atom(element = "H", coords = atom_to_protonate.coords + direction, id='HY'))

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


def coord_to_protein(protein):
    #Will update coords from the coord file in the protein passed
    #Needs to have a label and a coord file present in directory
    if not os.path.isfile("coord"):
        logger.error("No coord file found")
        raise FileNotFoundError("coord")

    if not os.path.isfile("label"):
        logger.error("No label file found")
        raise FileNotFoundError("label")

    protein = addH(protein)

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
            if protein.get_atom([chain, resNum, atomID]).element == atom[1]:
                protein.get_atom([chain, resNum, atomID]).coords = atom[0] / constants.A_TO_BOHR
                
            else:
                logger.error("Label file is out of order from coord file!")
                raise ValueError("Label file")

        except ValueError:
            logger.error("Cannot find atoms in the protein!")
            raise
    
    #Now we are tech. done!
    return protein

