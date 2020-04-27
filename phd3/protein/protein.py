#!/usr/bin/env python3
"""
Author  ==>> Matthew R. Hennefarth
Date    ==>> April 16, 2020
"""

#Standard Library Imports
import logging
import pkg_resources
import csv
import os
import numpy as np
from subprocess import Popen, PIPE

#PHD3 Imports
from ..utility import constants
from . import chain, residue

__all__=[
    'Protein'
]

class Protein:

    __slots__ = ['chains', '_logger', 'non_residues', 'metals', 'name', 'sub_chain', 'coords']

    def __init__(self, name: str, chains: [chain.Chain]):

        self._logger = logging.getLogger(__name__)

        self._logger.debug(f"Initializing vars for protein {name}")
        self.name = name
        self.chains = chains
        self.non_residues = []
        self.metals = []
        self.sub_chain = chain.Chain()
        self.coords = None

        self._logger.debug(f"Created protein {str(self)}")

    def reformat_protein(self, relabel_protein=True):
        # This is the BIG BIG BIG function that fixes EVERYTHING of a pdb for DMD
        # Don't question why it does things, it needs to
        res_renum = 1
        chain_let = 'A'

        chain_num = 0
        self._logger.debug("Looping over the chains")
        while chain_num < len(self.chains):
            atom_renum = 1
            residue_num = 0
            while residue_num < len(self.chains[chain_num].residues):
                # Case where we have a non-residue/metal, most likely substrate
                if self.chains[chain_num].residues[residue_num].name not in constants.AMINO_ACID_RESIDUES:
                    atom_num = 0
                    while atom_num < len(self.chains[chain_num].residues[residue_num].atoms):
                        # Remove metal from this, and make it its own residue essentially
                        if self.chains[chain_num].residues[residue_num].atoms[atom_num].element.lower() in constants.METALS:
                            self._logger.debug(f"Found a metal: {self.chains[chain_num].residues[residue_num].atoms[atom_num]} in residue {self.chains[chain_num].residues[residue_num]}")
                            
                            if self.chains[chain_num].residues[residue_num].atoms[atom_num].element.lower() == "zn":
                                if self.chains[chain_num].residues[residue_num].atoms[atom_num].id.lower() in constants.METALS:
                                    self.chains[chain_num].residues[residue_num].atoms[atom_num].element = self.chains[chain_num].residues[residue_num].atoms[atom_num].id.capitalize()

                            self.metals.append(self.chains[chain_num].residues[residue_num].atoms.pop(atom_num))
                            atom_num -= 1

                        atom_num += 1

                    self._logger.debug(f"Removing residue: {self.chains[chain_num].residues[residue_num]}")
                    if self.chains[chain_num].residues[residue_num].atoms:
                        self.non_residues.append(self.chains[chain_num].residues.pop(residue_num))

                    else:
                        del self.chains[chain_num].residues[residue_num]

                else:
                    # update residue, and atomic numbering for normal atoms
                    self._logger.debug(f"Renumbering residue: {self.chains[chain_num].residues[residue_num]} to {res_renum} and its atoms starting at {atom_renum}")
                    self.chains[chain_num].residues[residue_num].number = res_renum
                    self.chains[chain_num].residues[residue_num].inConstr_number = residue_num + 1
                    for atom in self.chains[chain_num].residues[residue_num].atoms:
                        atom.number = atom_renum
                        atom_renum += 1

                    res_renum += 1
                    residue_num += 1
                    self._logger.debug("Finished renumbering")

            if not self.chains[chain_num].residues:
                self._logger.debug(f"Removing chain: {self.chains[chain_num]}")
                del self.chains[chain_num]

            else:
                self.chains[chain_num].name = chain_let
                chain_let = chr(ord(chain_let) + 1)
                chain_num += 1

        #Now we add in the non-residues and the metals into a new chain!
        res_num = 1
        if self.non_residues or self.metals:
            self._logger.debug("Rearranging the non-residue substrates and metals")
            self._logger.debug("Creating a substrate/metal chain")
            atom_num = 1

            # in the event that we have only substrate
            if not len(self.chains):
                self.sub_chain.name = 'A'

            else:
                self.sub_chain.name = chr(ord(self.chains[-1].name) + 1)

            # Want to sort the metals!
            self.metals.sort(key=lambda metal: metal.element)

            cur_metal = ""
            metal_num = 1
            for metal in self.metals:
                metal.number = atom_num
                if cur_metal != metal.element:
                    metal_num = 1
                    cur_metal = metal.element

                if len(metal.element) == 1:
                    name = metal.element.upper() + f"{metal_num:02d}"

                elif len(metal.element) == 2:
                    name = metal.element.upper() + f"{metal_num:01d}"

                else:
                    self._logger.error(f"Encountered a metal with an unusual element ID: {metal}")
                    raise ValueError

                if len(name) != 3:
                    if metal_num > 35:
                        self._logger.error(f"The name for this is too long: {metal} with name: {name}")
                        raise ValueError

                    name = metal.element.upper() + f"{chr(metal_num-10 + ord('A'))}"

                # DMD does not know how to handle any other metal but zinc
                if metal.element.lower() != "zn":
                    metal.element = 'Zn'

                metal_residue = residue.Residue(name=name, number=res_num)
                self._logger.debug(f"Adding metal {metal_residue} as a residue to substrate chain")
                self.sub_chain.add_residue(metal_residue)
                metal_residue.add_atom(metal)

                res_num += 1
                metal_num += 1
                atom_num += 1

            # Issue with assigning chain to atom
            for res in self.non_residues:
                start = 100
                for atom in res.atoms:
                    atom.id = f"{atom.element.upper()[0]}{start}"
                    start += 1

                self._logger.debug(f"Adding residue {res} to substrate chain")
                self.sub_chain.add_residue(res)
                res.number = res_num
                res.inConstr_number = res_num
                for atom in res.atoms:
                    atom.number = atom_num
                    atom_num += 1

                res_num += 1

            self._logger.debug("Adding substrate chain to master chain")
            self.chains.append(self.sub_chain)

        if relabel_protein:
            self._logger.debug("Relabeling the protein")
            self.relabel()
        
        else: #relabel calls make_bond_table already
            self._logger.debug("Making the bond table for the protein")
            self.make_bond_table()

    def get_atom(self, identifier):
        for chain in self.chains:
            if chain.name == identifier[0]:
                for residue in chain.residues:
                    if residue.number == identifier[1]:
                        for atom in residue.atoms:
                            if atom.id == identifier[2]:
                                return atom

        self._logger.error(f"Could not find requested atom {identifier}")
        raise ValueError

    def get_residue(self, identifier):
        for chain in self.chains:
            if chain.name == identifier[0]:
                for residue in chain.residues:
                    if residue.number == identifier[1]:
                        return residue

        self._logger.error("Could not find requested residue")
        raise ValueError

    def get_chain(self, identifier):
        for chain in self.chains:
            if chain.name == identifier:
                return chain

        self._logger.error("Could not find the requested chain")
        raise ValueError

    def write_pdb(self, name=None):
        if name is None:
            name = self.name

        self._logger.debug(f"Writing out pdb: {name}")
        try:
            with open(name, 'w') as pdb:
                for chain in self.chains[:-1]:
                    for residue in chain.residues:
                        for atom in residue.atoms:
                            pdb.write(atom.pdb_line())
                    pdb.write('TER\n')

                if not self.sub_chain.residues:
                    for residue in self.chains[-1].residues:
                        for atom in residue.atoms:
                            pdb.write(atom.pdb_line())
                    pdb.write('TER\n')

                else:
                    for residue in self.sub_chain.residues:
                        for atom in residue.atoms:
                            pdb.write(atom.pdb_line())
                        pdb.write('TER\n')

                pdb.write('ENDMDL\n')

        except IOError:
            self._logger.exception(f"Error writing out to file {self.name}")
            raise

    def relabel(self, format: str="DMD"):

        #Need to make the bond table
        self.make_bond_table()

        atom_label_dict = {}
        with open(pkg_resources.resource_filename('phd3.resources', 'atom_label.csv')) as csvfile:
            csvreader = csv.reader(csvfile)
            schemenames = next(csvreader)[1:]
            try:
                newid = schemenames.index(format)

            except ValueError:
                raise ValueError("Format key not found in atom_label.csv")

            for row in csvreader:
                if row[0] in atom_label_dict:
                    atom_label_dict[row[0]].append(row[1:])

                else:
                    atom_label_dict[row[0]] = [row[1:]]

        self._logger.debug("Loaded in the atom_label.csv file")

        def rename_residue(residue, terminus: str = None):

            # Checks for any amino acids/molecules not in the csv file first!
            if residue.name not in atom_label_dict.keys():
                self._logger.warning(f"Residue: {residue.name}{residue.number} not in atom_label.csv!")
                return

            #This way we don't need so many naming schemes
            nterm_hydrogens = []
            cterm_oxygens = []
            if terminus == "NTERM":
                #Finding hydrogens attached to nitrogen
                try:
                    for a in residue.get_atom("N").bonds:
                        if a.element.upper() == "H":
                            nterm_hydrogens.append(a)
                
                except:
                    self._logger.warn("Protein does not have a nitrogen at n-terminus")

            if terminus == "CTERM":
                try:
                    for a in residue.get_atom("C").bonds:
                        if a.element.upper() == "O":
                            cterm_oxygens.append(a)

                except:
                    self._logger.warn("Protein does not have a carbonly at c-terminus")

            #Find the column that has the current naming scheme present
            self._logger.debug("Checking for the scheme id")
            for schemeid in range(len(schemenames)):
                scheme = []
                nterm_scheme = []
                for namelist in atom_label_dict[residue.name]:
                    self._logger.debug(f"Adding scheme name: {namelist[schemeid]}")
                    scheme.append(namelist[schemeid])

                if terminus is not None:
                    for namelist in atom_label_dict[terminus]:
                        scheme.append(namelist[schemeid])

                # Check to see if this is the naming scheme
                for atom in residue.atoms:
                    #Want to ignore and hydrogens attached to the nitrogen
                    if atom not in nterm_hydrogens and atom not in cterm_oxygens and atom.id not in scheme:
                        self._logger.debug(f"Atom: {atom.id} not in scheme: {scheme}")
                        break

                else:
                    break # This is the correct naming scheme!

            else:
                raise ValueError(f"Could not find the naming scheme for {residue.name}{residue.number} {atom}")

            #Loop over all of the atoms
            self._logger.debug("Found the scheme id")
            for atom in residue.atoms:
                if atom in nterm_hydrogens or atom in cterm_oxygens:
                    #We want to deal with these seperately
                    continue

                old_atomid = scheme.index(atom.id)
                if terminus is not None:
                    atom.id = (atom_label_dict[residue.name] + atom_label_dict[terminus])[old_atomid][newid]

                else:
                    atom.id = atom_label_dict[residue.name][old_atomid][newid]

            if nterm_hydrogens:
                for schemeid in range(len(schemenames)):
                    scheme = []
                    for namelist in atom_label_dict[residue.name]:
                        self._logger.debug(f"Adding scheme name: {namelist[schemeid]}")
                        scheme.append(namelist[schemeid])

                    if terminus is not None:
                        for namelist in atom_label_dict[terminus]:
                            scheme.append(namelist[schemeid])

                    # Check to see if this is the naming scheme
                    for atom in nterm_hydrogens:
                        if atom.id not in scheme:
                            self._logger.debug(f"Atom: {atom.id} not in scheme: {scheme}")
                            break

                    else:
                        break # This is the correct naming scheme!

                else:
                    raise ValueError(f"Could not find the naming scheme for {residue.name}{residue.number} {atom}")
            
                self._logger.debug("Found the scheme id for nterminus hydrogens")
                for atom in nterm_hydrogens:
                    old_atomid = scheme.index(atom.id)
                    if terminus is not None:
                        atom.id = (atom_label_dict[residue.name] + atom_label_dict[terminus])[old_atomid][newid]

                    else:
                        atom.id = atom_label_dict[residue.name][old_atomid][newid]
            
            if cterm_oxygens:
                for schemeid in range(len(schemenames)):
                    scheme = []
                    for namelist in atom_label_dict[residue.name]:
                        self._logger.debug(f"Adding scheme name: {namelist[schemeid]}")
                        scheme.append(namelist[schemeid])

                    if terminus is not None:
                        for namelist in atom_label_dict[terminus]:
                            scheme.append(namelist[schemeid])

                    # Check to see if this is the naming scheme
                    for atom in cterm_oxygens:
                        if atom.id not in scheme:
                            self._logger.debug(f"Atom: {atom.id} not in scheme: {scheme}")
                            break

                    else:
                        break # This is the correct naming scheme!

                else:
                    raise ValueError(f"Could not find the naming scheme for {residue.name}{residue.number} {atom}")
            
                self._logger.debug("Found the scheme id for nterminus hydrogens")
                for atom in cterm_oxygens:
                    old_atomid = scheme.index(atom.id)
                    if terminus is not None:
                        atom.id = (atom_label_dict[residue.name] + atom_label_dict[terminus])[old_atomid][newid]

                    else:
                        atom.id = atom_label_dict[residue.name][old_atomid][newid]


            # TODO add Jacks glycine hydrogen fixed so that naming convention is always the same with Chimera

        for chain in self.chains:
            index = 0
            rename_residue(chain.residues[0], "NTERM")
            rename_residue(chain.residues[-1], "CTERM")
            cterm = chain.residues[-1]

            for residue in chain.residues[1:-1]:
                self._logger.debug(f"Relabeling residue: {residue.name} {residue.number}")
                rename_residue(residue)

    #TODO, better way of determinging what is 'bound' to the metal my guess.
    def atoms_near_metal(self, metal, cutoff = 3.05):
        atom_list = [atom for c in self.chains[:-1]\
                            for r in c.residues\
                                for atom in r.atoms\
                                if atom.element.lower() in constants.HEAVY_ATOMS\
                                and np.linalg.norm(atom.coords - metal.coords) < cutoff]

        if not self.sub_chain.residues:
            atom_list.extend([atom for r in self.chains[-1].residues\
                                    for atom in r.atoms\
                                        if atom.element.lower() in constants.HEAVY_ATOMS\
                                        and np.linalg.norm(atom.coords - metal.coords) < cutoff])
       
        return atom_list

    def make_bond_table(self):
        self.write_pdb("bond.pdb")

        successful = False

        with Popen(f"babel bond.pdb bond.mol2", stdin=PIPE, stdout=PIPE, stderr=PIPE,
                   universal_newlines=True, shell=True, bufsize=1, env=os.environ) as shell:
            while shell.poll() is None:
                self._logger.debug(shell.stdout.readline().strip())
                output = shell.stderr.readline().strip()
                self._logger.debug(output)
                if "1 molecule converted" in output:
                    successful = True

        if not successful:
            self._logger.error("Could not create {residue.name} mol2 file!")
            raise OSError("mol2_file")

        #Clear the bond lists first:
        for chain in self.chains:
            for residue in chain.residues:
                for atom in residue.atoms:
                    atom.bonds.clear()

        atom_list = [atom for c in self.chains for r in c.residues for atom in r.atoms]

        with open("bond.mol2") as mol_file:
            bond_section = False
            for line in mol_file:
                if "BOND" in line:
                    bond_section = True

                elif bond_section:
                    line = line.split()
                    atom_list[int(line[1])-1].add_bond(atom_list[int(line[2]) - 1])

        self._logger.info("Succesfully created the bond lists for each atom")
        self._logger.debug("Cleaning up files created")
        os.remove("bond.pdb")
        os.remove("bond.mol2")

    def get_coords(self):
        if self.coords is None:
            self.coords = [atom.coords for chain in self.chains for residue in chain.residues for atom in residue.atoms]
        
        return self.coords
        
    def aa_rmsd(self, pro):
        
        if pro is self:
            return 0.0

        #Makes it easier to work with this sort of stuff

        this_coords = self.get_coords() - np.average(self.get_coords(), axis=0)
        other_coords = pro.get_coords() - np.average(pro.get_coords(), axis=0)

        #Now the coords have been centered, now we can start the rotation
        assert(len(this_coords) == len(other_coords))
        n_vec = np.shape(this_coords)[0]

        h = np.dot(np.transpose(this_coords), other_coords)
        v, s, w = np.linalg.svd(h)

        is_reflection = (np.linalg.det(v) * np.linalg.det(w)) < 0.0

        if is_reflection:
            s[-1] = -s[-1]

        other_coords = np.dot(other_coords, np.dot(v,w))
       
        diff = other_coords - this_coords

        return np.sqrt((diff*diff).sum()/ len(this_coords))




