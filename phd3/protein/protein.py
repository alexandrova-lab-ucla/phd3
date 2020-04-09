#!/usr/bin/env python3


import logging
import pkg_resources
import csv
import os
import numpy as np
from subprocess import Popen, PIPE

import phd3.utility as constants
import phd3.protein.chain as chain
import phd3.protein.residue as residue

__all__=[
    'Protein'
]

class Protein:

    __slots__ = ['chains', '_logger', 'non_residues', 'metals', 'name', 'sub_chain']

    def __init__(self, name: str, chains: [chain.Chain]):

        self._logger = logging.getLogger(__name__)

        self._logger.debug(f"Initializing vars for protein {name}")
        self.name = name
        self.chains = chains
        self.non_residues = []
        self.metals = []
        self.sub_chain = chain.Chain()

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
                        if self.chains[chain_num].residues[residue_num].atoms[atom_num].element in constants.METALS:
                            self._logger.debug(f"Found a metal: {self.chains[chain_num].residues[residue_num].atoms[atom_num]} in residue {self.chains[chain_num].residues[residue_num]}")
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
            for residue in self.non_residues:
                start = 100
                for atom in residue.atoms:
                    atom.id = f"{atom.element.upper()[0]}{start}"
                    start += 1

                self._logger.debug(f"Adding residue {residue} to substrate chain")
                self.sub_chain.add_residue(residue)
                residue.number = res_num
                residue.inConstr_number = res_num
                for atom in residue.atoms:
                    atom.number = atom_num
                    atom_num += 1

                res_num += 1

            self._logger.debug("Adding substrate chain to master chain")
            self.chains.append(self.sub_chain)

        if relabel_protein:
            self._logger.debug("Relabeling the protein")
            self.relabel()
        
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

            #Find the column that has the current naming scheme present
            self._logger.debug("Checking for the scheme id")
            for schemeid in range(len(schemenames)):
                scheme = []
                for namelist in atom_label_dict[residue.name]:
                    self._logger.debug(f"Adding scheme name: {namelist[schemeid]}")
                    scheme.append(namelist[schemeid])

                if terminus is not None:
                    for namelist in atom_label_dict[terminus]:
                        scheme.append(namelist[schemeid])

                # Check to see if this is the naming scheme
                for atom in residue.atoms:
                    if atom.id not in scheme:
                        self._logger.debug(f"Atom: {atom.id} not in scheme: {scheme}")
                        break

                else:
                    break # This is the correct naming scheme!

            else:
                raise ValueError(f"Could not find the naming scheme for {residue.name}{residue.number} {atom}")

            #Loop over all of the atoms
            self._logger.debug("Found the scheme id")
            for atom in residue.atoms:
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

    def atoms_near_metal(self, metal, cutoff = 3.05):
        atom_list = []
        for chain in self.chains[:-1]:
            for residue in chain.residues:
                for atom in residue.atoms:
                    if atom.element in constants.HEAVY_ATOMS:
                        if np.linalg.norm(atom.coords - metal.coords) < cutoff:
                            atom_list.append(atom)

        if not self.sub_chain.residues:
            for residue in self.chains[-1].residues:
                for atom in residue.atoms:
                    if atom.element in constants.HEAVY_ATOMS:
                        if np.linalg.norm(atom.coords - metal.coords) < cutoff:
                            atom_list.append(atom)

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

        atom_list = []
        for chain in self.chains:
            for residue in chain.residues:
                for atom in residue.atoms:
                    atom_list.append(atom)

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

    def center(self):
        cen = np.array([0.0, 0.0, 0.0])
        num_atoms = 0

        for chain in self.chains:
            for residue in chain.residues:
                for atom in residue.atoms:
                    cen += atom.coords
                    num_atoms += 1

        try:
            cen = cen / num_atoms
        
        except ZeroDivisionError:
            self._logger.error("No atoms, cannot get center")
            raise
        
        return cen

    def aa_rmsd(self, pro):
        
        #Makes it easier to work with this sort of stuff
        this_atoms = []
        for chain in self.chains:
            for residue in chain.residues:
                this_atoms.extend(residue.atoms)

        this_center = self.center()
        this_coords = [a.coords.copy() - this_center for a in this_atoms]
        this_coords = np.array(this_coords)


        other_atoms = []
        for chain in pro.chains:
            for residue in chain.residues:
                other_atoms.extend(residue.atoms)

        other_center = pro.center()
        other_coords = [a.coords.copy() - other_center for a in other_atoms]
        other_coords = np.array(other_coords)


        #Now the coords have been centered, now we can start the rotation
        assert(len(this_coords) == len(other_coords))
        n_vec = np.shape(this_coords)[0]

        h = np.dot(np.transpose(this_coords), other_coords)
        v, s, w = np.linalg.svd(h)

        is_reflection = (np.linalg.det(v) * np.linalg.det(w)) < 0.0

        if is_reflection:
            s[-1] = -s[-1]

        E0 = sum(sum(this_coords*this_coords)) + sum(sum(other_coords*other_coords))

        rmsd_sq = (E0-2.0*sum(s)) / float(n_vec)
        rmsd_sq = max([rmsd_sq, 0.0])

        rmsd = np.sqrt(rmsd_sq)
        return rmsd




