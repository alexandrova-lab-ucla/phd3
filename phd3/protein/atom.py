#!/usr/bin/env python3
"""
Author  ==>> Matthew R. Hennefarth
Date    ==>> April 16, 2020
"""

#Standard Library Imports
import numpy as np

#PHD3 Imports
import phd3.utility.constants as constants

__all__=[
    'Atom'
]

class Atom:

    __slots__ = ['element', 'coords', 'id', 'residue', 'chain', 'number', 'bonds', 'freeze']

    def __init__(self, line: str = None, element: str = None, coords: np.array = None, id=None, number=None):

        if line is None:
            self.element = element
            self.coords = coords
            self.id = id
            self.number = number

        elif line is not None:
            self.element = line[76:78].strip().lower()
            self.coords = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            self.id = line[12:16].strip().upper()
            self.number = int(line[6:11])
            #Fix the formatting
            self.element = self.element.capitalize()


        if self.element.lower() == 'eh':
            self.element = 'h'

        self.residue = None
        self.chain = None
        self.freeze = False
        self.bonds = []


    def write_inConstr(self):
        if self.residue.name in constants.AMINO_ACID_RESIDUES:
            return f"{ord(self.chain.name)-ord('A')+1}.{self.residue.inConstr_number}.{self.id.upper()}"

        else:
            return f"{ord(self.chain.name) - ord('A') + self.residue.inConstr_number}.1.{self.id.upper()}"

    def coord_line(self):
        if self.element.upper() == "ZN":
            actual_element = self.id

        else:
            actual_element = self.element

        return f"    {self.coords[0] * constants.A_TO_BOHR:.5f} {self.coords[1] * constants.A_TO_BOHR:.5f} {self.coords[2] * constants.A_TO_BOHR:.5f} {actual_element}{' f' if self.freeze else ''}\n"


    def pdb_line(self):
        if self.element.lower() in constants.METALS:
            return '{:<6}{:>5} {:<4} {} {}{:>4}    {:>8.3f}{:>8.3f}{:>8.3f}  1.00  0.00          {:>2}\n'.format(
              "HETATM",
              self.number, self.id, self.residue.name, self.residue.chain.name, self.residue.number,
              self.coords[0], self.coords[1], self.coords[2], self.element.capitalize())

        return '{:<6}{:>5} {:<4} {} {}{:>4}    {:>8.3f}{:>8.3f}{:>8.3f}  1.00  0.00          {:>2}\n'.format(
            'ATOM' if self.residue.name in constants.AMINO_ACID_RESIDUES else "HETATM",
            self.number, self.id if len(self.id) > 3 else f" {self.id}", self.residue.name, self.residue.chain.name, self.residue.number,
            self.coords[0], self.coords[1], self.coords[2], self.element.capitalize())

    def __str__(self):
        return f"{self.id} {self.coords} {self.element}"

    def label(self):
        return f"{self.chain.name}:{self.residue.number}:{self.id}"

    def add_bond(self, atom):
        self.bonds.append(atom)
        atom.bonds.append(self)
