#!/usr/bin/env python3
"""
Author  ==>> Matthew R. Hennefarth
Date    ==>> April 16, 2020
"""

import copy

#PHD3 Imports
from . import residue

__all__=[
    'Chain'
]

class Chain:

    __slots__ = ['name', 'residues']

    def __init__(self, name: str=''):
        self.name = name
        self.residues = []

    def add_residue(self, res: residue.Residue):
        res.set_chain(self)
        self.residues.append(res)

    def __str__(self):
        return self.name

    def write_inConstr(self):
        return f"{ord(self.chain.name) - ord('A') + 1}.*.*"

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result

        result.name = copy.deepcopy(self.name, memo)
        for res in self.residues:
            copy_res = copy.deepcopy(res)
            result.add_residue(copy_res)

        return result

