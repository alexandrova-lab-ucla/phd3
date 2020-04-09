#!/usr/bin/env python3

import phd3.protein.residue as Residue

__all__=[
    'Chain'
]

class Chain:

    __slots__ = ['name', 'residues']

    def __init__(self, name: str=''):
        self.name = name
        self.residues = []

    def add_residue(self, res: Residue):
        res.chain = self
        self.residues.append(res)

    def __str__(self):
        return self.name

    def write_inConstr(self):
        return f"{ord(self.chain.name) - ord('A') + 1}.*.*"
