#!/usr/bin/env python3
"""
Author  ==>> David J. Reilley
Date    ==>> April 16, 2020
"""

# Tabulated data used by Titratable-QM/DMD
# DICTIONARY: number of letters of the alphabet for chains
# There must be a built-in function for this in python, but I couldn't find it quickly
#given a character, you can get the ASCII number using ord()
#hence to get the number B to be 2, you would do ord('B')- ord('A') +1
# since the difference between acii number of A and B is 1...
alphabet_order = {}
alphabet_order['A'] = 1
alphabet_order['B'] = 2
alphabet_order['C'] = 3
alphabet_order['D'] = 4
alphabet_order['E'] = 5
alphabet_order['F'] = 6
alphabet_order['G'] = 7
alphabet_order['F'] = 8
alphabet_order['I'] = 9
alphabet_order['J'] = 10
# LIST: all amino acid 3 letter codes
amino_acids_3let = ['ALA', 'ARG', 'ASP', 'ASN', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SEC', 'SER', 'THR', 'TRP', 'TYR', 'VAL'] # proteinogenic amino acid 3 letter codes
amino_acids_3let += ['3CT'] # non-proteinogenic amino acid 3 letter codes
# LIST: 3 letter codes for amino acids which are potentially titratable near physiological conditions
titr_amino_acids_3let = ['ARG', 'ASP', 'CYS', 'GLU', 'HIS', 'LYS', 'TYR']
# LIST: protonation state 3 letter codes
titr_forms_3let = ['ARG', 'ASP', 'CYS', 'GLU', 'HIS', 'HIE', 'HID', 'HIP', 'LYS', 'LYZ', 'TYR', 'TYH']
titr_forms_3let = ['ASPH', 'GLUH']
# LIST: hydrogens implicated in protonation states by amino acid
titr_form_prots = ['ARG:HE', 'ARG:1HH1', 'ARG:2HH1', 'ARG:1HH2', 'ARG:2HH2', 'HIS:1HNE', 'HIS:HE2', 'HIS:HD1', 'LYS:HZ1', 'LYS:HZ2', 'LYS:HZ3', 'TYR:HO', 'TYR:HH'] # standard names just from the physiological forms of the amino acids
titr_form_prots += ['ARG:HH11', 'ARG:HH12', 'ARG:HH21', 'ARG:HH22', 'ASP:1HND', 'ASP:2HND', 'CYS:HG1', 'CYS:HG', 'GLU:1HNE', 'GLU:2HNE', 'HIE:1HNE', 'HID:HD1', 'HIP:1HNE', 'HIP:HD1', 'LYZ:HZ1', 'LYZ:HZ2', 'LYZ:HZ3'] # standard names from less common protonation states
titr_form_prots += ['ASP:HD1', 'ASP:HD2', 'GLU:HE1', 'GLU:HE2'] # For the alternative naming of protonated ASP and GLU
titr_form_prots += ['N+:H1', 'N+:H2', 'N+:H3', 'N+:1H', 'N+:2H', 'N+:HN', 'N+:1HN', 'N+:2HN', 'N+:HN1', 'N+:HN2', 'C-:HOXT1', 'C-:HOXT2', 'C-:HO1', 'C-:HO2'] # N-terminal and C-terminal protons #fix
# LIST: heretoatoms implicated in protonation states by amino acid
titr_form_heteroatoms = ['ARG:NE', 'ARG:NH1', 'ARG:NH2', 'ASP:OD1', 'ASP:OD2', 'CYS:SG', 'GLU:OE1', 'GLU:OE2', 'HIS:1NE', 'HIS:NE2', 'HIS:ND1', 'LYS:NZ', 'TYR:OH'] # standard names just from the physiological forms of the amino acids
titr_form_heteroatoms += ['ASPH:OD1', 'ASPH:OD2', 'GLUH:OE1', 'GLUH:OE2']
titr_form_heteroatoms += ['N+:N', 'C-:O', 'C-:OXT', 'C-:O1', 'C-:O2'] # N-terminal and C-terminal heteroatoms
# LIST: names for heteroatoms involved in covalent links
cov_link_heteroatoms = ['CYS:SG']
# DICTIONARY: definitions of protonation states based on hydrogen present
prots2titr_form = {}
prots2titr_form['ARG:1HH1:1HH2:2HH1:2HH2:HE:NE:NH1:NH2'] = 'ARG:+:1'
prots2titr_form['ARG:1HH1:1HH2:2HH1:2HH2:NE:NH1:NH2'] = 'ARG:-:1'
prots2titr_form['ARG:1HH1:1HH2:2HH1:HE:NE:NH1:NH2'] = 'ARG:-:2'
prots2titr_form['ARG:1HH1:1HH2:2HH2:HE:NE:NH1:NH2'] = 'ARG:-:3'
prots2titr_form['ARG:1HH1:2HH1:2HH2:HE:NE:NH1:NH2'] = 'ARG:-:4'
prots2titr_form['ARG:1HH2:2HH1:2HH2:HE:NE:NH1:NH2'] = 'ARG:-:5'
prots2titr_form['ARG:HE:HH11:HH12:HH21:HH22:NE:NH1:NH2'] = 'ARG:+:11'
prots2titr_form['ARG:HH11:HH12:HH21:HH22:NE:NH1:NH2'] = 'ARG:-:11'
prots2titr_form['ARG:HE:HH11:HH12:HH21:NE:NH1:NH2'] = 'ARG:-:21'
prots2titr_form['ARG:HE:HH11:HH12:HH22:NE:NH1:NH2'] = 'ARG:-:31'
prots2titr_form['ARG:HE:HH11:HH21:HH22:NE:NH1:NH2'] = 'ARG:-:41'
prots2titr_form['ARG:HE:HH12:HH21:HH22:NE:NH1:NH2'] = 'ARG:-:51'
prots2titr_form['ASP:OD1:OD2'] = 'ASP:-:1'
prots2titr_form['ASP:HD1:OD1:OD2'] = 'ASP:+:1'
prots2titr_form['ASP:HD2:OD1:OD2'] = 'ASP:+:2'
prots2titr_form['CYS:SG'] = 'CYS:-:1'
prots2titr_form['CYS:HG:SG'] = 'CYS:+:1'
prots2titr_form['GLU:OE1:OE2'] = 'GLU:-:1'
prots2titr_form['GLU:HE1:OE1:OE2'] = 'GLU:+:1'
prots2titr_form['GLU:HE2:OE1:OE2'] = 'GLU:+:2'
prots2titr_form['HIS:HD1:1NE:ND1'] = 'HIS:-:1'
prots2titr_form['HIS:1HNE:1NE:ND1'] = 'HIS:-:2'
prots2titr_form['HIS:HD1:ND1:NE2'] = 'HIS:-:11'
prots2titr_form['HIS:HE2:ND1:NE2'] = 'HIS:-:21'
prots2titr_form['HIS:1HNE:HD1:1NE:ND1'] = 'HIS:+:1'
prots2titr_form['HIS:HD1:HE2:ND1:NE2'] = 'HIS:+:11'
prots2titr_form['LYS:HZ1:HZ2:HZ3:NZ'] = 'LYS:+:1'
prots2titr_form['LYS:HZ1:HZ2:NZ'] = 'LYS:-:3'
prots2titr_form['LYS:HZ1:HZ3:NZ'] = 'LYS:-:2'
prots2titr_form['LYS:HZ2:HZ3:NZ'] = 'LYS:-:1'
prots2titr_form['TYR:HH:OH'] = 'TYR:+:1'
prots2titr_form['TYR:OH'] = 'TYR:-:1'
prots2titr_form['N+:H1:H2:H3:N'] = 'N+:+:1'
prots2titr_form['N+:H1:H2:N'] = 'N+:-:3'
prots2titr_form['N+:H1:H3:N'] = 'N+:-:2'
prots2titr_form['N+:H2:H3:N'] = 'N+:-:1'
prots2titr_form['N+:1H:2H:HN:N'] = 'N+:+:11'
prots2titr_form['N+:1H:2H:N'] = 'N+:-:31'
prots2titr_form['N+:1H:HN:N'] = 'N+:-:21'
prots2titr_form['N+:2H:HN:N'] = 'N+:-:11'
prots2titr_form['N+:1HN:2HN:HN:N'] = 'N+:+:12'
prots2titr_form['N+:1HN:2HN:N'] = 'N+:-:32'
prots2titr_form['N+:1HN:HN:N'] = 'N+:-:22'
prots2titr_form['N+:2HN:HN:N'] = 'N+:-:12'
prots2titr_form['N+:HN1:HN2:HN:N'] = 'N+:+:13'
prots2titr_form['N+:HN1:HN2:N'] = 'N+:-:33'
prots2titr_form['N+:HN1:HN:N'] = 'N+:-:23'
prots2titr_form['N+:HN2:HN:N'] = 'N+:-:13'
prots2titr_form['N+:HN:N'] = 'N+:-:43' # If there is only one N-terminal hydrogen, the second one is included implicitly
prots2titr_form['N+:N'] = 'N+:-:43' # If there are no N-terminal hydrogen, add them implicitly
#prots2titr_form['N+:HN1:HN2:HN3:N'] = 'N+:+:14'
#prots2titr_form['N+:HN1:HN2:N'] = 'N+:-:34'
#prots2titr_form['N+:HN1:HN3:N'] = 'N+:-:24'
#prots2titr_form['N+:HN2:HN3:N'] = 'N+:-:14'
#prots2titr_form['N+:HN1:N'] = 'N+:-:24'
#prots2titr_form['N+:HN2:N'] = 'N+:-:14'
#prots2titr_form['N+:HN3:N'] = 'N+:-:14'
prots2titr_form['C-:O1:O2'] = 'C-:-:1'
prots2titr_form['C-:HO1:O1:OXT2'] = 'C-:+:1' #fixprots2titr_form['C-:HOXT2:O:OXT'] = 'C-:+:2' #fix
prots2titr_form['C-:O'] = 'C-:-:11'
prots2titr_form['C-:HOXT1:O'] = 'C-:+:11' #fix
prots2titr_form['C-:HOXT2:O'] = 'C-:+:21' #fix
# DICTIONARY: new protonation state and changed atoms based on the old protonation state and whether a proton is added or removed
# Provides each possible state in a list
old_titr_form2new_titr_form = {}
old_titr_form2new_titr_form['ARG:Add:1'] = [[['+', 1],['HE']]]
old_titr_form2new_titr_form['ARG:Add:2'] = [[['+', 1],['2HH2']]]
old_titr_form2new_titr_form['ARG:Add:3'] = [[['+', 1],['2HH1']]]
old_titr_form2new_titr_form['ARG:Add:4'] = [[['+', 1],['1HH2']]]
old_titr_form2new_titr_form['ARG:Add:5'] = [[['+', 1],['1HH1']]]
old_titr_form2new_titr_form['ARG:Remove:1'] = [[['-', 4],['1HH2']]]
#old_titr_form2new_titr_form['ARG:Remove:1'] = [[['-', 1],['HE']], [['-', 2],['2HH2']], [['-', 3],['2HH1']], [['-', 4],['1HH2']], [['-', 5],['1HH1']]]
old_titr_form2new_titr_form['ARG:Add:11'] = [[['+', 11],['HE']]]
old_titr_form2new_titr_form['ARG:Add:21'] = [[['+', 11],['HH22']]]
old_titr_form2new_titr_form['ARG:Add:31'] = [[['+', 11],['HH21']]]
old_titr_form2new_titr_form['ARG:Add:41'] = [[['+', 11],['HH12']]]
old_titr_form2new_titr_form['ARG:Add:51'] = [[['+', 11],['HH11']]]
old_titr_form2new_titr_form['ARG:Remove:11'] = [[['-', 41],['HH12']]]
#old_titr_form2new_titr_form['ARG:Remove:11'] = [[['-', 11],['HE']], [['-', 21],['HH22']], [['-', 31],['HH21']], [['-', 41],['HH12']], [['-', 51],['HH11']]]
old_titr_form2new_titr_form['ASP:Add:1'] = [[['+', 2],['HD2']]]
#old_titr_form2new_titr_form['ASP:Add:1'] = [[['+', 1],'1HND'], [['+', 2],'2HND']] # Currently, only one ASP titratable proton is defined in DMD, so 2HND is unneccessary
old_titr_form2new_titr_form['ASP:Remove:1'] = [[['-', 1],['HD1']]]
old_titr_form2new_titr_form['ASP:Remove:2'] = [[['-', 1],['HD2']]]
old_titr_form2new_titr_form['CYS:Add:1'] = [[['+', 1],['HG']]]
old_titr_form2new_titr_form['CYS:Remove:1'] = [[['-', 1],['HG']]]
old_titr_form2new_titr_form['GLU:Add:1'] = [[['+', 2],['HE2']]]
#old_titr_form2new_titr_form['GLU:Add:1'] = [[['+', 1],'1HNE'], [['+', 2],'2HNE']] # Currently, only one GLU titratable proton is defined in DMD, so 2HNE is unneccessary
old_titr_form2new_titr_form['GLU:Remove:1'] = [[['-', 1],['HE1']]]
old_titr_form2new_titr_form['GLU:Remove:2'] = [[['-', 1],['HE2']]]
old_titr_form2new_titr_form['HIS:Add:1'] = [[['+', 1],['1NHE']]]
old_titr_form2new_titr_form['HIS:Add:2'] = [[['+', 1],['HD1']]]
old_titr_form2new_titr_form['HIS:Remove:1'] = [[['-', 1],['1NHE']], [['-', 2],['HD1']]]
old_titr_form2new_titr_form['HIS:Add:11'] = [[['+', 11],['HE2']]]
old_titr_form2new_titr_form['HIS:Add:21'] = [[['+', 11],['HD1']]]
old_titr_form2new_titr_form['HIS:Remove:11'] = [[['-', 11],['HE2']], [['-', 21],['HD1']]]
old_titr_form2new_titr_form['LYS:Add:1'] = [[['+', 1],['HZ1']]]
old_titr_form2new_titr_form['LYS:Add:2'] = [[['+', 1],['HZ2']]]
old_titr_form2new_titr_form['LYS:Add:3'] = [[['+', 1],['HZ3']]]
old_titr_form2new_titr_form['LYS:Remove:1'] = [[['-', 3],['HZ3']]]
#old_titr_form2new_titr_form['LYS:Remove:1'] = [[['-', 1],['HZ1']], [['-', 2],['HZ2']], [['-', 3],['HZ3']]]
old_titr_form2new_titr_form['TYR:Add:1'] = [[['+', 1],['HH']]]
old_titr_form2new_titr_form['TYR:Remove:1'] = [[['-', 1],['HH']]]
old_titr_form2new_titr_form['N+:Add:1'] = [[['+', 1],['H1']]]
old_titr_form2new_titr_form['N+:Add:2'] = [[['+', 1],['H2']]]
old_titr_form2new_titr_form['N+:Add:3'] = [[['+', 1],['H3']]]
old_titr_form2new_titr_form['N+:Remove:1'] = [[['-', 3],['H3']]]
#old_titr_form2new_titr_form['N+:Remove:1'] = [[['-', 1],['H1']], [['-', 2],['H2']], [['-', 3],['H3']]]
old_titr_form2new_titr_form['N+:Add:11'] = [[['+', 11],['1H']]]
old_titr_form2new_titr_form['N+:Add:21'] = [[['+', 11],['2H']]]
old_titr_form2new_titr_form['N+:Add:31'] = [[['+', 11],['HN']]]
old_titr_form2new_titr_form['N+:Remove:11'] = [[['-', 11],['1H']], [['-', 21],['2H']], [['-', 31],['HN']]]
old_titr_form2new_titr_form['N+:Add:12'] = [[['+', 12],['1HN']]]
old_titr_form2new_titr_form['N+:Add:22'] = [[['+', 12],['2HN']]]
old_titr_form2new_titr_form['N+:Add:32'] = [[['+', 12],['HN']]]
old_titr_form2new_titr_form['N+:Remove:12'] = [[['-', 12],['1HN']], [['-', 22],['2HN']], [['-', 32],['HN']]]
old_titr_form2new_titr_form['N+:Add:13'] = [[['+', 13],['HN1']]]
old_titr_form2new_titr_form['N+:Add:23'] = [[['+', 13],['HN2']]]
old_titr_form2new_titr_form['N+:Add:33'] = [[['+', 13],['HN']]]
old_titr_form2new_titr_form['N+:Add:43'] = [[['+', 13],['HN3','HN2']]]
old_titr_form2new_titr_form['N+:Remove:13'] = [[['-', 13],['HN1']], [['-', 23],['HN2']], [['-', 33],['HN']]]
#old_titr_form2new_titr_form['N+:Add:14'] = [[['+', 14],'HN1']]
#old_titr_form2new_titr_form['N+:Add:24'] = [[['+', 14],'HN2']]
#old_titr_form2new_titr_form['N+:Add:34'] = [[['+', 14],'HN3']]
#old_titr_form2new_titr_form['N+:Remove:14'] = [[['-', 14],'HN1'], [['-', 24],'HN2'], [['-', 34],'HN3']]
old_titr_form2new_titr_form['C-:Add:1'] = [[['+', 1],['HO1']]] #fix
#old_titr_form2new_titr_form['C-:Add:1'] = [[['+', 1],['HO1']], [['+', 2],['HO2']]] #fix
old_titr_form2new_titr_form['C-:Remove:1'] = [[['-', 1],['HO1']]] #fix
old_titr_form2new_titr_form['C-:Remove:2'] = [[['-', 1],['HO2']]] #fix
old_titr_form2new_titr_form['C-:Add:11'] = [[['+', 11],['HOXT1']], [['+', 21],['HOXT2']]] #fix
#old_titr_form2new_titr_form['C-:Add:11'] = [[['+', 11],['HOXT1']], [['+', 21],['HOXT2']]] #fix
old_titr_form2new_titr_form['C-:Remove:11'] = [[['-', 11],['HOXT1']]] #fix
old_titr_form2new_titr_form['C-:Remove:21'] = [[['-', 11],['HOXT2']]] #fix
# DICTIONARY: bound heteroatom data by added hydrogen
hydrogen2boundheteroatom = {}
hydrogen2boundheteroatom['ARG:HE:1'] = ['NE','CD','CZ']
hydrogen2boundheteroatom['ARG:2HH2:2'] = ['NH2','CZ','1HH2']
hydrogen2boundheteroatom['ARG:2HH1:3'] = ['NH2','CZ','1HH1']
hydrogen2boundheteroatom['ARG:1HH2:4'] = ['NH1','CZ','2HH2']
hydrogen2boundheteroatom['ARG:1HH1:5'] = ['NH1','CZ','2HH1']
hydrogen2boundheteroatom['ARG:HH22:21'] = ['NH2','CZ','HH21']
hydrogen2boundheteroatom['ARG:HH21:31'] = ['NH2','CZ','HH11']
hydrogen2boundheteroatom['ARG:HH12:41'] = ['NH1','CZ','HH22']
hydrogen2boundheteroatom['ARG:HH11:51'] = ['NH1','CZ','HH12']
hydrogen2boundheteroatom['ASP:HD1:1'] = ['OD1','CG','CB']
hydrogen2boundheteroatom['ASP:HD2:1'] = ['OD2','CG','CB']
hydrogen2boundheteroatom['CYS:HG:1'] = ['SG','CB','CA']
hydrogen2boundheteroatom['GLU:HE1:1'] = ['OE1','CD','CG']
hydrogen2boundheteroatom['GLU:HE2:1'] = ['OE2','CD','CG']
hydrogen2boundheteroatom['HIS:1NHE:1'] = ['1NE','CE1','CD2']
hydrogen2boundheteroatom['HIS:HD1:2'] = ['ND1','CE1','CG']
hydrogen2boundheteroatom['HIS:HE2:11'] = ['NE2','CE1','CD2']
hydrogen2boundheteroatom['HIS:HD1:21'] = ['ND1','CE1','CD2']
hydrogen2boundheteroatom['LYS:HZ1:1'] = ['NZ','CE','HZ2','HZ3']
hydrogen2boundheteroatom['LYS:HZ2:2'] = ['NZ','CE','HZ1','HZ3']
hydrogen2boundheteroatom['LYS:HZ3:3'] = ['NZ','CE','HZ1','HZ2']
hydrogen2boundheteroatom['TYR:HH:1'] = ['OH','CZ','CE1']
hydrogen2boundheteroatom['N+:H1:1'] = ['N','CA','H2','H3']
hydrogen2boundheteroatom['N+:H2:2'] = ['N','CA','H1','H3']
hydrogen2boundheteroatom['N+:H3:3'] = ['N','CA','H1','H2']
hydrogen2boundheteroatom['N+:1H:11'] = ['N','CA','2H','HN']
hydrogen2boundheteroatom['N+:2H:21'] = ['N','CA','1H','HN']
hydrogen2boundheteroatom['N+:HN:31'] = ['N','CA','1H','2H']
hydrogen2boundheteroatom['N+:1HN:12'] = ['N','CA','2HN','HN']
hydrogen2boundheteroatom['N+:2HN:22'] = ['N','CA','1HN','HN']
hydrogen2boundheteroatom['N+:HN:32'] = ['N','CA','1HN','2HN']
hydrogen2boundheteroatom['N+:HN3:13'] = ['N','CA','HN2','HN']
hydrogen2boundheteroatom['N+:HN2:23'] = ['N','CA','HN']
hydrogen2boundheteroatom['N+:HN:33'] = ['N','CA','HN1','HN2']
hydrogen2boundheteroatom['N+:HN3:43'] = ['N','CA','HN2','HN']
hydrogen2boundheteroatom['N+:HN2:43'] = ['N','CA','HN']
#hydrogen2boundheteroatom['N+:HN1:14'] = ['N','CA','HN2','HN3']
#hydrogen2boundheteroatom['N+:HN2:24'] = ['N','CA','HN1','HN3']
#hydrogen2boundheteroatom['N+:HN3:34'] = ['N','CA','HN1','HN2']
hydrogen2boundheteroatom['C-:HOXT1:1'] = ['O','C','CA'] # fix
hydrogen2boundheteroatom['C-:HOXT2:2'] = ['OXT','C','CA'] # fix
hydrogen2boundheteroatom['C-:HOXT1:11'] = ['O','C','CA'] # fix
hydrogen2boundheteroatom['C-:HOXT2:21'] = ['O','C','CA'] # fix
# DICTIONARY: Required protonation heteroatom names for constr files
#constr_heteroatom_name['1NE'] = 'NE2'
# DICTIONARY: Nearby atoms by added hydrogen, without titration form number
hydrogen2boundheteroatom_notitrform = {}
hydrogen2boundheteroatom_notitrform['ARG:HE'] = ['NE','CD','CZ']
hydrogen2boundheteroatom_notitrform['ARG:2HH2'] = ['NH2','CZ','1HH2']
hydrogen2boundheteroatom_notitrform['ARG:2HH1'] = ['NH2','CZ','1HH1']
hydrogen2boundheteroatom_notitrform['ARG:1HH2'] = ['NH1','CZ','2HH2']
hydrogen2boundheteroatom_notitrform['ARG:1HH1'] = ['NH1','CZ','2HH1']
hydrogen2boundheteroatom_notitrform['ASP:1HND'] = ['OD1','CG']
hydrogen2boundheteroatom_notitrform['ASP:2HND'] = ['OD2','CG']
hydrogen2boundheteroatom_notitrform['CYS:HG1'] = ['SG','CB']
hydrogen2boundheteroatom_notitrform['GLU:1HNE'] = ['OE1','CD']
hydrogen2boundheteroatom_notitrform['GLU:2HNE'] = ['OE2','CD']
hydrogen2boundheteroatom_notitrform['HIS:1NHE'] = ['NE2','CE1','CD2']
hydrogen2boundheteroatom_notitrform['HIS:HE2'] = ['NE2','CE1','CD2']
hydrogen2boundheteroatom_notitrform['HIS:HD1'] = ['ND1','CE1','CD2']
hydrogen2boundheteroatom_notitrform['LYS:HZ1'] = ['NZ','CE','HZ2','HZ3']
hydrogen2boundheteroatom_notitrform['LYS:HZ2'] = ['NZ','CE','HZ1','HZ3']
hydrogen2boundheteroatom_notitrform['LYS:HZ3'] = ['NZ','CE','HZ1','HZ2']
hydrogen2boundheteroatom_notitrform['TYR:HO'] = ['OH','CZ']
hydrogen2boundheteroatom_notitrform['N+:H1'] = ['N','CA','H2','H3']
hydrogen2boundheteroatom_notitrform['N+:H2'] = ['N','CA','H1','H3']
hydrogen2boundheteroatom_notitrform['N+:H3'] = ['N','CA','H1','H2']
hydrogen2boundheteroatom_notitrform['N+:1H'] = ['N','CA','2H','HN']
hydrogen2boundheteroatom_notitrform['N+:2H'] = ['N','CA','1H','HN']
hydrogen2boundheteroatom_notitrform['N+:HN'] = ['N','CA','1H','2H']
hydrogen2boundheteroatom_notitrform['N+:1HN'] = ['N','CA','2HN','HN']
hydrogen2boundheteroatom_notitrform['N+:2HN'] = ['N','CA','1HN','HN']
#hydrogen2boundheteroatom_notitrform['N+:HN1:14'] = ['N','CA','HN2','HN3']
#hydrogen2boundheteroatom_notitrform['N+:HN2:24'] = ['N','CA','HN1','HN3']
#hydrogen2boundheteroatom_notitrform['N+:HN3:34'] = ['N','CA','HN1','HN2']
hydrogen2boundheteroatom_notitrform['C-:HOXT1'] = ['O','C','CA'] # fix
hydrogen2boundheteroatom_notitrform['C-:HOXT2'] = ['OXT','C','CA'] # fix

