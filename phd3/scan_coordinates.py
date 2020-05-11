#!/usr/bin/env python3

import logging
import numpy as np
import json
import os
import shutil

from .setupjob import setupTMjob
from .qm_calculation import TMcalculation

__all__ = [
        'scan_coordinates',
        'scan_to_xyz'
        ]

logger = logging.getLogger(__name__)

def write_coords(coordinates):
    with open("coord", 'w+') as coord_file:
        coord_file.write("$coord\n")
        for line in coordinates:
            c = line[0]
            new_line = [str(c[0]), str(c[1]), str(c[2])]
            new_line += line[1:]
            coord_file.write(" ".join(new_line) + '\n')

        coord_file.write("$end")

def read_coords(coord_name="coord"):
    if not os.path.isfile("coord"):
        logger.error("No coord file specified")
        raise FileNotFoundError("coord")

    coords = []
    with open(coord_name) as coord_file:
        for line in coord_file:
            if line.startswith("$coord"):
                continue

            elif line.startswith("$"):
                break
            
            line = line.split()
            c = np.array([float(line[0]), float(line[1]), float(line[2])])
            c = [c]
            c.extend(line[3:])
            coords.append(c)

    return coords

def adjust_bond(coordinates, bond, distance):
    
    d = coordinates[bond[0]][0] - coordinates[bond[1]][0]
    d = d / np.linalg.norm(d)

    d = d * distance

    print(d)
    coordinates[bond[0]][0] = coordinates[bond[1]][0] + d
    return coordinates

def coord_to_xyz(coordinates):
    xyz_struct = []
    for line in coordinates:
        new_line = []
        new_line.append(line[1])
        ang_coords = line[0] * 0.529177
        new_line.append(ang_coords)
        xyz_struct.append(new_line)

    return xyz_struct

def scan_to_xyz(file_name="scan.xyz"):
    dirs = [d for d in os.listdir() if os.path.isdir(d)]
    dirs = sorted(dirs, key=lambda x: float(x.split("_")[-1]))

    structs = []
    for d in dirs:
        energy = TMcalculation.get_energy(energyfile=f"{d}/energy") 
        coords = read_coords(f"{d}/coord")

        structs.append([energy, coord_to_xyz(coords)])

    with open(file_name, "w") as output:
        for struct in structs:
            num_atoms=  len(struct[1])
            output.write(f"{num_atoms}\n")
            output.write(f" Energy =   {struct[0]}\n")
            for line in struct[1]:
                new_line = [line[0].capitalize(), str(line[1][0]), str(line[1][1]), str(line[1][2])]
                new_line = " ".join(new_line)
                output.write(f" {new_line}\n")

def scan_coordinates(submit=False, run=True):
    if not os.path.isfile("definput.json"):
        logger.error("No definput.json file!")
        raise FileNotFoundError("definput.json")

    with open("definput.json", 'r') as definput:
        parameters = json.load(definput)

    if "Scan" not in parameters["geometry"].keys():
        logger.error("No Scans specified")
        return
    
    coords = read_coords()

    if "Bond" in parameters["geometry"]["Scan"].keys():
        #scan the bond
        for bond in parameters["geometry"]["Scan"]["Bond"]:
            atoms = [bond[0]-1, bond[1]-1]

            steps = int(bond[2])
            step_size = bond[3]

            start_distance =  coords[atoms[0]][0] - coords[atoms[1]][0]
            start_distance =  np.linalg.norm(start_distance)

            parameters["geometry"]["idef"]["idef_on"] = True
            lbl = f"{atoms[0]+1},{atoms[1]+1}"
            if lbl not in parameters["geometry"]["idef"]["freeze_stretch"]:
                parameters["geometry"]["idef"]["freeze_stretch"].append(lbl)
            
            for i in range(steps):
                distance = start_distance + float(i)*float(step_size)
                coords = adjust_bond(coords, atoms, distance)

                dir_name = f"bond_{atoms[0]+1}_{atoms[1]+1}_{distance:.3f}"
                if not os.path.isdir(dir_name):
                    os.mkdir(dir_name)

                elif "GEO_OPT_CONVERGED" in os.listdir(dir_name):
                    continue
                
                else:
                    shutil.rmtree(dir_name)
                    os.mkdir(dir_name)

                os.chdir(dir_name)

                write_coords(coords)
                with open("definput.json", 'w') as param_file:
                    json.dump(parameters, param_file)
                

                if submit:
                    logger.warn("SUBMITTING OMG")

                elif run:
                    logger.warn("RUNNING")

                    tj = TMcalculation(cores=12)
                    #Good to just update now...
                    coords = read_coords()
                    
                else:
                    sj = setupTMjob()

                os.chdir("..")

            #Start over again
            coords = read_coords()

