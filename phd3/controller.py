#!/usr/bin/env python3

import logging
import iteration
import os
import json

class controller:

    def __init__(self):

        #Could create a list of functions that get called in sequential ordering of some sort? then save a numerical index of what step we are on maybe

        # We need to find the phd_control.json file and load that in
        # Then we need to see where we are in the calculations and set a variable to that step
        # Finally, we need to make sure that we delete any files that were created on a failed run if that is the step we are now performing

        # need to have it pass the parameters to the appropriate iteration when it so needs to!!!!!!
        # This will be the brains of the operation and spawn the iterations as needed and then abort when it is done
        pass


if __name__ == "__main__":


    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


    par = {"last pdb": os.path.abspath("./HG3.pdb"),
            "DMD CONVERGE": True,
            "MAX DMD STEPS": 5,
            "QM Chop":{
                "Residues": ["A:49", "A:126", "A:43-A:46", "A:264"]
                }
            }

    with open("dmdinput.json", 'r') as f:
        par["dmd params"]= json.load(f)

    with open("definput.json", 'r') as f:
        par["qm params"] = json.load(f)

    i = iteration.iteration( "iter_1", ".", par, 1)
    i.continue_calculation()
