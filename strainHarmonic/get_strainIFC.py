#
#  generate_straindisp.py
#
#  Generate supercells with atomic displacements
#  and other DFT inputs for strain-harmonic IFC coupling.
#

from alm import ALM
import numpy as np
from ase.io import read
from ase.io import write
from ase import Atoms
from ase.cell import Cell

import json
import copy
import os
import warnings
import shutil
import argparse

from model_with_strain import ModelWithStrain

parser = argparse.ArgumentParser()
parser.add_argument("-smag", "--strain_mag", help = "magnitude of the strain",
                    default = 0.005,
                    type = float)

args = parser.parse_args()

smag = args.strain_mag

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)
json_file = open(script_dir + '/strain_modes.json', 'r')
json_object = json.load(json_file)

print(json_object["strain_modes"][0]["id"])

supercell = read("original/VASP/POSCAR")

for item in json_object["strain_modes"]:
    print(item)
    strain_cell = ModelWithStrain(item["id"], smag*np.array(item["mode"]), supercell, None, args)

    strain_cell.get_IFCs()






