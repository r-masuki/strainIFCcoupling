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

parser.add_argument("-dmag", "--disp_mag", help = "magnitude of the atomic displacements [Ang]",
                    default = 0.01,
                    type = float)

parser.add_argument("--DFT", help = "the DFT engine",
                    choices = ["VASP"],
                    required = True)

parser.add_argument("--no_offset", help = "do not generate structure with strain but withou atomic displacements. can be used when the offsets of the force in the strained cells are zero.",
                    action = "store_true")

parser.add_argument("--copy_potcar", help = "whether to copy POTCAR in the script. used when --DFT = VASP.",
                    action = "store_true")

args = parser.parse_args()

smag = args.strain_mag
dmag = args.disp_mag

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)
json_file = open(script_dir + '/strain_modes.json', 'r')
json_object = json.load(json_file)

print(json_object["strain_modes"][0]["id"])

supercell = read("original/VASP/POSCAR")

# print(supercell.get_cell()[:])

# print(type(supercell))

# print(json_object["strain_modes"])
for item in json_object["strain_modes"]:
    print(item)
    strain_cell = ModelWithStrain(item["id"], smag*np.array(item["mode"]), supercell, args)

    # print(strain_cell._supercell.get_cell().cellpar())
    # print(strain_cell._supercell.get_cell()[:])
    # print(strain_cell._supercell.get_scaled_positions())
    # print(type(strain_cell._supercell.get_positions()))

    strain_cell.suggest_displacements()
    strain_cell.generate_disp_supercells(dmag)

    strain_cell.prepare_DFT_inputs()






