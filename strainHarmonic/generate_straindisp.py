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

from model_with_strain import ModelWithStrain

eta = 0.005

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
    strain_cell = ModelWithStrain(item["id"], eta*np.array(item["mode"]), supercell)

    # print(strain_cell._supercell.get_cell().cellpar())
    # print(strain_cell._supercell.get_cell()[:])
    # print(strain_cell._supercell.get_scaled_positions())
    # print(type(strain_cell._supercell.get_positions()))

    strain_cell.suggest_displacements()
    strain_cell.generate_disp_supercells(0.01)

    strain_cell.prepare_DFT_inputs()






