#
# generate_straindisp.py
#
# Script to make strained supercells and write DFT-supercells with atomic displacements.
#
# Copyright (c) 2023 Ryota Masuki
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#

from alm import ALM
import numpy as np
from ase.io import read, write
from ase.io.espresso import read_espresso_in, write_espresso_in, read_fortran_namelist
from ase import Atoms
from ase.cell import Cell

import json
import copy
import os
import warnings
import shutil
import argparse

from model_with_strain import ModelWithStrain
from read_input import *


parser = argparse.ArgumentParser()
parser.add_argument("-smag", "--strain_mag", help = "magnitude of the strain",
                    type = float)

parser.add_argument("-dmag", "--disp_mag", help = "magnitude of the atomic displacements [Ang]",
                    default = 0.01,
                    type = float)

parser.add_argument("--DFT", help = "the DFT engine",
                    choices = ["VASP", "QE"],
                    required = True)

parser.add_argument("--no_offset", help = "do not generate structure with strain but withou atomic displacements. can be used when the offsets of the force in the strained cells are zero.",
                    action = "store_true")

args = parser.parse_args()

smag = args.strain_mag
dmag = args.disp_mag


if(os.path.exists("./strain_modes.json")):
    with open("./strain_modes.json", "r") as json_file:
        json_object = json.load(json_file)

else:
    script_path = os.path.abspath(__file__)
    script_dir = os.path.dirname(script_path)
    with open(script_dir + "/strain_modes.json", "r") as json_file:
        json_object = json.load(json_file)

if args.DFT == "VASP":
    filename_in = "original/VASP/POSCAR"
    supercell = read(filename_in)
    dft_input = read_input(args, filename_in)

elif args.DFT == "QE":
    filename_in = "original/QE/pw.in"
    supercell = read_espresso_in(filename_in)
    dft_input = read_input(args, filename_in)



# print(supercell.get_cell()[:])

# print(type(supercell))

# print(json_object["strain_modes"])
for item in json_object["strain_modes"]:
    print(item)
    strain_cell = ModelWithStrain(item, supercell, dft_input, args)

    # print(strain_cell._supercell.get_cell().cellpar())
    # print(strain_cell._supercell.get_cell()[:])
    # print(strain_cell._supercell.get_scaled_positions())
    # print(type(strain_cell._supercell.get_positions()))

    strain_cell.suggest_displacements()
    strain_cell.generate_disp_supercells(dmag)

    strain_cell.prepare_DFT_inputs()






