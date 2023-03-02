#
#  generate_straindisp.py
#
#  Generate supercells with atomic displacements
#  and other DFT inputs for strain-harmonic IFC coupling.
#

from alm import ALM
import numpy as np
from ase.io import read
from ase import Atoms
from ase.cell import Cell


import json
import copy


class ModelWithStrain:
    """Calculate harmonic and anharmonic interatomic force constants.

    Attributes
    ----------
    id : int
        id of the model
    umn : ndarray
        Tensor u_{mu nu} that describes strain.
        shape=(3, 3), dtype='double'
    Fmn : ndarray
        Deformation gradient tensor F_{mu nu} = I_{mu nu} + u_{mu nu}
        shape=(3, 3), dtype='double'
    supercell : ase.atoms.Atoms
        Strained supercell.
    """
    
    def __init__(self, id, umn, supercell):
        self._id = id
        self._umn = umn
        self._Fmn = np.eye(3) + umn

        self._supercell = copy.deepcopy(supercell)
        # apply deformation
        lavec_tmp = np.dot(supercell.get_cell()[:], np.transpose(self._Fmn))
        self._supercell.set_cell(Cell.ascell(lavec_tmp))
        pos_tmp = np.dot(supercell.get_positions(), np.transpose(self._Fmn))
        self._supercell.set_positions(pos_tmp)

    
    @property
    def id(self):
        """Getter of id."""
        return self._id

    @property
    def umn(self):
        """Getter of umn."""
        return self._umn

    @property
    def Fmn(self):
        "Getter of Fmn (deformation gradient tensor)."
        return self._Fmn

    @property
    def supercell(self):
        "Getter of supercell."]
        return self._supercell



eta = 0.005

json_file = open('strain_modes.json', 'r')
json_object = json.load(json_file)

print(json_object["strain_modes"][0]["id"])

supercell = read("../temp/POSCAR_BTO222")

print(supercell.get_cell()[:])

print(type(supercell))

print(json_object["strain_modes"])
for item in json_object["strain_modes"]:
    print(item)
    strain_cell = ModelWithStrain(item["id"], eta*np.array(item["mode"]), supercell)

    print(strain_cell._supercell.get_cell().cellpar())
    print(strain_cell._supercell.get_cell()[:])
    print(strain_cell._supercell.get_scaled_positions())
    print(type(strain_cell._supercell.get_positions()))

print(type(strain_cell._supercell))