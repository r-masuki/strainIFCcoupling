#
#  model_with_strain.py
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
    
    def __init__(self, id, umn, supercell, args):
        self._id = id
        self._umn = umn
        self._Fmn = np.eye(3) + umn

        self._supercell = copy.deepcopy(supercell)
        # apply deformation
        lavec_tmp = np.dot(supercell.get_cell()[:], np.transpose(self._Fmn))
        self._supercell.set_cell(Cell.ascell(lavec_tmp))
        pos_tmp = np.dot(supercell.get_positions(), np.transpose(self._Fmn))
        self._supercell.set_positions(pos_tmp)

        self._ctrlargs = args

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
        "Getter of supercell."
        return self._supercell
    
    @property
    def n_disp(self):
        "Getter of the number of patterns of displacements."
        return self._n_disp

    def _gen_alm_crystal(self):
        BOHR_IN_AA = 0.52917721067
        lavec = self.supercell.get_cell() / BOHR_IN_AA
        xcoord = self.supercell.get_scaled_positions()
        numbers = self.supercell.get_atomic_numbers()
        crystal = (lavec, xcoord, numbers)

        return crystal
    
    def suggest_displacements(self):
        crystal = self._gen_alm_crystal()

        with ALM(*crystal) as alm:
            maxorder = 1
            num_elems = np.size(np.unique(self.supercell.get_atomic_numbers()))
            cutoff_radius = np.full([maxorder, num_elems, num_elems], -1, dtype=int)
            alm.define(maxorder, cutoff_radius, symmetrization_basis="Lattice")
            alm.set_verbosity(0)
            alm.suggest()

            self._disp_patterns = alm.get_displacement_patterns(alm.maxorder)

            self._n_disp = len(self._disp_patterns)

    def _gen_alm_dfset(self, fname):

        dfset = np.loadtxt(fname).reshape((-1, len(self._supercell), 6))
        force = dfset[:, :, 3:]
        disp = dfset[:, :, :3]

        return disp, force
        
    def get_IFCs(self):
        crystal = self._gen_alm_crystal()
        disp_harm, force_harm = self._gen_alm_dfset("DFSETS/DFSET_harmonic_"
                                            + "{:0>{}}".format(self._id, 3))

        with ALM(*crystal) as alm:
            maxorder = 1
            num_elems = np.size(np.unique(self.supercell.get_atomic_numbers()))
            cutoff_radius = np.full([maxorder, num_elems, num_elems], -1, dtype=int)
            alm.define(maxorder, cutoff_radius, symmetrization_basis="Lattice")
            alm.displacements = disp_harm
            alm.forces = force_harm

            optcontrol = {'linear_model': 1}
            alm.set_optimizer_control(optcontrol)
            alm.optimize()

            if not os.path.exists("results"):
                # If it doesn't exist, create it
                os.mkdir("results")
            alm.save_fc(filename="results/strain_" + "{:0>{}}".format(self._id, 3) + ".xml", 
                        format="alamode")


    def generate_disp_supercells(self, magnitude_disp):

        self._disp_supercells = []

        for disp_pattern in self._disp_patterns:
            disp_supercell = copy.deepcopy(self._supercell)
            pos_tmp = disp_supercell.get_positions()
            for atom_disp in disp_pattern:
                pos_tmp[atom_disp[0]] += atom_disp[1] * magnitude_disp
            
            disp_supercell.set_positions(pos_tmp)

            self._disp_supercells.append(disp_supercell)

    def prepare_DFT_inputs(self):

        # prepare working directory
        workdir_name = "strain_" + "{:0>{}}".format(self._id, 3)
        if not os.path.exists(workdir_name):
            os.mkdir(workdir_name)
        else:
            warnings.warn("The directory already exists.")

        # prepare DFT inputs
        num_width = len(str(self._n_disp)) + 1

        for i_disp in range(self._n_disp):

            dispdir_name = workdir_name + "/disp_" + "{:0>{}}".format(i_disp+1, num_width)

            if not os.path.exists(dispdir_name):
                # If it doesn't exist, create it
                os.mkdir(dispdir_name)
            else:
                warnings.warn("The directory already exists.")


            if(self._ctrlargs.DFT == "VASP"):
                write(dispdir_name + "/POSCAR", self._disp_supercells[i_disp])


        if(self._ctrlargs.DFT == "VASP"):
            write(workdir_name + "/POSCAR_no_disp", self._supercell)

        if(not self._ctrlargs.no_offset):
            nodispdir_name = workdir_name + "/nodisp"

            if not os.path.exists(nodispdir_name):
                # If it doesn't exist, create it
                os.mkdir(nodispdir_name)
            else:
                warnings.warn("The directory already exists.")

            if(self._ctrlargs.DFT == "VASP"):
                shutil.copy(workdir_name + "/POSCAR_no_disp", nodispdir_name + "/POSCAR")

            
        # prepare a jobscript
        if os.path.exists("original/job.sh"):
            shutil.copy("original/job.sh", workdir_name + "/job.sh")

        DFT_commands = []
        with open('original/DFT_command.sh') as f:
            for line in f:
                DFT_commands.append("  " + line)
            DFT_commands.append("\n")

        with open(workdir_name + "/job.sh", "a") as f:

            if (not self._ctrlargs.no_offset):
                f.write("\n\n")
                f.write("cd nodisp\n")

                for line in DFT_commands:
                    f.write(line)
                f.write("cd ..")

            f.write("\n\n")
            f.write("for i_disp in ")
            for i_disp in range(self._n_disp):
                f.write("{:0>{}} ".format(i_disp+1, num_width))
            f.write("\ndo\n\n")

            f.write("  cd disp_${i_disp}\n\n")
            for line in DFT_commands:
                f.write(line)

            f.write("  cd ..\n\n")
            f.write("done")

        # prepare a script to make DFSET
        if os.path.exists("original/extract.sh"):
            shutil.copy("original/extract.sh", workdir_name + "/extract.sh")

        with open(workdir_name + "/extract.sh", "a") as f:
            if(not self._ctrlargs.no_offset):
                f.write("\n\n")
                if(self._ctrlargs.DFT == "VASP"):
                    f.write("cp nodisp/vasprun.xml vasprun0.xml")

            f.write("\n\n")
            f.write("for i_disp in ")
            for i_disp in range(self._n_disp):
                f.write("{:0>{}} ".format(i_disp+1, num_width))
            f.write("\ndo\n\n")
            f.write("  cd disp_${i_disp}\n")
            if(self._ctrlargs.DFT == "VASP"):
                f.write("  cp vasprun.xml ../vasprun_${i_disp}.xml\n")
            f.write("  cd ..\n\n")

            f.write("done\n\n")

            if(self._ctrlargs.no_offset):
                if(self._ctrlargs.DFT == "VASP"):
                    f.write("python3 ${ALAMODE_TOOLS}/extract.py --VASP=POSCAR_no_disp vasprun_*.xml > DFSET_harmonic" + "\n\n")

            else:
                if(self._ctrlargs.DFT == "VASP"):
                    f.write("python3 ${ALAMODE_TOOLS}/extract.py --VASP=POSCAR_no_disp --offset vasprun0.xml vasprun_*.xml > DFSET_harmonic" + "\n\n")

            f.write("mkdir -p ../DFSETS\n")
            f.write("cp DFSET_harmonic ../DFSETS/DFSET_harmonic_" + "{:0>{}}\n\n".format(self._id, 3))

