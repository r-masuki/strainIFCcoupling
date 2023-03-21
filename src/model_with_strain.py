#
# model_with_strain.py
#
# Define the class to apply strain to a supercell and run ALM calculations.
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

from read_input import *

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
    
    def __init__(self, json_item, supercell, dft_input, args):
        self._id = json_item["id"]

        smag = 0.0
        if args.strain_mag:
            smag = args.strain_mag
        else:
            smag = json_item["strain_mag"]
        umn = smag*get_strain_mode(json_item["mode"])

        self._smag = smag
        self._umn = umn
        self._Fmn = np.eye(3) + umn

        self._supercell = copy.deepcopy(supercell)
        # apply deformation
        lavec_tmp = np.dot(supercell.get_cell()[:], np.transpose(self._Fmn))
        self._supercell.set_cell(Cell.ascell(lavec_tmp))
        pos_tmp = np.dot(supercell.get_positions(), np.transpose(self._Fmn))
        self._supercell.set_positions(pos_tmp)

        self._weight = json_item["weight"]
        self._mode = json_item["mode"]

        self._dft_input = dft_input
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

    @property
    def mode(self):
        "Getter of the mode of the strain."
        return self._mode

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
    
    def write_strain_harmonic_in(self):
        with open("results/strain_harmonic.in", "a") as f:
            f.write("{0:4s} {1:25.15f} {2:25.15f}".format(self._mode, self._smag, self._weight))
            f.write(" {:25s}\n".format("strain_" + "{:0>{}}".format(self._id, 3) + ".xml"))

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

            elif(self._ctrlargs.DFT == "QE"):
                with open(dispdir_name + "/pw.in", 'w') as f:
                    write_espresso_in(f, self._disp_supercells[i_disp], 
                                      **self._dft_input,
                                      crystal_coordinates=True)

        if(self._ctrlargs.DFT == "VASP"):
            write(workdir_name + "/POSCAR_no_disp", self._supercell)

        elif(self._ctrlargs.DFT == "QE"):
            with open(workdir_name + "/pw.no_disp.in", 'w') as f:
                write_espresso_in(f, self._supercell,
                                  **self._dft_input,
                                  crystal_coordinates=True)

        if(not self._ctrlargs.no_offset):
            nodispdir_name = workdir_name + "/nodisp"

            if not os.path.exists(nodispdir_name):
                # If it doesn't exist, create it
                os.mkdir(nodispdir_name)
            else:
                warnings.warn("The directory already exists.")

            if(self._ctrlargs.DFT == "VASP"):
                shutil.copy(workdir_name + "/POSCAR_no_disp", nodispdir_name + "/POSCAR")

            if(self._ctrlargs.DFT == "QE"):
                shutil.copy(workdir_name + "/pw.no_disp.in", nodispdir_name + "/pw.in")
            
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
                elif(self._ctrlargs.DFT == "QE"):
                    f.write("cp nodisp/pw.out pw.no_disp.out")

            f.write("\n\n")
            f.write("for i_disp in ")
            for i_disp in range(self._n_disp):
                f.write("{:0>{}} ".format(i_disp+1, num_width))
            f.write("\ndo\n\n")
            f.write("  cd disp_${i_disp}\n")

            if(self._ctrlargs.DFT == "VASP"):
                f.write("  cp vasprun.xml ../vasprun_${i_disp}.xml\n")
            elif(self._ctrlargs.DFT == "QE"):
                f.write("  cp pw.out ../pw.disp_${i_disp}.out\n")
            f.write("  cd ..\n\n")

            f.write("done\n\n")

            if(self._ctrlargs.no_offset):
                if(self._ctrlargs.DFT == "VASP"):
                    f.write("python3 ${ALAMODE_TOOLS}/extract.py --VASP=POSCAR_no_disp vasprun_*.xml > DFSET_harmonic" + "\n\n")
                elif(self._ctrlargs.DFT == "QE"):
                    f.write("python3 ${ALAMODE_TOOLS}/extract.py --QE=pw.no_disp.in pw.disp_*.out > DFSET_harmonic" + "\n\n")

            else:
                if(self._ctrlargs.DFT == "VASP"):
                    f.write("python3 ${ALAMODE_TOOLS}/extract.py --VASP=POSCAR_no_disp --offset vasprun0.xml vasprun_*.xml > DFSET_harmonic" + "\n\n")
                elif(self._ctrlargs.DFT == "QE"):
                    f.write("python3 ${ALAMODE_TOOLS}/extract.py --QE=pw.no_disp.in --offset pw.no_disp.out pw.disp_*.out > DFSET_harmonic" + "\n\n")

            f.write("mkdir -p ../DFSETS\n")
            f.write("cp DFSET_harmonic ../DFSETS/DFSET_harmonic_" + "{:0>{}}\n\n".format(self._id, 3))


    def prepare_DFT_primitive(self):

        # prepare working directory
        workdir_name = "strain_" + "{:0>{}}".format(self._id, 3)
        if not os.path.exists(workdir_name):
            os.mkdir(workdir_name)
        else:
            warnings.warn("The directory already exists.")

        # prepare DFT inputs

        primdir_name = workdir_name + "/primitive" 

        if not os.path.exists(primdir_name):
            # If it doesn't exist, create it
            os.mkdir(primdir_name)
        else:
            warnings.warn("The directory already exists.")

        if(self._ctrlargs.DFT == "VASP"):
            write(primdir_name + "/POSCAR", self._supercell)

        elif(self._ctrlargs.DFT == "QE"):
            with open(primdir_name + "/pw.in", 'w') as f:
                write_espresso_in(f, self._supercell,
                                  **self._dft_input,
                                  crystal_coordinates=True)

        # prepare a jobscript
        if os.path.exists("original/job.sh"):
            shutil.copy("original/job.sh", workdir_name + "/job_prim.sh")

        DFT_commands = []
        with open('original/DFT_primitive.sh') as f:
            for line in f:
                DFT_commands.append("  " + line)
            DFT_commands.append("\n")

        with open(workdir_name + "/job_prim.sh", "a") as f:

            f.write("\n\n")
            f.write("cd primitive\n")

            for line in DFT_commands:
                f.write(line)
            f.write("cd ..")

        # prepare a script to make DFSET
        if os.path.exists("original/extract.sh"):
            shutil.copy("original/extract.sh", workdir_name + "/extract_prim.sh")

        with open(workdir_name + "/extract_prim.sh", "a") as f:
            f.write("\n\n")
            f.write("cd primitive\n\n")
            # if(self._ctrlargs.DFT == "VASP"):
            #     f.write("cp nodisp/vasprun.xml vasprun.prim.xml")
            # elif(self._ctrlargs.DFT == "QE"):
            #     f.write("cp nodisp/pw.out pw.prim.out")

            # f.write("\n\n")
            # f.write("for i_disp in ")
            # for i_disp in range(self._n_disp):
            #     f.write("{:0>{}} ".format(i_disp+1, num_width))
            # f.write("\ndo\n\n")
            # f.write("  cd disp_${i_disp}\n")
# 
            # if(self._ctrlargs.DFT == "VASP"):
            #     f.write("  cp vasprun.xml ../vasprun_${i_disp}.xml\n")
            # elif(self._ctrlargs.DFT == "QE"):
            #     f.write("  cp pw.out ../pw.disp_${i_disp}.out\n")
            # f.write("  cd ..\n\n")

            # f.write("done\n\n")

            if(self._ctrlargs.DFT == "VASP"):
                f.write("python3 ${ALAMODE_TOOLS}/extract.py --VASP=POSCAR vasprun.xml > DFSET_primitive" + "\n\n")
            elif(self._ctrlargs.DFT == "QE"):
                f.write("python3 ${ALAMODE_TOOLS}/extract.py --QE=pw.in pw.out > DFSET_primitive" + "\n\n")

            f.write("mkdir -p ../../DFSETS_primitive\n")
            f.write("cp DFSET_primitive ../../DFSETS_primitive/DFSET_primitive_" + "{:0>{}}\n\n".format(self._id, 3))

    def write_strain_force_in(self):

        RYBOHR_TO_EVANG = 13.60569301/0.52917721067
        force = [0.0, 0.0, 0.0]

        if not os.path.exists("results"):
            # If it doesn't exist, create it
            os.mkdir("results")

        with open("results/strain_force.in", "a") as fout:
            fout.write("{0:4s} {1:25.15f} {2:25.15f}\n".format(self._mode, self._smag, self._weight))


            with open("DFSETS_primitive/DFSET_primitive_" + "{:0>{}}".format(self._id, 3)) as fin:
                line = fin.readline()

                line = fin.readline()
                while line:                
                    arr = line.rstrip().split()
                    for i in range(3):
                        force[i] = float(arr[i+3]) * RYBOHR_TO_EVANG
                        

                    fout.write("{0:25.15f} {1:25.15f} {2:25.15f}\n".format(force[0], force[1], force[2]))
                    line = fin.readline()

