# strainIFCcoupling

This is a script to calculate strain-IFC coupling by finite-displacement method with respect to the strain.
The obtained strain-IFC coupling constants can be used in structural optimization at finite temperatures implemented in the ALAMODE package[1-3].

In calculating the coupling between strain and harmonic IFCs

$\frac{\partial \Phi_{\mu_1 \mu_2}({0}\alpha_1, {R}\alpha_2)}{\partial u_{\mu \nu}}$,

we use the frozen-phonon calculation in strained supercells using the ALM package.

The strain-force coupling

$\frac{\partial \Phi_{\mu_1}({0}\alpha_1)}{\partial u_{\mu \nu}}$

are obtained by calculating the atomic forces in strained primitive cells using external DFT engines.

## Usage

We support interface with VASP and Quantum Espresso as external DFT packages, which is used to calculate atomic forces for strained supercells and primitive cells.

### strain-harmonic coupling (VASP)
---
#### Preparation

The directory structures and file names must be the same as those in `example/BaTiO3_VASP`. Please make `original` directory and prepare the following inputs. 

* `DFT_command.sh` : The shell command to run a VASP calculation in each directory. We recommend making INCAR, KPOINTS and POTCAR here because they are not generated in the python script.
* `job.sh` : The header of the shell script for the VASP calculations.
* `extract.sh` : The header of the shell script to extract DFSET. Please define ALAMODE_TOOLS.
* `VASP` : VASP input files without strain or atomic displacements. 

Note that `VASP/POTCAR` is left blank due to the license issues.


#### How to run

From the command line, please run
```
> python3 generate_straindisp.py --DFT=VASP -smag=0.005 -dmag=0.01
```
`-smag` is the magnitude of the strain. `-dmag` defines the magnitude of atomic displacements in Angstrom. You can also use `--no_offset` option if the atomic forces are zero when strain is applied, which condition is not checked in the script.

Here, six directories of different strain patterns will be generated.
If you run the VASP calculation in supercomputers or cluster computers, copy the whole working folder (`BaTiO3_VASP` in this case).

In each of `strain_001`, ..., `strain_006`, run VASP calculation and extract displacement-force data.

```
> qsub job.sh
> bash extract.bash
```

The displacement-force data will be stored in `DFSETS` directory. Copy `DFSETS` direcory if you continue with your local environment.

Run the python script to get the input of the strain-harmonic coupling.
```
> python3 get_strainIFC.py --DFT=VASP -smag=0.005
```
Note that the value of `-smag` need to be the same as in `generate_straindisp.py`, which is not checked in the script.

### strain-harmonic coupling (Quantum Espresso)
---

#### Preparation

The directory structures and file names must be the same as those in `example/ZnO_QE`. Please make `original` directory and prepare the following inputs. 

* `DFT_command.sh` : The shell command to run a QE calculation in each directory. Note that the name of the output file must be `pw.out`.
* `job.sh` : The header of the shell script for the QE calculations.
* `extract.sh` : The header of the shell script to extract DFSET. Please define ALAMODE_TOOLS.
* `QE` : QE input files without strain or atomic displacements. The name of the input file need to be `pw.in` and `ibrav` need to be set as zero.

#### How to run

The flow of the calculation is the same as in the previous section.
Please replace `--DFT=VASP` by `--DFT=QE` when running the python script.

### strain-force coupling (VASP)
---
#### Preparation

The directory structures and file names must be the same as those in `example/ZnO_VASP`. Please make `original` directory and prepare the following inputs. 

* `DFT_primitive.sh` : The shell command to run a VASP calculation in each directory.
* `job.sh` : The header of the shell script for the VASP calculations.
* `extract.sh` : The header of the shell script to extract DFSET. Please define ALAMODE_TOOLS.
* `VASP_primitive` : VASP input files without strain or atomic displacements. 

Note that `VASP_primitive/POTCAR` is left blank due to the license issues.

#### How to run

From the command line, please run
```
> python3 generate_strainprim.py --DFT=VASP -smag=0.005 
```
`-smag` is the magnitude of the strain.

Here, six directories of different strain patterns will be generated.
If you run the VASP calculation in supercomputers or cluster computers, copy the whole working folder (`ZnO_VASP` in this case).

In each of `strain_001`, ..., `strain_006`, run VASP calculation and extract displacement-force data.

```
> qsub job_prim.sh
> bash extract_prim.bash
```

The displacement-force data will be stored in `DFSETS_primitive` directory. Copy `DFSETS_primitive` direcory if you continue with your local environment.

Run the python script to get the input of the strain-harmonic coupling.
```
> python3 get_strainForce.py --DFT=VASP -smag=0.005
```
Note that the value of `-smag` need to be the same as in `generate_strainprim.py`, which is not checked in the script. The resultant strain-force coupling is written in `results/strain_force.in`.

### strain-force coupling (Quantum Espresso)
---

#### Preparation

The directory structures and file names must be the same as those in `example/ZnO_QE`. Please make `original` directory and prepare the following inputs. 

* `DFT_primitive.sh` : The shell command to run a QE calculation in each directory. The name of the output file must be `pw.out`
* `job.sh` : The header of the shell script for the QE calculations.
* `extract.sh` : The header of the shell script to extract DFSET. Please define ALAMODE_TOOLS.
* `QE_primitive` : VASP input files without strain or atomic displacements. The name of the input file need to be `pw.in` and `ibrav` need to be set as zero.

#### How to run

The flow of the calculation is the same as in the previous section.
Please replace `--DFT=VASP` by `--DFT=QE` when running the python script.

## Dependencies

- ALM [[link](https://github.com/ttadano/ALM)]
- ase
- numpy

Since we use the python interface of ALM, a user need to get ALM [[link](https://github.com/ttadano/ALM)], which is provided separately from the ALAMODE package. Please build ALM using conda, which is explained in [[link](https://alm.readthedocs.io/en/develop/compile-with-conda-packages.html#building-alm-using-conda)].

## References

[1] R. Masuki, T. Nomoto, R. Arita, T. Tadano. "Ab initio structural optimization at finite temperatures based on anharmonic phonon theory: Application to the structural phase transitions of BaTiO3" [Physical Review B 106, 224104 (2022)](https://doi.org/10.1103/PhysRevB.106.224104)

[2] R. Masuki, T. Nomoto, R. Arita, T. Tadano. "Full optimization of quasiharmonic free energy with anharmonic lattice model: Application to thermal expansion and pyroelectricity of wurtzite GaN and ZnO" [arXiv:2302.04537 (2023)](https://doi.org/10.1103/PhysRevB.106.224104)

[3] T. Tadano, Y. Gohda, S. Tsuneyuki, "Anharmonic force constants extracted from first-principles molecular dynamics: applications to heat transfer simulations" [J. Phys.: Condens. Matter 26, 225402 (2014)](http://iopscience.iop.org/0953-8984/26/22/225402/)