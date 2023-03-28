#!/bin/sh
#PBS -l ...
#PBS -q ...

VASPDIR=...
VASP=${VASPDIR}/vasp_std
NUMPRO=16

cd $PBS_O_WORKDIR

RUN_DFT_CALCULATION