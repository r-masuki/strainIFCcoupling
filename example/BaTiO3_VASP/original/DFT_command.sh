cp ../../original/VASP/INCAR INCAR
cp ../../original/VASP/KPOINTS KPOINTS
ln -s ../../original/VASP/POTCAR POTCAR

mpirun -n ${NUMPRO} $VASP > vasp.out