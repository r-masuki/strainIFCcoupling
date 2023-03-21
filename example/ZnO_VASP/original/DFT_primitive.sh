cp ../../original/VASP_primitive/INCAR INCAR
cp ../../original/VASP_primitive/KPOINTS KPOINTS
ln -s ../../original/VASP_primitive/POTCAR POTCAR

mpirun -n ${NUMPRO} $VASP > vasp.out