from ase.io.espresso import read_espresso_in, write_espresso_in, read_fortran_namelist
import sys
from collections import OrderedDict

def read_input(args, filename):

    if args.DFT == "VASP":
        return [{}, {}, {}, False]

    if args.DFT == "QE":
        namelist = {}
        otherlist = []
        with open(filename, 'r') as f:
            obj = read_fortran_namelist(f)
            for entry in obj:
                if isinstance(entry, OrderedDict):
                    namelist = entry
                else:
                    otherlist.extend(entry)
        # Needed to flattern the nested dictionary into 1d
        namelist_flattern = {}
        for key in namelist.keys():
            for key2, val2 in namelist[key].items():
                namelist_flattern[key2] = val2

        # get pseudopotential data
        ntype = namelist_flattern['ntyp']
        pseudo_dic = {}
        kpts = [0, 0, 0]
        koffset = False

        for i, entry in enumerate(otherlist):
            if "ATOMIC_SPECIES" in entry:
                try:
                    for j in range(ntype):
                        line = otherlist[i + j + 1].split()
                        pseudo_dic[line[0]] = line[2]
                except:
                    raise RuntimeError

            if "K_POINTS" in entry:
                try:
                    if(not ("automatic" in entry)):
                        raise RuntimeError("K_POINTS must be used with automatic.")

                    line = otherlist[i+1].split()
                    for j in range(3):
                        kpts[j] = int(line[j])
                    
                    if (int(line[3]) == int(line[4]) == int(line[5]) == 1):
                        koffset = True
                    elif (int(line[3]) == int(line[4]) == int(line[5]) == 0):
                        koffset = False
                    else:
                        raise RuntimeError("Offset in KPOINTS must be 0 0 0 or 1 1 1.")
                except:
                    raise RuntimeError

            

        if len(pseudo_dic) != ntype:
            raise RuntimeError

        return [namelist_flattern, pseudo_dic, kpts, koffset]
