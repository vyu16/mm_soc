import numpy as np

def find_org_atom(n_atom):
    is_org = np.zeros(n_atom,dtype=np.int8)
    i_atom = -1

    with open("geometry.in","r") as f:
        for line in f:
            t = line.split()

            if len(t) > 0:
                if t[0] == "atom" or t[0] == "atom_frac":
                    i_atom += 1

                    if t[4] != "Pb" and t[4] != "I":
                        is_org[i_atom] = 1

    return is_org

def parse_mulliken(dim_dict):
    atom_is_org = find_org_atom(dim_dict["n_atom"])
    n_kpt = dim_dict["n_kpt"]
    i_max = dim_dict["i_max"]
    i_min = dim_dict["i_min"]
    n_state = i_max-i_min+1
    i_atom = -1

    kwt = np.zeros((n_kpt))
    val = np.zeros((n_kpt,n_state))
    tmp = np.zeros((n_kpt,n_state))
    is_org = np.zeros((n_kpt,n_state),dtype=np.int8)
    is_occ = np.zeros((n_kpt,n_state),dtype=np.int8)

    with open("Mulliken.out","r") as f:
        for line in f:
            t = line.split()

            if len(t) > 0:
                if t[0] == "Atom":
                    i_atom += 1
                    i_kpt = -1
                elif t[0] == "k":
                    i_kpt += 1

                    if i_atom == 0:
                        kwt[i_kpt] = float(t[10])
                elif t[0] == "Spin":
                    i_kpt = -1
                elif t[0] == "State":
                    continue
                elif t[0] == "#":
                    continue
                elif i_min <= int(t[0]) <= i_max:
                    if atom_is_org[i_atom]:
                        tmp[i_kpt,int(t[0])-i_min] += float(t[3])

                    if i_atom == 0:
                        val[i_kpt,int(t[0])-i_min] = float(t[1])

                        if float(t[2]) >= 0.5:
                            is_occ[i_kpt,int(t[0])-i_min] = 1

    for i_kpt in range(n_kpt):
        for i_state in range(n_state):
            if tmp[i_kpt,i_state] > 0.8:
                is_org[i_kpt,i_state] = 1

    return is_org,is_occ,val,kwt
