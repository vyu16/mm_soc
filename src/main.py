#!/usr/bin/env python3

import numpy as np
from dimension import find_dim
from mulliken import parse_mulliken
from smear import gauss
from read_elsi import read_elsi_to_den
from plot import plot_data

dim_dict = find_dim()
n_ks = dim_dict["i_max"]-dim_dict["i_min"]+1

print("Parsing Mulliken.out",flush=True)

state_val,max_occ,state_is_org,kwt = parse_mulliken(dim_dict)

wid = 0.1
nw = 1000
dw = 0.01
w0 = 0.01
w1 = w0+(nw-1)*dw
xw = np.linspace(w0,w1,num=nw)

for i_dir in ["x","y","z"]:
    yw = np.zeros((nw,4))

    for i_kpt in range(dim_dict["n_kpt"]):
        filename = "moment_soc_"+i_dir+"_kpt_"+"{:06}".format(i_kpt+1)+".csc"

        print("Processing file "+filename,flush=True)

        mommat = read_elsi_to_den(filename)

        de = [[],[],[]]
        m2 = [[],[],[]]
        w = w0-dw

        for i1 in range(max_occ[i_kpt]):
            for i2 in range(max_occ[i_kpt],n_ks):
                tmp = state_val[i_kpt,i2]-state_val[i_kpt,i1]

                if tmp < w1+5.0:
                    mm = abs(mommat[i1,i2])

                    if state_is_org[i_kpt,i1] and state_is_org[i_kpt,i2]:
                        m2[0].append(mm*mm)
                        de[0].append(tmp)
                    elif state_is_org[i_kpt,i1] or state_is_org[i_kpt,i2]:
                        m2[1].append(mm*mm)
                        de[1].append(tmp)
                    else:
                        m2[2].append(mm*mm)
                        de[2].append(tmp)

        for iw in range(nw):
            w += dw

            for i1 in range(3):
                yw[iw,i1] = sum([gauss(w,de[i1][i2],wid)*m2[i1][i2] for i2 in range(len(de[i1]))])
                yw[iw,i1] *= kwt[i_kpt]

        yw[:,3] = [sum(yw[iw,:3]) for iw in range(nw)]

    plot_data(i_dir,xw,yw)
