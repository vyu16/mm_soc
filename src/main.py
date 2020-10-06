#!/usr/bin/env python3

import sys
import numpy as np
from dimension import find_dim
from mulliken import parse_mulliken
from smear import gauss
from read_elsi import read_elsi_to_den
from plot import plot_all

dim_dict = find_dim()
n_ks = dim_dict["i_max"]-dim_dict["i_min"]+1

print("Parsing Mulliken.out")

state_val,max_occ,state_is_org,kwt = parse_mulliken(dim_dict)

wid = 0.1
nw = 2000
dw = 0.005
w0 = 0.001
w1 = w0+(nw-1)*dw
xw = np.linspace(w0,w1,num=nw)
yw = np.zeros((nw,4))

#for i_kpt in range(dim_dict["n_kpt"]):
for i_kpt in range(1):
    filename = "moment_soc_x_kpt_"+"{:06}".format(i_kpt+1)+".csc"

    print("Processing file "+filename)

    mommat = read_elsi_to_den(filename)

    de = []
    m2 = []
    tt = []
    w = w0

    for i1 in range(max_occ[i_kpt]):
        for i2 in range(max_occ[i_kpt],n_ks):
            if state_is_org[i_kpt,i1] and state_is_org[i_kpt,i2]:
                t123 = 0
            elif not (state_is_org[i_kpt,i1]) and (not state_is_org[i_kpt,i2]):
                t123 = 2
            else:
                t123 = 1

            mm = abs(mommat[i1,i2])
            m2.append(mm*mm)
            de.append(state_val[i_kpt,i2]-state_val[i_kpt,i1])
            tt.append(t123)

    for iw in range(nw):
        for i1 in range(len(de)):
            yw[iw,tt[i1]] += kwt[i_kpt]*gauss(w,de[i1],wid)*m2[i1]
            yw[iw,3] += kwt[i_kpt]*gauss(w,de[i1],wid)*m2[i1]

        w += dw

plot_all(xw,yw)
