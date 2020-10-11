#!/usr/bin/env python3

import numpy as np
from mulliken import parse_mulliken
from tool import find_dim,gauss,save_and_plot
from read_elsi import read_elsi_to_den

dim_dict = find_dim()
n_ks = dim_dict["i_max"]-dim_dict["i_min"]+1

print("Parsing Mulliken.out",flush=True)

state_val,n_occ,state_is_org,kwt = parse_mulliken(dim_dict)

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

        # [[org-org], [org-inorg], [inorg-inorg]]
        de = [[],[],[]]
        m2 = [[],[],[]]
        w = w0-dw

        for i1 in range(n_occ[i_kpt]):
            for i2 in range(n_occ[i_kpt],n_ks):
                tmp = state_val[i_kpt,i2]-state_val[i_kpt,i1]

                if tmp < w1+wid*10:
                    mm = abs(mommat[i1,i2])

                    # org-org
                    if state_is_org[i_kpt,i1] and state_is_org[i_kpt,i2]:
                        de[0].append(tmp)
                        m2[0].append((mm/tmp)**2)
                    # org-inorg
                    elif state_is_org[i_kpt,i1] or state_is_org[i_kpt,i2]:
                        de[1].append(tmp)
                        m2[1].append((mm/tmp)**2)
                    # inorg-inorg
                    else:
                        de[2].append(tmp)
                        m2[2].append((mm/tmp)**2)

        for iw in range(nw):
            w += dw

            for i1 in range(3):
                tmp = 0.0

                for i2 in range(len(de[i1])):
                    if abs(w-de[i1][i2]) < wid*10:
                        tmp += gauss(w,de[i1][i2],wid)*m2[i1][i2]

                yw[iw,i1] += tmp*kwt[i_kpt]

        yw[:,3] = [sum(yw[iw,:3]) for iw in range(nw)]

    save_and_plot(i_dir,xw,yw)
