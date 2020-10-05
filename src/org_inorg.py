#!/usr/bin/env python

import numpy as np

def find_org_atom(n_atom,filename):
    is_org = np.zeros(n_atom,dtype=int)
    i_atom = -1

    with open(filename) as f:
        for line in f:
            t = line.split()

            if len(t) > 0:
                if t[0] == "atom" or t[0] == "atom_frac":
                    i_atom += 1

                    if t[4] != "Pb" and t[4] != "I":
                        is_org[i_atom] = 1

    return is_org

def find_org_state(i_min,i_max,n_kpt,atom_is_org,filename):
    n_state = i_max-i_min+1
    tmp = np.zeros((n_kpt,n_state))
    i_atom = -1

    with open(filename) as f:
        for line in f:
            t = line.split()

            if len(t) > 0:
                if t[0] == "Atom":
                    i_atom += 1
                    i_kpt = -1
                elif t[0] == "k":
                    i_kpt += 1
                elif t[0] == "Spin":
                    i_kpt = -1
                elif t[0] == "State":
                    continue
                elif t[0] == "#":
                    continue
                elif i_min <= int(t[0]) <= i_max:
                    if atom_is_org[i_atom]:
                        tmp[i_kpt,int(t[0])-i_min] += float(t[3])

    is_org = np.zeros((n_kpt,n_state),dtype=int)

    for i_kpt in range(n_kpt):
        for i_state in range(n_state):
            if tmp[i_kpt,i_state] > 0.8:
                is_org[i_kpt,i_state] = 1

    return is_org
