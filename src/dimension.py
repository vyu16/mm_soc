import numpy as np

def find_dim():
    dim_dict = {}

    with open("aims.out","r") as f:
        for line in f:
            if line.find("| Number of atoms") != -1:
                t = line.split()
                dim_dict["n_atom"] = int(t[5])

            if line.find("| Number of k-points") != -1:
                t = line.split()
                dim_dict["n_kpt"] = int(t[5])

            if line.find("Index of first state to include in dielectric") != -1:
                t = line.split()
                dim_dict["i_min"] = int(t[10])

            if line.find("Index of last state to include in dielectric") != -1:
                t = line.split()
                dim_dict["i_max"] = int(t[10])

    return dim_dict
