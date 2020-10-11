import numpy as np
import matplotlib.pyplot as plt

# Parse FHI-aims output to find a few numbers
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

# Gaussian (omitting norm as it doesn't matter)
def gauss(x,mu,sigma):
#    pre = 1.0/(sigma*np.sqrt(2*np.pi))
    tmp = (x-mu)/sigma
    tmp = -0.5*tmp*tmp

    return np.exp(tmp)

# Save and plot data
def save_and_plot(i_dir,x,y):
    with open("decomp_"+i_dir+".dat","w") as f:
        for i in range(len(x)):
            f.write(("{:12.6f}"*5).format(x[i],y[i,0],y[i,1],y[i,2],y[i,3])+"\n")

    plt.plot(x,y[:,0],"r-",label="org-org")
    plt.plot(x,y[:,1],"b-",label="org-inorg")
    plt.plot(x,y[:,2],"y-",label="inorg-inorg")
    plt.plot(x,y[:,3],"k--",label="all")
    plt.xlim(0.0,np.ceil(x[-1]))
    plt.legend()
    plt.savefig("decomp_"+i_dir+".png",dpi=300)
    plt.clf()
