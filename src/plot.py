import matplotlib.pyplot as plt

def plot_data(i_dir,x,y):
    plt.plot(x,y[:,0],"r-",label="org-org")
    plt.plot(x,y[:,1],"b-",label="org-inorg")
    plt.plot(x,y[:,2],"y-",label="inorg-inorg")
    plt.plot(x,y[:,3],"k--",label="all")
    plt.xlim(0.0,x[-1]+0.1)
    plt.legend()
    plt.savefig("decomp_"+i_dir+".png",dpi=300)
    plt.clf()
