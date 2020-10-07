import numpy as np

def gauss(x,mu,sigma):
#    pre = 1.0/(sigma*np.sqrt(2*np.pi))
    tmp = (x-mu)/sigma
    tmp = -0.5*tmp*tmp

    return np.exp(tmp)
