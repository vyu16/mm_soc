import matplotlib.pyplot as plt

def plot_all(x,y):
    plt.plot(x,y)
    plt.xlim(0.0,x[-1]+0.1)
    plt.show()
