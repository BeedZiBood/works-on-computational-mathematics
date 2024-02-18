import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

if __name__ == '__main__':
    data = pd.read_csv("test.csv")
    x = np.array(data['x'])
    f = np.array(data['f'])
    lagrange = sp.interpolate.lagrange(x, f)
    print(lagrange)
    print(lagrange.coef)
    x_new = np.arange(-1, 0.9, 0.1)
    plt.plot(x_new, lagrange(x_new), label='Lagrange')
    plt.show()
