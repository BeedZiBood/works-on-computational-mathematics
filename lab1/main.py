import pandas as pd
import numpy as np
import scipy as sp
import QUANC8
import matplotlib.pyplot as plt

if __name__ == '__main__':
    data = pd.read_csv("test.csv")
    x = np.array(data['x'])
    f = np.array(data['f'])
    lagrange = sp.interpolate.lagrange(x, f)
    spline = sp.interpolate.CubicSpline(x, f, bc_type='natural')
    x_new = np.arange(-1, 0.9, 0.1)
    plt.scatter(x, f, label='data')
    plt.plot(x_new, lagrange(x_new), label='Lagrange')
    plt.plot(x_new, spline(x_new), label='Spline')
    plt.legend()
    plt.show()
