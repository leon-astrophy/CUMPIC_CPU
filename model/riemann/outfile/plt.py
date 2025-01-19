import h5py
import numpy as np
import matplotlib.pyplot as plt

data = h5py.File("star_weno_data_5.hdf5", "r")
prim = data['prim'][()].T
xf = data['xf'][()]
rho = prim[0,:,0,0]
plt.plot(xf[0:1000], rho)
plt.show()
