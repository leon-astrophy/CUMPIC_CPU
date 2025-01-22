import h5py
import numpy as np
import matplotlib.pyplot as plt

data = h5py.File("star_weno_data_0.hdf5", "r")

prim = data['prim'][()].T
r = data['xf'][()]
th = data['yf'][()]
phi = data['zf'][()]

rnew = np.zeros(len(r)-1)
thnew = np.zeros(len(th)-1)

for i in range (0, len(rnew)):
  rnew[i] = 0.5*(r[i+1]+r[i])
 
for i in range (0, len(thnew)):
  thnew[i] = 0.5*(th[i+1]+th[i])
   
R, Th = np.meshgrid(rnew, thnew, indexing="ij")
X = R*np.sin(Th); Y = R*np.cos(Th)

rho = prim[0,:,:,0]

plt.pcolormesh(X, Y, np.log10(rho))
plt.show()