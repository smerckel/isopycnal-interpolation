import numpy as np
import struct
import netCDF4


data = netCDF4.Dataset("out.nc")


O2 = data.variables['O2'][0,...]
salt = data.variables['salt'][0,...]
temp = data.variables['temp'][0,...]

h = data.variables['h'][...]
s_rho = data.variables['s_rho'][...]
zeta = data.variables['zeta'][...]
mask_rho = data.variables['mask_rho'][...]

# file structure:
# nx unsigned int (I)
# ny unsigned int (I)
# nx*ny double (d) for data (H) for mask
# Row major order

nz, ny, nx = O2.shape
#2D fields
with open("mask_rho.dat", "wb") as fd:
    fd.write(struct.pack("I"*2, ny, nx))
    fd.write(struct.pack("d"*ny*nx, *mask_rho.astype(np.double).ravel()))

    
with open("zeta.dat", "wb") as fd:
    fd.write(struct.pack("I"*2, ny, nx))
    fd.write(struct.pack("d"*ny*nx, *zeta.astype(np.double).ravel()))

with open("h.dat", "wb") as fd:
    fd.write(struct.pack("I"*2, ny, nx))
    fd.write(struct.pack("d"*ny*nx, *h.astype(np.double).ravel()))
# 1D
with open("s_rho.dat", "wb") as fd:
    fd.write(struct.pack("I", nz))
    fd.write(struct.pack("d"*nz, *s_rho.astype(np.double).ravel()))

#3D fields
    
with open("O2.dat", "wb") as fd:
    fd.write(struct.pack("I"*3, nz, ny, nx))
    fd.write(struct.pack("d"*nz*ny*nx, *O2.astype(np.double).ravel()))

with open("salt.dat", "wb") as fd:
    fd.write(struct.pack("I"*3, nz, ny, nx))
    fd.write(struct.pack("d"*nz*ny*nx, *salt.astype(np.double).ravel()))

with open("temp.dat", "wb") as fd:
    fd.write(struct.pack("I"*3, nz, ny, nx))
    fd.write(struct.pack("d"*nz*ny*nx, *temp.astype(np.double).ravel()))

