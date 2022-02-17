import numpy as np
from scipy.io import loadmat
import struct


data = loadmat("data_for_lucas.mat", squeeze_me=True)

zr = data['zr']
O2mean = data['O2mean']
h271mean = data['h271mean']

mask = np.isnan(h271mean).astype(int)

# file structure:
# nx unsigned int (I)
# ny unsigned int (I)
# nx*ny double (d) for data (H) for mask
# Row major order

nz, ny, nx = zr.shape

with open("mask.dat", "wb") as fd:
    fd.write(struct.pack("I"*2, ny, nx))
    fd.write(struct.pack("H"*ny*nx, *mask.astype(np.uint16).ravel()))

with open("s.dat", "wb") as fd:
    fd.write(struct.pack("I"*2, ny, nx))
    fd.write(struct.pack("d"*ny*nx, *h271mean.astype(np.double).ravel()))


with open("z.dat", "wb") as fd:
    fd.write(struct.pack("I"*3, nz, ny, nx))
    fd.write(struct.pack("d"*nz*ny*nx, *zr.astype(np.double).ravel()))

with open("v.dat", "wb") as fd:
    fd.write(struct.pack("I"*3, nz, ny, nx))
    fd.write(struct.pack("d"*nz*ny*nx, *O2mean.astype(np.double).ravel()))
k=30
j=20
i=427
print(h271mean[j, i])
m = j*nx + i
print(h271mean.ravel()[m])

print(O2mean[k, j, i])
m = k*nx*ny + j*nx + i
print(O2mean.ravel()[m])

mm = 16247
print(zr.ravel()[mm])

with open("vi.dat", "rb") as fd:
    nbytes = struct.calcsize("I"*2)
    b = fd.read(nbytes)
    ny, nx = struct.unpack("I"*2, b)
    nbytes = struct.calcsize("d"*ny*nx)
    b = fd.read(nbytes)
    vi = struct.unpack("d"*ny*nx, b)
vi = np.array(vi).reshape(ny, nx)

data2 = loadmat("O2slices.mat", squeeze_me=True)

O2meaniso=data2['O2meaniso']
dO2 = (O2meaniso - vi).ravel()
dO2s = dO2.compress(np.isfinite(dO2))
