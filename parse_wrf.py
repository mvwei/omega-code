from netCDF4 import Dataset
import numpy as np

from quasigeostrophic_omega import run_qg

f = Dataset("data/wrfout_d01_2013-11-18_06_00_00")

# for now, let's just get the data for eta=0.75 or so.
#
# We want to store dx and dy right now--that is, the distance between
# points on the grid.
#
#
# Imagine a grid as follows:
#
# a --- b --- c
# |     |     |
# d --- e --- f
# |     |     |
# g --- h --- i
#
# dx will be the distance used at each point to calculate the gradient.
# For interior points, this will be the surrounding points.
#
# d(b,a) --- d(c,a) --- d(c,b)
#   |           |         |
# d(e,d) --- d(f,d) --- d(f,e)
#   |           |         |
# d(h,g) --- d(i,g) --- d(i,h)
#
#  dy will be:
#
# (d-a) --- (e-b) --- (f-c)
#   |         |         |
# (e-d) --- (f-d) --- (f-e)
#   |         |         |
# (h-g) --- (i-g) --- (i-h)

# Let's code as if we have all the data. Because fuck this noise.

# For now, let's do only the first timestep.
timestep = 0

lon = f.variables['XLONG'][timestep]
lat = f.variables['XLAT'][timestep]

# expect this to be a X by Y by 5 array. The middle section (3) is
# the one we want to look at; the rest are going to be boundary conditions.
# Yay.
u = f.variables['U'][timestep][0:5]
v = f.variables['V'][timestep][0:5]
temp = f.variables['T'][timestep][0:5]

# bullshit specs for now
p = [900, 875, 850, 825, 800]

run_qg(u=u, v=v, temp=temp, p=p, lat=lat, lon=lon)
