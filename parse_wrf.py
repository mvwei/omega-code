import xarray as xr
import numpy as np

from quasigeostrophic_omega import run_qg

f = xr.open_dataset('data/2013.11.18_regrid.nc')

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

# For now, let's do only the first timestep. Eventually, this will
# be a for loop.
timestep = 0

lon = f.lon
lat = f.lat
pressure_levels = f.levs

# expect this to be a X by Y by 5 array. The middle section (3) is
# the one we want to look at; the rest are going to be boundary conditions.
# Yay.
u = f.U
v = f.V
w = f.W
temp = f.TK
q = f.Q

# surface variables
hfx = f.HFX
lh = f.LH
u10 = f.U10
v10 = f.V10
q2 = f.Q2
t2 = f.T2


run_qg(u=u, v=v, temp=temp, pressure_levels=pressure_levels, lat=lat, lon=lon)
