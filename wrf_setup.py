import xarray as xr
import numpy as np

from quasigeostrophic_omega import (
    run_qg,
    run_continuity
)

from helpers import w_to_omega

f = xr.open_dataset('../data/2013.11.18_regrid.nc')

# we're getting only the 850Pa values.
f = f.isel(lev=slice(1, 4), lat=slice(100, 110), lon=slice(100, 110))

lat = f.lat
lon = f.lon
levels = f.lev

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

omega = w_to_omega(w, q, temp, levels)

results = run_continuity(u, v, omega, lat, lon, levels)

print("-----------------OMEGA-----------------")
print(omega[1, 1, 4, :].values)

print("-----------------CALC-----------------")
print(results[1, 4, :])
