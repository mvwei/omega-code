import xarray as xr
import numpy as np

from quasigeostrophic_omega import (
    run_qg,
    run_continuity
)

from helpers import (
    w_to_omega,
    heating_rate_from_surface_flux
)

f = xr.open_dataset('../data/2013.11.18_regrid.nc')

# we're getting only the 850Pa values.
f = f.sel(lev=[880, 850, 820]).isel(lat=slice(100, 140), lon=slice(100, 140))

lat = f.lat
lon = f.lon
# convert to pascals
levels = f.lev * 100

u = f.U
v = f.V
w = f.W
temp = f.TK
q = f.Q
z = f.Z

# surface variables
hfx = f.HFX
lh = f.LH
u10 = f.U10
v10 = f.V10
q2 = f.Q2
t2 = f.T2

omega = w_to_omega(w, q, temp, levels)

results = run_qg(u, v, omega, temp, q, lh, hfx, z, lat, lon, levels)

# results = run_continuity(u, v, omega, lat, lon, levels, w)

print("-----------------OMEGA-----------------")
print(omega[1, 1, 5, :].values)

print("-----------------CALC-----------------")
print(results[1, 5, :])

# Things we want in the output nc file (we want output file so we can reaccess
# the data without having to rerun the code):
#
# lat
# lon
# pressure levels -- 850, 500
# omega_q
# omega_sufrace_latent
# omega_surface_sensible
# omega_latent_horizontal
# omega_latent_phrase_change
#
# lat-lon of low
