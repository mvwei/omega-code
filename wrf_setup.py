import xarray as xr
import numpy as np
import pandas as pd

from quasigeostrophic_omega import (
    run_qg,
    run_continuity
)

from helpers import (
    w_to_omega,
    heating_rate_from_surface_flux
)

f = xr.open_dataset('../data/2013.11.18_regrid.nc')
low_file = pd.read_csv("../data/PL_2013-11-18-0600_20-4km_stormtrack1.csv")

nb = 2
points_around_low=15

# parse the file to get only the layers we want.
if nb == 2:
    f = f.isel(lev=slice(0, 5))
else:
    f = f.sel(lev=[880, 850, 820])

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

low_data = {
    "lat": low_file["Latitude"].values,
    "lon": low_file["Longitude"].values,
    "slp": low_file["SLP"].values,
}

omega = w_to_omega(w, q, temp, levels)

run_qg(u, v, omega, temp, q, lh, hfx, z, lat, lon, levels, low_data, nb=nb, points_around_low=points_around_low)

# results = run_continuity(u, v, omega, lat, lon, levels, w)

# mid_layer = int((len(levels) - 1)/2)

# test_slice = (1, mid_layer, 5, slice(None))

# print("-----------------OMEGA-----------------")
# print(omega[test_slice].values)

# print("-----------------CALC-----------------")
# print(results[test_slice[1:]])

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
