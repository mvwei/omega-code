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

nb = 1
points_around_low = 5

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
ust = f.UST
temp = f.TK
q = f.Q
q_c = f.Qcl
z = f.Z

# surface variables
hfx = f.HFX
lh = f.LH
u10 = f.U10
v10 = f.V10
q2 = f.Q2
t2 = f.T2
psfc = f.PSFC

low_data = {
    "lat": low_file["Latitude"].values,
    "lon": low_file["Longitude"].values,
    "slp": low_file["SLP"].values,
}

omega = w_to_omega(w, q, temp, levels)

run_qg(u, v, ust, omega, temp, q, q_c, lh, hfx, z, u10, v10, t2, q2, psfc, lat, lon, levels, low_data, nb=nb, points_around_low=points_around_low)

