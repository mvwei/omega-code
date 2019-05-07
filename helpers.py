import numpy as np
import pandas

"""
For model-specific functions. Calculating things like distances,
heat flux, etc.
"""

def great_circle_distance(coor1, coor2):
    """
    Given two sets of coordinate arrays [lat, lon], calculate the great circle distance
    between them.

    We're given things of the format:
    [
        [ [lat1, lon1], [lat2, lon2], [lat3, lon3], ... ],
        [ [lat4, lon4], [lat5, lon5], [lat6, lon6], ... ],
    ]

    And to maximize utility of numpy, we want to split them into arrays of strictly
    lat and strictly lon.

    Formula taken from
    https://stackoverflow.com/questions/29545704/fast-haversine-approximation-python-pandas
    """
    rad_coor1 = np.radians(coor1)
    rad_coor2 = np.radians(coor2)

    nx, ny, ncoor = np.shape(coor1)

    lat1, lon1 = np.dsplit(rad_coor1, 2)
    lat2, lon2 = np.dsplit(rad_coor2, 2)

    lat1 = np.reshape(lat1, (nx, ny))
    lon1 = np.reshape(lon1, (nx, ny))

    lat2 = np.reshape(lat2, (nx, ny))
    lon2 = np.reshape(lon2, (nx, ny))

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2

    c = 2 * np.arcsin(np.sqrt(a))
    km = 6367 * c
    return km

def w_to_omega(w, qv, temp, p):
    """
    w: vertical velocity
    qv: water vapor mixing ratio
    temp: temperature (K)
    p: pressure levels

    Taken from omgcode in wrf_rip_phys_routines.f and adapted from Python.
    Rewritten because it's faster than trying to get wrf-python
    to cooperate with me.

    omega = -rho * g * w

    The bulk of this is trying to estimate the density.
    """
    g = 9.81        # m/s**2
    rgas = 287.04   # J/K/kg
    eps = 0.622

    # get ready to divide this guy along the pressure axis, which is
    # either axis 0 or 1 depending on whether or not we're passing in time.
    plvl = p[:, np.newaxis, np.newaxis]

    # P/RT * (1 + qv) / (1 + 1/eps * qv), except with the right part multiplied
    # by eps/eps
    density = (plvl / (rgas * T)) * (eps * (1 + qv)) / (eps + qv)

    return -g * w * density

def visualize(data):
    """
    I don't know, what the fuck do I know?
    """
    return

def get_2d_matrix_row(i, j, k, nx, ny):
    return k*nx*ny + i*nx + j
