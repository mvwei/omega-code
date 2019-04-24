import numpy as np
from classes import (
    ScalarField,
    VectorField
)

def run_qg(u, v, temp, pressure_levels, lat, lon):
    """
    The main equation. Takes in all pre-processed data (so that it's pointing
    in the right directions, correct vertical coordinates, etc) and outputs
    the calculated vertical velocity.

    Shape: a -- b -- c -- d

    nx: number of points along a constant latitude + constant vertical level.
    ny: number of points along a constant longitude + constant vertical level.
    nz: number of points, vertical stacked, along a constant lon-lat pair.

    lon: 1-D xarray of shape (nx)
    lat: 1-D xarray of shape (ny)

    u: numpy of shape (nz, ny, nx)
    v: numpy of shape (nz, ny, nx)
    temp: numpy of shape (nz, ny, nx)
    """

    return
