import numpy as np
from classes import (
    ScalarField,
    VectorField
)

def run_divergence(u, v, w):
    # Calculate du/dx.
    #
    # Calculate dv/dy.
    #
    # Create the empty shell of a massive 3 x nx x ny omega coefficient array.
    # It is all 1s.
    #
    # Create the empty shell of a 3 x nx x ny omega array.
    #
    # Create the w matrix at 850, slot it in.
    #
    # Create the constants matrix at 850, slot it in.
    #
    # Flatten the two arrays.
    #
    # Solve.
    return

def run_qg(u, v, w, temp, pressure_levels, lat, lon):
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

    # First thing's first, get all our objects. Get derivatives working.

    # Caclulate Q.
    #
    # Calculate heat.
    #
    # Create the empty shell of a massive 5 x nx x ny omega coefficient array.
    # All zeros, for now.
    #
    # Create the empty shell of a 5 x nx x ny omega array.
    #
    # Create the omega matrix at 850, slot it in.
    #
    # Create the constants matrix at 850, slot it in.
    #
    # Flatten the two arrays.
    #
    # Solve.

    return
