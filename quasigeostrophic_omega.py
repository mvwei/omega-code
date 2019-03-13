import numpy as np
from classes import (
    ScalarField,
    VectorField
)

def run_qg(u, v, temp, p, lat, lon):
    """
    The main equation. Takes in all pre-processed data (so that it's pointing
    in the right directions, correct vertical coordinates, etc) and outputs
    the calculated vertical velocity.

    Shape: a -- b -- c -- d

    lon: numpy of shape (ny, nx)
    lat: numpy of shape (ny, nx)

    u: numpy of shape (nvert, ny, nx)
    v: numpy of shape (nvert, ny, nx)
    temp: numpy of shape (nvert, ny, nx)
    """

    return
