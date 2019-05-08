import numpy as np
from classes import (
    ScalarField,
    VectorField,
    CoordinateField,
)

from functions import (
    get_coefficient_matrix,
    get_rhs,
)

from helpers import (
    great_circle_distance,
    w_to_omega,
    omega_to_w,
)

def run_continuity(u, v, omega, lat, lon, levels):
    def continuity_helper(time_slice, omega, dp):
        print("Getting RHS")
        b = get_rhs(omega, time_slice, 1)

        print(b.size)

        print("Getting coefficient matrix")
        A = get_coefficient_matrix(
            nb=1,
            shape=time_slice.shape,
            point_coeff=(-2/dp**2),
            dp_coeff=(1/dp**2),
        )
        print(A)

        print("Solver")
        return np.linalg.solve(A, b)

    # import pdb
    # pdb.set_trace()

    # Generate the coordinate field.
    c = CoordinateField(lat, lon, great_circle_distance, levels)

    # Calculate du/dx.
    u_field = ScalarField(u, c.dx, c.dy)

    # Calculate dv/dy.
    v_field = ScalarField(v, c.dx, c.dy)

    # dw/dp = -du/dx - dv/dy
    calculated_values = -u_field.get_ddx() - v_field.get_ddy()

    # results = np.empty((calculated_values.shape[0]))

    # let's just do the first time splice right now
    results = continuity_helper(calculated_values[1], omega[1], c.dp)

    return np.reshape(results, (levels.size, lat.size, lon.size))

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
    #
    def qg_omega_helper_func():
        return

    return
