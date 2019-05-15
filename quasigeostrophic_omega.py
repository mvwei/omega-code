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

def run_qg(u, v, w, temp, qv, lhf, shf, height, lat, lon, levels):
    def qg_helper(time_slice, omega, dp):
        print("Getting RHS")
        b = get_rhs(omega, time_slice, 1)

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

    # Let's get some constants and useful things here.
    plvl = levels[:, np.newaxis, np.newaxis] # this is for easy divisibility.

    sigma = 0.000002      # units are m2Pa-2s-2
    fo = 0.0001           # units are s-1
    R = 287.05            # units are J/kgK
    cp = 1004.            # units are J/kgK
    Lv = 2500000.         # units are kg/m3
    ts = 21600.           # units are seconds

    # Generate the coordinate field.
    c = CoordinateField(lat, lon, great_circle_distance, levels)

    # Let's get some fields.
    u_f = ScalarField(u, c.dx, c.dy)
    v_f = ScalarField(v, c.dx, c.dy)
    t_f = ScalarField(temp, c.dx, c.dy)
    qv_f = ScalarField(qv, c.dx, c.dy)

    # Calculate Q.
    q1 = (-R / plvl) * (u_f.get_ddx() * t_f.get_ddx() + v_f.get_ddx() * t_f.get_ddy())
    q2 = (-R / plvl) * (u_f.get_ddy() * t_f.get_ddx() + v_f.get_ddy() * t_f.get_ddy())

    Q_vector = VectorField(q1, q2, c.dx, c.dy)

    # Calculate heat contribution form surface sensible heat flux.
    H_surface_sensible_flux = None

    # Calculate heat contribution from surface latent heat flux.
    H_surface_latent_flux = None

    # Calculate latent heat contribution from change in vapor over time.
    H_latent_horizontal_flux = None

    # Calculate horizontal latent heat flux (advection of moisture).
    H_latent_phase_change_flux = None

    calculated_values = -2 * Q_vector.divergence() - (R / (plvl * cp)) * (
        H_surface_sensible_flux.get_laplacian() +
        H_surface_latent_flux.get_laplacian() +
        H_latent_horizontal_flux.get_laplacian() +
        H_latent_phase_change_flux.get_laplacian()
    )

    return
