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
    heating_rate_from_surface_flux,
    heating_rate_from_moisture_change,
    heating_rate_from_horizontal_flux,
)

def run_continuity(u, v, omega, lat, lon, levels, w):
    def continuity_helper(time_slice, omega, dp):
        print("Getting RHS")
        b = get_rhs(omega, time_slice, 1)

        print("Getting coefficient matrix")
        # this is just wrong as a depiction of the continuity equation,
        # but we're leaving it here mostly to prove that orders of magnitude
        # are okay.
        A = get_coefficient_matrix(
            nb=1,
            shape=time_slice.shape,
            point_coeff=(-2/dp),
            dp_coeff=(1/dp),
        )
        # print(A)

        print("Solver")
        return np.linalg.solve(A, b)

    # Generate the coordinate field.
    c = CoordinateField(lat, lon, great_circle_distance, levels)

    # Calculate du/dx.
    u_field = ScalarField(u, c.dx, c.dy)

    # Calculate dv/dy.
    v_field = ScalarField(v, c.dx, c.dy)

    # dw/dp = -du/dx - dv/dy
    calculated_values = -u_field.get_ddx() - v_field.get_ddy()

    print("-----------------C FIELD-----------------")
    print(c.dx)
    print("-----------------U FIELD-----------------")
    print(u[1, 1])
    print("-----------------DU/DX-----------------")
    print(u_field.get_ddx()[1, 1])

    results = np.empty((calculated_values.shape[0]))

    # let's just do the first time slice right now
    results = continuity_helper(calculated_values[1], omega[1], c.dp)

    return np.reshape(results, (levels.size, lat.size, lon.size))

def run_qg(u, v, omega, temp, qv, lhf, shf, height, lat, lon, levels, solved_for_plvl):
    def qg_helper(time_slice, omega, dx, dy, dp):
        sigma = 0.000002      # units are m2Pa-2s-2
        fo = 0.0001           # units are s-1

        print("Getting RHS")
        b = get_rhs(omega, time_slice, 1)

        print("Getting coefficient matrix")
        A = get_coefficient_matrix(
            nb=1,
            shape=time_slice.shape,
            point_coeff=None,
            dp_coeff=None,
            d2p_coeff=None,
            dx_coeff=None,
            d2x_coeff=None,
            dy_coeff=None,
            d2y_coeff=None,
        )
        print(A)

        print("Solver")
        return np.linalg.solve(A, b)

    shape = u.shape

    # Let's get some constants and useful things here.
    plvl = levels.values[:, np.newaxis, np.newaxis] # this is for easy divisibility.

    dt = 3600             # time in seconds between timesteps

    R = 287.05            # units are J/kgK
    cp = 1004.            # units are J/kgK

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

    # HEAT TERM TIME
    # Calculate heat contribution form surface sensible heat flux.
    H_surface_sensible_flux = np.zeros(shape)
    H_surface_sensible_flux[:, 2] = heating_rate_from_surface_flux(shf, height, levels)

    # # Calculate heat contribution from surface latent heat flux.
    H_surface_latent_flux = np.zeros(shape)
    H_surface_latent_flux[:, 2] = heating_rate_from_surface_flux(lhf, height, levels)

    # Calculate horizontal latent heat flux (advection of moisture).
    H_latent_horizontal_flux = heating_rate_from_horizontal_flux(u_f, v_f, qv_f, temp, levels, height, dt)

    # Calculate latent heat contribution from change in vapor over time.
    # It will be zero for the first time step, because we don't have any previous
    # moisture data. This is fine because we're throwing out the first data,
    # anyhow.
    H_latent_phase_change_flux = np.zeros(shape)
    H_latent_phase_change_flux[1:] = heating_rate_from_moisture_change(qv)

    # all our forcing values! Here we are, at long last.
    Q_vector = VectorField(q1, q2, c.dx, c.dy)

    H_sflx_f = ScalarField(H_surface_sensible_flux, c.dx, c.dy)
    H_lflx_f = ScalarField(H_surface_latent_flux, c.dx, c.dy)
    H_lhor_f = ScalarField(H_latent_horizontal_flux, c.dx, c.dy)
    H_lphase_f = ScalarField(H_latent_phase_change_flux, c.dx, c.dy)

    calculated_values = -2 * Q_vector.divergence() - (R / (plvl * cp)) * (
        H_sflx_f.get_laplacian() +
        H_lflx_f.get_laplacian() +
        H_lhor_f.get_laplacian() +
        H_lphase_f.get_laplacian()
    )

    # let's just do the first time slice right now
    results = continuity_helper(calculated_values[1], omega[1], c.dx, c.dy, c.dp)

    return np.reshape(results, (levels.size, lat.size, lon.size))
