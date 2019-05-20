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
    w_friction,
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
            center_coeff=(-1/(2*dp)),
            dp_coeff=[0, (1/(2*dp))],
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

    results = np.empty((calculated_values.shape[0]))

    # let's just do the first time slice right now
    results = continuity_helper(calculated_values[1], omega[1], c.dp)

    return np.reshape(results, (levels.size, lat.size, lon.size))

def run_qg(u, v, omega, temp, qv, lhf, shf, height, u10, v10, t2, lat, lon, levels, low_data, nb=2, points_around_low=5):
    def qg_helper(calc_values, omega, dx, dy, dp):
        """
        This is the function that actually implements the solver. It is broken
        into two parts: the coefficient matrix generation, and the result matrix
        generation. All of these functions are for the purpose of
        generating these functions.
        """
        sigma = 0.000002      # units are m2Pa-2s-2
        f0 = 0.0001           # units are s-1

        def get_d2x_coeff(z, y, x):
            """
            Only used if nb==2.
            """
            x_dist = dx[y]

            return (-1./12.) * sigma / x_dist**2

        def get_dx_coeff(z, y, x):
            x_dist = dx[y]

            if nb == 2:
                return (4./3.) * sigma / x_dist**2

            elif nb == 1:
                return sigma / x_dist**2

        def get_center_coeff(z, y, x):
            x_dist = dx[y]

            if nb == 2:
                x_coeff = -(5./2.) * sigma / x_dist**2
                y_coeff = -(5./2.) * sigma / dy**2
                p_coeff = -(5./2.) * f0**2 / dp**2

            elif nb == 1:
                x_coeff = -2 * sigma / x_dist**2
                y_coeff = -2 * sigma / dy**2
                p_coeff = -2 * f0**2 / dp**2

            return x_coeff + y_coeff + p_coeff

        b = get_rhs(boundary_values=omega, calculated_values=calc_values, nb=nb)

        if nb == 2:
            A = get_coefficient_matrix(
                nb=nb,
                shape=calc_values.shape,
                center_coeff=get_center_coeff,
                dx_coeff=get_dx_coeff,
                d2x_coeff=get_d2x_coeff,
                dy_coeff=(4./3.) * (sigma / dy**2),
                d2y_coeff=(-1./12.) * (sigma / dy**2),
                dp_coeff=(4./3.) * (f0**2 / dp**2),
                d2p_coeff=(-1./12.) * (f0**2 / dp**2),
            )
        elif nb == 1:
            A = get_coefficient_matrix(
                nb=nb,
                shape=calc_values.shape,
                center_coeff=get_center_coeff,
                dx_coeff=get_dx_coeff,
                dy_coeff=(sigma/dy**2),
                dp_coeff=(f0**2/dp**2),
            )

        return np.linalg.solve(A, b)

    if nb != 2 and nb != 1:
        raise ValueError("nb should either be 1 or 2.")

    shape = u.shape

    # Let's get some constants and useful things here.
    plvl = levels.values[:, np.newaxis, np.newaxis] # this is for easy divisibility.

    dt = 3600             # time in seconds between timesteps
    R = 287.05            # units are J/kgK
    cp = 1004.            # units are J/kgK

    # # Generate the coordinate field.
    c = CoordinateField(lat, lon, great_circle_distance, levels)

    # Let's get some fields.
    u_f = ScalarField(u, c.dx, c.dy)
    v_f = ScalarField(v, c.dx, c.dy)
    t_f = ScalarField(temp, c.dx, c.dy)
    qv_f = ScalarField(qv, c.dx, c.dy)

    # Calculate Q.
    print("Calculating Q components")
    q1 = (-R / plvl) * (u_f.get_ddx() * t_f.get_ddx() + v_f.get_ddx() * t_f.get_ddy())
    q2 = (-R / plvl) * (u_f.get_ddy() * t_f.get_ddx() + v_f.get_ddy() * t_f.get_ddy())

    print("Calculating surface sensible heat flux")
    # Calculate heat contribution form surface sensible heat flux.
    H_surface_sensible_flux = ScalarField(
        heating_rate_from_surface_flux(shf, height, levels),
        c.dx,
        c.dy
    )

    print("Calculating surface latent heat flux")
    # # Calculate heat contribution from surface latent heat flux.
    H_surface_latent_flux = ScalarField(
        heating_rate_from_surface_flux(lhf, height, levels),
        c.dx,
        c.dy
    )

    # Latent phase change!
    print("Calculating phase change latent heat flux")
    H_latent_phase_change_flux = ScalarField(
        heating_rate_from_moisture_change(qv),
        c.dx,
        c.dy
    )

    # all our forcing values! Here we are, at long last.
    Q_vector = VectorField(q1, q2, c.dx, c.dy)

    # Time to calculate the variables in our equation! Exciting!
    print("Calculating Q Divergence")
    q_terms = - 2 * Q_vector.divergence()

    print("Calculating heat laplacians")
    sfc_sensible = - (R / (plvl * cp)) * H_surface_sensible_flux.get_laplacian()
    sfc_latent = - (R / (plvl * cp)) * H_surface_latent_flux.get_laplacian()
    latent_phase = - (R / (plvl * cp)) * H_latent_phase_change_flux.get_laplacian()

    heat_terms = sfc_sensible + sfc_latent + latent_phase

    total_values = q_terms + heat_terms

    # get boundaries, if at all.
    boundaries = omega

    if nb == 1:
        print("Calculating w friction")
        w_fric = w_friction(u10, v10, t2, c.dx, c.dy)
        boundaries[:, 0] = w_fric[:]

    # some indexing stuff
    vertical_midpoint = int((levels.size-1)/2)
    meridional_midpoint = int((lat.size-1)/2)
    zonal_midpoint = int((lon.size-1)/2)

    low_lat_data = low_data["lat"]
    low_lon_data = low_data["lon"]
    low_slp_data = low_data["slp"]

    # this is where the solving begins. For each time step:
    # 1. Get the area around the low
    # 2. Solve for the area around the low for multiple "b" in Ax=b
    # 3. Print out.
    with open('output.txt', 'w') as f:

        f.write("Time,SLP,lat,lon,closest lat,closest lon,WRF omega,Calculated omega,Q omega,Surface Sensible omega,Surface Latent omega,Latent phase omega\n")

        for time_idx in range(1, u.XTIME.size):
            time = u.XTIME[time_idx].values

            print("Calculating for %s (%s out of %s)" % (time, time_idx, u.XTIME.size - 1))

            # get the calculated arrays. We're just going to slice to get the
            # things we want, we're not going to run different boundaries each time.
            # That'd suck.
            #
            # The vast majority of the code below is for getting the values around
            # the low.
            low_lat = low_lat_data[time_idx]
            low_lon = low_lon_data[time_idx]
            low_slp = low_slp_data[time_idx]

            if np.isnan(low_lat) or np.isnan(low_lon) or np.isnan(low_slp):
                continue

            closest_lat = lat.sel(lat=low_lat, method="nearest")
            closest_lon = lon.sel(lon=low_lon, method="nearest")

            closest_lat_idx = np.where(lat==closest_lat)[0][0]
            closest_lon_idx = np.where(lon==closest_lon)[0][0]

            # sanity check to make sure that we don't go over
            if (
                closest_lat_idx - points_around_low < 0 or
                closest_lat_idx + points_around_low + 1 >= shape[-2] or
                closest_lon_idx - points_around_low < 0 or
                closest_lon_idx + points_around_low + 1 >= shape[-1]
            ):
                raise Exception("We have gone out of bounds! Considering reducing points_around_low.")

            calc_slice = (
                time_idx,
                slice(None),
                slice(closest_lat_idx - points_around_low, closest_lat_idx + points_around_low + 1),
                slice(closest_lon_idx - points_around_low, closest_lon_idx + points_around_low + 1),
            )

            calc_arr = [
                total_values[calc_slice],
                q_terms[calc_slice],
                sfc_sensible[calc_slice],
                sfc_latent[calc_slice],
                latent_phase[calc_slice]
            ]

            boundary = boundaries[calc_slice]

            results = []

            for arr in calc_arr:
                calculated = qg_helper(arr, boundary, c.dx, c.dy, c.dp)
                reshaped = np.reshape(
                    calculated,
                    (levels.size, points_around_low * 2 + 1, points_around_low * 2 + 1)
                )
                result = reshaped[vertical_midpoint, points_around_low, points_around_low]

                results.append(str(result))

            # this is the omega we're basing off of!
            wrf_omega = omega[time_idx, vertical_midpoint, closest_lat_idx, closest_lon_idx]

            float_formatter = lambda x: "%.3f" % x

            context_data = [
                str(time),
                float_formatter(low_slp),
                float_formatter(low_lat),
                float_formatter(low_lon),
                float_formatter(closest_lat.values),
                float_formatter(closest_lon.values),
                str(wrf_omega.values),
            ]

            row = context_data + results
            f.write(",".join(row) + "\n")
