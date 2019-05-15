import numpy as np
import pandas

"""
For model-specific functions. Calculating things like distances,
heat flux, etc.
"""

def great_circle_distance(coor1, coor2, unit='m'):
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

    if unit == 'km':
        return km
    elif unit == 'm':
        return km * 1000
    else:
        print("Unit % is not supported, returning in km." % unit)
        return km

def moist_air_density(qv, temp, p):
    rgas = 287.04   # J/K/kg
    eps = 0.622

    # get ready to have this guy as a constant along the pressure axis, which is
    # either axis 0 or 1 depending on whether or not we're passing in time.
    try:
        # if it's a numpy array, this will succeed.
        # If it's a data array, this will fail.
        plvl = p[:, np.newaxis, np.newaxis]
    except IndexError:
        plvl = p.values[:, np.newaxis, np.newaxis]

    # P/RT * (1 + qv) / (1 + 1/eps * qv), except with the right part multiplied
    # by eps/eps
    density = (plvl / (rgas * temp)) * (eps * (1 + qv)) / (eps + qv)

    return density

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
    g = 9.81      # m/s**2

    density = moist_air_density(qv, temp, p)

    return -g * w * density

def omega_to_w(omega, qv, temp, p):
    density = moist_air_density(qv, temp, p)

    return omega / (-g * density)

def heating_rate_from_surface_flux(sflx, height, plvl, density=1.225):
    # hardcoding for now. We only want surface flux to be
    # on the lower level calculations. Should be 0 for upper level (~500 hPa)
    # calculations.
    if plvl[-1] < 500:
        return np.zeros(sflx.shape)

    z = height[:, -1, :, :] - height[:, 0, :, :]

    return sflx / (density * z)

def heating_rate_from_moisture_change(qv, dt=3600):
    Lv = 2500000.   # units are J/kg
    rho = 1.0687    # units are kg/m^3

    t2 = qv[1:]
    t1 = qv[:-1]

    return (Lv * rho * (t2 - t1)) / dt

def heating_rate_from_horizontal_flux(u_f, v_f, qv_f, temp, p, z, dt=3600):
    # hardcoding for the middle layer.
    density = moist_air_density(qv_f.values, temp, p)

    # hardcoding this guy. The target pressure we want is at index 2.
    box_height = np.zeros(temp.shape)

    box_height[:, 2] = z[:, 3] - z[:, 1]
    dx = u_f.dx[:, np.newaxis]
    volume = box_height * (2*dx) * (2*u_f.dy)

    heat_rate = density * volume * dt * (
        (qv_f.values * (u_f.get_ddx() + v_f.get_ddy())) +
        (u_f.values * qv_f.get_ddx() + v_f.values * qv_f.get_ddy())
    )

    return heat_rate

