import numpy as np

# gradient function (we can either try numpy or otherwise)

def get_shifted_by_one(array):
    # the x[i-1] array. The roll function shifts values n steps in the increasing
    # direction. Thus, what once was at n-1 will now be at n, which is what we want!
    xm1 = np.roll(array, 1, axis=1)
    xp1 = np.roll(array, -1, axis=1)
    ym1 = np.roll(array, 1, axis=0)
    yp1 = np.roll(array, -1, axis=0)

    # now set the boundary conditions
    xm1[:, 0] = array[:, 0]
    xp1[:, -1] = array[:, -1]

    ym1[0] = array[0]
    yp1[-1] = array[-1]

    return xm1, xp1, ym1, yp1

def get_shifted_by_two(array):
    return

def horizontal_gradient(values):
    return

def vertical_gradient(values):
    return

def laplacian(values):
    return

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
