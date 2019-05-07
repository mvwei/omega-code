import numpy as np

def first_derivative(f, axis, distance):
    """
    We aren't getting fancy here. This is the numpy gradient method, which is
    second order central finite difference in the middle and first order (forward
    difference on the left, backward difference on the right) on the edges.

    The variable distance can either be a scalar or an array. Keep in mind
    that the array's usage is different from numpy's gradient array usage, and is customized
    to our use case.
    """
    def gradient_helper(row):
        """
        Given a row of format [dx, ...data], calculate the gradient of the data
        given dx. Not a function I'm proud of, but until apply_along_axis passes in
        the index of the row, this is what we're doing.
        """
        dx = row[0]
        actual_data = np.delete(row, 0)

        return np.gradient(actual_data, dx)

    # if scalar value, then just do the gradient.
    if not np.ndims(distance):
        return np.gradient(f, distance, axis=axis)

    elif np.ndims(distance) == 1 and distance.size == f.shape[axis - 1]:
        # this is cringeworthy. I will take suggestions about how to make this less
        # bad, but from what I can tell it is impossible to concurrently iterate
        # over two np arrays of different dimensions simultaneously. So we're just
        # going to throw dx in as the first variable, and then pop it out in the
        # function.
        temp_arr = np.insert(f, 0, distance, axis=-1)

        return np.apply_along_axis(gradient_helper, axis, temp_arr)

    else:
        print("What the hell did you pass in?")
        return

def second_derivative(f, axis, distance):
    """
    Oh god. Help me.

    Fourth order central difference in range [2, n-3].
    Second order central difference for 1 and n-2.
    First order forward difference for 0.
    First order backward difference for n-1.

    Likely copied from numpy gradient code with some tweaks.
    """
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

def generate_coefficient_matrix(dx_coeff, dy_coeff, dp_coeff, d2x_coeff, d2y_coeff, d2p_coeff):
    """
    Generate the coefficient matrix.

    If we have any d2 coefficients, we assume that we have boundary conditions 2 above,
    2 below, and 2 along the sides.

    If we have only d1 coefficients, we assume boundary conditions 1 above, 1 below,
    and 1 along the sides.

    Boundary conditions will appear because their rows will just be identity.
    """
    return

def result_matrix():
    """
    The result matrix. Will involve flattening an (nz x ny x nx) array into a (1 x (nx*ny*nx))
    array. Well, this will be fun.
    """
    return
