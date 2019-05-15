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
    if not np.ndim(distance):
        return np.gradient(f, distance, axis=axis)

    elif np.ndim(distance) == 1 and distance.size == f.shape[axis - 1]:
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

def second_derivative(f, axis, dx):
    """
    Similar to the first derivative.
    """
    def gradient_helper(row):
        dx = row[0]
        actual_data = np.delete(row, 0)

        return calculate_second_derivative(actual_data, axis=0, dx=dx)

    # if scalar value, then just do the gradient.
    if not np.ndim(distance):
        return calculate_second_derivative(f, axis=axis, dx=dx)

    elif np.ndim(distance) == 1 and distance.size == f.shape[axis - 1]:
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

def calculate_second_derivative(f, axis=0, dx=1.):
    """
    Fourth order central difference in range [2, n-3].
    Second order central difference for 1 and n-2.
    First order forward difference for 0.
    First order backward difference for n-1.

    Brute forcing. We're not trying to be subtle here.
    Axis is required.

    Based off the numpy gradient function to be found in
    numpy/lib/function_base.py, except with some modified coefficients
    to account for the second derivative and some reduced functionality.
    """
    N = f.ndim

    out = np.empty_like(f, dtype=f.dtype)

    slice1 = [slice(None)]*N
    slice2 = [slice(None)]*N
    slice3 = [slice(None)]*N
    slice4 = [slice(None)]*N
    slice5 = [slice(None)]*N
    slice6 = [slice(None)]*N

    # Numerical differentiation: fourth order, interior
    slice1[axis] = slice(2, -2)

    slice2[axis] = slice(1, -3)
    slice3[axis] = slice(3, -1)
    slice4[axis] = slice(None, -4)
    slice5[axis] = slice(4, None)

    out[tuple(slice1)] = (-1 * (f[tuple(slice4)] + f[tuple(slice5)]) + 16 * (f[tuple(slice2)] + f[tuple(slice3)]) -30 * f[tuple(slice1)]) / (12 * dx**2)

    # Numerical differentiation: second order, inner boundary
    slice1[axis] = 1
    slice2[axis] = 0
    slice3[axis] = 2

    out[tuple(slice1)] = (f[tuple(slice2)] - 2*f[tuple(slice1)] + f[tuple(slice3)]) / dx**2

    slice1[axis] = -2
    slice2[axis] = -1
    slice3[axis] = -3

    out[tuple(slice1)] = (f[tuple(slice2)] - 2*f[tuple(slice1)] + f[tuple(slice3)]) / dx**2

    # Numerical differentiation: first order, outer boundary. Forward / backward difference.
    slice1[axis] = 0
    slice2[axis] = 1
    slice3[axis] = 2

    out[tuple(slice1)] = (f[tuple(slice1)] - 2*f[tuple(slice2)] + f[tuple(slice3)]) / dx**2

    slice1[axis] = -1
    slice2[axis] = -2
    slice3[axis] = -3

    out[tuple(slice1)] = (f[tuple(slice1)] - 2*f[tuple(slice2)] + f[tuple(slice3)]) / dx**2

    return out

def get_coefficient_matrix(nb, shape, point_coeff=0, x_coeff=0, dx_coeff=0, dy_coeff=0, dp_coeff=0, d2x_coeff=0, d2y_coeff=0, d2p_coeff=0):
    """
    Generate the coefficient matrix.

    If we have any d2 coefficients, we assume that we have boundary conditions 2 above,
    2 below, and 2 along the sides.

    If we have only d1 coefficients, we assume boundary conditions 1 above, 1 below,
    and 1 along the sides.

    Boundary conditions will appear because their rows will just be identity.

    nb is the number of points along each "side" that will be boundaries.

    Collapsing in C-memory order.
    """
    (nz, ny, nx) = shape

    def get_row_from_coordinate(p_index, lat_index, lon_index):
        n_horizontal_points = ny*nx

        return p_index * n_horizontal_points + lat_index * nx + lon_index

    ntotal = nz * ny * nx

    full_arr = np.identity(ntotal)

    # I don't think there's any value in trying to be clever here, seeing as
    # we're placing things all sorts of willy-nilly.

    for z in range(nb, nz-nb):
        for y in range(nb, ny-nb):
            for x in range(nb, nx-nb):
                arr = np.zeros((ntotal))

                point = get_row_from_coordinate(z, y, x)
                arr[point] = point_coeff

                if d2x_coeff:
                    xm2 = get_row_from_coordinate(z, y, x-2)
                    xp2 = get_row_from_coordinate(z, y, x+2)

                    arr[xm2] = d2x_coeff
                    arr[xp2] = d2x_coeff

                if d2y_coeff:
                    ym2 = get_row_from_coordinate(z, y-2, x)
                    yp2 = get_row_from_coordinate(z, y+2, x)

                    arr[ym2] = d2y_coeff
                    arr[yp2] = d2y_coeff

                if d2p_coeff:
                    zm2 = get_row_from_coordinate(z-2, y, x)
                    zp2 = get_row_from_coordinate(z+2, y, x)

                    arr[zm2] = d2p_coeff
                    arr[zp2] = d2p_coeff

                if dx_coeff:
                    xm1 = get_row_from_coordinate(z, y, x-1)
                    xp1 = get_row_from_coordinate(z, y, x+1)

                    arr[xm1] = dx_coeff
                    arr[xp1] = dx_coeff

                if dy_coeff:
                    ym1 = get_row_from_coordinate(z, y-1, x)
                    yp1 = get_row_from_coordinate(z, y+1, x)

                    arr[ym1] = dy_coeff
                    arr[yp1] = dy_coeff

                if dp_coeff:
                    zm1 = get_row_from_coordinate(z-1, y, x)
                    zp1 = get_row_from_coordinate(z+1, y, x)

                    arr[zm1] = dp_coeff
                    arr[zp1] = dp_coeff

                full_arr[point, :] = arr

    return full_arr

def get_rhs(boundary_values, calculated_values, nb):
    """
    Generate the RHS of an equation.

    Will flatten in C-memory order, aka the furthest "right" index
    is the one that increases first. None of that Fortran nonsense.

    boundary_values: the values of the boundaries. In the case of qg-omega, it's omega.

    calculated_values: the right hand side values that we've calculated. In this case,
    divergence of Q and heat terms.

    nb: number of boundary 'layers' needed in order our equation to work.
    Because we're doing fourth order second derivative, we need two points "above"
    and two points "below".
    """
    min_datapoints = 2 * nb + 1

    if np.ndim(boundary_values) > 3:
        print("Do not pass time in as a dimension.")
        return

    if any(npoints < min_datapoints for npoints in np.shape(boundary_values)):
        print("Insufficient boundary layers provided, please provide at least %s." % min_datapoints)
        return

    result_arr = np.copy(boundary_values)

    result_arr[nb:-nb, nb:-nb, nb:-nb] = calculated_values[nb:-nb, nb:-nb, nb:-nb]

    # reshape this into one long monstrosity.
    return np.reshape(result_arr, (-1))
