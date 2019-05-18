import numpy as np
from numpy.testing import assert_array_equal

# get from the directory above.
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from functions import (
    get_rhs,
    get_coefficient_matrix,
)

def test_get_coefficient_matrix():
    result = get_coefficient_matrix(
        nb=1,
        shape=(3, 3, 3),
        center_coeff=2,
        dx_coeff=1,
        dy_coeff=3,
        dp_coeff=4,
    )

    expected = np.identity(27)

    row = np.zeros(27)

    mid = 13

    # center
    row[mid] = 2
    # x (+/- 1)
    row[mid-1] = 1
    row[mid+1] = 1
    # y (+/- 3)
    row[mid-3] = 3
    row[mid+3] = 3
    # z (+/- 9)
    row[mid-9] = 4
    row[mid+9] = 4

    expected[mid] = row

    assert_array_equal(result, expected)

def test_get_coefficient_matrix_funcs_as_args():
    result = get_coefficient_matrix(
        nb=2,
        shape=(5, 5, 5),
        center_coeff=lambda z, y, x: 100,
        dx_coeff=lambda z, y, x: 1,
        d2x_coeff=lambda z, y, x: 10,
        dy_coeff=lambda z, y, x: 2,
        d2y_coeff=lambda z, y, x: 20,
        dp_coeff=lambda z, y, x: 3,
        d2p_coeff=lambda z, y, x: 30,
    )

    expected = np.identity(125)

    row = np.zeros(125)

    mid = 62

    # center
    row[mid] = 100

    # x
    row[mid-1] = 1
    row[mid+1] = 1

    # d2x
    row[mid-2] = 10
    row[mid+2] = 10

    # y
    row[mid-5] = 2
    row[mid+5] = 2

    # d2y
    row[mid-10] = 20
    row[mid+10] = 20

    # p
    row[mid-25] = 3
    row[mid+25] = 3

    # d2p
    row[mid-50] = 30
    row[mid+50] = 30

    expected[mid] = row

    assert_array_equal(result, expected)

def test_get_rhs():
    shape = (3, 3, 3)

    boundaries = np.ones(shape)
    calculated_values = np.ones(shape) * 6

    result = get_rhs(boundaries, calculated_values, 1)

    unique, counts = np.unique(result, return_counts=True)
    counter = dict(zip(unique, counts))

    assert(counter[1] == 26)
    assert(counter[6] == 1)

    assert(result[13] == 6)
