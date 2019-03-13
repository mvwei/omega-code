import numpy as np

from functions import (
    horizontal_gradient,
    vertical_gradient,
    laplacian,
    great_circle_distance,
)

# Given a three-dimensional NumPy array, make a class around this. This is
# so that we can easily access derivative information, so that we can figure
# out how to write some readable equations. We definitely need some readable
# equations.
class ScalarField:
    def __init__(self, np_array):
        self.values = values

    def ddx(self):
        return

    def ddy(self):
        return

    def ddp(self):
        return

    def gradient(self):
        return

    def laplacian(self):
        return

# Values should be of form (a, b), corresponding to the vector ai + bj
class VectorField:
    def __init__(self, values):
        self.values = values

    def divergence(self):
        return

class CoordinateField:
    def __init__(self, lat, lon):
        self.coors = np.dstack(lat, lon)

    def dx(self):
        """
        Should return a field that represents dx at that specific point
        while calculating a normal x gradient. For interior points, we calculate
        dx as the distance between x+1 and x-1. On the left edge, we calculate
        dx as the difference between x and x+1. On the right edge, we calculate
        dx as x and x-1.
        """
        if self.dx:
            return self.dx

        # the x[i-1] array. The roll function shifts values n steps in the increasing
        # direction. Thus, what once was at n-1 will now be at n, which is what we want!
        xm1 = np.roll(self.coors, 1, axis=1)
        xm1[:, 0] = self.coors[:, 0]

        # the x[i+1] array. Similarly, by shifting everything over by -1, what once was at
        # n+1 will be at n.
        xp1 =  np.roll(self.coors, -1, axis=1)
        xp1[:, -1] = self.coors[:, -1]

        dx = great_circle_distance(xm1, xp1)

        self.dx = dx

        return dx

    def dy(self):
        """
        Should return a field that represents dx at that specific point
        while calculating a normal y gradient. For interior points, we calculate
        dy as the distance between y+1 and y-1. On the "north" edge, we calculate
        dy as the difference between y and y+1. On the  edge, we calculate
        dy as y and y-1.
        """
        if self.dy:
            return self.dy

        ym1 = np.roll(self.coors, 1, axis=0)
        ym1[0] = self.coors[0]

        yp1 = np.roll(self.coors, -1, axis=0)
        yp1[-1] = self.coors[-1]

        dy = great_circle_distance(ym1, yp1)

        return dy

