#!/usr/bin/env python

import math
import random
import unittest

import numpy as np
np.seterr(all="raise")

from pylal import sphericalutils as su
from pylal.lalconstants import *

n_tests = 1000

#
# Utility functions
#

def random_polaz():
    return math.acos(2 * random.uniform(0, 1) - 1), \
        random.uniform(0, 2 * LAL_PI)

def random_euler_angle():
    return random.uniform(0, LAL_TWOPI), random.uniform(0, LAL_PI), \
        random.uniform(0, 2 * LAL_PI)

def random_unit_xyz():
    theta, phi = random_polaz()
    return np.array([np.sin(theta) * np.cos(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(theta)])

def norm(xyz_vec):
    return np.sqrt(np.dot(xyz_vec, xyz_vec))

#
# Unit tests
#

class test_euler(unittest.TestCase):
    def test_identity(self):
        """
        Test that Euler angles of zero are the identity operation.
        """
        for i in xrange(n_tests):
            vector = np.array([random_polaz()])
            err = su.angle_between_points(vector,
                su.rotate_euler(vector, 0., 0., 0.))
            self.assertAlmostEqual(err, 0.)

    def test_inverse(self):
        """
        Test that rotation * inverse = identity.
        """
        for i in xrange(n_tests):
            vector = np.array([random_polaz()])
            a, b, c = random_euler_angle()
            err = su.angle_between_points(vector, \
                su.rotate_euler(su.rotate_euler(vector, a, b, c), -c, -b, -a))
            self.assertAlmostEqual(err, 0.)

    def test_angle_preserved(self):
        """
        Test that the angle between two vectors is preserved under rotation.
        """
        for i in xrange(n_tests):
            vectors_0 = np.array((random_polaz(), random_polaz()),
                dtype=float)
            d_0 = su.angle_between_points(*vectors_0)

            a, b, c = random_euler_angle()
            vectors_1 = su.rotate_euler(vectors_0, a, b, c)
            d_1 = su.angle_between_points(*vectors_1)
            self.assertAlmostEqual(d_0, d_1)

    def test_north_pole(self):
        """
        Test that new_z_to_euler() gives Euler angles that will rotate the
        original z axis to the requested z axis.
        """
        north_pole = np.array([[0., 0.]])
        for i in xrange(n_tests):
            new_dir = random_polaz()
            a, b = su.new_z_to_euler(new_dir)
            new_north = su.rotate_euler(north_pole, a, b, 0)
            self.assertAlmostEqual(0., \
                su.angle_between_points(new_north, np.array([new_dir])))

class test_axis(unittest.TestCase):
    def test_norm(self):
        for i in xrange(n_tests):
            x = random_unit_xyz()
            phi = random.uniform(0, LAL_TWOPI)
            axis = random_unit_xyz()
            self.assertAlmostEqual(norm(su.rotate_about_axis(x, axis, phi)), 1)

    def test_identity(self):
        # spinning a vector about itself should always return the original
        for i in xrange(n_tests):
            x = random_unit_xyz()
            phi = random.uniform(0, LAL_TWOPI)
            self.assertTrue(np.allclose(x, su.rotate_about_axis(x, x, phi)))

    def test_right_handed(self):
        x = np.array([1, 0, 0])
        y = np.array([0, 1, 0])
        z = np.array([0, 0, 1])
        self.assertTrue(np.allclose(su.rotate_about_axis(x, z, LAL_PI / 2), y))
        self.assertTrue(np.allclose(su.rotate_about_axis(y, x, LAL_PI / 2), z))
        self.assertTrue(np.allclose(su.rotate_about_axis(z, y, LAL_PI / 2), x))

# construct and run the test suite
suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_euler))
suite.addTest(unittest.makeSuite(test_axis))
unittest.TextTestRunner(verbosity=2).run(suite)
