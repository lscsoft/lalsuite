#!/usr/bin/env python

import itertools
import unittest

from pylal import galaxyutils

class test_galaxyutils(unittest.TestCase):
    def test_load_CBCGC(self):
        """
        Test that we can load a CBC galaxy catalog without error.
        """
        galaxyutils.CBCGC.from_file(open("test_galaxyutils_CBCGC.txt"))

    def test_load_GWGC(self):
        """
        Test that we can load a gravitational-wave galaxy catalog without error.
        """
        galaxyutils.GWGC.from_file(open("test_galaxyutils_GWGC.txt"))

    def test_load_deprecated_GalaxyCatalog(self):
        """
        Test that we can use the deprecated interface to the CBC galaxy
        catalog without error.
        """
        galaxyutils.GalaxyCatalog.from_file(open("test_galaxyutils_CBCGC.txt"))

    def test_grid_in_polygon(self):
        """
        Test that is_within_polygon returns the correct result for a few
        simple test cases with known-correct results.
        """
        # integer grid
        points = [(i, j) for i in xrange(10) for j in xrange(10)]

        # triangle
        vertices = [(1.5, 1.7), (4.1, 4.1), (4.1, 1.7)]
        self.assert_(sum(galaxyutils.is_inside_polygon(p, vertices) for \
            p in points) == 6)

        # square
        vertices = [(1.5, 1.7), (1.5, 4.1), (4.1, 4.1), (4.1, 1.7)]
        self.assert_(sum(galaxyutils.is_inside_polygon(p, vertices) for \
            p in points) == 9)

    def test_random_in_polygon(self):
        """
        Put random points inside a polygon and output a PNG for the user
        to see by eye that the result is correct.
        """
        import numpy as np
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg \
            as FigureCanvas

        points = np.random.uniform(low=0.0, high=5.0, size=(5000,2))
        polygon = np.array([[2, 2], [2, 3], [3, 2]], dtype=float)
        found_list = [galaxyutils.is_inside_polygon(point, polygon) for point \
            in points]
        found_points = np.array([point for point,found in \
            itertools.izip(points, found_list) if found])
        missed_points = np.array([point for point,found in \
            itertools.izip(points, found_list) if not found])

        fig = Figure()
        canvas = FigureCanvas(fig)
        ax = fig.add_subplot(111)
        ax.plot(found_points[:, 0], found_points[:, 1], "bx", label="found")
        ax.plot(missed_points[:, 0], missed_points[:, 1], "r.", label="missed")

        closed_poly = np.vstack((polygon, polygon[0]))
        ax.plot(closed_poly[:, 0], closed_poly[:, 1], "k-", label="_nolabel_")
        canvas.print_figure("test_random_in_polygon.png")
        print "Check test_random_in_polygon.png manually."

# construct and run the test suite.
suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_galaxyutils))
unittest.TextTestRunner(verbosity=2).run(suite)
