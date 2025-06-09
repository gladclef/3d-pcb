import unittest

import numpy as np
import numpy.testing as npt

from geometry import *

class TestGeometry(unittest.TestCase):
    def test_distance_along_line(self):
        sq2 = np.sqrt(1/2)
        eq = npt.assert_array_almost_equal
        arr = np.array
        eq( actual=arr(distance_along_line((1,       0),     1)), desired=arr((sq2, sq2)),    decimal=4 )
        eq( actual=arr(distance_along_line((-1,      0),     1)), desired=arr((sq2, -sq2)),   decimal=4 )
        eq( actual=arr(distance_along_line((1/10,    0),     1)), desired=arr((.995, .0995)), decimal=4 )
        eq( actual=arr(distance_along_line((10,      0),     1)), desired=arr((.0995, .995)), decimal=4 )
        eq( actual=arr(distance_along_line((np.inf,  0),     1)), desired=arr((0, 1)),        decimal=4 )
        eq( actual=arr(distance_along_line((-np.inf, 0),     1)), desired=arr((0, -1)),       decimal=4 )
        eq( actual=arr(distance_along_line((0,       0),     1)), desired=arr((1, 0)),        decimal=4 )
        eq( actual=arr(distance_along_line((-0,      0),     1)), desired=arr((1, 0)),        decimal=4 )
        eq( actual=arr(distance_along_line(((0,0), (-10,0)), 1)), desired=arr((-1, 0)),       decimal=4 )
    
    def test_get_tangent_line(self):
        eq = npt.assert_array_almost_equal
        arr = np.array
        eq( actual=arr(get_tangent_line( (1, 0),       (0, 0)     )),  desired=arr((-1, 0)),             decimal=4 )
        eq( actual=arr(get_tangent_line( (-1, 0),      (0, 0)     )),  desired=arr((1, 0)),              decimal=4 )
        eq( actual=arr(get_tangent_line( (1/10, 0),    (0, 0)     )),  desired=arr((-10, 0)),            decimal=4 )
        eq( actual=arr(get_tangent_line( (10, 0),      (0, 0)     )),  desired=arr((-1/10, 0)),          decimal=4 )
        eq( actual=arr(get_tangent_line( (np.inf, 0),  (0, 0)     )),  desired=arr((0, 0)),              decimal=4 )
        eq( actual=arr(get_tangent_line( (-np.inf, 0), (0, 0)     )),  desired=arr((0, 0)),              decimal=4 )
        eq( actual=arr(get_tangent_line( (0, 0),       (0, 0)     )),  desired=arr((-np.inf, 0)),        decimal=4 )
        eq( actual=arr(get_tangent_line( (-0, 0),      (0, 0)     )),  desired=arr((-np.inf, 0)),        decimal=4 )
        eq( actual=arr(get_tangent_line( ((0,0), (-10,0)), (0, 0) )),  desired=arr(((0, 0), (0, 10))),   decimal=4 )
        eq( actual=arr(get_tangent_line( ((0,0), (5,5)),   (0, 0) )),  desired=arr(((0, 0), (10, -10))), decimal=4 )
        eq( actual=arr(get_tangent_line( ((0,0), (-10,0)), (-5, 0) )), desired=arr(((-5, 0), (-5, 10))), decimal=4 )
        eq( actual=arr(get_tangent_line( ((0,0), (5,5)),   (2, 2) )),  desired=arr(((2, 2), (12, -8))),  decimal=4 )

    def test_line_two_points_to_slope_intercept(self):
        eq = npt.assert_array_almost_equal
        arr = np.array
        eq( actual=arr(line_two_points_to_slope_intercept( ((0, 0),     (10, 0))    )), desired=arr((0, 0)),       decimal=4 )
        eq( actual=arr(line_two_points_to_slope_intercept( ((0, 0),     (-10, 0))   )), desired=arr((0, 0)),       decimal=4 )
        eq( actual=arr(line_two_points_to_slope_intercept( ((0, 0),     (0, 10))    )), desired=arr((np.inf, 0)),  decimal=4 )
        eq( actual=arr(line_two_points_to_slope_intercept( ((0, 0),     (0, -10))   )), desired=arr((-np.inf, 0)), decimal=4 )
        eq( actual=arr(line_two_points_to_slope_intercept( ((0, 0),     (10, 10))   )), desired=arr((1, 0)),       decimal=4 )
        eq( actual=arr(line_two_points_to_slope_intercept( ((0, 0),     (10, -10))  )), desired=arr((-1, 0)),      decimal=4 )
        eq( actual=arr(line_two_points_to_slope_intercept( ((0, 0),     (-10, 10))  )), desired=arr((-1, 0)),      decimal=4 )
        eq( actual=arr(line_two_points_to_slope_intercept( ((0, 0),     (-10, -10)) )), desired=arr((1, 0)),       decimal=4 )
        eq( actual=arr(line_two_points_to_slope_intercept( ((10, 10),   (0, 0))     )), desired=arr((1, 0)),       decimal=4 )
        eq( actual=arr(line_two_points_to_slope_intercept( ((10, -10),  (0, 0))     )), desired=arr((-1, 0)),      decimal=4 )
        eq( actual=arr(line_two_points_to_slope_intercept( ((-10, 10),  (0, 0))     )), desired=arr((-1, 0)),      decimal=4 )
        eq( actual=arr(line_two_points_to_slope_intercept( ((-10, -10), (0, 0))     )), desired=arr((1, 0)),       decimal=4 )

if __name__ == '__main__':
    unittest.main()
