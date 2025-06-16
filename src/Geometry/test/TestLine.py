import random
import unittest

import numpy as np
import numpy.testing as npt

from Geometry.Line import Line

class TestLine(unittest.TestCase):
    def test_is_point_on_right(self):
        self.assertTrue(Line( 1,  0).is_point_on_right(( 0, -1)))
        self.assertTrue(Line( 1,  1).is_point_on_right(( 1,  0)))
        self.assertTrue(Line( 0,  1).is_point_on_right(( 1,  1)))
        self.assertTrue(Line(-1,  1).is_point_on_right(( 0,  1)))
        self.assertTrue(Line(-1,  0).is_point_on_right((-1,  1)))
        self.assertTrue(Line(-1, -1).is_point_on_right((-1,  0)))
        self.assertTrue(Line( 0, -1).is_point_on_right((-1, -1)))
        self.assertTrue(Line( 1, -1).is_point_on_right(( 0, -1)))
        
        self.assertTrue( Line( 1,  0, y_intercept=1).is_point_on_right((  0,  .5)))
        self.assertFalse(Line( 1,  0, y_intercept=1).is_point_on_right((  0, 1.5)))
        self.assertTrue( Line( 0,  1, x_intercept=1).is_point_on_right((1.5,   1)))
        self.assertFalse(Line( 0,  1, x_intercept=1).is_point_on_right(( .5,   1)))
        self.assertTrue( Line(-1,  0, y_intercept=1).is_point_on_right(( -1, 1.5)))
        self.assertFalse(Line(-1,  0, y_intercept=1).is_point_on_right(( -1,  .5)))
        self.assertTrue( Line( 0, -1, x_intercept=1).is_point_on_right(( .5,  -1)))
        self.assertFalse(Line( 0, -1, x_intercept=1).is_point_on_right((1.5,  -1)))

    def test_intersection(self):
        arr = np.array
        ap = Line.from_angle_point
        eq = npt.assert_array_almost_equal
        π = np.pi

        # 45 degree angles
        eq( actual=arr(ap(  π/4, (1,2)).intersection(ap(3*π/4, (1,4)))), desired=arr((2,3)), decimal=4 )
        eq( actual=arr(ap(3*π/4, (1,2)).intersection(ap(  π/4, (1,4)))), desired=arr((0,3)), decimal=4 )
        eq( actual=arr(ap(5*π/4, (1,2)).intersection(ap(3*π/4, (1,4)))), desired=arr((2,3)), decimal=4 )
        eq( actual=arr(ap(7*π/4, (1,2)).intersection(ap(  π/4, (1,4)))), desired=arr((0,3)), decimal=4 )

        # 90 degree angles
        eq( actual=arr(ap(    0, (1,2)).intersection(ap(  π/2, (1,4)))), desired=arr((1,2)), decimal=4 )
        eq( actual=arr(ap(  π/2, (1,2)).intersection(ap(    0, (1,4)))), desired=arr((1,4)), decimal=4 )
        eq( actual=arr(ap(    π, (1,2)).intersection(ap(  π/2, (1,4)))), desired=arr((1,2)), decimal=4 )
        eq( actual=arr(ap(3*π/2, (1,2)).intersection(ap(    0, (1,4)))), desired=arr((1,4)), decimal=4 )

        # 45 degree angle and 90 degree angle
        eq( actual=arr(ap(  π/4, (1,2)).intersection(ap(    0, (1,4)))), desired=arr((3,4)), decimal=4 )
        eq( actual=arr(ap(3*π/4, (1,2)).intersection(ap(    0, (1,4)))), desired=arr((-1,4)), decimal=4 )
        eq( actual=arr(ap(5*π/4, (1,2)).intersection(ap(    0, (1,4)))), desired=arr((3,4)), decimal=4 )
        eq( actual=arr(ap(7*π/4, (1,2)).intersection(ap(    0, (1,4)))), desired=arr((-1,4)), decimal=4 )
        eq( actual=arr(ap(  π/4, (1,2)).intersection(ap(  π/2, (1,4)))), desired=arr((1,2)), decimal=4 )
        eq( actual=arr(ap(3*π/4, (1,2)).intersection(ap(  π/2, (1,4)))), desired=arr((1,2)), decimal=4 )
        eq( actual=arr(ap(5*π/4, (1,2)).intersection(ap(  π/2, (1,4)))), desired=arr((1,2)), decimal=4 )
        eq( actual=arr(ap(7*π/4, (1,2)).intersection(ap(  π/2, (1,4)))), desired=arr((1,2)), decimal=4 )

    def test_distance_along_line(self):
        sq2 = np.sqrt(1/2)
        eq = npt.assert_array_almost_equal
        arr = np.array
        si = Line.from_slope_intercept
        tp = Line.from_two_points

        # no offset
        eq( actual=arr(si(1,       0).distance_along_line(1)), desired=arr((sq2, sq2)),    decimal=4 )
        eq( actual=arr(si(-1,      0).distance_along_line(1)), desired=arr((sq2, -sq2)),   decimal=4 )
        eq( actual=arr(si(1/10,    0).distance_along_line(1)), desired=arr((.995, .0995)), decimal=4 )
        eq( actual=arr(si(10,      0).distance_along_line(1)), desired=arr((.0995, .995)), decimal=4 )
        eq( actual=arr(si(np.inf,  0).distance_along_line(1)), desired=arr((0, 1)),        decimal=4 )
        eq( actual=arr(si(-np.inf, 0).distance_along_line(1)), desired=arr((0, -1)),       decimal=4 )
        eq( actual=arr(si(0,       0).distance_along_line(1)), desired=arr((1, 0)),        decimal=4 )
        eq( actual=arr(si(-0,      0).distance_along_line(1)), desired=arr((1, 0)),        decimal=4 )
        eq( actual=arr(tp((0,0), (-10,0)).distance_along_line(1)), desired=arr((-1, 0)),   decimal=4 )
        
        # x or y intercept outside of zero
        eq( actual=arr(si(1,       2).distance_along_line(1)),     desired=arr((sq2, sq2+2)),   decimal=4 )
        eq( actual=arr(si(-1,      2).distance_along_line(1)),     desired=arr((sq2, -sq2+2)),  decimal=4 )
        eq( actual=arr(si(1/10,    2).distance_along_line(1)),     desired=arr((.995, 2.0995)), decimal=4 )
        eq( actual=arr(si(10,      2).distance_along_line(1)),     desired=arr((.0995, 2.995)), decimal=4 )
        eq( actual=arr(tp((2,0), (2,10)).distance_along_line(1)),  desired=arr((2, 1)),         decimal=4 )
        eq( actual=arr(tp((2,10), (2,0)).distance_along_line(1)),  desired=arr((2, -1)),        decimal=4 )
        eq( actual=arr(si(0,       2).distance_along_line(1)),     desired=arr((1, 2)),         decimal=4 )
        eq( actual=arr(si(-0,      2).distance_along_line(1)),     desired=arr((1, 2)),         decimal=4 )
        eq( actual=arr(tp((0,2), (-10,2)).distance_along_line(1)), desired=arr((-1, 2)),        decimal=4 )
    
    def test_get_tangent_line(self):
        eq = npt.assert_array_almost_equal
        arr = np.array
        yrr = lambda v: np.array((v.slope, v.y_intercept))
        xrr = lambda v: np.array((v.slope, v.x_intercept))
        trr = lambda v: np.array(((v.x1, v.y1), (v.x2, v.y2)))
        si = Line.from_slope_intercept
        tp = Line.from_two_points
        sq2 = np.sqrt(2)

        eq( actual=yrr( si(1, 0).get_tangent_line((0, 0)) ),            desired=arr((-1, 0)),             decimal=4 )
        eq( actual=yrr( si(-1, 0).get_tangent_line((0, 0)) ),           desired=arr((1, 0)),              decimal=4 )
        eq( actual=yrr( si(1/10, 0).get_tangent_line((0, 0)) ),         desired=arr((-10, 0)),            decimal=4 )
        eq( actual=yrr( si(10, 0).get_tangent_line((0, 0)) ),           desired=arr((-1/10, 0)),          decimal=4 )
        eq( actual=yrr( si(np.inf, 0).get_tangent_line((0, 0)) ),       desired=arr((0, 0)),              decimal=4 )
        eq( actual=yrr( si(-np.inf, 0).get_tangent_line((0, 0)) ),      desired=arr((0, 0)),              decimal=4 )
        eq( actual=xrr( si(0, 0).get_tangent_line((0, 0)) ),            desired=arr((-np.inf, 0)),        decimal=4 )
        eq( actual=xrr( si(-0, 0).get_tangent_line((0, 0)) ),           desired=arr((-np.inf, 0)),        decimal=4 )
        eq( actual=trr( tp((0,0), (-10,0)).get_tangent_line((0, 0)) ),  desired=arr(((0, 0), (0, 1))),   decimal=4 )
        eq( actual=trr( tp((0,0), (5,5)).get_tangent_line((0, 0)) ),    desired=arr(((0, 0), (1, -1))), decimal=4 )
        eq( actual=trr( tp((0,0), (-10,0)).get_tangent_line((-5, 0)) ), desired=arr(((-5, 0), (-5, 1))), decimal=4 )
        eq( actual=trr( tp((0,0), (5,5)).get_tangent_line((2, 2)) ),    desired=arr(((0, 4), (1, 3))),  decimal=4 )

    def test_line_two_points_to_slope_intercept(self):
        eq = npt.assert_array_almost_equal
        arr = np.array
        tp = Line.from_two_points
        to_si = lambda l: (l.slope, l.y_intercept)
        tx_si = lambda l: (l.slope, l.x_intercept)

        eq( actual=arr( to_si(tp((0, 0),     (10, 0)))    ), desired=arr((0, 0)),       decimal=4 )
        eq( actual=arr( to_si(tp((0, 0),     (-10, 0)))   ), desired=arr((0, 0)),       decimal=4 )
        eq( actual=arr( tx_si(tp((0, 0),     (0, 10)))    ), desired=arr((np.inf, 0)),  decimal=4 )
        eq( actual=arr( tx_si(tp((0, 0),     (0, -10)))   ), desired=arr((-np.inf, 0)), decimal=4 )
        eq( actual=arr( to_si(tp((0, 0),     (10, 10)))   ), desired=arr((1, 0)),       decimal=4 )
        eq( actual=arr( to_si(tp((0, 0),     (10, -10)))  ), desired=arr((-1, 0)),      decimal=4 )
        eq( actual=arr( to_si(tp((0, 0),     (-10, 10)))  ), desired=arr((-1, 0)),      decimal=4 )
        eq( actual=arr( to_si(tp((0, 0),     (-10, -10))) ), desired=arr((1, 0)),       decimal=4 )
        eq( actual=arr( to_si(tp((10, 10),   (0, 0)))     ), desired=arr((1, 0)),       decimal=4 )
        eq( actual=arr( to_si(tp((10, -10),  (0, 0)))     ), desired=arr((-1, 0)),      decimal=4 )
        eq( actual=arr( to_si(tp((-10, 10),  (0, 0)))     ), desired=arr((-1, 0)),      decimal=4 )
        eq( actual=arr( to_si(tp((-10, -10), (0, 0)))     ), desired=arr((1, 0)),       decimal=4 )

    def test_from_two_poitns(self):
        for i in range(100):
            x1 = random.random()*100
            y1 = random.random()*100
            x2, y2 = x1, y1
            while abs(x2 - x1) < 0.1:
                x2 = random.random()*100
            while abs(y2 - y1) < 0.1:
                y2 = random.random()*100
            
            l1 = Line.from_two_points((x1, y1), (x2, y2))
            l2 = Line.from_two_points((l1.x1, l1.y1), (l1.x2, l1.y2))

            self.assertAlmostEqual(l1.angle, l2.angle, places=5)
            self.assertAlmostEqual(l1.y_intercept, l2.y_intercept, places=5)
    
    def test_is_parallel_to(self):
        si = Line.from_slope_intercept
        ap = Line.from_angle_point
        π = np.pi

        self.assertFalse( si(0, 0).is_parallel_to(         si(np.inf, 0)    ) )
        self.assertFalse( si(0, 2.3).is_parallel_to(       si(1, 4.9)       ) )
        self.assertFalse( si(1, 2.3).is_parallel_to(       si(-1, 4.9)      ) )
        self.assertFalse( si(-1, 2.3).is_parallel_to(      si(1, 4.9)       ) )

        self.assertTrue(  si(0, 2.3).is_parallel_to(       si(0, 4.9)       ) )
        self.assertTrue(  si(1, 2.3).is_parallel_to(       si(1, 4.9)       ) )
        self.assertTrue(  si(-1, 2.3).is_parallel_to(      si(-1, 4.9)      ) )
        self.assertTrue(  si(np.inf, 0).is_parallel_to(    si(-np.inf, 0)   ) )
        self.assertTrue(  si(np.inf, 10).is_parallel_to(   si(np.inf, -5)   ) )
        self.assertTrue(  si(-np.inf, 2.3).is_parallel_to( si(np.inf, 4.9)  ) )
        self.assertTrue(  ap(π/4, (0,0)).is_parallel_to(   ap(π/4, (0,0))   ) )
        self.assertTrue(  ap(π/4, (0,0)).is_parallel_to(   ap(5*π/4, (0,0)) ) )
        self.assertTrue(  ap(5*π/4, (0,0)).is_parallel_to( ap(π/4, (0,0))   ) )

if __name__ == '__main__':
    unittest.main()
