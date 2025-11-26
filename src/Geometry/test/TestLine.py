import math
import random
import unittest

import numpy as np
import numpy.testing as npt

from Geometry.Line import Line
from Geometry.Point import Point as P2

class TestLine(unittest.TestCase):
    def test_xy_property(self):
        # Create a line with points (1, 2) and slope 2
        line = Line(P2(0, 0), P2(1, 2))

        self.assertAlmostEqual(line.xy.x, math.cos(np.atan2(2, 1)))
        self.assertAlmostEqual(line.xy.y, math.sin(np.atan2(2, 1)))

        # Create a horizontal line
        line = Line(P2(0, 5), P2(10, 5))

        self.assertAlmostEqual(line.xy.x, 1)
        self.assertAlmostEqual(line.xy.y, 0)

    def test_x_intercept(self):
        # Vertical line
        line = Line(P2(3, 0), P2(3, 10))
        self.assertEqual(line.x_intercept, 3)

        # Horizontal line
        line = Line(P2(0, 3), P2(10, 3))
        self.assertEqual(line.x_intercept, 0)

        # Ascending line (y = 2x + 1)
        line = Line(P2(10, 21), P2(20, 41))
        self.assertAlmostEqual(line.x_intercept, -0.5)

        # Ascending line (y = 2x + 1)
        line = Line(P2(-20, -39), P2(-10, -19))
        self.assertAlmostEqual(line.x_intercept, -0.5)

        # Descending line (y = -1/2x - 10)
        line = Line(P2(-30, 5), P2(-20, 0))
        self.assertAlmostEqual(line.x_intercept, -20)

        # Descending line (y = -1/2x - 10)
        line = Line(P2(-10, -5), P2(10, -15))
        self.assertAlmostEqual(line.x_intercept, -20)

    def test_y_intercept(self):
        # Vertical line
        line = Line(P2(3, 0), P2(3, 10))
        self.assertEqual(line.y_intercept, 0)

        # Horizontal line
        line = Line(P2(0, 3), P2(10, 3))
        self.assertEqual(line.y_intercept, 3)

        # Ascending line (y = 2x + 1)
        line = Line(P2(10, 21), P2(20, 41))
        self.assertAlmostEqual(line.y_intercept, 1)

        # Ascending line (y = 2x + 1)
        line = Line(P2(-20, -39), P2(-10, -19))
        self.assertAlmostEqual(line.y_intercept, 1)

        # Descending line (y = -1/2x - 10)
        line = Line(P2(-30, 5), P2(-20, 0))
        self.assertAlmostEqual(line.y_intercept, -10)

        # Descending line (y = -1/2x - 10)
        line = Line(P2(-10, -5), P2(10, -15))
        self.assertAlmostEqual(line.y_intercept, -10)

    def test_is_point_on_right(self):
        self.assertTrue(Line.from_vector(P2( 1,  0)).is_point_on_right(P2( 0, -1)))
        self.assertTrue(Line.from_vector(P2( 1,  1)).is_point_on_right(P2( 1,  0)))
        self.assertTrue(Line.from_vector(P2( 0,  1)).is_point_on_right(P2( 1,  1)))
        self.assertTrue(Line.from_vector(P2(-1,  1)).is_point_on_right(P2( 0,  1)))
        self.assertTrue(Line.from_vector(P2(-1,  0)).is_point_on_right(P2(-1,  1)))
        self.assertTrue(Line.from_vector(P2(-1, -1)).is_point_on_right(P2(-1,  0)))
        self.assertTrue(Line.from_vector(P2( 0, -1)).is_point_on_right(P2(-1, -1)))
        self.assertTrue(Line.from_vector(P2( 1, -1)).is_point_on_right(P2( 0, -1)))
        
        self.assertTrue( Line.from_vector(P2( 1,  0), y_intercept=1).is_point_on_right(P2(  0,  .5)))
        self.assertFalse(Line.from_vector(P2( 1,  0), y_intercept=1).is_point_on_right(P2(  0, 1.5)))
        self.assertTrue( Line.from_vector(P2( 0,  1), x_intercept=1).is_point_on_right(P2(1.5,   1)))
        self.assertFalse(Line.from_vector(P2( 0,  1), x_intercept=1).is_point_on_right(P2( .5,   1)))
        self.assertTrue( Line.from_vector(P2(-1,  0), y_intercept=1).is_point_on_right(P2( -1, 1.5)))
        self.assertFalse(Line.from_vector(P2(-1,  0), y_intercept=1).is_point_on_right(P2( -1,  .5)))
        self.assertTrue( Line.from_vector(P2( 0, -1), x_intercept=1).is_point_on_right(P2( .5,  -1)))
        self.assertFalse(Line.from_vector(P2( 0, -1), x_intercept=1).is_point_on_right(P2(1.5,  -1)))

    def test_intersection(self):
        arr = np.array
        ap = Line.from_angle_point
        eq = npt.assert_array_almost_equal
        π = np.pi
        intersection = lambda a1, p1, a2, p2: ap(a1, P2(*p1)).intersection(ap(a2, P2(*p2))).as_tuple()

        # 45 degree angles
        eq( actual=arr(intersection(  π/4, (1,2), 3*π/4, (1,4))), desired=arr((2,3)),  decimal=4 )
        eq( actual=arr(intersection(3*π/4, (1,2),   π/4, (1,4))), desired=arr((0,3)),  decimal=4 )
        eq( actual=arr(intersection(5*π/4, (1,2), 3*π/4, (1,4))), desired=arr((2,3)),  decimal=4 )
        eq( actual=arr(intersection(7*π/4, (1,2),   π/4, (1,4))), desired=arr((0,3)),  decimal=4 )

        # 90 degree angles
        eq( actual=arr(intersection(    0, (1,2),   π/2, (1,4))), desired=arr((1,2)),  decimal=4 )
        eq( actual=arr(intersection(  π/2, (1,2),     0, (1,4))), desired=arr((1,4)),  decimal=4 )
        eq( actual=arr(intersection(    π, (1,2),   π/2, (1,4))), desired=arr((1,2)),  decimal=4 )
        eq( actual=arr(intersection(3*π/2, (1,2),     0, (1,4))), desired=arr((1,4)),  decimal=4 )

        # 45 degree angle and 90 degree angle
        eq( actual=arr(intersection(  π/4, (1,2),     0, (1,4))), desired=arr((3,4)),  decimal=4 )
        eq( actual=arr(intersection(3*π/4, (1,2),     0, (1,4))), desired=arr((-1,4)), decimal=4 )
        eq( actual=arr(intersection(5*π/4, (1,2),     0, (1,4))), desired=arr((3,4)),  decimal=4 )
        eq( actual=arr(intersection(7*π/4, (1,2),     0, (1,4))), desired=arr((-1,4)), decimal=4 )
        eq( actual=arr(intersection(  π/4, (1,2),   π/2, (1,4))), desired=arr((1,2)),  decimal=4 )
        eq( actual=arr(intersection(3*π/4, (1,2),   π/2, (1,4))), desired=arr((1,2)),  decimal=4 )
        eq( actual=arr(intersection(5*π/4, (1,2),   π/2, (1,4))), desired=arr((1,2)),  decimal=4 )
        eq( actual=arr(intersection(7*π/4, (1,2),   π/2, (1,4))), desired=arr((1,2)),  decimal=4 )

    def test_distance_along_line(self):
        sq2 = np.sqrt(1/2)
        eq = npt.assert_array_almost_equal
        arr = np.array
        si = Line.from_slope_intercept
        tp = lambda p1, p2: Line.from_two_points(P2(*p1), P2(*p2))
        at_dist = lambda l: l.distance_along_line(1).as_tuple()

        # no offset
        eq( actual=arr(at_dist(si(      1, 0      ))), desired=arr((  sq2, sq2  )), decimal=4 )
        eq( actual=arr(at_dist(si(     -1, 0      ))), desired=arr((  sq2, -sq2 )), decimal=4 )
        eq( actual=arr(at_dist(si(   1/10, 0      ))), desired=arr(( .995, .0995)), decimal=4 )
        eq( actual=arr(at_dist(si(     10, 0      ))), desired=arr((.0995, .995 )), decimal=4 )
        eq( actual=arr(at_dist(si( np.inf, 0      ))), desired=arr((    0, 1    )), decimal=4 )
        eq( actual=arr(at_dist(si(-np.inf, 0      ))), desired=arr((    0, -1   )), decimal=4 )
        eq( actual=arr(at_dist(si(      0, 0      ))), desired=arr((    1, 0    )), decimal=4 )
        eq( actual=arr(at_dist(si(     -0, 0      ))), desired=arr((    1, 0    )), decimal=4 )
        eq( actual=arr(at_dist(tp(  (0,0), (-10,0)))), desired=arr((   -1, 0    )), decimal=4 )
        
        # x or y intercept outside of zero
        eq( actual=arr(at_dist(si(     1, 2      ))), desired=arr((  sq2, sq2+2 )), decimal=4 )
        eq( actual=arr(at_dist(si(    -1, 2      ))), desired=arr((  sq2, -sq2+2)), decimal=4 )
        eq( actual=arr(at_dist(si(  1/10, 2      ))), desired=arr(( .995, 2.0995)), decimal=4 )
        eq( actual=arr(at_dist(si(    10, 2      ))), desired=arr((.0995, 2.995 )), decimal=4 )
        eq( actual=arr(at_dist(tp( (2,0), (2,10) ))), desired=arr((    2, 1     )), decimal=4 )
        eq( actual=arr(at_dist(tp((2,10), (2,0)  ))), desired=arr((    2, -1    )), decimal=4 )
        eq( actual=arr(at_dist(si(     0, 2      ))), desired=arr((    1, 2     )), decimal=4 )
        eq( actual=arr(at_dist(si(    -0, 2      ))), desired=arr((    1, 2     )), decimal=4 )
        eq( actual=arr(at_dist(tp( (0,2), (-10,2)))), desired=arr((   -1, 2     )), decimal=4 )
    
    def test_get_tangent_line(self):
        eq = npt.assert_array_almost_equal
        arr = np.array
        yrr = lambda v: np.array((v.slope, v.y_intercept)) # slope + x intercept array
        xrr = lambda v: np.array((v.slope, v.x_intercept)) # slope + y intercept array
        trr = lambda v: np.array((v.xy0.as_tuple(), v.xy1.as_tuple())) # tuple array
        si = Line.from_slope_intercept
        tp = lambda t1, t2: Line.from_two_points(P2(*t1), P2(*t2))

        eq( actual=yrr( si(      1, 0).get_tangent_line(     P2(0, 0))),  desired=arr((-1, 0)),              decimal=4 )
        eq( actual=yrr( si(     -1, 0).get_tangent_line(     P2(0, 0))),  desired=arr((1, 0)),               decimal=4 )
        eq( actual=yrr( si(   1/10, 0).get_tangent_line(     P2(0, 0))),  desired=arr((-10, 0)),             decimal=4 )
        eq( actual=yrr( si(     10, 0).get_tangent_line(     P2(0, 0))),  desired=arr((-1/10, 0)),           decimal=4 )
        eq( actual=yrr( si( np.inf, 0).get_tangent_line(     P2(0, 0))),  desired=arr((0, 0)),               decimal=4 )
        eq( actual=yrr( si(-np.inf, 0).get_tangent_line(     P2(0, 0))),  desired=arr((0, 0)),               decimal=4 )
        eq( actual=xrr( si(      0, 0).get_tangent_line(     P2(0, 0))),  desired=arr((-np.inf, 0)),         decimal=4 )
        eq( actual=xrr( si(     -0, 0).get_tangent_line(     P2(0, 0))),  desired=arr((-np.inf, 0)),         decimal=4 )
        eq( actual=trr( tp((0,0), (-10,0)).get_tangent_line( P2(0, 0))),  desired=arr(( (0, 0), (0, 10)  )), decimal=4 )
        eq( actual=trr( tp((0,0), (5,5)  ).get_tangent_line( P2(0, 0))),  desired=arr(( (0, 0), (10, -10))), decimal=4 )
        eq( actual=trr( tp((0,0), (-10,0)).get_tangent_line( P2(-5, 0))), desired=arr(((-5, 0), (-5, 10) )), decimal=4 )
        eq( actual=trr( tp((0,0), (5,5)  ).get_tangent_line( P2(2, 2))),  desired=arr(( (0, 4), (10, -6) )), decimal=4 )

    def test_line_two_points_to_slope_intercept(self):
        eq = npt.assert_array_almost_equal
        arr = np.array
        tp = lambda t1, t2: Line.from_two_points(P2(*t1), P2(*t2))
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

    def test_from_two_points(self):
        for i in range(100):
            x1 = random.random()*100
            y1 = random.random()*100
            x2, y2 = x1, y1
            while abs(x2 - x1) < 0.1:
                x2 = random.random()*100
            while abs(y2 - y1) < 0.1:
                y2 = random.random()*100
            
            l1 = Line.from_two_points(P2(x1, y1), P2(x2, y2))
            l2 = Line.from_two_points(l1.xy0, l1.xy1)

            self.assertAlmostEqual(l1.angle, l2.angle, places=5)
            self.assertAlmostEqual(l1.y_intercept, l2.y_intercept, places=5)
    
    def test_is_parallel_to(self):
        si = Line.from_slope_intercept
        ap = lambda a, t: Line.from_angle_point(a, P2(*t))
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

    def test_intersection(self):
        line1 = Line.from_slope_intercept(1, 0)
        line2 = Line.from_slope_intercept(-1, 5)

        intersection_point = line1.intersection(line2)
        self.assertEqual(intersection_point.x, 2.5)
        self.assertEqual(intersection_point.y, 2.5)

        # Test parallel lines
        line3 = Line.from_slope_intercept(1, 2)
        self.assertIsNone(line1.intersection(line3))

    def test_distance_along_line(self):
        line = Line.from_slope_intercept(1.0, 0.0)

        point = line.distance_along_line(5)
        expected_from_0 = P2(np.sqrt(25/2), np.sqrt(25/2))
        self.assertAlmostEqual(point.x, expected_from_0.x)
        self.assertAlmostEqual(point.y, expected_from_0.y)

        point = line.distance_along_line(5, expected_from_0)
        expected_from_5 = P2(np.sqrt(25/2)*2, np.sqrt(25/2)*2)
        self.assertAlmostEqual(point.x, expected_from_5.x)
        self.assertAlmostEqual(point.y, expected_from_5.y)

        point = line.distance_along_line(5, P2(-5, 5))
        expected_from_offset = P2(np.sqrt(25/2)-5, np.sqrt(25/2)+5)
        self.assertAlmostEqual(point.x, expected_from_offset.x)
        self.assertAlmostEqual(point.y, expected_from_offset.y)

    def test_get_tangent_line(self):
        line = Line.from_slope_intercept(1.0, 0.0)

        tangent = line.get_tangent_line(P2(3, 3))
        self.assertEqual(tangent.slope, -1)
        self.assertAlmostEqual(tangent.y_intercept, 6)

    def test_reversed(self):
        # line with positive slope
        line = Line.from_slope_intercept(2, 0)
        reversed_line = line.reversed()
        self.assertAlmostEqual(reversed_line.run, -line.run)
        self.assertAlmostEqual(reversed_line.rise, -line.rise)
        
        # line with negative slope
        line = Line.from_slope_intercept(-2, 0)
        reversed_line = line.reversed()
        self.assertAlmostEqual(reversed_line.run, -line.run)
        self.assertAlmostEqual(reversed_line.rise, -line.rise)
        
        # horizontal line
        line = Line.from_two_points(P2(0, 0), P2(1, 0))
        reversed_line = line.reversed()
        self.assertAlmostEqual(reversed_line.run, -line.run)
        self.assertAlmostEqual(reversed_line.rise, -line.rise)
        
        # vertical line
        line = Line.from_two_points(P2(0, 0), P2(0, 1))
        reversed_line = line.reversed()
        self.assertAlmostEqual(reversed_line.run, -line.run)
        self.assertAlmostEqual(reversed_line.rise, -line.rise)

    def test_angle_between(self):
        # 90-degree angle for 45-degree lines
        line1 = Line.from_slope_intercept(1, 0)
        line2 = Line.from_slope_intercept(-1, 0)
        angle = line1.angle_between(line2)
        self.assertEqual(angle, np.pi/2)

        # 90-degree angle for horizontal and vertical lines
        line3 = Line.from_angle_point(0, P2(0, 0))
        line4 = Line.from_angle_point(np.pi/2, P2(0, 0))
        angle = line3.angle_between(line4)
        self.assertEqual(angle, np.pi/2)

        # 0-degree angle for parallel lines
        line5 = Line.from_angle_point(0, P2(0, 0))
        line6 = Line.from_angle_point(0, P2(1, 1))
        angle = line5.angle_between(line6)
        self.assertEqual(angle, 0)

if __name__ == '__main__':
    unittest.main()
