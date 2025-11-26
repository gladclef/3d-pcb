import unittest

import numpy as np
import numpy.testing as npt

from Geometry.LineSegment import LineSegment
from Geometry.Point import Point as P2

class TestLineSegment(unittest.TestCase):
    def setUp(self):
        self.line_segment = LineSegment(P2(0, 0), P2(3, 3))

    def test_length(self):
        """Test the length of a line segment."""
        self.assertAlmostEqual(self.line_segment.length, np.sqrt(9+9))

    def test_intersection_with(self):
        """Test intersection between two line segments."""
        another_line = LineSegment(P2(-3, 3), P2(3, -3))
        expected_intersection = P2(0, 0)
        npt.assert_array_almost_equal(
            np.array(self.line_segment.intersection(another_line).as_tuple()),
            np.array(expected_intersection.as_tuple()))
        
        another_line = LineSegment(P2(-3, 5), P2(3, -1))
        expected_intersection = P2(1, 1)
        npt.assert_array_almost_equal(
            np.array(self.line_segment.intersection(another_line).as_tuple()),
            np.array(expected_intersection.as_tuple()))
        
        another_line = LineSegment(P2(-3, 4), P2(1, 1))
        expected_intersection = P2(1, 1)
        npt.assert_array_almost_equal(
            np.array(self.line_segment.intersection(another_line).as_tuple()),
            np.array(expected_intersection.as_tuple()))
        
        another_line = LineSegment(P2(0, 1), P2(2, 1))
        expected_intersection = P2(1, 1)
        npt.assert_array_almost_equal(
            np.array(self.line_segment.intersection(another_line).as_tuple()),
            np.array(expected_intersection.as_tuple()))
        
        another_line = LineSegment(P2(-3, -3), P2(3, 3))
        self.assertIsNone(self.line_segment.intersection(another_line))
        
        another_line = LineSegment(P2(-3, 4), P2(.5, 1.5))
        self.assertIsNone(self.line_segment.intersection(another_line))

if __name__ == '__main__':
    unittest.main()