import unittest

import numpy as np
import numpy.testing as npt

from Geometry.LineSegment import LineSegment

class TestLineSegment(unittest.TestCase):
    def setUp(self):
        self.line_segment = LineSegment((0, 0), (3, 3))

    def test_length(self):
        """Test the length of a line segment."""
        self.assertAlmostEqual(self.line_segment.length, np.sqrt(9+9))

    def test_intersection_with(self):
        """Test intersection between two line segments."""
        another_line = LineSegment((-3, 3), (3, -3))
        expected_intersection = (0, 0)
        npt.assert_array_almost_equal(
            np.array(self.line_segment.intersection(another_line)),
            np.array(expected_intersection))
        
        another_line = LineSegment((-3, 5), (3, -1))
        expected_intersection = (1, 1)
        npt.assert_array_almost_equal(
            np.array(self.line_segment.intersection(another_line)),
            np.array(expected_intersection))
        
        another_line = LineSegment((-3, 4), (1, 1))
        expected_intersection = (1, 1)
        npt.assert_array_almost_equal(
            np.array(self.line_segment.intersection(another_line)),
            np.array(expected_intersection))
        
        another_line = LineSegment((1, 1), (1, 1))
        expected_intersection = (1, 1)
        npt.assert_array_almost_equal(
            np.array(self.line_segment.intersection(another_line)),
            np.array(expected_intersection))
        
        another_line = LineSegment((-3, -3), (3, 3))
        self.assertIsNone(self.line_segment.intersection(another_line))
        
        another_line = LineSegment((-3, 4), (.5, 1.5))
        self.assertIsNone(self.line_segment.intersection(another_line))

if __name__ == '__main__':
    unittest.main()