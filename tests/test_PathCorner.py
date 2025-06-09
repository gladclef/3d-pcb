import copy
import json
import os
import unittest

import numpy as np
from scipy.spatial.transform import Rotation

from PipeShape import PipeBasicBox
from SingleTrace import SingleTrace
from PathCorner import PathCorner
from units import *

class TestPathCorner(unittest.TestCase):
    def _tst_arc_center(self, xy_points: list[tuple[float, float]], expected_centers: list[tuple[float, float]], debug=False):
        xy_points_original = copy.deepcopy(xy_points)
        expected_centers_original = copy.deepcopy(expected_centers)

        # Helper function to rotate a point (x,y) by angle (radians) about the origin.
        def rotate_point(point: list[float], angle: float) -> list[float]:
            x, y = tuple(point)
            r = Rotation.from_euler('z', angle)
            xy_rotated: np.ndarray = r.apply(np.array([x, y, 0]))
            return xy_rotated[0], xy_rotated[1]
        
        # Helper function for more readable printing of float arrays
        def print_pnts(points: list[float] | list[list]) -> str:
            if isinstance(points[0], tuple) or isinstance(points[0], list):
                return "[" + ",".join([print_pnts(pnts) for pnts in points]) + "]"
            else:
                return "[" + ",".join([f"{p:0.3f}" for p in points]) + "]"

        for inverted in [False, True]:
            xy_points_inverted = copy.deepcopy(xy_points_original)
            if inverted:
                xy_points_inverted = [
                    [xy_points_inverted[2][0], xy_points_inverted[2][1]],
                    [xy_points_inverted[1][0], xy_points_inverted[1][1]],
                    [xy_points_inverted[0][0], xy_points_inverted[0][1]],
                ] 

            for rotation in range(7):
                angle = rotation * np.pi / 4
                xy_points = [rotate_point(pnt, angle) for pnt in xy_points_inverted]
                expected_centers = [rotate_point(exp_cent, angle) for exp_cent in expected_centers_original]

                segments = [
                    [0, 1],
                    [1, 2]
                ]
                shape = PipeBasicBox(awg2mm(26))
                trace = SingleTrace(xy_points, segments, shape)
                corner = PathCorner(trace, trace.segments, 1)
                corner.debug = debug

                angle_diff, mid_angle, arc_length, arc_center = corner.get_arc_properties()
                # print(f"{angle_diff=}\n{mid_angle=}\n{arc_length=}\n{arc_center=}")

                # find one of the expected center options that matches, if any
                delta = 0.1
                for exp_cent in expected_centers:
                    if abs(exp_cent[0] - arc_center[0]) < delta and abs(exp_cent[1] - arc_center[1]) < 0.1:
                        break

                self.assertAlmostEqual(arc_center[0], exp_cent[0], msg=f"\n\texpected {print_pnts(exp_cent)}\n\tactual {print_pnts(arc_center)}\n\t({angle=}, {inverted=}, xy_points={print_pnts(xy_points)})", delta=delta)
                self.assertAlmostEqual(arc_center[1], exp_cent[1], msg=f"\n\texpected {print_pnts(exp_cent)}\n\tactual {print_pnts(arc_center)}\n\t({angle=}, {inverted=}, xy_points={print_pnts(xy_points)})", delta=delta)
    
    def test_arc_center_straight(self):
        xy_points = [
            [0, 0],
            [5, 0],
            [10, 0]
        ]
        expected_centers = ((5, -1), (5, 1))
        self._tst_arc_center(xy_points, expected_centers)

    def test_arc_center_45_degree(self):
        xy_points = [
            [0, 0],
            [5, 0],
            [0, -5]
        ]
        expected_center = [(2.586, -1)]
        self._tst_arc_center(xy_points, expected_center)
        
    def test_arc_center_90_degree(self):
        xy_points = [
            [0, 0],
            [5, 0],
            [5, -5]
        ]
        expected_center = [(4, -1)]
        self._tst_arc_center(xy_points, expected_center)

    def test_arc_center_135_degree(self):
        xy_points = [
            [0, 0],
            [5, 0],
            [10, -5]
        ]
        expected_center = [(4.586, -1)]
        self._tst_arc_center(xy_points, expected_center)
        
    def test_arc_center_225_degree(self):
        xy_points = [
            [0, 0],
            [5, 0],
            [10, 5]
        ]
        expected_center = [(4.586, 1)]
        self._tst_arc_center(xy_points, expected_center)
        
    def test_arc_center_270_degree(self):
        xy_points = [
            [0, 0],
            [5, 0],
            [5, 5]
        ]
        expected_center = [(4, 1)]
        self._tst_arc_center(xy_points, expected_center)
        
    def test_arc_center_315_degree(self):
        xy_points = [
            [0, 0],
            [5, 0],
            [0, 5]
        ]
        expected_center = [(2.586, 1)]
        self._tst_arc_center(xy_points, expected_center)

    def test_center_line_points_regression(self):
        """ Check that the values we got before are the same as the values now. """
        with open(os.path.join(os.path.dirname(__file__), "test_PathCorners_regression_vals.json"), "r") as fin:
            regression_values = json.load(fin)

        for reg_val in regression_values:
            deg100: float = reg_val["deg100"]
            deg: float = reg_val["deg"]
            rad: float = reg_val["rad"]
            angle_diff: float = reg_val["angle_diff"]
            mid_angle: float = reg_val["mid_angle"]
            arc_length: float = reg_val["arc_length"]
            arc_center: tuple[float, float] = reg_val["arc_center"]
            center_line_points: list[tuple[ tuple[float, float], float ]] = reg_val["center_line_points"]

            xy_points = [
                [-3.5, 0],
                [0, 0],
                [0+np.cos(rad)*3.5, 0+np.sin(rad)*3.5]
            ]
            segments = [
                [0, 1],
                [1, 2]
            ]
            shape = PipeBasicBox(awg2mm(26))
            trace = SingleTrace(xy_points, segments, shape)
            corner = PathCorner(trace, trace.segments, 1)

            for i, clp in enumerate(corner.get_center_line_points()):
                (act_x, act_y), act_rad = clp
                (exp_x, exp_y), exp_rad = center_line_points[i]
                self.assertAlmostEqual(act_x, exp_x, delta=0.01, msg=f"Actual center line point {clp} does not match expected {center_line_points[i]} for degrees {deg}")
                self.assertAlmostEqual(act_y, exp_y, delta=0.01, msg=f"Actual center line point {clp} does not match expected {center_line_points[i]} for degrees {deg}")
                self.assertAlmostEqual(act_rad, exp_rad, delta=0.01, msg=f"Actual center line point {clp} does not match expected {center_line_points[i]} for degrees {deg}")

if __name__ == '__main__':
    unittest.main()
