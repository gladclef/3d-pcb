import copy
import json
import os
import unittest

import numpy as np
from scipy.spatial.transform import Rotation

from FileIO.Line import Line as FLine
from Geometry.PathCorner import PathCorner
from Geometry.Point import Point as P2
from Trace.PipeShape import PipeBasicCircle
from Trace.SingleTrace import SingleTrace, _TraceLine
from tool.units import *

class TestPathCorner(unittest.TestCase):
    def _tst_arc_center(self, xy_points: list[P2], expected_centers: list[P2], debug=False):
        xy_points_original = copy.deepcopy(xy_points)
        expected_centers_original = copy.deepcopy(expected_centers)

        # Helper function to rotate a point (x,y) by angle (radians) about the origin.
        def rotate_point(point: P2, angle: float) -> list[float]:
            r = Rotation.from_euler('z', angle)
            xy_rotated: np.ndarray = r.apply(np.array([point.x, point.y, 0]))
            return P2(xy_rotated[0], xy_rotated[1])
        
        for inverted in [False, True]:
            xy_points_inverted = copy.copy(xy_points_original)
            if inverted:
                xy_points_inverted = list(reversed(xy_points_inverted))

            for rotation in range(7):
                angle = rotation * np.pi / 4
                source_lines = [FLine("N/A", i, "N/A") for i in range(len(xy_points_inverted)-1)]
                xy_points = [rotate_point(pnt, angle) for pnt in xy_points_inverted]
                expected_centers = [rotate_point(exp_cent, angle) for exp_cent in expected_centers_original]

                segments = [
                    _TraceLine(source_lines[0], 0, 1, is_end=True),
                    _TraceLine(source_lines[1], 1, 2, is_end=True)
                ]
                shape = PipeBasicCircle(awg2mm(26))
                trace = SingleTrace(source_lines, "TOP", xy_points, segments, shape)
                corner = PathCorner(trace, trace.segments, 1)
                corner.debug = debug

                angle_diff, mid_angle, arc_length, arc_center = corner.get_arc_properties()
                # print(f"{angle_diff=}\n{mid_angle=}\n{arc_length=}\n{arc_center=}")

                # find one of the expected center options that matches, if any
                delta = 0.1
                for exp_cent in expected_centers:
                    if abs(exp_cent.x - arc_center.x) < delta and abs(exp_cent.y - arc_center.y) < 0.1:
                        break

                self.assertAlmostEqual(arc_center.x, exp_cent.x, msg=f"\n\texpected {exp_cent}\n\tactual {arc_center}\n\t({angle=}, {inverted=}, xy_points={xy_points})", delta=delta)
                self.assertAlmostEqual(arc_center.y, exp_cent.y, msg=f"\n\texpected {exp_cent}\n\tactual {arc_center}\n\t({angle=}, {inverted=}, xy_points={xy_points})", delta=delta)
    
    def test_arc_center_straight(self):
        xy_points = [
            P2(0, 0),
            P2(5, 0),
            P2(10, 0)
        ]
        expected_centers = [P2(5, -1), P2(5, 1)]
        self._tst_arc_center(xy_points, expected_centers)

    def test_arc_center_45_degree(self):
        xy_points = [
            P2(0, 0),
            P2(5, 0),
            P2(0, -5)
        ]
        expected_center = [P2(2.586, -1)]
        self._tst_arc_center(xy_points, expected_center)
        
    def test_arc_center_90_degree(self):
        xy_points = [
            P2(0, 0),
            P2(5, 0),
            P2(5, -5)
        ]
        expected_center = [P2(4, -1)]
        self._tst_arc_center(xy_points, expected_center)

    def test_arc_center_135_degree(self):
        xy_points = [
            P2(0, 0),
            P2(5, 0),
            P2(10, -5)
        ]
        expected_center = [P2(4.586, -1)]
        self._tst_arc_center(xy_points, expected_center)
        
    def test_arc_center_225_degree(self):
        xy_points = [
            P2(0, 0),
            P2(5, 0),
            P2(10, 5)
        ]
        expected_center = [P2(4.586, 1)]
        self._tst_arc_center(xy_points, expected_center)
        
    def test_arc_center_270_degree(self):
        xy_points = [
            P2(0, 0),
            P2(5, 0),
            P2(5, 5)
        ]
        expected_center = [P2(4, 1)]
        self._tst_arc_center(xy_points, expected_center)
        
    def test_arc_center_315_degree(self):
        xy_points = [
            P2(0, 0),
            P2(5, 0),
            P2(0, 5)
        ]
        expected_center = [P2(2.586, 1)]
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
            arc_center: P2 = reg_val["arc_center"]
            center_line_points: list[tuple[ tuple[float, float], float ]] = reg_val["center_line_points"]
            center_line_points = [(P2(*clp[0]), clp[1]) for clp in center_line_points]

            xy_points = [
                P2(-3.5, 0),
                P2(0, 0),
                P2(0+np.cos(rad)*3.5, 0+np.sin(rad)*3.5)
            ]
            source_lines = [FLine("N/A", i, "N/A") for i in range(len(xy_points)-1)]
            segments = [
                _TraceLine(source_lines[0], 0, 1, is_end=True),
                _TraceLine(source_lines[1], 1, 2, is_end=True)
            ]
            shape = PipeBasicCircle(awg2mm(26))
            trace = SingleTrace(source_lines, "TOP", xy_points, segments, shape)
            corner = PathCorner(trace, trace.segments, 1)

            for i, clp in enumerate(corner.get_center_line_points()):
                act_pnt, act_rad = clp
                exp_pnt, exp_rad = center_line_points[i]
                self.assertAlmostEqual(act_pnt.x, exp_pnt.x, delta=0.01, msg=f"Actual center line point {clp} does not match expected {center_line_points[i]} for degrees {deg}")
                self.assertAlmostEqual(act_pnt.y, exp_pnt.y, delta=0.01, msg=f"Actual center line point {clp} does not match expected {center_line_points[i]} for degrees {deg}")
                self.assertAlmostEqual(act_rad, exp_rad, delta=0.01, msg=f"Actual center line point {clp} does not match expected {center_line_points[i]} for degrees {deg}")

if __name__ == '__main__':
    unittest.main()
