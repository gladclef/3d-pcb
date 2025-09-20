import copy
import re

import matplotlib.axis as maxis
import matplotlib.patches as mpatches
import numpy as np

from FileIO.Line import Line as FLine
from Geometry.LineSegment import LineSegment
from tool.units import *
import Geometry.geometry_tools as geo

class Arc:
    """ Represents an arc used in the outline of a component shape. """

    def __init__(self, start_point: tuple[float, float], end_point: tuple[float, float], center_point: tuple[float, float]):
        self.start_point = start_point
        self.end_point = end_point
        self.center_point = center_point

    @property
    def radius(self) -> float:
        return np.sqrt((self.start_point[0] - self.center_point[0])**2 + (self.start_point[1] - self.center_point[1])**2)

    @property
    def start_angle(self) -> float:
        seg = LineSegment(self.center_point, self.start_point)
        return seg.angle
    
    @property
    def end_angle(self) -> float:
        seg = LineSegment(self.center_point, self.end_point)
        return seg.angle

    def apply_translation_rotation_layer(self, translation: tuple[float, float], rotation: float, is_bottom: bool) -> "Arc":
        """
        Apply translation and rotation to the arc, with optional flipping.

        Parameters
        ----------
        translation : tuple[float, float]
            Translation vector (x, y).
        rotation : float
            Rotation angle in radians.
        is_bottom : bool
            If True, flip the coordinates along the x-axis after applying translation and rotation.

        Returns
        -------
        Arc
            A new Arc instance with the applied transformations.
        """
        start_point = geo.apply_translation_rotation_flip(self.start_point, translation, rotation, is_bottom)
        end_point = geo.apply_translation_rotation_flip(self.end_point, translation, rotation, is_bottom)
        center_point = geo.apply_translation_rotation_flip(self.center_point, translation, rotation, is_bottom)

        ret = copy.deepcopy(self)
        ret.start_point, ret.end_point, ret.center_point = start_point, end_point, center_point

        return ret

    @classmethod
    def from_cad_file(cls, lines: list[FLine]) -> tuple[list["Arc"], list[FLine]]:
        """
        Create an Arc instance by parsing lines from a CAD file.

        Parameters
        ----------
        lines : list[FLine]
            Lines of text representing the Arc in a CAD file.

        Returns
        -------
        tuple[Arc, list[FLine]]
            A tuple containing the created Arc and any remaining unprocessed lines.

        Raises
        ------
        RuntimeError
            If the Arc line does not match the expected pattern.
        """
        # example Arc line:
        #     ARC -0.0507874 -0.0608201 0.167717 0 0.05 0

        arc_line = None
        for line_idx, line in enumerate(lines):
            if line.v.startswith("ARC"):
                arc_line = line
                ret_lines = lines[:line_idx] + lines[line_idx+1:]
                break
        if arc_line is None:
            return [], lines

        # regex explanation:           start_x    start_y    end_x      end_y      center_x   center_y
        arc_pattern = re.compile(r"ARC ([-\d\.]+) ([-\d\.]+) ([-\d\.]+) ([-\d\.]+) ([-\d\.]+) ([-\d\.]+)")

        # break out the various parts of the arc
        match = arc_pattern.match(arc_line.v)
        if match is None:
            raise RuntimeError("Error in Arc.from_cad_file(): failed to match arc_pattern to line:\n\t" + arc_line.v)

        start_x, start_y, end_x, end_y, center_x, center_y = match.groups()
        start_x, start_y = in2mm(float(start_x)), in2mm(float(start_y))
        end_x, end_y = in2mm(float(end_x)), in2mm(float(end_y))
        center_x, center_y = in2mm(float(center_x)), in2mm(float(center_y))

        arc = cls((start_x, start_y), (end_x, end_y), (center_x, center_y))
        return [arc], ret_lines

    def draw(self, ax: maxis.Axis):
        start_angle = np.rad2deg(self.start_angle)
        end_angle = np.rad2deg(self.end_angle)
        ax.add_patch(mpatches.Arc(self.center_point,
                                  self.radius*2,
                                  self.radius*2,
                                  angle=start_angle,
                                  theta1=0.0,
                                  theta2=end_angle-start_angle,
                                  edgecolor="grey",
                                  facecolor="None"))