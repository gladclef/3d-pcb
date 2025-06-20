import copy
import re

import matplotlib.axis as maxis
import matplotlib.pyplot as plt
import numpy as np

from tool.units import *
import Geometry.geometry_tools as geo

class Circle:
    """ Represents a circle used in the outline of a component shape. """

    def __init__(self, center_point: tuple[float, float], radius: float):
        self.center_point = center_point
        self.radius = radius

    def apply_translation_rotation_layer(self, translation: tuple[float, float], rotation: float, is_bottom: bool) -> "Circle":
        """
        Apply translation and rotation to the circle, with optional flipping.

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
        Circle
            A new Circle instance with the applied transformations.
        """
        center_point = geo.apply_translation_rotation_flip(self.center_point, translation, rotation, is_bottom)

        ret = copy.deepcopy(self)
        ret.center_point = center_point

        return ret

    @classmethod
    def from_cad_file(cls, lines: list[str]) -> tuple["Circle", list[str]]:
        """
        Create a Circle instance by parsing lines from a CAD file.

        Parameters
        ----------
        lines : list[str]
            Lines of text representing the Circle in a CAD file.

        Returns
        -------
        tuple[Circle, list[str]]
            A tuple containing the created Circle and any remaining unprocessed lines.

        Raises
        ------
        RuntimeError
            If the Circle line does not match the expected pattern.
        """
        # example Circle line:
        #     CIRCLE 0.05 0 0.0984252

        circle_line = None
        for line_idx, line in enumerate(lines):
            if line.startswith("CIRCLE"):
                circle_line = line
                ret_lines = lines[:line_idx] + lines[line_idx+1:]
                break
        if circle_line is None:
            return None, lines

        # regex explanation:                 center_x   center_y   radius
        circle_pattern = re.compile(r"CIRCLE ([-\d\.]+) ([-\d\.]+) ([-\d\.]+)")

        # break out the various parts of the circle
        match = circle_pattern.match(circle_line)
        if match is None:
            raise RuntimeError("Error in Circle.from_cad_file(): failed to match circle_pattern to line:\n\t" + circle_line)

        center_x, center_y, radius = match.groups()
        center_x, center_y = in2mm(float(center_x)), in2mm(float(center_y))
        radius = in2mm(float(radius))

        circle = cls((center_x, center_y), radius)
        return circle, ret_lines

    def draw(self, ax: maxis.Axis):
        ax.add_patch(plt.Circle(self.center_point, self.radius, edgecolor="grey", facecolor="None"))
