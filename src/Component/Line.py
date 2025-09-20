import copy
import re

import matplotlib.axis as maxis
import numpy as np

from tool.units import *
from FileIO.Line import Line as FLine
import Geometry.geometry_tools as geo

class Line:
    """ Represents a line used in the outline of a component shape. """

    def __init__(self, xy1: tuple[float, float], xy2: tuple[float, float]):
        self.xy1 = xy1
        self.xy2 = xy2

    def apply_translation_rotation_layer(self, translation: tuple[float, float], rotation: float, is_bottom: bool) -> "Line":
        """
        Apply translation and rotation to the line, with optional flipping.

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
        Line
            A new Line instance with the applied transformations.
        """
        xy1 = geo.apply_translation_rotation_flip(self.xy1, translation, rotation, is_bottom)
        xy2 = geo.apply_translation_rotation_flip(self.xy2, translation, rotation, is_bottom)

        ret = copy.deepcopy(self)
        ret.xy1, ret.xy2 = xy1, xy2

        return ret

    @classmethod
    def from_cad_file(cls, lines: list[FLine]) -> tuple[list["Line"], list[FLine]]:
        """
        Create a Line instance by parsing lines from a CAD file.

        Parameters
        ----------
        lines : list[FLine]
            Lines of text representing the Line in a CAD file.

        Returns
        -------
        tuple[Line, list[FLine]]
            A tuple containing the created Line and any remaining unprocessed lines.

        Raises
        ------
        RuntimeError
            If the Line line does not match the expected pattern.

        """
        # example Line line:
        #     LINE 0.00511811 -0.133858 0.00511811 -0.0472441

        line_line = None
        for line_idx, line in enumerate(lines):
            if line.v.startswith("LINE"):
                line_line = line
                ret_lines = lines[:line_idx] + lines[line_idx+1:]
                break
        if line_line is None:
            return [], lines
        
        # regex explanation:             x1         y1         x2         y2
        line_pattern = re.compile(r"LINE ([-\d\.]+) ([-\d\.]+) ([-\d\.]+) ([-\d\.]+)")

        # break out the various parts of the line
        match = line_pattern.match(line_line.v)
        if match is None:
            raise RuntimeError("Error in Line.from_cad_file(): failed to match line_pattern to line:\n\t" + line_line.v)

        x1, y1, x2, y2 = match.groups()
        x1 = in2mm(float(x1))
        y1 = in2mm(float(y1))
        x2 = in2mm(float(x2))
        y2 = in2mm(float(y2))

        line = cls((x1, y1), (x2, y2))
        return [line], ret_lines

    def draw(self, ax: maxis.Axis):
        x = (self.xy1[0], self.xy2[0])
        y = (self.xy1[1], self.xy2[1])
        ax.plot(x, y, color="grey")