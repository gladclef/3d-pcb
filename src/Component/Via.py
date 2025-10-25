import copy
import re

import matplotlib.axis as maxis
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation
import vtk

from Component.DrillHole import DrillHole
from Geometry.Point import Point
from tool.units import *
from FileIO.Line import Line as FLine
import Geometry.geometry_tools as geo
from Trace.VtkPointGroup import VtkPointGroup
from tool.globals import board_parameters as g
import tool.vtk_tools as vt

class Via(DrillHole):
    """Represents a single via for a trace."""

    def __init__(self, location: Point, name: str="?", source_lines: list[FLine] = None):
        """
        Initialize the Via instance.

        Parameters
        ----------
        location : Point
            The location of the via relative to the board origin.

        """
        super().__init__(location, source_lines=source_lines, is_via=True)

        self.name = name

    def apply_translation_rotation_layer(self, translation: tuple[float, float] | Point, rotation: float, is_bottom: bool) -> "Via":
        """
        Apply translation and rotation to the via, with optional flipping.

        Parameters
        ----------
        translation : tuple[float, float] | Point
            Translation vector (x, y).
        rotation : float
            Rotation angle in radians.
        is_bottom : bool
            If True, flip the coordinates along the x-axis after applying translation and rotation.

        Returns
        -------
        Via
            A new Via instance with the applied transformations.

        """
        ret = copy.deepcopy(self)

        super(self.__class__, ret).apply_translation_rotation_layer(translation, rotation, is_bottom)

        return ret

    @classmethod
    def from_cad_file(cls, cad_lines: list[FLine]) -> tuple[list["Via"], list[FLine]]:
        """
        Create a Via instance by parsing lines from a CAD file.

        Parameters
        ----------
        lines : list[FLine]
            Lines of text representing the via in a CAD file.

        Returns
        -------
        tuple[list[Via], list[FLine]]
            A tuple containing the created Via(s) and any remaining unprocessed lines.

        Raises
        ------
        RuntimeError
            If the via line does not match the expected pattern.

        """
        # find the next via line, if any
        pre_via, via, post_via = [], None, []
        for i, line in enumerate(cad_lines):
            if line.v.startswith("VIA "):
                pre_via = cad_lines[:i]
                via = line
                post_via = cad_lines[i+1:]
                break
        if via is None:
            return [], cad_lines

        # Example via line:
        #     ?          ?       ?       x     y      ?   ?         name
        # VIA VIA1200000.1000000.1010101 2.125 -6.275 ALL 0.0393701 via3
        #                              ?            ?            ?           x            y            ?        ?            name
        via_pattern = re.compile(r"VIA ([0-9A-Z]+)\.([0-9A-Z]+)\.([0-9A-Z]+) (-?[0-9\.]+) (-?[0-9\.]+) ([A-Z]+) (-?[0-9\.]+) (.*)")
        m = via_pattern.match(via.v.strip())
        if m is None:
            raise RuntimeError(f"Error in Via.from_cad_file: via line {via.lineno} from doesn't match expected pattern!\n\tLine: {via.v.strip()}\n\tSource file: {via.sourcefile}")
        x = in2mm(float(m.groups()[3]))
        y = in2mm(float(m.groups()[4]))
        name = m.groups()[7]
        ret = Via(Point(x, y), name, [via])

        print(f"Parsed via on line {via.lineno}")

        return [ret], pre_via + post_via