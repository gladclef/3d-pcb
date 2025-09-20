import copy
import re

import matplotlib.axis as maxis
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation
import vtk

from Component.DrillHole import DrillHole
from tool.units import *
from FileIO.Line import Line as FLine
import Geometry.geometry_tools as geo
from Trace.VtkPointGroup import VtkPointGroup
from tool.globals import board_parameters as g
import tool.vtk_tools as vt

class Pin(DrillHole):
    """Represents a single through hole for a component shape."""

    def __init__(self, parent: "Shape", pad_name: str, x_offset: float, y_offset: float, layer: str, is_pad=False):
        """
        Initialize the Pin instance.

        Parameters
        ----------
        parent : Shape
            The shape that contains this pin.
        pad_name : str
            The name of the pad stack for this pin.
        x_offset : float
            The x offset of the pin from its origin.
        y_offset : float
            The y offset of the pin from its origin.
        layer : str
            The layer on which the pin is placed. Example "TOP" or "BOTTOM"
        is_pad : bool
            True if this instance is a pad, or False if it's a through-hole.

        """
        super().__init__(x_offset, y_offset, is_through_hole=True)

        self.parent = parent
        """The shape that contains this pin."""
        self.pad_name = pad_name
        """The name of the pin."""
        self.layer = layer
        """The layer on which the pin is placed."""
        self.is_pad = is_pad
        """True if this instance is a pad, or False if it's a through-hole."""

    def apply_translation_rotation_layer(self, translation: tuple[float, float], rotation: float, is_bottom: bool) -> "Pin":
        """
        Apply translation and rotation to the pin, with optional flipping.

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
        Pin
            A new Pin instance with the applied transformations.

        """
        ret = copy.deepcopy(self)

        super(self.__class__, ret).apply_translation_rotation_layer(translation, rotation, is_bottom)

        return ret

    @classmethod
    def from_cad_file(cls, lines: list[FLine]) -> tuple[list["Pin"], list[FLine]]:
        """
        Create a Pin instance by parsing lines from a CAD file.

        Parameters
        ----------
        lines : list[FLine]
            Lines of text representing the pin in a CAD file.

        Returns
        -------
        tuple[Pin, list[FLine]]
            A tuple containing the created Pin and any remaining unprocessed lines.

        Raises
        ------
        RuntimeError
            If the pin line does not match the expected pattern.

        """
        # example pin line:
        # PIN "2" PAD1 0.5 0 TOP 0 0

        pin_line = None
        for line_idx, line in enumerate(lines):
            if line.v.startswith("PIN"):
                pin_line = line
                ret_lines = lines[:line_idx] + lines[line_idx+1:]
                break
        if pin_line is None:
            return [], lines
        
        # regex explanation:                   pad num  x offset   y offset   layer   ???        ???
        pin_pattern = re.compile(r"PIN \"\d+\" (PAD\d+) ([-\d\.]+) ([-\d\.]+) ([^ ]+) ([-\d\.]+) ([-\d\.]+)")

        # break out the various parts of the pin line
        # print(pin_line.v.strip())
        match = pin_pattern.match(pin_line.v.strip())
        if match is None:
            raise RuntimeError("Error in Pin.from_cad_file(): failed to match pin_pattern to line:\n\t" + pin_line.v)

        pad_name, x_offset, y_offset, layer, _, _ = match.groups()
        x_offset = in2mm(float(x_offset))
        y_offset = in2mm(float(y_offset))

        parent = None
        pin = cls(parent, pad_name, x_offset, y_offset, layer)
        return [pin], ret_lines