import copy
import re

import matplotlib.axis as maxis
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation
import vtk

from Component.DrillHole import DrillHole
from tool.units import *
import Geometry.geometry_tools as geo
from Trace.VtkPointGroup import VtkPointGroup
from tool.globals import board_parameters as g
import tool.vtk_tools as vt

class Via(DrillHole):
    """Represents a single via for a trace."""

    def __init__(self, x_offset: float, y_offset: float):
        """
        Initialize the Via instance.

        Parameters
        ----------
        x_offset : float
            The x offset of the via from the board origin.
        y_offset : float
            The y offset of the via from the board origin.

        """
        super().__init__(x_offset, y_offset, is_via=True)

    def apply_translation_rotation_layer(self, translation: tuple[float, float], rotation: float, is_bottom: bool) -> "Via":
        """
        Apply translation and rotation to the via, with optional flipping.

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
        Via
            A new Via instance with the applied transformations.

        """
        ret = copy.deepcopy(self)

        super(self.__class__, ret).apply_translation_rotation_layer(translation, rotation, is_bottom)

        return ret

    @classmethod
    def from_cad_file(cls, lines: list[str]) -> tuple[list["Via"], list[str]]:
        """
        Create a Via instance by parsing lines from a CAD file.

        Parameters
        ----------
        lines : list[str]
            Lines of text representing the via in a CAD file.

        Returns
        -------
        tuple[Via, list[str]]
            A tuple containing the created Via and any remaining unprocessed lines.

        Raises
        ------
        RuntimeError
            If the via line does not match the expected pattern.

        """
        # TODO
        print("TODO: import vias from gencad files")
        return [], lines