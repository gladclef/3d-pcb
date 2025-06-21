import copy
import re

import matplotlib.axis as maxis
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation
import vtk

from tool.units import *
import Geometry.geometry_tools as geo
from Trace.VtkPointGroup import VtkPointGroup
import tool.vtk_tools as vt

class Pin:
    """Represents a single pad or via for a component shape."""

    def __init__(self, pad_name: str, x_offset: float, y_offset: float, layer: str):
        """
        Initialize the Pin instance.

        Parameters
        ----------
        pad_name : str
            The name of the pad stack for this pin.
        x_offset : float
            The x offset of the pin from its origin.
        y_offset : float
            The y offset of the pin from its origin.
        layer : str
            The layer on which the pin is placed. Example "TOP" or "BOTTOM"

        """
        self.pad_name = pad_name
        """The name of the pin."""
        self.x_offset = x_offset
        """The x offset of the pin from its origin."""
        self.y_offset = y_offset
        """The y offset of the pin from its origin."""
        self.layer = layer
        """The layer on which the pin is placed."""

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
        x, y = self.x_offset, self.y_offset
        x, y = geo.apply_translation_rotation_flip((x, y), translation, rotation, is_bottom)

        ret = copy.deepcopy(self)
        ret.x_offset, ret.y_offset = x, y

        return ret

    @classmethod
    def from_cad_file(cls, lines: list[str]) -> tuple["Pin", list[str]]:
        """
        Create a Pin instance by parsing lines from a CAD file.

        Parameters
        ----------
        lines : list[str]
            Lines of text representing the pin in a CAD file.

        Returns
        -------
        tuple[Pin, list[str]]
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
            if line.startswith("PIN"):
                pin_line = line
                ret_lines = lines[:line_idx] + lines[line_idx+1:]
                break
        if pin_line is None:
            return None, lines
        
        # regex explanation:                   pad num  x offset   y offset   layer   ???        ???
        pin_pattern = re.compile(r"PIN \"\d+\" (PAD\d+) ([-\d\.]+) ([-\d\.]+) ([^ ]+) ([-\d\.]+) ([-\d\.]+)")

        # break out the various parts of the pin line
        match = pin_pattern.match(pin_line)
        if match is None:
            raise RuntimeError("Error in Pin.from_cad_file(): failed to match pin_pattern to line:\n\t" + pin_line)

        pad_name, x_offset, y_offset, layer, _, _ = match.groups()
        x_offset = in2mm(float(x_offset))
        y_offset = in2mm(float(y_offset))

        pin = cls(pad_name, x_offset, y_offset, layer)
        return pin, ret_lines

    def to_vtk(self, polydata: vtk.vtkPolyData) -> vtk.vtkPolyData:
        vtk_points: vtk.vtkPoints = polydata.GetPoints()
        vtk_cells: vtk.vtkCellData = polydata.GetCellData()
        radius = 0.4

        # build the cylinder
        cylinder = vtk.vtkCylinderSource()
        cylinder.SetRadius(radius)
        cylinder.SetHeight(2)
        cylinder.SetResolution(32)
        cylinder.CappingOn()
        cylinder.Update()

        # join the cylinder with the input polydata
        cylinder_polydata = cylinder.GetOutput()
        vt.rotate(cylinder_polydata, Rotation.from_euler('x', np.pi/2))
        vt.translate(cylinder_polydata, (self.x_offset, self.y_offset, -1))
        vt.join(polydata, cylinder_polydata)
        
        return polydata

    def draw(self, ax: maxis.Axis):
        center = (self.x_offset, self.y_offset)
        ax.add_patch(plt.Circle(center, .3, color="tab:orange"))