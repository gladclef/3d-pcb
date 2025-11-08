from abc import ABC
import copy
import re

import matplotlib.axis as maxis
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation
import vtk

from Geometry.Point import Point
from tool.units import *
from FileIO.Line import Line as FLine
import Geometry.geometry_tools as geo
from Trace.VtkPointGroup import VtkPointGroup
from tool.globals import board_parameters as g
import tool.vtk_tools as vt

class DrillHole(ABC):
    """ Represents a drill hole in a circuit board. """

    def __init__(self, location: Point, is_via=False, is_through_hole=False, diameter: float=None, source_lines: list[FLine]=None):
        """
        Initialize the DrillHole instance.

        Parameters
        ----------
        location : Point
            The location of the instance relative to parent's origin.
        diameter : float
            The diameter used for manually specified drill holes.

        """
        self.location = location
        """ The location of the instance relative to parent's origin. """
        self.is_via = is_via
        self.is_through_hole = is_through_hole

        self._diameter = diameter
        """ The diameter used for manually specified drill holes. """
        self.source_lines = [] if source_lines is None else source_lines

    @property
    def diameter(self) -> float:
        if self.is_via:
            return self.via_diameter()
        elif self.is_through_hole:
            return self.through_hole_diameter()
        else:
            return self._diameter
    
    @staticmethod
    def through_hole_diameter() -> float:
        return g.THROUGH_HOLE_DIAMETER + g.THROUGH_HOLE_DIAMETER_CLEARANCE

    def apply_translation_rotation_layer(self, translation: tuple[float, float] | Point, rotation: float, is_bottom: bool):
        """
        Apply translation and rotation directly to this drill hole, with optional flipping.

        Parameters
        ----------
        translation : tuple[float, float] | Point
            Translation vector (x, y).
        rotation : float
            Rotation angle in radians.
        is_bottom : bool
            If True, flip the coordinates along the x-axis after applying translation and rotation.
        """
        xy = self.location
        xy = geo.apply_translation_rotation_flip(xy, translation, rotation, is_bottom)
        self.location = xy

    def to_vtk(self, polydata: vtk.vtkPolyData) -> vtk.vtkPolyData:
        radius = self.diameter / 2
        height = g.BOARD_THICKNESS + 1.0

        # build the cylinder
        cylinder = vtk.vtkCylinderSource()
        cylinder.SetRadius(radius)
        cylinder.SetHeight(height)
        cylinder.SetResolution(g.CIRCLE_RESOLUTION)
        cylinder.CappingOn()
        cylinder.Update()

        # join the cylinder with the input polydata
        cylinder_polydata = cylinder.GetOutput()
        vt.rotate(cylinder_polydata, Rotation.from_euler('x', np.pi/2))
        vt.translate(cylinder_polydata, (self.location.x, self.location.y, -height/2))
        vt.join(polydata, cylinder_polydata)
        
        return polydata
    
    def _draw(self, ax: maxis.Axis, diameter: float = None, color: str = "tab:olive"):
        diameter = diameter or self.through_hole_diameter()
        center = (self.location.x, self.location.y)
        ax.add_patch(plt.Circle(center, diameter, color=color))

    def draw(self, ax: maxis.Axis):
        self._draw(ax)