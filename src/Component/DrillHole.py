from abc import ABC
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
from tool.globals import board_parameters as g
import tool.vtk_tools as vt

class DrillHole(ABC):
    """ Represents a drill hole in a circuit board. """

    def __init__(self, x_offset: float, y_offset: float, is_via=False, is_through_hole=False, diameter: float=None):
        """
        Initialize the DrillHole instance.

        Parameters
        ----------
        x_offset : float
            The x offset of the instance from the parent's origin.
        y_offset : float
            The y offset of the instance from the parent's origin.
        diameter : float
            The diameter used for manually specified drill holes.

        """
        self.x_offset = x_offset
        """ The x offset of the instance from the parent's origin. """
        self.y_offset = y_offset
        """ The y offset of the instance from the parent's origin. """
        self.is_via = is_via
        self.is_through_hole = is_through_hole

        self._diameter = diameter
        """ The diameter used for manually specified drill holes. """

    @property
    def diameter(self) -> float:
        if self.is_via:
            return self.via_diameter()
        elif self.is_through_hole:
            return self.through_hole_diameter()
        else:
            return self._diameter
    
    @staticmethod
    def via_diameter() -> float:
        return g.VIA_DIAMETER + g.VIA_DIAMETER_CLEARANCE
    
    @staticmethod
    def through_hole_diameter() -> float:
        return g.THROUGH_HOLE_DIAMETER + g.THROUGH_HOLE_DIAMETER_CLEARANCE

    def apply_translation_rotation_layer(self, translation: tuple[float, float], rotation: float, is_bottom: bool):
        """
        Apply translation and rotation directly to this drill hole, with optional flipping.

        Parameters
        ----------
        translation : tuple[float, float]
            Translation vector (x, y).
        rotation : float
            Rotation angle in radians.
        is_bottom : bool
            If True, flip the coordinates along the x-axis after applying translation and rotation.
        """
        x, y = self.x_offset, self.y_offset
        x, y = geo.apply_translation_rotation_flip((x, y), translation, rotation, is_bottom)
        self.x_offset, self.y_offset = x, y

    def to_vtk(self, polydata: vtk.vtkPolyData) -> vtk.vtkPolyData:
        radius = self.diameter / 2
        height = g.BOARD_THICKNESS + 1.0

        # build the cylinder
        cylinder = vtk.vtkCylinderSource()
        cylinder.SetRadius(radius)
        cylinder.SetHeight(height)
        cylinder.SetResolution(32)
        cylinder.CappingOn()
        cylinder.Update()

        # join the cylinder with the input polydata
        cylinder_polydata = cylinder.GetOutput()
        vt.rotate(cylinder_polydata, Rotation.from_euler('x', np.pi/2))
        vt.translate(cylinder_polydata, (self.x_offset, self.y_offset, -height/2))
        vt.join(polydata, cylinder_polydata)
        
        return polydata

    def draw(self, ax: maxis.Axis):
        center = (self.x_offset, self.y_offset)
        ax.add_patch(plt.Circle(center, .3, color="tab:orange"))