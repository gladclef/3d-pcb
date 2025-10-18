from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import vtk

from Geometry.Point import Point
from Trace.AbstractVtkPointTracker import AbstractVtkPointTracker as PntInc
from FileIO.Line import Line as FLine
from Geometry.LineSegment import LineSegment
from Geometry.Path import Path
from Trace.PipeShape import PipeShape, DEFAULT_PIPE_SHAPE
from tool.globals import board_parameters as g
from tool.units import *

if TYPE_CHECKING:
    from Trace.SingleTrace import _TraceLine

class AbstractTrace(Path, PntInc, ABC):
    """
    Represents a wire trace on a PCB.
    """
    def __init__(self,
                 source_lines: list[FLine],
                 xy_points: list[Point],
                 segments: list["_TraceLine"] | list[LineSegment],
                 shape: PipeShape=None):
        Path.__init__(self, source_lines, xy_points, segments)

        # set some defaults
        if shape is None:
            shape = DEFAULT_PIPE_SHAPE(g.DEFAULT_WIRE_DIAMETER)

        self.shape = shape

    @abstractmethod
    def to_vtk(self) -> vtk.vtkPolyData:
        raise NotImplementedError