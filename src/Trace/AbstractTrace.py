from abc import ABC, abstractmethod

import vtk

from Geometry.Point import Point
from Trace.AbstractVtkPointTracker import AbstractVtkPointTracker as PntInc
from FileIO.Line import Line as FLine
from Geometry.Path import Path
from Trace.PipeShape import PipeShape, DEFAULT_PIPE_SHAPE
from Trace.TraceLine import TraceLine
from tool.globals import board_parameters as g
from tool.units import *

class AbstractTrace(Path, PntInc, ABC):
    """
    Represents a wire trace on a PCB.
    """
    def __init__(self,
                 source_lines: list[FLine],
                 segments: list[TraceLine],
                 shape: PipeShape=None):
        Path.__init__(self, source_lines, segments)

        # set some defaults
        if shape is None:
            shape = DEFAULT_PIPE_SHAPE(g.DEFAULT_WIRE_DIAMETER)

        self.shape = shape

    @abstractmethod
    def to_vtk(self) -> vtk.vtkPolyData:
        raise NotImplementedError