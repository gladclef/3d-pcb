from abc import ABC, abstractmethod

import vtk

from Trace.AbstractVtkPointTracker import AbstractVtkPointTracker as PntInc
from Geometry.LineSegment import LineSegment
from Geometry.Path import Path
from Trace.PipeShape import PipeShape, PipeBasicCircle
from tool.units import *

class AbstractTrace(Path, PntInc, ABC):
    """
    Represents a wire trace on a PCB.
    """
    def __init__(self, xy_points: list[tuple[float, float]], segments: list[tuple[int, int]] | list[LineSegment], shape: PipeShape=None):
        Path.__init__(self, xy_points, segments)

        # set some defaults
        if shape is None:
            shape = PipeBasicCircle(awg2mm(26))

        self.shape = shape

    @abstractmethod
    def to_vtk(self) -> vtk.vtkPolyData:
        raise NotImplementedError