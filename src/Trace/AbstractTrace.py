from abc import ABC, abstractmethod

import vtk

from Trace.AbstractVtkPointTracker import AbstractVtkPointTracker as PntInc
from Geometry.Path import Path
from Trace.PipeShape import PipeShape
from Geometry.Segment import Segment

class AbstractTrace(Path, PntInc, ABC):
    """
    Represents a wire trace on a PCB.
    """
    def __init__(self, xy_points: list[tuple[float, float]], segments: list[tuple[int, int]] | list[Segment], shape: PipeShape):
        Path.__init__(self, xy_points, segments)
        self.shape = shape

    @abstractmethod
    def to_vtk(self) -> vtk.vtkPolyData:
        raise NotImplementedError