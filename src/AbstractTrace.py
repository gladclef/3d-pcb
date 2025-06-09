from abc import ABC, abstractmethod

import vtk

from AbstractVtkPointTracker import AbstractVtkPointTracker as PntInc
from AbstractPath import AbstractPath
from PipeShape import PipeShape
from Segment import Segment

class AbstractTrace(AbstractPath, PntInc, ABC):
    """
    Represents a wire trace on a PCB.
    """
    def __init__(self, xy_points: list[tuple[float, float]], segments: list[tuple[int, int]] | list[Segment], shape: PipeShape):
        AbstractPath.__init__(self, xy_points, segments)
        self.shape = shape

    @abstractmethod
    def to_vtk(self) -> vtk.vtkPolyData:
        raise NotImplementedError