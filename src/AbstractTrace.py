from abc import ABC, abstractmethod

import vtk

from AbstractVtkPointTracker import AbstractVtkPointTracker as PntInc
from PipeShape import PipeShape

class AbstractTrace(PntInc, ABC):
    """
    Represents a wire trace on a PCB.
    """
    def __init__(self, xy_points: list[tuple[float, float]], shape: PipeShape):
        self.xy_points = xy_points
        self.shape = shape

    @abstractmethod
    def to_vtk(self) -> vtk.vtkPolyData:
        raise NotImplementedError