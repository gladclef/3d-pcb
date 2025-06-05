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
    
    @abstractmethod
    def get_setment_angle(self, segment: tuple[int, int]) -> float:
        """ Get the angle that the given segment runs at, from 0 to 2pi """
        raise NotImplementedError