from abc import ABC, abstractmethod

class AbstractVtkPointTracker(ABC):
    @abstractmethod
    def inc_vtk_indicies(self, start: int, cnt: int):
        raise NotImplementedError()