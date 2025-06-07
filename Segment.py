import math

import numpy as np
import geometry as geo
from AbstractTrace import AbstractTrace

class Segment:
    def __init__(self, parent: AbstractTrace, xy_point_indicies: tuple[int, int]):
        self.parent = parent
        self.xy_point_indicies = tuple(xy_point_indicies)

    @property
    def xy0(self) -> tuple[float, float]:
        return self.parent.xy_points[self.xy_point_indicies[0]]
    
    @property
    def xy1(self) -> tuple[float, float]:
        return self.parent.xy_points[self.xy_point_indicies[1]]
    
    @property
    def x0(self) -> float:
        return self.xy0[0]
    
    @property
    def y0(self) -> float:
        return self.xy0[1]
    
    @property
    def x1(self) -> float:
        return self.xy1[0]
    
    @property
    def y1(self) -> float:
        return self.xy1[1]
    
    @property
    def rise(self) -> float:
        return self.y1 - self.y0
    
    @property
    def run(self) -> float:
        return self.x1 - self.x0

    @property
    def angle(self) -> float:
        return geo.line_segment_angle(self.xy0, self.xy1)
    
    @property
    def length(self) -> float:
        return math.sqrt(self.rise**2 + self.run**2)
    
    def __repr__(self):
        return f"Seg<{self.xy0}:{self.xy1}>"