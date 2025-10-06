import math
from typing import Union
from Geometry.Line import Line

from Geometry.Point import Point
import Geometry.geometry_tools as geo
from FileIO.Line import Line as FLine


class LineSegment(Line):
    def __init__(self, pnt1: Point, pnt2: Point, source_line: FLine = None):
        l1 = Line.from_two_points(pnt1, pnt2)
        super().__init__(l1.xy, l1.y_intercept, l1.x_intercept)

        self.pnt1 = pnt1
        self.pnt2 = pnt2
        self.source_line = source_line

    @classmethod
    def from_vector(cls, x: float, y: float, y_intercept: float=None, x_intercept: float=None, length: float=1) -> "LineSegment":
        l1 = Line(x, y, y_intercept, x_intercept)
        pnt1 = l1.x1, l1.y1
        pnt2 = l1.distance_along_line(length, pnt1)
        cls(pnt1, pnt2)

    @classmethod
    def from_slope_intercept(cls, slope: float, y_intercept: float, length: float=1) -> "LineSegment":
        l1 = Line.from_slope_intercept(slope, y_intercept)
        pnt1 = l1.x1, l1.x2
        pnt2 = l1.distance_along_line(length, pnt1)
        return cls(pnt1, pnt2)
    
    @classmethod
    def from_angle_point(cls, angle: float, point: Point, length: float=1) -> "LineSegment":
        l1 = Line.from_angle_point(angle, point)
        pnt1 = l1.xy1
        pnt2 = l1.distance_along_line(length, pnt1)
        return cls(pnt1, pnt2)
    
    @classmethod
    def from_two_points(cls, pnt1: Point, pnt2: Point) -> "LineSegment":
        cls(pnt1, pnt2)

    @property
    def x1(self) -> float:
        """ The x component of the start of this segment. """
        return self.pnt1.x
    
    @property
    def x2(self) -> float:
        """ The x component of the end of this segment. """
        return self.pnt2.x
    
    @property
    def y1(self) -> float:
        """ The y component of the start of this segment. """
        return self.pnt1.y
    
    @property
    def y2(self) -> float:
        """ The y component of the end of this segment. """
        return self.pnt2.y

    @property
    def rise(self) -> float:
        return self.y2 - self.y1
    
    @property
    def run(self) -> float:
        return self.x2 - self.x1

    @property
    def length(self) -> float:
        return math.sqrt(self.rise**2 + self.run**2)
    
    def intersection(self, other: Union["Line", "LineSegment"]) -> Point | None:
        """ Get the intersection point of two line (segments).

        If the line (segments) don't have an intersection point, then return None.
        """
        xing = super().intersection(other)
        if xing is None:
            return xing

        # Check that the intersection point is within the bounding boxes of the two line segments
        x0 = min(self.x1, self.x2)
        y0 = min(self.y1, self.y2)
        x1 = max(self.x1, self.x2)
        y1 = max(self.y1, self.y2)
        for v1, v2 in [(x0, xing.x), (xing.x, x1), (y0, xing.y), (xing.y, y1)]:
            if v1 > v2 and v1 - v2 > geo.ZERO_THRESH:
                return None

        if isinstance(other, LineSegment):
            x2 = min(other.x1, other.x2)
            y2 = min(other.y1, other.y2)
            x3 = max(other.x1, other.x2)
            y3 = max(other.y1, other.y2)
            for v1, v2 in [(x2, xing.x), (xing.x, x3), (y2, xing.y), (xing.y, y3)]:
                if v1 > v2 and v1 - v2 > geo.ZERO_THRESH:
                    return None
        
        return xing

    def distance_along_line(self, distance: float, from_point: Point = None, limit_range=False) -> Point:
        """
        Get a point along the line that is a distance away from the given from_point.

        Parameters
        ----------
        distance : float
            The distance along the line.
        from_point : Point | None
            A point along the reference line to find the distant points at,
            or None to just use (0, y-intercept) as the from_point.
        limit_range : bool
            True if the returned value should be limited such that it is on the segment,
            False to allow the returned point to be beyond the segment's end.
        """
        x, y = super().distance_along_line(distance, from_point)

        if limit_range:
            x0 = min(self.x1, self.x2)
            y0 = min(self.y1, self.y2)
            x1 = max(self.x1, self.x2)
            y1 = max(self.y1, self.y2)

            x = min(max(x, x0), x1)
            y = min(max(y, y0), y1)
        
        return x, y

    def reversed(self) -> "LineSegment":
        # override parent method to return a line segment
        return LineSegment(self.xy2, self.xy1)

    def __repr__(self):
        return f"LineSeg<{self.xy1}:{self.xy2}>"