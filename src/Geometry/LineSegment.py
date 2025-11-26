import math
from typing import Union
from Geometry.Line import Line

from Geometry.Point import Point
import Geometry.geometry_tools as geo
from FileIO.Line import Line as FLine


class LineSegment(Line):
    def __init__(self, xy0: Point, xy1: Point, fline: FLine = None):
        super().__init__(xy0, xy1)

        self.fline = fline

    @classmethod
    def from_vector(cls, x: float, y: float, y_intercept: float=None, x_intercept: float=None, length: float=1, fline: FLine = None) -> "LineSegment":
        l0 = Line(x, y, y_intercept, x_intercept)
        xy0 = l0.xy0
        xy1 = l0.distance_along_line(length, xy0)
        cls(xy0, xy1, fline)

    @classmethod
    def from_slope_intercept(cls, slope: float, y_intercept: float, length: float=1, fline: FLine = None) -> "LineSegment":
        l0 = Line.from_slope_intercept(slope, y_intercept)
        xy0 = l0.xy0
        xy1 = l0.distance_along_line(length, xy0)
        return cls(xy0, xy1, fline)
    
    @classmethod
    def from_angle_point(cls, angle: float, point: Point, length: float=1, fline: FLine = None) -> "LineSegment":
        l0 = Line.from_angle_point(angle, point)
        xy0 = l0.xy0
        xy1 = l0.distance_along_line(length, xy0)
        return cls(xy0, xy1, fline)
    
    @classmethod
    def from_two_points(cls, xy1: Point, pnt2: Point, fline: FLine = None) -> "LineSegment":
        cls(xy1, pnt2, fline)

    @property
    def rise(self) -> float:
        return self.y1 - self.y0
    
    @property
    def run(self) -> float:
        return self.x1 - self.x0

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
        x0 = min(self.x0, self.x1)
        y0 = min(self.y0, self.y1)
        x1 = max(self.x0, self.x1)
        y1 = max(self.y0, self.y1)
        for v0, v1 in [(x0, xing.x), (xing.x, x1), (y0, xing.y), (xing.y, y1)]:
            if v0 > v1 and v0 - v1 > geo.ZERO_THRESH:
                return None

        if isinstance(other, LineSegment):
            x2 = min(other.x0, other.x1)
            y2 = min(other.y0, other.y1)
            x3 = max(other.x0, other.x1)
            y3 = max(other.y0, other.y1)
            for v0, v1 in [(x2, xing.x), (xing.x, x3), (y2, xing.y), (xing.y, y3)]:
                if v0 > v1 and v0 - v1 > geo.ZERO_THRESH:
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
        xy = super().distance_along_line(distance, from_point)

        if limit_range:
            x0 = min(self.x0, self.x1)
            y0 = min(self.y0, self.y1)
            x1 = max(self.x0, self.x1)
            y1 = max(self.y0, self.y1)

            x = min(max(xy.x, x0), x1)
            y = min(max(xy.y, y0), y1)
            xy = Point(x, y)
        
        return xy

    def reversed(self) -> "LineSegment":
        # override parent method to return a line segment
        return self.__class__(self.xy1, self.xy0, self.fline)

    def __repr__(self):
        return f"LineSeg<{self.xy0}:{self.xy1}>"