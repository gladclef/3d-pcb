import numpy as np

from Geometry.Point import Point
import Geometry.geometry_tools as geo


class Line:
    def __init__(self, xy: Point, y_intercept: float=None, x_intercept: float=None):
        """
        Parameters
        ----------
        xy : Point
            The x and y components of the line slope.
        y_intercept : float, optional
            Where the line crosses the y axis, by default 0
        x_intercept : float, optional
            Where the line crosses the x axis, by default 0.
            Only necessary when the line is vertical.
        """
        # get default values
        slope = self.xy_to_slope(xy.x, xy.y)
        if slope == 0:
            # x intercept is undefined
            y_intercept = y_intercept if y_intercept is not None else 0
            if x_intercept is not None:
                raise ValueError("2D line x intercept is undefined for a horizontal line")
        elif abs(slope) == np.inf:
            # y intercept is undefined
            x_intercept = x_intercept if x_intercept is not None else 0
            if y_intercept is not None:
                raise ValueError("2D line y intercept is undefined for a horizontal line")
        else:
            # assume/calculate x and y intercepts
            if x_intercept is not None:
                calc_y_intercept = -slope*x_intercept
                if y_intercept is not None and abs(y_intercept - calc_y_intercept) > geo.ZERO_THRESH:
                    raise ValueError(f"2D line has both an x and y intercept, but provided y-intercept {y_intercept} != calculated y-intercept {calc_y_intercept}")
                y_intercept = calc_y_intercept
            else:
                if y_intercept is not None:
                    x_intercept = -y_intercept/slope
                else:
                    x_intercept = 0
                    y_intercept = 0

        self.xy = xy
        self.y_intercept = y_intercept
        self.x_intercept = x_intercept

    @staticmethod
    def xy_to_slope(x: float, y: float) -> float:
        if abs(x) <= geo.ZERO_THRESH:
            if abs(y) <= geo.ZERO_THRESH:
                raise RuntimeError("2D line with zero x slope component and zero y slope component is undefined.")
            else:
                y_dir = 1 if y > 0 else -1
                return y_dir * np.inf
        else: # x != 0
            if abs(y) <= geo.ZERO_THRESH:
                return 0
            else: # x != 0 and y != 0
                if abs(x) >= geo.INF_THRESH:
                    if abs(y) >= geo.INF_THRESH:
                        raise RuntimeError("2D line with infinite x slope component and infinite y slope component is undefined.")
                    else: # y < inf
                        return 0
                else: # x < inf
                    if abs(y) >= geo.INF_THRESH:
                        y_dir = 1 if y > 0 else -1
                        return y_dir * np.inf
                    else:
                        return y / x
    
    @classmethod
    def from_slope_intercept(cls, slope: float, y_intercept: float) -> "Line":
        if abs(slope) >= geo.INF_THRESH:
            if slope > 0:
                return cls(Point(0, np.inf))
            else:
                return cls(Point(0, -np.inf))
        elif abs(slope) <= geo.ZERO_THRESH:
            return cls(Point(10, 0), y_intercept)
        else:
            return cls(Point(10, slope*10), y_intercept)
    
    @classmethod
    def from_angle_point(cls, angle: float, point: Point) -> "Line":
        angle = geo.normalize_angle(angle)

        # check for vertical or horizontal lines
        if abs(np.sin(angle)*geo.INF_THRESH) >= geo.INF_THRESH-1:
            return cls(Point(np.cos(angle), np.sin(angle)), x_intercept=point.x)
        elif abs(np.cos(angle)*geo.INF_THRESH) >= geo.INF_THRESH-1:
            return cls(Point(np.cos(angle), np.sin(angle)), y_intercept=point.y)
        
        slope = np.sin(angle) / np.cos(angle)
        y_intercept = point.y - slope*point.x
        x = np.cos(angle) * np.sign(angle)
        y = np.sin(angle) * np.sign(angle)
        return cls(Point(x, y), y_intercept=y_intercept)
    
    @classmethod
    def from_two_points(cls, pnt0: Point, pnt1: Point) -> "Line":
        diff = pnt1 - pnt0
        angle = np.atan2(diff.y, diff.x)
        return cls.from_angle_point(angle, pnt0)

    @property
    def is_vertical(self) -> bool:
        return (abs(self.xy.y) > geo.INF_THRESH) or \
               (abs(self.xy.x) < geo.ZERO_THRESH) or \
               (abs(self.xy.y / self.xy.x) > geo.INF_THRESH)
    
    @property
    def is_horizontal(self) -> bool:
        return ( (abs(self.xy.x) > geo.INF_THRESH) or \
                 (abs(self.xy.y) < geo.ZERO_THRESH) or \
                 (abs(self.xy.y / self.xy.x) < geo.ZERO_THRESH) ) and \
               ( not self.is_vertical )

    @property
    def x0(self) -> float:
        """
        The first x value that can be used to define the line from two points.
        If this line is vertical, then this value will be the x intercept.
        This value is usually going to be 0.
        """
        if self.is_vertical:
            return self.x_intercept
        return 0

    @property
    def x1(self) -> float:
        """
        The second x value that can be used to define the line from two points.
        If this line is vertical, then this value will be the x intercept.
        This value is usually going to be x1+1.
        """
        if self.is_vertical:
            return self.x_intercept
        elif self.is_horizontal:
            return self.x0+1 if self.xy.x > 0 else self.x0-1
        sign = 1 if (self.xy.x > 0) else -1
        return self.x0+sign

    @property
    def y0(self) -> float:
        """
        The first y value that can be used to define the line from two points.
        If this line is vertical, then this value will be 0.
        This value is usually going to be the y-intercept.
        """
        if self.is_vertical:
            return 0
        return self.y_intercept

    @property
    def y1(self) -> float:
        """
        The first y value that can be used to define the line from two points.
        If this line is vertical, then this value will be +/- 10.
        This value is usually going to be slope*x2.
        """
        if self.is_vertical:
            return 1 if self.xy.y > 0 else -1
        return self.slope*self.x1 + self.y_intercept

    @property
    def xy0(self) -> Point:
        return Point(self.x0, self.y0)
    
    @property
    def xy1(self) -> Point:
        return Point(self.x1, self.y1)
    
    @property
    def xy_points(self) -> tuple[Point, Point]:
        return self.xy0, self.xy1
    
    @property
    def angle(self) -> float:
        """ Angle of this line, in the range 0-2pi """
        ang = np.atan2(self.y1 - self.y0, self.x1 - self.x0)
        return geo.normalize_angle(ang)
    
    @property
    def slope(self) -> float:
        if self.is_vertical:
            return np.inf if self.xy.y > 0 else -np.inf
        elif self.is_horizontal:
            return 0
        return self.xy.y / self.xy.x
    
    @property
    def rise(self) -> float:
        return np.sin(self.angle)
    
    @property
    def run(self) -> float:
        return np.cos(self.angle)

    @staticmethod
    def angle_to_slope(angle: float) -> float:
        """ Returns the slope (rise / run) for the given angle. """
        return np.sin(angle) / np.cos(angle)

    def is_point_on_right(self, test_point: Point) -> bool:
        if self.is_vertical:
            a, b = Point(self.x_intercept, 0), Point(self.xy.x+self.x_intercept, self.xy.y)
        elif self.is_horizontal or True:
            a, b = Point(0, self.y_intercept), Point(self.xy.x, self.xy.y+self.y_intercept)
        c = test_point

        is_on_left = (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x) > 0
        return not is_on_left

    def is_parallel_to(self, other: "Line") -> bool:
        ang1, ang2 = self.angle, other.angle
        if abs(ang1 - ang2) < 1/1e9 or abs(np.pi - abs(ang1 - ang2)) < 1/1e9:
            return True
        return False

    def intersection(self, other: "Line") -> Point | None:
        """ Get the intersection point of two lines.
        
        Returns
        -------
        Point
            The x,y intersection point, or None if the lines are parallel
        """
        # Don't check for intersection between two lines that are parallel
        if self.is_parallel_to(other):
            return None

        elif self.is_vertical:
            x = self.x_intercept
            y = other.slope*x + other.y_intercept
        elif other.is_vertical:
            x = other.x_intercept
            y = self.slope*x + self.y_intercept
        elif self.is_horizontal:
            y = self.y_intercept
            x = (y-other.y_intercept)/other.slope
        elif other.is_horizontal:
            y = other.y_intercept
            x = (y-self.y_intercept)/self.slope
        else:
            # Get the intersection point for two infinite lines
            # assume the y values are equal at the intersection point
            #   y1 = slope_a*x + y_int_a
            #   y2 = slope_b*x + y_int_b
            #   slope_a*x + y_int_a = slope_b*x + y_int_b
            # solve for x
            #   slope_a*x - slope_b*x = y_int_b - y_int_a
            #   x = (y_int_b - y_int_a) / (slope_a - slope_b)
            # get the y value
            #   y = slope_a*x + y_int_a
            slope_a, y_int_a = self.slope, self.y_intercept
            slope_b, y_int_b = other.slope, other.y_intercept
            x = (y_int_b - y_int_a) / (slope_a - slope_b)
            y = slope_a*x + y_int_a

        return Point(x, y)

    def distance_along_line(self, distance: float, from_point: Point = None) -> Point:
        """
        Get a point along the line that is a distance away from the given from_point.

        Parameters
        ----------
        distance : float
            The distance along the line.
        from_point : Point | None
            A point along the reference line to find the distant points at,
            or None to just use (0, y-intercept) as the from_point.
        """
        from_x, from_y = (from_point.x, from_point.y) if (from_point is not None) else (None, None)

        # Check for an infinite or 0 slope.
        if self.is_vertical:
            from_x = from_x if from_x is not None else self.x_intercept
            from_y = from_y if from_y is not None else 0
            if self.xy.y > 0:
                return Point(from_x, from_y + distance)
            else:
                return Point(from_x, from_y - distance)
        if self.is_horizontal:
            from_x = from_x if from_x is not None else 0
            from_y = from_y if from_y is not None else self.y_intercept
            if self.xy.x > 0:
                return Point(from_x + distance, from_y)
            else:
                return Point(from_x - distance, from_y)
        
        # Use some defaults.
        from_x = from_x if from_x is not None else 0
        from_y = from_y if from_y is not None else self.y_intercept

        # Get the distant point.
        # Here 'x' and 'y' are the x and y values for the
        # point along the line relative to the from_point.
        #
        #     y = slope*x + 0                         equation 1
        #
        #     d = sqrt(x^2 + y^2)
        #     y^2 = d^2 - x^2
        #     y = sqrt(d^2 - x^2)                equation 2
        #
        #     sqrt(d^2 - x^2) = slope*x
        #     d^2 - x^2 = (slope^2)*(x^2)
        #     (1 + slope^2)*(x^2) = d^2
        #     x = sqrt(d^2 / (1 + slope^2))      equation 3
        #
        x = np.sqrt(distance**2 / (1 + self.slope**2))
        if self.xy.x < 0:
            x = -x
        dx = from_x + x
        dy = self.slope*x + from_y
        
        return Point(dx, dy)

    def get_tangent_line(self, from_point: Point) -> "Line":
        """ Get the tangent line to the given reference_line.

        Parameters
        ----------
        from_point: Point
            The x,y point along the reference line that the
            tangent line should start at.
        """
        
        # Set some defaults.
        if from_point is None:
            from_point = (0, self.y_intercept)

        # Check for vertical or horizontal lines.
        if self.is_vertical:
            if self.angle < np.pi:
                return self.__class__.from_angle_point(0, from_point)
            else:
                return self.__class__.from_angle_point(np.pi, from_point)
        if self.is_horizontal:
            if self.angle < np.pi*1/2 or self.angle > np.pi*3/2:
                return self.__class__.from_angle_point(np.pi*3/2, from_point)
            else:
                return self.__class__.from_angle_point(np.pi*1/2, from_point)
            
        # Get the tangent slope.
        tan_slope = -1 / self.slope

        # Get the tangent y-intercept
        tan_y_intercept = -tan_slope*from_point.x + from_point.y

        return self.__class__.from_slope_intercept(tan_slope, tan_y_intercept)
    
    def reversed(self) -> "Line":
        return self.__class__.from_two_points(self.xy1, self.xy0)
    
    def angle_between(self, other: "Line"):
        """ Returns the angle between this line and the other line, between [0-pi) radians. """
        rel_angle = abs(self.angle - other.angle)
        if rel_angle > np.pi:
            rel_angle = 2*np.pi - rel_angle
        return rel_angle

    def __repr__(self) -> str:
        xi = "N/A" if self.x_intercept is None else f"{self.x_intercept:.3f}"
        yi = "N/A" if self.y_intercept is None else f"{self.y_intercept:.3f}"
        return f"Line<x:{self.run:.3f},y:{self.rise:.3f},xi:{xi},yi:{yi}>"