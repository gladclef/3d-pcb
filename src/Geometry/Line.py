import numpy as np

import Geometry.geometry_tools as geo


class Line:
    def __init__(self, x: float, y: float, y_intercept: float=None, x_intercept: float=None):
        """
        Parameters
        ----------
        x : float
            The x component of the line slope.
        y : float
            The y component of the line slope.
        y_intercept : float, optional
            Where the line crosses the y axis, by default 0
        x_intercept : float, optional
            Where the line crosses the x axis, by default 0.
            Only necessary when the line is vertical.
        """
        # get default values
        slope = self.xy_to_slope(x, y)
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

        self.x = x
        self.y = y
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
                return cls(0, np.inf)
            else:
                return cls(0, -np.inf)
        elif abs(slope) <= geo.ZERO_THRESH:
            return cls(10, 0, y_intercept)
        else:
            return cls(10, slope*10, y_intercept)
    
    @classmethod
    def from_angle_point(cls, angle: float, point: tuple[float, float]) -> "Line":
        angle = geo.normalize_angle(angle)
        x, y = point

        # check for vertical or horizontal lines
        if abs(np.sin(angle)*geo.INF_THRESH) >= geo.INF_THRESH-1:
            return Line(np.cos(angle), np.sin(angle), x_intercept=x)
        elif abs(np.cos(angle)*geo.INF_THRESH) >= geo.INF_THRESH-1:
            return Line(np.cos(angle), np.sin(angle), y_intercept=y)
        
        slope = np.sin(angle) / np.cos(angle)
        y_intercept = y - slope*x
        return Line(np.cos(angle), np.sin(angle), y_intercept=y_intercept)
    
    @classmethod
    def from_two_points(cls, pnt1: tuple[float, float], pnt2: tuple[float, float]) -> "Line":
        x, y = pnt2[0] - pnt1[0], pnt2[1] - pnt1[1]
        angle = np.atan2(y, x)
        return cls.from_angle_point(angle, pnt1)

    @property
    def is_vertical(self) -> bool:
        return (abs(self.y) > geo.INF_THRESH) or \
               (abs(self.x) < geo.ZERO_THRESH) or \
               (abs(self.y / self.x) > geo.INF_THRESH)
    
    @property
    def is_horizontal(self) -> bool:
        return ( (abs(self.x) > geo.INF_THRESH) or \
                 (abs(self.y) < geo.ZERO_THRESH) or \
                 (abs(self.y / self.x) < geo.ZERO_THRESH) ) and \
               ( not self.is_vertical )

    @property
    def x1(self) -> float:
        """
        The first x value that can be used to define the line from two points.
        If this line is vertical, then this value will be the x intercept.
        This value is usually going to be 0.
        """
        if self.is_vertical:
            return self.x_intercept
        return 0

    @property
    def x2(self) -> float:
        """
        The second x value that can be used to define the line from two points.
        If this line is vertical, then this value will be the x intercept.
        This value is usually going to be x1+1.
        """
        if self.is_vertical:
            return self.x_intercept
        elif self.is_horizontal:
            return self.x1+1 if self.x > 0 else self.x1-1
        return self.x1+1

    @property
    def y1(self) -> float:
        """
        The first y value that can be used to define the line from two points.
        If this line is vertical, then this value will be 0.
        This value is usually going to be the y-intercept.
        """
        if self.is_vertical:
            return 0
        return self.y_intercept

    @property
    def y2(self) -> float:
        """
        The first y value that can be used to define the line from two points.
        If this line is vertical, then this value will be +/- 10.
        This value is usually going to be slope*x2.
        """
        if self.is_vertical:
            return 1 if self.y > 0 else -1
        return self.slope*self.x2 + self.y_intercept
    
    @property
    def xy1(self) -> tuple[float, float]:
        return self.x1, self.y1
    
    @property
    def xy2(self) -> tuple[float, float]:
        return self.x2, self.y2
    
    @property
    def angle(self) -> float:
        """ Angle of this line, in the range 0-2pi """
        ang = np.atan2(self.y2 - self.y1, self.x2 - self.x1)
        return geo.normalize_angle(ang)
    
    @property
    def slope(self) -> float:
        if self.is_vertical:
            return np.inf if self.y > 0 else -np.inf
        elif self.is_horizontal:
            return 0
        return self.y / self.x
    
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

    def is_point_on_right(self, test_point: tuple[float, float]) -> bool:
        if self.is_vertical:
            a, b = (self.x_intercept, 0), (self.x+self.x_intercept, self.y)
        elif self.is_horizontal or True:
            a, b = (0, self.y_intercept), (self.x, self.y+self.y_intercept)
        c = test_point

        is_on_left = (b[0] - a[0])*(c[1] - a[1]) - (b[1] - a[1])*(c[0] - a[0]) > 0
        return not is_on_left

    def is_parallel_to(self, other: "Line") -> bool:
        ang1, ang2 = self.angle, other.angle
        if abs(ang1 - ang2) < 1/1e9 or abs(np.pi - abs(ang1 - ang2)) < 1/1e9:
            return True
        return False

    def intersection(self, other: "Line") -> tuple[float, float] | None:
        """ Get the intersection point of two lines.
        
        Returns
        -------
        tuple[float, float]
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

        return x, y

    def distance_along_line(self, distance: float, from_point: tuple[float, float] = None) -> tuple[float, float]:
        """
        Get a point along the line that is a distance away from the given from_point.

        Parameters
        ----------
        distance : float
            The distance along the line.
        from_point : tuple[float, float] | None
            A point along the reference line to find the distant points at,
            or None to just use (0, y-intercept) as the from_point.
        """
        from_x, from_y = from_point if from_point is not None else (None, None)

        # Check for an infinite or 0 slope.
        if self.is_vertical:
            from_x = from_x if from_x is not None else self.x_intercept
            from_y = from_y if from_y is not None else 0
            if self.y > 0:
                return (from_x, from_y + distance)
            else:
                return (from_x, from_y - distance)
        if self.is_horizontal:
            from_x = from_x if from_x is not None else 0
            from_y = from_y if from_y is not None else self.y_intercept
            if self.x > 0:
                return (from_x + distance, from_y)
            else:
                return (from_x - distance, from_y)
        
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
        if self.x < 0:
            x = -x
        dx = from_x + x
        dy = self.slope*x + from_y
        
        return (dx, dy)

    def get_tangent_line(self, from_point: tuple[float, float]) -> "Line":
        """ Get the tangent line to the given reference_line.

        Parameters
        ----------
        from_point: tuple[float, float]
            The x,y point along the reference line that the
            tangent line should start at.
        """
        
        # Set some defaults.
        if from_point is None:
            from_point = (0, self.y_intercept)
        from_x, from_y = from_point

        # Check for vertical or horizontal lines.
        if self.is_vertical:
            if self.angle < np.pi:
                return Line.from_angle_point(0, from_point)
            else:
                return Line.from_angle_point(np.pi, from_point)
        if self.is_horizontal:
            if self.angle < np.pi*1/2 or self.angle > np.pi*3/2:
                return Line.from_angle_point(np.pi*3/2, from_point)
            else:
                return Line.from_angle_point(np.pi*1/2, from_point)
            
        # Get the tangent slope.
        tan_slope = -1 / self.slope

        # Get the tangent y-intercept
        tan_y_intercept = -tan_slope*from_x + from_y

        return self.__class__.from_slope_intercept(tan_slope, tan_y_intercept)
    
    def reversed(self) -> "Line":
        pnt1 = self.x1, self.y1
        pnt2 = self.x2, self.y2
        return self.__class__.from_two_points(pnt2, pnt1)
    
    def __repr__(self) -> str:
        xi = "N/A" if self.x_intercept is None else f"{self.x_intercept:.3f}"
        yi = "N/A" if self.y_intercept is None else f"{self.y_intercept:.3f}"
        return f"Line<x:{self.x:.3f},y:{self.y:.3f},xi:{xi},yi:{yi}>"