import math
from typing import overload

import numpy as np


def normalize_angle(angle: float) -> float:
    """ Returns an angle between 0-2pi """
    while angle < 0:
        angle += 2*np.pi
    while angle > 2*np.pi:
        angle -= 2*np.pi
    return angle

def angle_to_slope(angle: float) -> float:
    """ Returns the slope (rise / run) for the given angle. """
    return math.sin(angle) / math.cos(angle)

def is_point_on_right(line_points: tuple[tuple[float, float], tuple[float, float]], test_point: tuple[float, float]) -> bool:
    a, b, c = line_points[0], line_points[1], test_point
    is_on_left = (b[0] - a[0])*(c[1] - a[1]) - (b[1] - a[1])*(c[0] - a[0]) > 0
    return not is_on_left

def line_segment_angle(p0: tuple[float, float], p1: tuple[float, float]) -> float:
    """
    Parameters
    ----------
    p0 : tuple[float, float]
        The x,y values of one end of the line segment.
    p1 : tuple[float, float]
        The x,y values of the other end of the line segment.

    Returns
    -------
    float
        The angle from p0 to p1, in the range 0-2pi
    """
    ang = math.atan2(p1[1] - p0[1], p1[0] - p0[0])
    if ang < 0:
        ang = np.pi * 2 + ang
    return ang

def line_segment_slope(p0: tuple[float, float], p1: tuple[float, float]) -> float:
    """
    Returns
    -------
    float
        The slope of the line segment, as in 'a' in 'y = ax + b'
    """
    rise = p1[1] - p0[1]
    run = p1[0] - p0[0]
    return rise / run

def line_segment_yintercept(p0: tuple[float, float], p1: tuple[float, float]) -> float:
    """
    Returns
    -------
    float
        The y-intercept of the infinite line defined by the two end points.
    """
    slope = line_segment_slope(p0, p1)
    y_intercept = p0[1] - slope * p0[0]
    return y_intercept

def lines_intersection(line_a: tuple[tuple[float, float], tuple[float, float]],
                       line_b: tuple[tuple[float, float], tuple[float, float]]) -> tuple[float, float] | None:
    """ Get the intersection point of two lines.

    Parameters
    ----------
    p0 : tuple[float, float]
        The x,y values of one end of the first line segment.
    p1 : tuple[float, float]
        The x,y values of the other end of the first line segment.
    p2 : tuple[float, float]
        The x,y values of one end of the second line segment.
    p3 : tuple[float, float]
        The x,y values of the other end of the second line segment.
    Returns
    -------
    tuple[float, float]
        The intersection point, or None if the lines are parallel
    """
    # Don't check for intersection between two lines that are parallel
    ang1, ang2 = line_segment_angle(*line_a), line_segment_angle(*line_b)
    if ang1 - ang2 < 0.001 or 2*np.pi - (ang1 - ang2) < 0.001:
        return None

    # Get the intersection point for two infinite lines
    # assume the y values are equal at the intersection point
    #   y1 = a1*x + b1
    #   y2 = a2*x + b2
    #   a1*x + b1 = a2*x + b2
    # solve for x
    #   a1*x - a2*x = b2 - b1
    #   x = (b2 - b1) / (a1 - a2)
    # get the y value
    #   y = a1*x + b1
    a1, a2 = line_segment_slope(*line_a), line_segment_slope(*line_b)
    b1, b2 = line_segment_yintercept(*line_a), line_segment_yintercept(*line_b)
    x = (b2 - b1) / (a1 - a2)
    y = a1*x + b1

    return x, y

def line_segments_intersection(line_a: tuple[tuple[float, float], tuple[float, float]],
                               line_b: tuple[tuple[float, float], tuple[float, float]]) -> tuple[float, float] | None:
    """ Get the intersection point of two line segments, if it exists.
    
    See lines_intersection() for a more detailed description.

    Returns
    -------
    tuple[float, float] | None
        The intersection point, or None if the segments don't intersect
    """
    x_y = lines_intersection(line_a, line_b)
    if x_y is None:
        return None
    x, y = x_y

    # Check that the intersection point is within the bounding boxes of the two line segments
    x0 = min(line_a[0][0], line_a[1][0])
    y0 = min(line_a[0][1], line_a[1][1])
    x1 = max(line_a[0][0], line_a[1][0])
    y1 = max(line_a[0][1], line_a[1][1])

    x2 = min(line_b[0][0], line_b[1][0])
    y2 = min(line_b[0][1], line_b[1][1])
    x3 = max(line_b[0][0], line_b[1][0])
    y3 = max(line_b[0][1], line_b[1][1])

    if (x < x0 or x > x1 or y < y0 or y > y1):
        return None

    if (x < x2 or x > x3 or y < y2 or y > y3):
        return None

    return (x, y)

def line_slope_intercept_to_two_points(line: tuple[float, float]) -> tuple[tuple[float, float], tuple[float, float]]:
    slope, y_intercept = line

    if abs(slope) == np.inf:
        if slope > 0:
            return (0, 0), (0, 10)
        else:
            return (0, 0), (0, -10)
    elif slope == 0:
        return (0, y_intercept), (10, y_intercept)
    else:
        x2, y2 = 10, slope*10 + y_intercept
        return (0, y_intercept), (x2, y2)

def line_two_points_to_slope_intercept(line: tuple[tuple[float, float], tuple[float, float]]) -> tuple[float, float]:
    (x1, y1), (x2, y2) = line
    
    if x1 != x2:
        # the line is not perfectly vertical
        slope = (y2 - y1) / (x2 - x1)
        y_intercept = -slope*x1 + y1
    else:
        # the line is vertical
        slope = np.inf if y2 > y1 else -np.inf
        y_intercept = 0
    
    return slope, y_intercept

@overload
def distance_along_line(line: tuple[float, float], distance: float, from_point: tuple[float, float]) -> tuple[float, float]:
    """
    Get a point along the line that is a distance away from the given from_point.

    Parameters
    ----------
    reference_line : tuple[float, float]
        The slope and y-intercept of the line to get the point along.
    distance : float
        The distance along the line.
    from_point : tuple[float, float] | None
        A point along the reference line to find the distant points at,
        or None to just use (0, y-intercept) as the from_point.
    """
    pass
@overload
def distance_along_line(line: tuple[tuple[float, float], tuple[float, float]], distance: float, from_point: tuple[float, float]) -> tuple[float, float]:
    """
    Get a point along the line that is a distance away from the given from_point.

    Parameters
    ----------
    reference_line : tuple[tuple[float, float], tuple[float, float]]
        Two points that represent the line.
    distance : float
        The distance along the line.
    from_point : tuple[float, float] | None
        A point along the reference line to find the distant points at,
        or None to just use (0, y-intercept) as the from_point.
    """
    pass
def distance_along_line(line: tuple[float, float] | tuple[tuple, tuple], distance: float, from_point: tuple[float, float]=None) -> tuple[float, float]:
    # Standardize the line input values
    try:
        (x1, y1), (x2, y2) = line
        slope, y_intercept = line_two_points_to_slope_intercept(line)
    except:
        slope, y_intercept = line
        (x1, y1), (x2, y2) = line_slope_intercept_to_two_points(line)
    
    # Set some defaults.
    if from_point is None:
        from_point = (0, y_intercept)
    from_x, from_y = from_point

    # Check for an infinite or 0 slope.
    if abs(slope) == np.inf:
        if slope > 0:
            return (from_x, from_y + distance)
        else:
            return (from_x, from_y - distance)
    if slope == 0:
        if x2 > x1:
            return (from_x + distance, from_y)
        else:
            return (from_x - distance, from_y)

    # Get the distant point.
    # Here 'x' and 'y' are the x and y values for the
    # point along the line relative to the from_point.
    #
    #     y = slope*x + 0                         equation 1
    #
    #     d = math.sqrt(x^2 + y^2)
    #     y^2 = d^2 - x^2
    #     y = math.sqrt(d^2 - x^2)                equation 2
    #
    #     math.sqrt(d^2 - x^2) = slope*x
    #     d^2 - x^2 = (slope^2)*(x^2)
    #     (1 + slope^2)*(x^2) = d^2
    #     x = math.sqrt(d^2 / (1 + slope^2))      equation 3
    #
    x = math.sqrt(distance**2 / (1 + slope**2))
    if x2 < x1:
        x = -x
    dx = from_x + x
    dy = slope*x + from_y
    
    return (dx, dy)

@overload
def get_tangent_line(reference_line: tuple[float, float], from_point: tuple[float, float]) -> tuple[float, float]:
    """ Get the tangent line to the given reference_line.

    Parameters
    ----------
    reference_line: tuple[float, float]
        The reference line as (slope, y-intercept).
    from_point: tuple[float, float]
        The x,y point along the reference line that the
        tangent line should start at.
    """
    pass
@overload
def get_tangent_line(reference_line: tuple[tuple[float, float], tuple[float, float]], from_point: tuple[float, float]) -> tuple[tuple[float, float], tuple[float, float]]:
    """ Get the tangent line to the given reference_line.

    Parameters
    ----------
    reference_line: tuple[tuple[float, float], tuple[float, float]]
        The reference line as defined by two points along the line.
    from_point: tuple[float, float]
        The x,y point along the reference line that the
        tangent line should start at.
    """
    pass
def get_tangent_line(reference_line: tuple[float, float] | tuple[tuple, tuple], from_point: tuple[float, float]=None) -> tuple[float, float] | tuple[tuple, tuple]:
    # Standardize the line input values
    try:
        (x1, y1), (x2, y2) = reference_line
        slope, y_intercept = line_two_points_to_slope_intercept(reference_line)
        ret_type = "two points"
    except:
        slope, y_intercept = reference_line
        (x1, y1), (x2, y2) = line_slope_intercept_to_two_points(reference_line)
        ret_type = "slope intercept"
    
    # Set some defaults.
    if from_point is None:
        from_point = (0, y_intercept)
    from_x, from_y = from_point

    # Check for an infinite or 0 slope.
    if ret_type == "two points":
        if abs(slope) == np.inf:
            xdiff = 10 if y2 > y1 else -10
            return (from_x, from_y), (from_x + xdiff, from_y)
        if slope == 0:
            ydiff = -10 if x2 > x1 else 10
            return (from_x, from_y), (from_x, from_y + ydiff)
    else:
        if abs(slope) == np.inf:
            return (0, from_y)
        if slope == 0:
            neg = -1 if x2 > x1 else 1
            return (neg*np.inf, from_x)
        
    # Get the tangent slope.
    tan_slope = -1 / slope

    # Get the tangent y-intercept
    tan_y_intercept = -tan_slope*from_x + from_y

    if ret_type == "two points":
        return ((from_x, from_y), (from_x+10, 10*tan_slope + from_y))
    else:
        return (tan_slope, tan_y_intercept)

if __name__ == "__main__":
    # testing distance_along_line
    print(f"{distance_along_line((1, 0), 1)=}")
    print(f"{distance_along_line((-1, 0), 1)=}")
    print(f"{distance_along_line((1/10, 0), 1)=}")
    print(f"{distance_along_line((10, 0), 1)=}")
    print(f"{distance_along_line((np.inf, 0), 1)=}")
    print(f"{distance_along_line((-np.inf, 0), 1)=}")
    print(f"{distance_along_line((0, 0), 1)=}")
    print(f"{distance_along_line((-0, 0), 1)=}")
    print(f"{distance_along_line(((0,0), (-10,0)), 1)=}")
    
    # testing get_tangent_line
    print(f"{get_tangent_line((1, 0), (0, 0))=}")
    print(f"{get_tangent_line((-1, 0), (0, 0))=}")
    print(f"{get_tangent_line((1/10, 0), (0, 0))=}")
    print(f"{get_tangent_line((10, 0), (0, 0))=}")
    print(f"{get_tangent_line((np.inf, 0), (0, 0))=}")
    print(f"{get_tangent_line((-np.inf, 0), (0, 0))=}")
    print(f"{get_tangent_line((0, 0), (0, 0))=}")
    print(f"{get_tangent_line((-0, 0), (0, 0))=}")
    print(f"{get_tangent_line(((0,0), (-10,0)), (0, 0))=}")
    print(f"{get_tangent_line(((0,0), (5,5)), (0, 0))=}")

    # testing line_two_points_to_slope_intercept
    print(f"{line_two_points_to_slope_intercept(((0, 0), (10, 0)))=}")
    print(f"{line_two_points_to_slope_intercept(((0, 0), (-10, 0)))=}")
    print(f"{line_two_points_to_slope_intercept(((0, 0), (0, 10)))=}")
    print(f"{line_two_points_to_slope_intercept(((0, 0), (0, -10)))=}")
    print(f"{line_two_points_to_slope_intercept(((0, 0), (10, 10)))=}")
    print(f"{line_two_points_to_slope_intercept(((0, 0), (10, -10)))=}")
    print(f"{line_two_points_to_slope_intercept(((0, 0), (-10, 10)))=}")
    print(f"{line_two_points_to_slope_intercept(((0, 0), (-10, -10)))=}")
    print(f"{line_two_points_to_slope_intercept(((10, 10), (0, 0)))=}")
    print(f"{line_two_points_to_slope_intercept(((10, -10), (0, 0)))=}")
    print(f"{line_two_points_to_slope_intercept(((-10, 10), (0, 0)))=}")
    print(f"{line_two_points_to_slope_intercept(((-10, -10), (0, 0)))=}")