import math
from typing import overload

import numpy as np

INF_THRESH = 1e6
""" The threshold at which values are considered functionally equivalent to infinity. """
ZERO_THRESH = 1/INF_THRESH
""" The threshold at which values are considered functionally equivalent to zero. """


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
