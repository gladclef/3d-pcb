import math
from typing import overload

import numpy as np
from scipy.spatial.transform import Rotation
import tool.array_tools as at

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

def apply_translation_rotation_flip(xy_or_xyz: tuple[float, float] | tuple[float, float, float] | np.ndarray,
                                    translation: tuple[float, float],
                                    rotation: float,
                                    flip: bool) -> np.ndarray | tuple:
    """
    Apply translation, rotation, and optional flipping to a point or vector.

    Parameters
    ----------
    xy_or_xyz : Union[tuple[float, float], tuple[float, float, float], np.ndarray]
        A tuple or numpy array representing the input coordinates.
        Can be 2D (x, y) or 3D (x, y, z).
    translation : tuple[float, float]
        Translation vector (tx, ty).
    rotation : float
        Rotation angle in radians around the z-axis.
    flip : bool
        If True, flip the coordinates along the x-axis and invert the rotation.

    Returns
    -------
    xy_or_xyz_transformed: np.ndarray | tuple
        The transformed coordinates. The return type matches the input type.

    """
    (x, y, z), input_len, input_type = at.tuple_from_array(xy_or_xyz, [2,3], 3)

    # apply flip property
    if flip:
        x = -x
        rotation = -rotation

    # apply rotation
    r = Rotation.from_euler('z', rotation)
    vec = np.array([x, y, z])
    rotated = r.apply(vec)
    x, y, z = tuple(rotated.tolist())

    # apply translation
    x += translation[0]
    y += translation[1]

    # convert to the desired return type
    ret = at.retval_from_tuple((x, y, z), input_len, input_type)

    return ret