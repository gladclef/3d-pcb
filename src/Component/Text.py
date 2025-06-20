import copy
import re

import matplotlib.axis as maxis
import numpy as np

import Geometry.geometry_tools as geo
from tool.units import *


class Text:
    def __init__(self,
                 value: str,
                 font_size: float,
                 x_offset: float,
                 y_offset: float,
                 width: float,
                 height: float,
                 rotation: float,
                 layer: str):
        self.value = value
        self.font_size = font_size
        self.x_offset = x_offset
        self.y_offset = y_offset
        self.width = width
        self.height = height
        self.rotation = rotation
        self.layer = layer
    
    def apply_translation_rotation_layer(self, translation: tuple[float, float], rotation: float, is_bottom: bool) -> "Text":
        """
        Apply translation and rotation to the text, with optional flipping.

        Parameters
        ----------
        translation : tuple[float, float]
            Translation vector (x, y).
        rotation : float
            Rotation angle in radians.
        is_bottom : bool
            If True, flip the coordinates along the x-axis after applying translation and rotation.

        Returns
        -------
        Text
            A new Text instance with the applied transformations.

        """
        x, y = geo.apply_translation_rotation_flip((self.x_offset, self.y_offset), translation, rotation, is_bottom)

        ret = copy.deepcopy(self)
        ret.x_offset, ret.y_offset = x, y
        # ret.rotation += rotation

        return ret

    @classmethod
    def from_cad_file(cls, lines: list[str]) -> tuple["Text", list[str]]:
        # example text line:
        #     TEXT 0.25 0.107087 0.0393701 180 0 SILKSCREEN_TOP "R1" 0 0 0.0867079 0.0667913
        text_line = None
        for line_idx, line in enumerate(lines):
            if line.startswith("TEXT"):
                text_line = line
                ret_lines = lines[:line_idx] + lines[line_idx+1:]
        if text_line is None:
            return None, lines
        
        # regex explanation:             x offset  y offset   font sz?   rotation   ???        layer     value      ???        ???        text width text height
        text_pattern = re.compile(r"TEXT ([\d\.]+) ([-\d\.]+) ([-\d\.]+) ([-\d\.]+) ([-\d\.]+) ([^ ]+) \"([^\"]+)\" ([-\d\.]+) ([-\d\.]+) ([-\d\.]+) ([-\d\.]+)")

        # break out the various parts of the text line
        match = text_pattern.match(text_line)
        if match is None:
            raise RuntimeError("Error in Text.from_cad_file(): failed to match text_pattern to line:\n\t" + text_line)
        x_offset, y_offset, font_size, rotation, _, layer, value, _, _, width, height = match.groups()
        x_offset_mm = in2mm(float(x_offset))
        y_offset_mm = in2mm(float(y_offset))
        rotation_rad = geo.normalize_angle(np.deg2rad(float(rotation)))
        width_mm = in2mm(float(width))
        height_mm = in2mm(float(height))

        text = Text(value, font_size, x_offset_mm, y_offset_mm, width_mm, height_mm, rotation_rad, layer)
        return text, ret_lines

    def draw(self, ax: maxis.Axis):
        text_loc = self.x_offset, self.y_offset
        rotation = np.rad2deg(self.rotation)
        if rotation > 170:
            rotation -= 180
        ax.text(text_loc[0], text_loc[1], self.value, color="black", fontsize=10, rotation=rotation, horizontalalignment='center', verticalalignment='center')