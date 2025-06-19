import re

import numpy as np

from FileIO.CadFileHelper import CadFileHelper
from tool.units import *


class Text:
    def __init__(self,
                 value: str,
                 font_size: float,
                 group_x_offset: float,
                 group_y_offset: float,
                 val_x_offset: float,
                 val_y_offset: float,
                 rotation: float,
                 layer: str):
        self.value = value
        self.font_size = font_size
        self.group_x_offset = group_x_offset
        self.group_y_offset = group_y_offset
        self.val_x_offset = val_x_offset
        self.val_y_offset = val_y_offset
        self.rotation = rotation
        self.layer = layer

    @property
    def x_offset(self) -> float:
        return self.group_x_offset + self.val_x_offset
    
    @property
    def y_offset(self) -> float:
        return self.group_y_offset + self.val_y_offset

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
        
        # regex explanation:             font sz   group x    group y    rotation   ???        layer     value      ???        ???        x offset   y offset
        text_pattern = re.compile(r"TEXT ([\d\.]+) ([-\d\.]+) ([-\d\.]+) ([-\d\.]+) ([-\d\.]+) ([^ ]+) \"([^\"]+)\" ([-\d\.]+) ([-\d\.]+) ([-\d\.]+) ([-\d\.]+)")

        # break out the various parts of the text line
        match = text_pattern.match(text_line)
        if match is None:
            raise RuntimeError("Error in Text.from_cad_file(): failed to match text_pattern to line:\n\t" + text_line)
        font_size, group_x_offset, group_y_offset, rotation, _, layer, value, _, _, val_x_offset, val_y_offset = match.groups()
        group_x_offset = in2mm(float(group_x_offset))
        group_y_offset = in2mm(float(group_y_offset))
        rotation = np.deg2rad(float(rotation))
        val_x_offset = in2mm(float(val_x_offset))
        val_y_offset = in2mm(float(val_y_offset))

        text = Text(value, font_size, group_x_offset, group_y_offset, val_x_offset, val_y_offset, rotation, layer)
        return text, ret_lines