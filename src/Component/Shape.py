import copy
import re

from Component.Pin import Pin
from FileIO.CadFileHelper import CadFileHelper

class Shape:
    """
    A class representing an electronic component's shape.
    """

    def __init__(self, short_name: str, type_name: str, pins: list[Pin]):
        """
        Parameters
        ----------
        short_name : str
            Short name of the shape.
        type_name : str
            Type name of the shape.
        pins : list[Pin]
            List of `Pin` objects associated with this shape. Can be pads or vias.
        """
        self.short_name = short_name
        """Short name of the shape."""
        self.type_name = type_name
        """Type name of the shape."""
        self.pins = pins
        """List of `Pin` objects associated with this shape. Can be pads or vias."""

    @property
    def full_name(self) -> str:
        """
        Returns the full name of the shape. Used by the component class to assign shapes.
        """
        return self.short_name + ":" + self.type_name

    def apply_translation_rotation_layer(self, translation: tuple[float, float], rotation: float, is_bottom: bool) -> "Shape":
        """
        Applies a transformation to this `Shape` by translating and rotating its pins.

        This method creates a new instance of the shape with transformed pins.

        Parameters
        ----------
        translation : tuple[float, float]
            Translation vector as (x, y).
        rotation : float
            Rotation angle in radians.
        is_bottom : bool
            Flag indicating that this shape is located on the bottom layer.

        Returns
        -------
        Shape
            A new instance of `Shape` with transformed pins.
        """
        new_pins: list[Pin] = []
        for pin in self.pins:
            pin = copy.deepcopy(pin)
            pin.apply_translation_rotation_layer(translation, rotation, is_bottom)
            new_pins.append(pin)

        ret = copy.deepcopy(self)
        ret.pins = new_pins

    @classmethod
    def from_cad_file(cls, lines: list[str]) -> tuple["Shape", list[str]]:
        """
        Creates a `Shape` instance by parsing lines from a CAD file.

        Parameters
        ----------
        lines : list[str]
            Lines of text representing the shape in a CAD file.

        Returns
        -------
        tuple[Shape, list[str]]
            A tuple containing the created `Shape` and any remaining unprocessed lines.
        """
        helper = CadFileHelper(re.compile(r"SHAPE.*"), re.compile(r"^(\$ENDSHAPES)?$"))
        pre_lines, shape_lines, post_lines = helper.get_next_region(lines)
        if len(shape_lines) == 0:
            return None, lines

        # example "SHAPE" line:
        # SHAPE "Resistor_THT:R_Axial_DIN0309_L9.0mm_D3.2mm_P12.70mm_Horizontal"
        shape_header_line = shape_lines[0]
        assert shape_header_line.startswith("SHAPE ")
        full_name = shape_header_line.split('"', 1)[1].rstrip().rstrip('"')
        short_name, type_name = tuple(full_name.split(":", 1))

        pins = []
        pin, unmatched_lines = Pin.from_cad_file(lines)
        while pin is not None:
            pins.append(pin)
            pin, unmatched_lines = Pin.from_cad_file(unmatched_lines)

        shape = cls(short_name, type_name, pins)
        return shape, pre_lines + post_lines
