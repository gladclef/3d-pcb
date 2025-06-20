import copy
import re

import matplotlib.axis as maxis

from Component.Arc import Arc
from Component.Circle import Circle
from Component.Line import Line
from Component.Pin import Pin
from FileIO.CadFileHelper import CadFileHelper

class Shape:
    """
    A class representing an electronic component's shape.
    """

    def __init__(self, short_name: str, type_name: str, pins: list[Pin], lines: list[Line], arcs: list[Arc], circles: list[Circle]):
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
        self.lines = lines
        """List of `Line` objects for outlining this shape."""
        self.arcs = arcs
        """List of `Arc` objects for outlining this shape."""
        self.circles = circles
        """List of `Circle` objects for outlining this shape."""

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
        def apply_to_children(children: list) -> list:
            new_children: list[Pin] = []
            
            for child in children:
                new_child = child.apply_translation_rotation_layer(translation, rotation, is_bottom)
                new_children.append(new_child)

            return new_children

        new_pins = apply_to_children(self.pins)
        new_lines = apply_to_children(self.lines)
        new_arcs = apply_to_children(self.arcs)
        new_circles = apply_to_children(self.circles)

        ret = copy.deepcopy(self)
        ret.pins = new_pins
        ret.lines = new_lines
        ret.arcs = new_arcs
        ret.circles = new_circles

        return ret

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
        helper = CadFileHelper(re.compile(r"SHAPE.*"), re.compile(r"^(\$ENDSHAPES)?$"), ignore_false_endings=True)
        pre_lines, shape_lines, post_lines = helper.get_next_region(lines)
        if len(shape_lines) == 0:
            return None, lines
        
        # debugging
        # print("In Shape.from_cad_file(), shape_lines are:\n\t" + "\t".join(shape_lines))

        # example "SHAPE" line:
        # SHAPE "Resistor_THT:R_Axial_DIN0309_L9.0mm_D3.2mm_P12.70mm_Horizontal"
        shape_header_line = shape_lines[0]
        assert shape_header_line.startswith("SHAPE ")
        full_name = shape_header_line.split('"', 1)[1].rstrip().rstrip('"')
        short_name, type_name = tuple(full_name.split(":", 1))

        def children_from_cad_file(cls, unmatched_lines: list[str]) -> tuple[list, list[str]]:
            children = []

            child, unmatched_lines = cls.from_cad_file(unmatched_lines)
            while child is not None:
                children.append(child)
                child, unmatched_lines = cls.from_cad_file(unmatched_lines)
            
            return children, unmatched_lines
        
        unmatched_lines = shape_lines[1:]
        pins, unmatched_lines = children_from_cad_file(Pin, unmatched_lines)
        lines, unmatched_lines = children_from_cad_file(Line, unmatched_lines)
        arcs, unmatched_lines = children_from_cad_file(Arc, unmatched_lines)
        circles, unmatched_lines = children_from_cad_file(Circle, unmatched_lines)

        shape = cls(short_name, type_name, pins, lines, arcs, circles)
        return shape, pre_lines + post_lines
    
    def draw(self, ax: maxis.Axis):
        for line in self.lines:
            line.draw(ax)
        
        for arc in self.arcs:
            arc.draw(ax)
        
        for circle in self.circles:
            circle.draw(ax)
        
        for pin in self.pins:
            pin.draw(ax)
