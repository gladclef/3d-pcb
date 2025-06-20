import copy
import re

import matplotlib.axis as maxis
import numpy as np

from Component.Shape import Shape
from Component.Text import Text
from FileIO.CadFileHelper import CadFileHelper
from tool.units import *

class Component:
    """
    A class representing an electronic component from a CAD file.
    """

    def __init__(self,
                 name: str,
                 device_name: str,
                 translation: tuple[float, float],
                 layer: str,
                 rotation: float,
                 shape_name: str,
                 texts: list[Text],
                 sheet: str):
        """
        Parameters
        ----------
        name : str
            The name of the component.
        device_name : str
            The device name associated with the component.
        translation : tuple[float, float]
            The placement coordinates as (x, y).
        layer : str
            The layer on which this component is placed.
        rotation : float
            Rotation angle in radians.
        shape_name : str
            Name of the shape associated with this component.
        texts : list[Text]
            List of text elements related to this component.
        sheet : str
            Additional information about the component from the CAD file.
        """
        self.name = name
        """ The name of the component. """
        self.device_name = device_name
        """ The device name associated with the component. """
        self.translation = translation
        """ The placement coordinates as (x, y). """
        self.layer = layer
        """ The layer on which this component is placed. """
        self.rotation = rotation
        """ Rotation angle in radians. """
        self.shape_name = shape_name
        """ Name of the shape associated with this component. """
        self.texts = texts
        """ List of text elements related to this component. """
        self.sheet = sheet
        """ Additional information about the component from the CAD file. """

        self.shape: Shape = None
        """ The assigned shape object for this component if any. """

    def assign_shape(self, available_shapes: list[Shape]) -> bool:
        """
        Assign a shape to the component from a list of available shapes.

        Parameters
        ----------
        available_shapes : list[Shape]
            List of available shapes to match against.

        Returns
        -------
        bool
            True if a matching shape was found and assigned, otherwise False.
        """
        for shape in available_shapes:
            if shape.full_name == self.shape_name:
                self.shape = shape
                return True
        return False

    def get_transformed(self) -> "Component":
        """ Get a copy of this component with translation, rotation, and layer flip applied. """
        if self.shape is None:
            raise RuntimeError("ERror in Component.get_transformed(): " + "no assigned shape, call 'assign_shape()' first!")

        is_bottom = self.layer.lower() == "bottom"
        transformed_shape = self.shape.apply_translation_rotation_layer(self.translation, self.rotation, is_bottom)

        new_texts: list[Text] = []
        for text in self.texts:
            new_text = text.apply_translation_rotation_layer(self.translation, self.rotation, is_bottom)
            new_texts.append(new_text)

        ret = copy.deepcopy(self)
        ret.shape = transformed_shape
        ret.texts = new_texts
        
        return ret

    @classmethod
    def from_cad_file(cls, lines: list[str], shapes: list[Shape]=None) -> tuple["Component", list[str]]:
        """
        Creates a Component instance by parsing lines from a CAD file.

        Parameters
        ----------
        lines : list[str]
            Lines of text representing the component in a CAD file.

        Returns
        -------
        tuple[Component, list[str]]
            A tuple containing the created Component and any remaining unprocessed lines.
        """
        # example component lines:
        #     COMPONENT "R1"
        #     DEVICE "DEV_Resistor_THT:R_Axial_DIN0309_L9.0mm_D3.2mm_P12.70mm_Horizontal"
        #     PLACE 2.40157 -1.53543
        #     LAYER TOP
        #     ROTATION 180
        #     SHAPE "Resistor_THT:R_Axial_DIN0309_L9.0mm_D3.2mm_P12.70mm_Horizontal" 0 0
        #     TEXT 0.25 0.107087 0.0393701 180 0 SILKSCREEN_TOP "R1" 0 0 0.0867079 0.0667913
        #     TEXT 0.25 -0.107087 0.0393701 180 0 SILKSCREEN_TOP "250" 0 0 0.122328 0.0667913
        #     SHEET "RefDes: R1, Value: 250"

        helper = CadFileHelper(re.compile(r"COMPONENT.*"), re.compile(r"^(\$ENDCOMPONENTS)?$"), ignore_false_endings=True)
        pre_lines, component_lines, post_lines = helper.get_next_region(lines)
        if len(component_lines) == 0:
            return None, lines

        component_name: str = None
        device_name: str = None
        place: tuple[float, float] = None
        layer: str = None
        rotation: float = None
        shape_name: str = None
        texts: list[Text] = []
        sheet: str = None

        for line in component_lines:
            if line.startswith("COMPONENT"):
                component_name = line.split('"', 1)[1].rstrip().rstrip('"')
            elif line.startswith("DEVICE"):
                device_name = line.split('"', 1)[1].rstrip().rstrip('"')
            elif line.startswith("PLACE"):
                place_coords = re.findall(r"[-\d\.]+", line)
                place = in2mm(float(place_coords[0])), in2mm(float(place_coords[1]))
            elif line.startswith("LAYER"):
                layer = line.split(maxsplit=1)[1].strip()
            elif line.startswith("ROTATION"):
                rotation = np.deg2rad(float(line.split()[1]))
            elif line.startswith("SHAPE"):
                shape_name = re.findall(r'"(.*?)"', line)
                if shape_name:
                    shape_name = shape_name[0]
            elif line.startswith("TEXT"):
                text, _ = Text.from_cad_file([line])
                assert text is not None
                texts.append(text)
            elif line.startswith("SHEET"):
                sheet = line.split('"', 1)[1].rstrip().rstrip('"')

        component = Component(
            component_name,
            device_name,
            place,
            layer,
            rotation,
            shape_name,
            texts,
            sheet
        )
        if shapes is not None:
            component.assign_shape(shapes)

        return component, pre_lines + post_lines
    
    def draw(self, ax: maxis.Axis):
        self.shape.draw(ax)

        for text in self.texts:
            text.draw(ax)

    def __repr__(self) -> str:
        if len(self.texts) > 0:
            return f"{self.texts[0].value} ({self.texts[1].value})"
        else:
            return self.device_name