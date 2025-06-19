from typing import TypeVar

import pyvista
from scipy.spatial.transform import Rotation
import vtk

from Component.Component import Component
from Component.Shape import Shape
from FileIO.CadFileHelper import CadFileHelper
from Trace.SingleTrace import SingleTrace
import tool.vtk_tools as vt

T = TypeVar('T')


class Board:
    def __init__(self,
                 gencad_file: str,
                 traces: list[SingleTrace],
                 shapes: list[Shape],
                 components: list[Component]):
        self.gencad_file = gencad_file
        self.traces = traces
        self.shapes = shapes
        self.components = components

    @classmethod
    def from_cad_file(cls, gencad_file: str):
        with open(gencad_file, "r") as fin:
            lines = fin.readlines()

        shapes_helper = CadFileHelper("$SHAPES", "$ENDSHAPES")
        components_helper = CadFileHelper("$COMPONENTS", "$ENDCOMPONENTS")

        def get_instances(cls: T, lines: list[str]) -> tuple[list[T], list[str]]:
            ret: list[T] = []

            instance, unmatched_lines = cls.from_cad_file(lines)
            while instance is not None:
                ret.append(instance)
                instance, unmatched_lines = cls.from_cad_file(unmatched_lines)
            
            return ret, unmatched_lines

        traces, lines = get_instances(SingleTrace, lines)
        pre_lines, shapes_lines, post_lines = shapes_helper.get_next_region(lines)
        lines = pre_lines + post_lines
        shapes: list[Shape] = get_instances(Shape, shapes_lines)[0]
        pre_lines, components_lines, post_lines = components_helper.get_next_region(lines)
        lines = pre_lines + post_lines
        components: list[Component] = get_instances(Component, components_lines)[0]

        for component in components:
            component.assign_shape(shapes)
        
        board = cls(gencad_file, traces, shapes, components)
        return board

    def to_vtk(self, polydata: vtk.vtkPolyData) -> vtk.vtkPolyData:
        for trace in self.traces:
            trace.to_vtk(polydata)
        return polydata
    
    def draw_board(self):
        """
        For debugging: draw simplified versions of the traces and vias that
        will be generated for this board.
        """
        import matplotlib.pyplot as plt
        
        # create the plot
        fig, ax = plt.subplots(figsize=(10,10))

        # draw the segments
        for trace in self.traces:
            trace.draw(ax)

        # draw the components
        for component in self.components:
            component = component.get_transformed()
            component.draw(ax)

        # show the plot
        plt.show(block=True)

        
if __name__ == "__main__":
    board = Board.from_cad_file("../test schematics/hello_light/exports/hello_light.cad")
    board.draw_board()

    polydata = vt.new_polydata()
    polydata = board.to_vtk(polydata)

    print(f"{polydata.GetNumberOfPoints()=}")
    test_trace_mesh = pyvista.PolyData(polydata)
    pyvista.global_theme.allow_empty_mesh = True
    test_trace_mesh.plot(show_edges=True, opacity=1, show_vertices=True)
    # pyvista.PolyDataFilters.plot_normals(test_trace_mesh, mag=0.5, flip=False, faces=False, show_edges=True, opacity=0.95, show_verticies=True)