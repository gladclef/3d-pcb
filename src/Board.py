import copy
import os
from typing import TypeVar

import pyvista
from scipy.spatial.transform import Rotation
import vtk

from Component.Component import Component
from Component.Shape import Shape
from FileIO.Line import Line
from FileIO.CadFileHelper import CadFileHelper
from Trace.SingleTrace import SingleTrace
import tool.vtk_tools as vt
from tool.globals import board_parameters as g

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

    def cleanup(self):
        """ Fixes and/or detects problems with the board that would make it difficult to boolean. """
        for trace in self.traces:
            trace.cleanup()

    @classmethod
    def from_cad_file(cls, gencad_file: str, limit_layers: list[str] = None) -> "Board":
        lines = Line.from_file(gencad_file)

        shapes_helper = CadFileHelper("$SHAPES", "$ENDSHAPES")
        components_helper = CadFileHelper("$COMPONENTS", "$ENDCOMPONENTS")

        def get_instances(cls: T, lines: list[Line]) -> tuple[list[T], list[Line]]:
            ret: list[T] = []

            instances, unmatched_lines = cls.from_cad_file(lines)
            while len(instances) > 0:
                ret += instances
                instances, unmatched_lines = cls.from_cad_file(unmatched_lines)
            
            return ret, unmatched_lines

        traces, lines = get_instances(SingleTrace, lines)
        traces: list[SingleTrace] = traces
        if limit_layers is not None:
            traces = list(filter(lambda t: t.layer in limit_layers, traces))
        pre_lines, shapes_lines, post_lines = shapes_helper.get_next_region(lines)
        lines = pre_lines + post_lines
        shapes: list[Shape] = get_instances(Shape, shapes_lines)[0]
        pre_lines, components_lines, post_lines = components_helper.get_next_region(lines)
        lines = pre_lines + post_lines
        components: list[Component] = get_instances(Component, components_lines)[0]

        for component in components:
            component.assign_shape(shapes)
        for trace in traces:
            trace.add_trace_end_pins(components)
        
        board = cls(gencad_file, traces, shapes, components)
        return board

    def to_vtk(self) -> tuple[vtk.vtkPolyData, vtk.vtkPolyData, vtk.vtkPolyData]:
        traces_polydata = vt.new_polydata()
        vias_polydata = vt.new_polydata()
        component_polydata = vt.new_polydata()

        for trace in self.traces:
            trace = copy.deepcopy(trace)
            trace.to_vtk(traces_polydata, vias_polydata)

        for component in self.components:
            component = component.get_transformed()
            component.to_vtk(component_polydata)

        return traces_polydata, vias_polydata, component_polydata
    
    def draw_board(self):
        """
        For debugging: draw simplified versions of the traces and vias that
        will be generated for this board.
        """
        import matplotlib.pyplot as plt
        
        # create the plot
        fig, ax = plt.subplots(figsize=(10,10))

        # draw the components
        for component in self.components:
            component = component.get_transformed()
            component.draw(ax)

        # draw the segments
        for trace in self.traces:
            trace.draw(ax)

        # show the plot
        plt.axis('equal')
        plt.show(block=True)

        
if __name__ == "__main__":
    g.TRACE_CORNER_RADIUS = 0.6

    example_name = "deej"
    example_dir = os.path.join(os.path.dirname(__file__), "..", "examples", example_name)
    limit_layers, layer_name = ["TOP"], "_top"

    board = Board.from_cad_file(os.path.join(example_dir, "exports", f"{example_name}.cad"), limit_layers=limit_layers)
    board.cleanup()
    board.draw_board()

    traces_pd, vias_pd, component_pd = board.to_vtk()
    pyvista.PolyData(traces_pd).save(os.path.join(example_dir, "stls", f"{example_name}_traces{layer_name}.stl"))
    pyvista.PolyData(vias_pd).save(os.path.join(example_dir, "stls", f"{example_name}_vias{layer_name}.stl"))
    pyvista.PolyData(component_pd).save(os.path.join(example_dir, "stls", f"{example_name}_components{layer_name}.stl"))

    polydata = vt.new_polydata()
    vt.join(polydata, traces_pd)
    vt.join(polydata, vias_pd)
    vt.join(polydata, component_pd)

    print(f"{polydata.GetNumberOfPoints()=}")
    mesh = pyvista.PolyData(polydata)
    pyvista.global_theme.allow_empty_mesh = True
    mesh.plot(show_edges=True, opacity=1, show_vertices=True)
    # pyvista.PolyDataFilters.plot_normals(mesh, mag=0.5, flip=False, faces=False, show_edges=True, opacity=0.95, show_verticies=True)
    mesh.save(os.path.join(example_dir, "stls", f"{example_name}_full{layer_name}.stl"))