import math
import re
from typing import Union

import numpy as np
import pyvista
from scipy.spatial.transform import Rotation
import vtk
from vtk.util import numpy_support # type: ignore
from vtkbool import vtkBool

from FileIO.CadFileHelper import CadFileHelper
import Geometry.geometry_tools as geo
from Geometry.LineSegment import LineSegment
from Trace.AbstractTrace import AbstractTrace
from Trace.AbstractVtkPointTracker import AbstractVtkPointTracker as PntInc
from Trace.PipeShape import PipeShape
from Trace.SingleTrace import SingleTrace
from Trace.TraceCorner import TraceCorner
from Trace.VtkPointGroup import VtkPointGroup
import tool.vtk_tools as vt


class Board:
    def __init__(self, gencad_file: str, traces: list[SingleTrace]):
        self.gencad_file = gencad_file
        self.traces = traces

    @classmethod
    def from_cad_file(cls, gencad_file: str):
        traces = []

        with open(gencad_file, "r") as fin:
            lines = fin.readlines()

        trace, unmatched_lines = SingleTrace.from_cad_file(lines)
        while trace is not None:
            traces.append(trace)
            trace, unmatched_lines = SingleTrace.from_cad_file(unmatched_lines)
        
        board = cls(gencad_file, traces)
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
            for seg in trace.segments:
                ax.arrow(seg.x1, seg.y1, seg.x2-seg.x1, seg.y2-seg.y1, color="teal", head_width=.3)

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