import re
from typing import Union

import matplotlib.axis as maxis
import numpy as np
import pyvista
import vtk

from FileIO.CadFileHelper import CadFileHelper
import Geometry.geometry_tools as geo
from Geometry.LineSegment import LineSegment
from Trace.AbstractTrace import AbstractTrace
from Trace.PipeShape import PipeShape
from Trace.TraceCorner import TraceCorner
from Trace.VtkPointGroup import VtkPointGroup
import tool.vtk_tools as vt
from tool.units import *


class SingleTrace(AbstractTrace):
    """
    Represents a singular wire trace on a PCB.

    This is limited to a single continuous line. If there are
    junctions they should be built up of multiple SingleTrace objects
    using the more general Trace class.
    """
    def __init__(self, xy_points: list[tuple[float, float]], segments: list[tuple[int, int]] | list[LineSegment], shape: PipeShape=None, bend_radius: float=1, allow_overlap=False):
        super().__init__(xy_points, segments, shape)

        # set some defualts
        if bend_radius is None:
            bend_radius = 1

        self._xypnt_vtk_verticies: dict[int, VtkPointGroup] = {}
        """ Dictionary from xy point index to vtk points. """
        self.xypnt_trace_corners: dict[int, TraceCorner] = {}
        """ Dictionary from xy point index to trace corners. """
        self.bend_radius = bend_radius

        self.check_segment_duplicates(self.segments)
        if not allow_overlap:
            self.check_segments_overlap()
    
    def check_segment_duplicates(self, segments: list[LineSegment]):
        segment_tuples = [(s.xy1, s.xy2) for s in segments]
        segments_set = set(segment_tuples)
        if len(segment_tuples) != len(segments_set):
            raise ValueError(f"Error in SingleTrace.check_segment_duplicates(): there are {len(segments_set)} non-duplicate segments out of {len(segments)} total segments.")

    def check_segments_overlap(self):
        for s1idx, segment1 in enumerate(self.segments):
            for s2idx, segment2 in enumerate(self.segments):
                if s1idx == s2idx:
                    # don't check for self-segment intersections
                    continue

                if segment1.xy1 == segment2.xy1 or \
                    segment1.xy1 == segment2.xy2 or \
                    segment1.xy2 == segment2.xy1 or \
                    segment1.xy2 == segment2.xy2:
                    # don't check for intersections with segments that share one of the end points
                    continue

                # check if these segments are parallel
                angle_diff = abs(segment1.angle - segment2.angle)
                if angle_diff < geo.ZERO_THRESH or 2*np.pi - angle_diff < geo.ZERO_THRESH:
                    continue

                intersection = segment1.intersection(segment2)
                if intersection is not None:
                    raise ValueError(f"Error in SingleTrace.check_segments_overlap(): segments {s1idx} ({segment1}) and {s2idx} ({segment2}) overlap at [{intersection}].")

    def get_trace_corner(self, xy_idx: int) -> TraceCorner | None:
        """
        Get the trace corner for the given xy index.
        Builds the trace corner as necessary.

        Note: returns None for the first and last xy indicies,
        since there is no corner at the ends of the trace.

        Parameters
        ----------
        xy_idx : int
            The xy point index, for self.xy_points.

        Returns
        -------
        TraceCorner
            The corner for the given xy point, or None.
        """
        if xy_idx == 0 or xy_idx == len(self.xy_points)-1:
            return None
        
        if xy_idx not in self.xypnt_trace_corners:
            seg_idx = xy_idx-1
            segment_a = self.segments[seg_idx]
            segment_b = self.segments[seg_idx+1]
            self.xypnt_trace_corners[xy_idx] = TraceCorner(self, (segment_a, segment_b), self.bend_radius)
        
        return self.xypnt_trace_corners[xy_idx]

    def get_xypnt_vtk_verticies(self, xy_idx: int, segment: LineSegment) -> VtkPointGroup:
        """ Get the verticies for the given end-point of a trace segment.

        Parameters
        ----------
        xy_idx : int
            Either end of an trace segment, from self.xy_points.
        angle : float
            The angle of the trace segment in radians, as determined by its two end points.
        """
        if xy_idx not in self._xypnt_vtk_verticies:
            xy_point = self.xy_points[xy_idx]
            corner = self.get_trace_corner(xy_idx)

            if corner is None:
                angle = segment.angle
                xyz_points = np.array(self.shape.oriented_points(angle, xy_point))
                self._xypnt_vtk_verticies[xy_idx] = VtkPointGroup(xyz_points)
            
            else:
                if segment.xy1 == xy_point:
                    return corner.get_vtk_group(corner.n_points-1)
                elif segment.xy2 == xy_point:
                    return corner.get_vtk_group(0)
                else:
                    raise RuntimeError("Error in SingleTrace.get_xypnt_vtk_verticies(): " + f"expected the xy_point to be at the beggining or end of the given segment, " + f"but {xy_point=} and {segment.xy1=} and {segment.xy2=}!")

        return self._xypnt_vtk_verticies[xy_idx]
    
    def segment_to_vtk(self, polydata: vtk.vtkPolyData, segment: LineSegment):
        # get core segment values
        xy_a, xy_b = self.segment_xypnts(segment)
        xy_idx_a, xy_idx_b = self.xy_points.index(xy_a), self.xy_points.index(xy_b)

        # get the verticies at either end of the trace segment
        va = self.get_xypnt_vtk_verticies(xy_idx_a, segment)
        vb = self.get_xypnt_vtk_verticies(xy_idx_b, segment)

        # Add to vtk points
        va.add_missing_vtk_points(polydata)
        vb.add_missing_vtk_points(polydata)
        
        # build the sides along the length of the trace
        vt.adjoin_with_quads(polydata, va.vtk_idx_0, vb.vtk_idx_0, len(va))
        
        # Assign point normals
        start_idx = min(va.vtk_idx_0, vb.vtk_idx_0)
        end_idx = max(max(va.vtk_indices), max(vb.vtk_indices)) + 1
        vt.calculate_point_normals(polydata, start_idx, end_idx)
    
    def inc_vtk_indicies(self, start: int, cnt: int):
        for segment_points in self._xypnt_vtk_verticies.values():
            segment_points.inc_vtk_indicies(start, cnt)

    @classmethod
    def from_cad_file(cls, cad_lines: list[str], shape: PipeShape=None, bend_radius: float=None) -> tuple[Union["SingleTrace",None], list[str]]:
        routes_helper = CadFileHelper("$ROUTES", "$ENDROUTES")
        route_helper = CadFileHelper(re.compile(r"ROUTE.*"), re.compile(r"(ROUTE.*|\$ENDROUTES)"))

        # get the lines from the cad file for the next route
        pre_routes, routes, post_routes = routes_helper.get_next_region(cad_lines)
        if len(routes) == 0:
            return None, cad_lines
        
        pre_route, route, post_route = route_helper.get_next_region(routes)
        if len(route) == 0:
            return None, cad_lines
        if not route[-1].startswith("$ENDROUTES"):
            post_route.insert(0, route[-1])
        route = route[:-1]
        
        # parse the lines for this route
        segment_lines = list(filter(lambda l: l.startswith("LINE "), route))
        xy_points_orig: list[tuple[float, float]] = []
        edges: list[tuple[int, int]] = []
        for segment_line in segment_lines:

            x1, y1, x2, y2 = tuple(map(float, segment_line.strip()[5:].split(" ")))
            xy1, xy2 = (x1, y1), (x2, y2)
            if xy1 not in xy_points_orig:
                xy_points_orig.append(xy1)
            if xy2 not in xy_points_orig:
                xy_points_orig.append(xy2)
            
            xy1_idx, xy2_idx = xy_points_orig.index(xy1), xy_points_orig.index(xy2)
            edges.append((xy1_idx, xy2_idx))
        xy_points: list[tuple[float, float]] = [(in2mm(x), in2mm(y)) for x, y in xy_points_orig]
        
        if len(edges) > 1:
            # orient adjacent edges
            for edge_idx, (xy1_idx, xy2_idx) in enumerate(edges[:-1]):
                next_edge = edges[edge_idx]
                if xy1_idx in next_edge:
                    edges[edge_idx] = (xy2_idx, xy1_idx)
                    
            # orient the last edge
            xy1_idx, xy2_idx = edges[-1]
            prev_edge = edges[-2]
            if xy2_idx in prev_edge:
                edges[-1] = (xy2_idx, xy1_idx)

        # debugging
        for edge in edges:
            print(f"LINE   {xy_points_orig[edge[0]]}   {xy_points_orig[edge[1]]}")

        ret = cls(xy_points, edges, shape, bend_radius)
        return ret, pre_routes + pre_route + post_route + post_routes

    def to_vtk(self, polydata: vtk.vtkPolyData):
        vtk_points: vtk.vtkPoints = polydata.GetPoints()
        vtk_cells: vtk.vtkCellData = polydata.GetCellData()

        # build the segments verts and cells
        for segment in self.segments:
            self.segment_to_vtk(polydata, segment)

        # build the corners
        for xy_idx in range(len(self.xy_points)):
            corner = self.get_trace_corner(xy_idx)
            if corner is not None:
                corner.to_vtk(polydata)
            
        # Close the ends
        xy_idx_a = 0
        xy_idx_b = max(list(self._xypnt_vtk_verticies.keys()))
        for xy_idx in [xy_idx_a, xy_idx_b]:
            verticies = self._xypnt_vtk_verticies[xy_idx]
            point_id_0 = verticies.vtk_indices[0]
            new_vtk_points = self.shape.add_vtk_cells(polydata, point_id_0)
            self.inc_vtk_indicies(point_id_0 + len(verticies), len(new_vtk_points))
            for new_vtk_id in new_vtk_points:
                xyz_point = list(vtk_points.GetPoint(new_vtk_id))
                verticies.xyz_points = np.concat((verticies.xyz_points, np.array([xyz_point])), axis=0)
                verticies.vtk_indices.append(new_vtk_id)
        
        return polydata
    
    def draw(self, ax: maxis.Axis):
        for seg in self.segments:
            ax.arrow(seg.x1, seg.y1, seg.x2-seg.x1, seg.y2-seg.y1, color="teal", head_width=.3)


def _tst_trace_from_points():
    # simple trace
    trace_points = [
        (0, 0),
        (1, 0),
        (2, 1),
    ]
    trace_edges = [
        (0, 1),
        (1, 2)
    ]

    # hello trace
    trace_points = [
        # h
        (0, 0),
        (1, 10),
        (0, 10),
        (0, 5),
        (5, 5),
        (5, 0),
        # e
        (12, 4),
        (12, 5),
        (7, 5),
        (7, 0),
        # l
        (14, 0),
        (15, 10),
        (14, 10),
        (15, 0),
        # l
        (21, 0),
        (22, 10),
        (21, 10),
        (22, 0),
        # o
        (33, 0),
        (33, 5),
        (28, 5),
        (28, 0)
    ]
    trace_edges = [(i, i+1) for i in range(len(trace_points)-1)]

    test_trace = SingleTrace(trace_points, trace_edges, PipeBasicBox(awg2mm(26)), allow_overlap=True)
    return test_trace


if __name__ == "__main__":
    from Trace.PipeShape import PipeBasicBox
    from tool.units import *

    test_trace = _tst_trace_from_points()

    vtk_test_trace = test_trace.to_vtk(vt.new_polydata())
    print(f"{vtk_test_trace.GetNumberOfPoints()=}")
    vt.save_to_vtk(vtk_test_trace, "test_trace.vtk")

    test_trace_mesh = pyvista.PolyData(vtk_test_trace)
    pyvista.global_theme.allow_empty_mesh = True
    test_trace_mesh.plot(show_edges=True, opacity=1, show_vertices=True)
    # pyvista.PolyDataFilters.plot_normals(test_trace_mesh, mag=0.5, flip=False, faces=False, show_edges=True, opacity=0.95, show_verticies=True)