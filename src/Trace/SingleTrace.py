import copy
import re
from typing import Union

import matplotlib.axis as maxis
import numpy as np
import pyvista
import vtk

from Component.Component import Component
from Component.Pin import Pin
from FileIO.CadFileHelper import CadFileHelper
from FileIO.Line import Line as FLine
import Geometry.geometry_tools as geo
from Geometry.LineSegment import LineSegment
from Trace.AbstractTrace import AbstractTrace
from Trace.PipeShape import PipeShape, DEFAULT_PIPE_SHAPE
from Trace.TraceCorner import TraceCorner
from Trace.VtkPointGroup import VtkPointGroup
from tool.globals import board_parameters as g
import tool.vtk_tools as vt
from tool.units import *


ALLOW_MULTIPLE_TRACES_PER_ROUTE = True


class SingleTrace(AbstractTrace):
    """
    Represents a singular wire trace on a PCB.

    This is limited to a single continuous line. If there are
    junctions they should be built up of multiple SingleTrace objects
    using the more general Trace class.
    """

    def __init__(self,
                 source_lines: list[FLine],
                 layer: str,
                 xy_points: list[tuple[float, float]],
                 segments: list[tuple[int, int] | tuple[int, int, FLine]] | list[LineSegment],
                 shape: PipeShape=None,
                 bend_radius: float=1,
                 allow_overlap=False):
        """
        Initializes a SingleTrace object.

        Parameters
        ----------
        xy_points : list[tuple[float, float]]
            List of XY coordinates defining the trace points.
        segments : Union[list[tuple[int, int]], list[LineSegment]]
            List of line segments. Each segment can be defined as either
            a tuple of xy_point indices or a LineSegment object.
        shape : PipeShape, optional
            Shape to extrude along the trace's path, or None to use the
            DEFAULT_PIPE_SHAPE. By default None.
        bend_radius : float, optional
            Radius for rounding out trace bends in millimeters. Trace segment
            intersections will be rounded out using TraceCorners to match this
            bend radius. By default TRACE_CORNER_RADIUS.
        allow_overlap : bool, optional
            If True, overlapping segments are allowed, by default False.
        """
        linenos = [l.lineno for l in source_lines]
        print(f"Creating {self.__class__.__name__} instance on layer {layer} from lines ({min(linenos)+1}-{max(linenos)+1})")# +
        #      ":\n\t" + "\n\t".join([l.v.rstrip() for l in source_lines]))
        super().__init__(source_lines, xy_points, segments, shape)

        # set some defaults
        if bend_radius is None:
            bend_radius = g.TRACE_CORNER_RADIUS

        self.layer = layer
        """ Which layer of the board this trace is on """
        self._xypnt_vtk_verticies: dict[tuple[float, float], VtkPointGroup] = {}
        """ Dictionary from xy point index to vtk points. """
        self.xypnt_trace_corners: dict[tuple[float, float], TraceCorner] = {}
        """ Dictionary from xy point index to trace corners. """
        self.bend_radius = bend_radius
        """
        Radius of the allowed trace bends. Trace segment intersections will
        be rounded out using TraceCorners to match this bend radius.
        """
        self.pins: dict[str, Pin] = { "a": None, "b": None }
        """ The pins that indicate where to put vias at either end of this trace. """

        self.check_segment_duplicates(self.segments)
        if not allow_overlap:
            self.check_segments_overlap()

    def check_segment_duplicates(self, segments: list[LineSegment]):
        """
        Checks for duplicate segments in the given list.

        Parameters
        ----------
        segments : list[LineSegment]
            List of line segments to be checked for duplicates.
        """
        segment_tuples = [(s.xy1, s.xy2) for s in segments]
        segments_set = set(segment_tuples)
        if len(segment_tuples) != len(segments_set):
            raise ValueError(f"Error in SingleTrace.check_segment_duplicates(): there are {len(segments_set)} non-duplicate segments out of {len(segments)} total segments.")

    def check_segments_overlap(self):
        """
        Checks for overlapping segments within the trace.
        """
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

    def get_trace_corner(self, xy_pnt: tuple[float, float], segment_idx: int, segment: LineSegment) -> TraceCorner | None:
        """
        Get the trace corner for the given xy point.
        Builds the trace corner as necessary.

        Note: returns None for the first and last xy points,
        since there is no corner at the ends of the trace.

        Parameters
        ----------
        xy_pnt : tuple[float, float]
            The xy point (probably from self.xy_points).
        segment_idx : int
            The index of the segment (probably from self.segments).
        segment : LineSegment
            The segment that xy_pnt belongs to.

        Returns
        -------
        TraceCorner | None
            The corner for the given xy point, or None.
        """
        if segment_idx == 0 and xy_pnt == segment.xy1:
            return None
        if segment_idx == len(self.segments)-1 and xy_pnt == segment.xy2:
            return None

        if xy_pnt not in self.xypnt_trace_corners:
            if xy_pnt == segment.xy2:
                segment_a, segment_b = segment, self.segments[segment_idx+1]
            elif xy_pnt == segment.xy1:
                segment_a, segment_b = self.segments[segment_idx-1], segment
            else:
                raise RuntimeError("Error in SingleTrace.get_trace_corner(): " + "xy_pnt != segment.xy1 and xy_pnt != segment.xy2")

            self.xypnt_trace_corners[xy_pnt] = TraceCorner(self, (segment_a, segment_b), self.bend_radius)

        return self.xypnt_trace_corners[xy_pnt]

    def get_xypnt_vtk_verticies(self, xy_pnt: tuple[float, float], segment_idx: int, segment: LineSegment) -> VtkPointGroup:
        """
        Get the vertices for the given end-point of a trace segment.

        Parameters
        ----------
        xy_pnt : tuple[float, float]
            The xy point (probably from self.xy_points).
        segment_idx : int
            The index of the segment (probably from self.segments).
        segment : LineSegment
            The segment that xy_pnt belongs to.

        Returns
        -------
        VtkPointGroup
            Vertices at the given point as a VtkPointGroup object.
        """
        if xy_pnt not in self._xypnt_vtk_verticies:
            corner = self.get_trace_corner(xy_pnt, segment_idx, segment)

            if corner is None:
                angle = segment.angle
                xyz_points = np.array(self.shape.oriented_points(angle, xy_pnt))
                self._xypnt_vtk_verticies[xy_pnt] = VtkPointGroup(xyz_points)

            else:
                if segment.xy1 == xy_pnt:
                    return corner.get_vtk_group(corner.n_points-1)
                elif segment.xy2 == xy_pnt:
                    return corner.get_vtk_group(0)
                else:
                    raise RuntimeError("Error in SingleTrace.get_xypnt_vtk_verticies(): " + f"expected the xy_point to be at the beginning or end of the given segment, " + f"but {xy_pnt=} and {segment.xy1=} and {segment.xy2=}!")

        return self._xypnt_vtk_verticies[xy_pnt]

    def add_trace_end_pins(self, components: list[Component]):
        """
        Find the component pins matching this trace and adds those
        pins to the end points of this trace.

        Parameters
        ----------
        components : list[Component]
            List of components to search for closest through-hole pins.
        """
        # get the ends of the first and last segments
        xy_locs = { "a": self.segments[0].xy1, "b": self.segments[-1].xy2 }

        # find the component through holes that most closely match this instance
        closest_pins: dict[str, Pin] = { "a": None, "b": None }
        for component in components:
            component = component.get_transformed()
            for pin in component.shape.pins:
                for end in ["a", "b"]:
                    closest_pin, xy_loc = closest_pins[end], xy_locs[end]

                    if abs(pin.x_offset - xy_loc[0]) < Pin.through_hole_diameter():
                        if abs(pin.y_offset - xy_loc[1]) < Pin.through_hole_diameter():
                            if closest_pin is None:
                                closest_pin = pin
                            else:
                                dist_pin = np.sqrt((pin.x_offset - xy_loc[0])**2 + (pin.y_offset - xy_loc[1])**2)
                                dist_closest = np.sqrt((closest_pin.x_offset - xy_loc[0])**2 + (closest_pin.y_offset - xy_loc[1])**2)
                                if dist_pin < dist_closest:
                                    closest_pin = pin

                    if closest_pin is not None:
                        closest_pins[end] = closest_pin

        # add the vias, as necessary
        for end in ["a", "b"]:
            closest_pin = closest_pins[end]
            if closest_pin is None:
                # don't need to add a via if the trace doesn't end at a pad
                continue
            if closest_pin.is_pad:
                # don't need to add a via if it's not a through-hole component
                continue

            self.pins[end] = closest_pin

    def segment_to_vtk(self, polydata: vtk.vtkPolyData, segment_idx: int, segment: LineSegment):
        """
        Inserts a segment of the trace as points and cells into the polydata.

        Parameters
        ----------
        polydata : vtk.vtkPolyData
            The VTK PolyData object where the segment will be added.
        segment_idx : int
            The index of the segment (probably from self.segments).
        segment : LineSegment
            The Segment to insert.
        """
        # get core segment values
        xy_a, xy_b = segment.xy1, segment.xy2

        # get the vertices at either end of the trace segment
        va = self.get_xypnt_vtk_verticies(xy_a, segment_idx, segment)
        vb = self.get_xypnt_vtk_verticies(xy_b, segment_idx, segment)

        # Add to vtk points
        va.add_missing_vtk_points(polydata)
        vb.add_missing_vtk_points(polydata)

        # build the sides along the length of the trace
        vt.adjoin_with_quads(polydata, va.vtk_idx_0, vb.vtk_idx_0, len(va))

        # Assign point normals
        start_idx = min(va.vtk_idx_0, vb.vtk_idx_0)
        mid_idx = va.vtk_idx_0 if start_idx == vb.vtk_idx_0 else vb.vtk_idx_0
        end_idx = max(max(va.vtk_indices), max(vb.vtk_indices)) + 1
        vt.calculate_point_normals(polydata, start_idx, mid_idx)
        vt.calculate_point_normals(polydata, mid_idx, end_idx)

    def inc_vtk_indicies(self, start: int, cnt: int):
        """
        Increments the VTK indices starting from a given index.

        Parameters
        ----------
        start : int
            Starting index to increment.
        cnt : int
            Number to increment the vertices by.
        """
        for segment_points in self._xypnt_vtk_verticies.values():
            segment_points.inc_vtk_indicies(start, cnt)

    @classmethod
    def get_lines_for_next_trace(cls, cad_lines: list[FLine]):
        routes_helper = CadFileHelper("$ROUTES", "$ENDROUTES")
        route_helper = CadFileHelper(re.compile(r"ROUTE .*"), re.compile(r"(ROUTE.*|\$ENDROUTES)"))
        layer_helper = CadFileHelper(re.compile(r"(LAYER|LINE) .*"), re.compile(r"(ROUTE.*|\$ENDROUTES|LAYER .*)"), allow_multiple_starts=True)

        # get the lines from the cad file for the next route
        pre_routes, routes, post_routes = routes_helper.get_next_region(cad_lines)
        if len(routes) == 0:
            return cad_lines, [], []
        if len(routes) == 1:
            # just the end matcher "$ENDROUTES"
            assert routes[0].v.strip() == "$ENDROUTES"
            return pre_routes, [], post_routes
        if len(routes) == 2:
            if routes[0].v.strip() == "$ROUTES" and routes[1].v.strip() == "$ENDROUTES":
                return pre_routes, [], post_routes
        assert routes[-1].v.strip() in ["ROUTE", "$ENDROUTES"]

        # break up on route
        pre_route, route, post_route = route_helper.get_next_region(routes)
        if len(route) == 0:
            return cad_lines, [], []
        if len(route) == 1:
            # just the end matcher "ROUTE" or "$ENDROUTES"
            assert route[0].v.strip() in ["ROUTE", "$ENDROUTES"]
            # check if there are any more routes
            pre_trace, trace, post_trace = cls.get_lines_for_next_trace(pre_routes + pre_route + post_route + post_routes)
            if len(trace) != 0:
                return pre_trace, trace, post_trace
            else:
                return pre_routes + pre_route, [], post_route + post_routes
        assert route[-1].v.strip().startswith("ROUTE") or route[-1].v.strip() == "$ENDROUTES"
        post_route.insert(0, route.pop())

        # break up on layers
        pre_layer, layer, post_layer = layer_helper.get_next_region(route[1:])
        if len(layer) <= 1:
            if len(layer) == 1:
                # just the end matcher "ROUTE", "LAYER", or "$ENDROUTES"
                assert layer[0].v.strip() != "LAYER"
                assert layer[0].v.strip().startswith("ROUTE") or layer[0].v.strip() == "$ENDROUTES"
            # check if there are any more routes
            assert all([l.v.strip().startswith("TRACK") for l in pre_layer])
            assert all([l.v.strip().startswith("TRACK") for l in post_layer])
            # pre_trace, trace, post_trace = cls.get_lines_for_next_trace(pre_routes + pre_route + pre_layer + post_layer + post_route + post_routes)
            pre_trace, trace, post_trace = cls.get_lines_for_next_trace(pre_routes + pre_route + post_route + post_routes)
            if len(trace) != 0:
                return pre_trace, trace, post_trace
            else:
                return pre_routes + pre_route + pre_layer, [], post_layer + post_route + post_routes
        else:
            # multiple layers for this route
            assert route[0].v.startswith("ROUTE ")
            pre_layer.insert(0, route[0])
            if layer_helper.end_matches(layer[-1]):
                post_layer.insert(0, layer.pop())

        return pre_routes + pre_route + pre_layer, layer, post_layer + post_route + post_routes

    @classmethod
    def from_cad_file(cls, cad_lines: list[FLine], shape: PipeShape=None, bend_radius: float=None) -> tuple[list["SingleTrace"], list[FLine]]:
        """
        Creates a SingleTrace object from CAD file lines.

        Parameters
        ----------
        cad_lines : list[FLine]
            Lines from the CAD file defining routes.
        shape : PipeShape, optional
            Shape to extrude along the trace's path, or None to use the
            DEFAULT_PIPE_SHAPE. By default None.
        bend_radius : float, optional
            Radius for rounding out trace bends in millimeters,
            by default TRACE_CORNER_RADIUS.

        Returns
        -------
        tuple[Union["SingleTrace",None], list[FLine]]
            A tuple containing a SingleTrace instance (if found) and remaining lines from the CAD file.
        """
        global ALLOW_MULTIPLE_TRACES_PER_ROUTE
        
        # get the lines
        pre_trace, trace, post_trace = cls.get_lines_for_next_trace(cad_lines)
        if len(trace) == 0:
            return [], pre_trace + post_trace
        
        # get the layer name
        name_lines = list(filter(lambda l: l.v.strip().startswith("LAYER "), trace))
        assert len(name_lines) <= 1
        if len(name_lines) == 0:
            layer_name = "TOP"
        else:
            layer_name = name_lines[0].v.split("LAYER ")[1].strip()

        # parse the lines for this route+layer
        segment_lines = list(filter(lambda l: l.v.startswith("LINE "), trace))
        xy_points_orig: list[tuple[float, float]] = []
        edges: list[tuple[int, int, FLine]] = []
        for segment_line in segment_lines:

            x1, y1, x2, y2 = tuple(map(float, segment_line.v.strip()[5:].split(" ")))
            xy1, xy2 = (x1, y1), (x2, y2)
            if xy1 not in xy_points_orig:
                xy_points_orig.append(xy1)
            if xy2 not in xy_points_orig:
                xy_points_orig.append(xy2)

            xy1_idx, xy2_idx = xy_points_orig.index(xy1), xy_points_orig.index(xy2)
            edges.append((xy1_idx, xy2_idx, segment_line))
        xy_points: list[tuple[float, float]] = [(in2mm(x), in2mm(y)) for x, y in xy_points_orig]

        if len(edges) == 1:
            edge_groups = [edges]
        else:
            # Group edges that constitute a single trace (there can be
            # multiple traces per route+layer). To do this we look for all
            # edges that join to each other.
            edge_groups: list[list[tuple[int, int, FLine]]] = []
            solo_edge_groups: list[list[tuple[int, int, FLine]]] = []
            if ALLOW_MULTIPLE_TRACES_PER_ROUTE:
                # find all end edges
                solo_edges: list[tuple[int, int, FLine]] = []
                end_edges: list[tuple[int, int, FLine]] = []
                inner_edges: list[tuple[int, int, FLine]] = []
                for edge_idx, edge in enumerate(edges):
                    n_matches = 0
                    matching_edges = ""
                    for other_edge in filter(lambda e: e != edge, edges):
                        if edge[0] in other_edge[:2] or edge[1] in other_edge[:2]:
                            n_matches += 1
                            matching_edges += str(other_edge) + "\n\t\t"
                    assert n_matches <= 2, f"Found more than 2 matching edges for !\n\tSource edge:\n\t\t{edge}\n\tMatching edges:\n\t\t{matching_edges}\nIs it possible that you have traces that fork?"
                    if n_matches == 0:
                        solo_edges.append(edge)
                    elif n_matches == 1:
                        end_edges.append(edge)
                    elif n_matches == 2:
                        inner_edges.append(edge)
                assert len(end_edges) % 2 == 0

                # all solo edges are, by definition, an edge group
                for edge in solo_edges:
                    solo_edge_groups.append([edge])
                solo_edges.clear()
                
                # get the edges for each group
                if len(end_edges) > 0:
                    while True:
                        end_a = end_edges.pop()
                        edge_group = [end_a]
                        end_b: tuple[int, int, FLine] = None

                        # find the group inner edges
                        group_xy_points = [end_a[0], end_a[1]]
                        found_inner_edge = True
                        while found_inner_edge:
                            found_inner_edge = False
                            for edge in copy.copy(inner_edges):
                                if edge[0] in group_xy_points or edge[1] in group_xy_points:
                                    edge_group.append(edge)
                                    group_xy_points.append(edge[0])
                                    group_xy_points.append(edge[1])
                                    inner_edges.remove(edge)
                                    found_inner_edge = True
                                    break

                        # find the other end
                        for edge in end_edges:
                            if edge[0] in group_xy_points or edge[1] in group_xy_points:
                                assert end_b is None
                                end_b = edge
                                group_xy_points.append(edge[0])
                                group_xy_points.append(edge[1])
                        end_edges.remove(end_b)
                        edge_group.append(end_b)

                        # add the group
                        assert len(set(group_xy_points)) == len(edge_group)+1
                        edge_groups.append(edge_group)

                        # check if we've found all the edge groups
                        if len(sum(edge_groups+solo_edge_groups, start=[])) >= len(edges):
                            break
                
                # sanity check
                assert len(sum(edge_groups+solo_edge_groups, start=[])) == len(edges)
                assert len(solo_edges) == 0
                assert len(end_edges) == 0
                assert len(inner_edges) == 0
            else:
                edge_groups.append(edges)

            # order the edges in the edge groups
            for edge_group_idx, edge_group in enumerate(copy.copy(edge_groups)):
                end_a, end_b = edge_group[0], edge_group[-1]
                group_inner_edges = edge_group[1:-1]
                new_edge_group = [end_a]
                while len(group_inner_edges) > 0:
                    for edge in copy.copy(group_inner_edges):
                        if edge[0] in new_edge_group[-1][:2] or edge[1] in new_edge_group[-1][:2]:
                            new_edge_group.append(edge)
                            group_inner_edges.remove(edge)
                            break
                new_edge_group.append(end_b)
                edge_groups[edge_group_idx] = new_edge_group
            
            # orient adjacent edges
            for edge_group in edge_groups:
                for edge_idx, edge in enumerate(edge_group[:-1]):
                    next_edge = edge_group[edge_idx]
                    if edge[0] in next_edge[:2]:
                        edge_group[edge_idx] = (edge[1], edge[0], *edge[2:])

                # orient the last edge
                end_b = edge_group[-1]
                prev_edge = edge_group[-2]
                if end_b[1] in prev_edge[:2]:
                    edge_group[-1] = (end_b[1], end_b[0], *end_b[2:])
            
            # include the solo edges
            edge_groups += solo_edge_groups

        # # debugging
        # for edge in edges:
        #     print(f"LINE   {xy_points_orig[edge[0]]}   {xy_points_orig[edge[1]]}")

        ret = []
        for edge_group in edge_groups:
            source_lines = [e[2] for e in edge_group]
            instance = cls(source_lines, layer_name, xy_points, edge_group, shape, bend_radius)
            ret.append( instance )

        return ret, pre_trace + post_trace

    def _get_segments_with_through_holes(self) -> list[LineSegment]:
        """
        Gets the segments of the trace, with the first and last segment
        adjusted to match the pins on either end.

        Returns
        -------
        list[LineSegment]
            Updated segments including adjusted through-holes.
        """
        segments = copy.copy(self.segments)
        pins = self._get_pins_ajusted()

        # Recalculate the starting and ending segments depending on
        # if we should be adding through-holes for the traces.
        if pins["a"] is not None:
            x, y = pins["a"].x_offset, pins["a"].y_offset
            segments[0] = LineSegment((x, y), segments[0].xy2)

        if pins["b"] is not None:
            x, y = pins["b"].x_offset, pins["b"].y_offset
            segments[-1] = LineSegment(segments[-1].xy1, (x, y))

        return segments

    def _get_pins_ajusted(self) -> dict[str, Pin|None]:
        """
        Adjusts end point pin locations to be offset from the component pins.

        The distance of the offset is based on through-hole diameter and via diameter.

        Returns
        -------
        dict[str, Union[Pin, None]]
            Adjusted pins for the trace ends.
        """
        ret = copy.copy(self.pins)

        # Recalculate the pin locations to be offset from the component through holes
        for end in ["a", "b"]:
            pin = self.pins[end]
            if pin is None:
                continue

            x, y = pin.x_offset, pin.y_offset
            y += Pin.through_hole_diameter() / 2 + Pin.via_diameter() / 2
            adjusted_pin = Pin(None, "", "", x, y, pin.layer, pin.is_pad)
            ret[end] = adjusted_pin

        return ret

    def to_vtk(self, trace_polydata: vtk.vtkPolyData, via_polydata: vtk.vtkPolyData):
        """
        Inserts the entire trace and vias into the polydata objects.

        Parameters
        ----------
        trace_polydata : vtk.vtkPolyData
            The VTK PolyData object for trace.
        via_polydata : vtk.vtkPolyData
            The VTK PolyData object for trace vias.
        """
        vtk_points: vtk.vtkPoints = trace_polydata.GetPoints()
        vtk_cells: vtk.vtkCellData = trace_polydata.GetCellData()

        # Recalculate the starting and ending segments depending on
        # if we should be adding through-holes for the traces.
        segments = self._get_segments_with_through_holes()

        # build the list of xy points
        xy_points = [s.xy1 for s in segments]
        xy_points.append(segments[-1].xy2)

        # build the segments verts and cells
        for segment_idx, segment in enumerate(segments):
            self.segment_to_vtk(trace_polydata, segment_idx, segment)

        # build the corners
        for xy_pnt in xy_points:
            segment = list(filter(lambda s: xy_pnt in [s.xy1, s.xy2], segments))[0]
            segment_idx = segments.index(segment)
            corner = self.get_trace_corner(xy_pnt, segment_idx, segment)
            if corner is not None:
                corner.to_vtk(trace_polydata)

        # Close the ends
        for xy_pnt in [xy_points[0], xy_points[-1]]:
            verticies = self._xypnt_vtk_verticies[xy_pnt]
            point_id_0 = verticies.vtk_indices[0]
            new_vtk_points = self.shape.add_vtk_cells(trace_polydata, point_id_0)
            self.inc_vtk_indicies(point_id_0 + len(verticies), len(new_vtk_points))
            for new_vtk_id in new_vtk_points:
                xyz_point = list(vtk_points.GetPoint(new_vtk_id))
                verticies.xyz_points = np.concat((verticies.xyz_points, np.array([xyz_point])), axis=0)
                verticies.vtk_indices.append(new_vtk_id)

        # add the through holes
        for end, pin in self._get_pins_ajusted().items():
            if pin is not None:
                pin.to_vtk(via_polydata)

    def draw(self, ax: maxis.Axis):
        """
        Draws the trace for 2d visual verification.

        Parameters
        ----------
        ax : maxis.Axis
            Axis to draw the trace on.
        """
        # Recalculate the starting and ending segments depending on
        # if we should be adding through-holes for the traces.
        segments = self._get_segments_with_through_holes()

        for end, pin in self._get_pins_ajusted().items():
            if pin is not None:
                pin.draw(ax)

        for seg in segments:
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

    test_trace = SingleTrace(trace_points, trace_edges, DEFAULT_PIPE_SHAPE(awg2mm(26)), allow_overlap=True)
    return test_trace


if __name__ == "__main__":
    from Trace.PipeShape import DEFAULT_PIPE_SHAPE
    from tool.units import *

    test_trace = _tst_trace_from_points()

    vtk_test_trace = test_trace.to_vtk(vt.new_polydata())
    print(f"{vtk_test_trace.GetNumberOfPoints()=}")
    vt.save_to_vtk(vtk_test_trace, "test_trace.vtk")

    test_trace_mesh = pyvista.PolyData(vtk_test_trace)
    pyvista.global_theme.allow_empty_mesh = True
    test_trace_mesh.plot(show_edges=True, opacity=1, show_vertices=True)
    # pyvista.PolyDataFilters.plot_normals(test_trace_mesh, mag=0.5, flip=False, faces=False, show_edges=True, opacity=0.95, show_verticies=True)