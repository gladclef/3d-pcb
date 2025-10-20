import copy
import dataclasses
import re
from typing import Union

import matplotlib.axis as maxis
import numpy as np
import pyvista
import vtk

from Component.Component import Component
from Component.Pin import Pin
from Component.Via import Via
from FileIO.CadFileHelper import CadFileHelper
from FileIO.Line import Line as FLine
from Geometry.Point import Point
import Geometry.geometry_tools as geo
from Geometry.Line import Line
from Geometry.LineSegment import LineSegment
from Geometry.Path import Path
from Trace.AbstractTrace import AbstractTrace
from Trace.PipeShape import PipeShape, DEFAULT_PIPE_SHAPE
from Trace.TraceCorner import TraceCorner
from Trace.VtkPointGroup import VtkPointGroup
from tool.globals import board_parameters as g
import tool.vtk_tools as vt
from tool.units import *


ALLOW_MULTIPLE_TRACES_PER_ROUTE = True


@dataclasses.dataclass
class _TraceLine:
    fline: FLine
    xy0_idx: int
    xy1_idx: int
    is_end: bool = False
    """ True if this edge only connects to one other edge.
    Mutually exclusive with is_solo and is_inner_edge. """
    is_solo: bool = False
    """ True if this edge doesn't connect to any other edges.
    Mutually exclusive with is_end and is_inner_edge. """
    is_inner_edge: bool = False
    """ True if this edge connects to another edge on both ends.
    Mutually exclusive with is_end and is_solo. """
    joined_ends: list[int] = dataclasses.field(default_factory=list)
    """ If this contains a 0, then self.xy0_idx connects to another edge.
    If this contains a 1, then self.xy1_idx connects to another edge. """

    @property
    def xy_idxs(self) -> tuple[int, int]:
        return (self.xy0_idx, self.xy1_idx)

    @property
    def joined_idxs(self) -> list[int]:
        return [self.xy_idxs[je] for je in self.joined_ends]

    def joined_points(self, xy_points: list[Point]) -> list[Point]:
        """ Return the points that correspond to any indexes in self.joined_ends """
        return self._points(xy_points, self.joined_idxs)

    @property
    def unjoined_ends(self) -> list[int]:
        """ If this contains a 0, then self.xy0_idx does not connect to another edge.
        If this contains a 1, then self.xy1_idx does not connect to another edge. """
        return list(filter(lambda idx: idx not in self.joined_ends, [0, 1]))

    @property
    def unjoined_idxs(self) -> list[int]:
        return [self.xy_idxs[uje] for uje in self.unjoined_ends]

    def unjoined_points(self, xy_points: list[Point]) -> list[Point]:
        """ Return the points that correspond to any indexes in self.unjoined_ends """
        return self._points(xy_points, self.unjoined_idxs)

    def _points(self, xy_points: list[Point], idxs: list[int]) -> list[Point]:
        return [xy_points[idx] for idx in idxs]

    def points(self, xy_points: list[Point]) -> list[Point]:
        """ Return the points that correspond to self.xy0_idx and self.xy1_idx """
        return self._points(xy_points, self.xy_idxs)

    @classmethod
    def uncategorized(cls, trace_lines: list["_TraceLine"]) -> list["_TraceLine"]:
        return list(filter(lambda tl: not tl.is_end and not tl.is_solo and not tl.is_inner_edge, trace_lines))

    @classmethod
    def only_ends(cls, trace_lines: list["_TraceLine"]) -> list["_TraceLine"]:
        return list(filter(lambda tl: tl.is_end, trace_lines))

    @classmethod
    def only_solos(cls, trace_lines: list["_TraceLine"]) -> list["_TraceLine"]:
        return list(filter(lambda tl: tl.is_solo, trace_lines))

    @classmethod
    def only_inners(cls, trace_lines: list["_TraceLine"]) -> list["_TraceLine"]:
        return list(filter(lambda tl: tl.is_inner_edge, trace_lines))

    @classmethod
    def not_solo(cls, trace_lines: list["_TraceLine"]) -> list["_TraceLine"]:
        return list(filter(lambda tl: not tl.is_solo, trace_lines))

    def __repr__(self):
        return f"<{self.xy0_idx}, {self.xy1_idx}, {self.fline}>"


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
                 xy_points: list[Point],
                 segments: list["_TraceLine"] | list[LineSegment],
                 shape: PipeShape=None,
                 bend_radius: float=1,
                 allow_overlap=False,
                 inner_vias: list[Via]=None):
        """
        Initializes a SingleTrace object.

        Parameters
        ----------
        source_lines : list[FLine]
            The lines that were parsed in order to create this instance.
        layer : str
            The name of the layer that this instance came from (eg TOP, BOTTOM).
        xy_points : list[Point]
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
        inner_vias : list, optional
            Vias defined as part of a route.
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
        self._xypnt_vtk_verticies: dict[Point, VtkPointGroup] = {}
        """ Dictionary from xy point index to vtk points. """
        self.xypnt_trace_corners: dict[Point, TraceCorner] = {}
        """ Dictionary from xy point index to trace corners. """
        self.bend_radius = bend_radius
        """
        Radius of the allowed trace bends. Trace segment intersections will
        be rounded out using TraceCorners to match this bend radius.
        """
        self.pins: dict[str, Pin] = { "a": None, "b": None }
        """ The pins that indicate where to put vias at either end of this trace. """
        self.inner_vias: dict[Point, Via] = {}
        """ Vias defined as part of a route. """

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
        segment_tuples = [(s.xy0, s.xy1) for s in segments]
        segments_set = set(segment_tuples)
        if len(segment_tuples) != len(segments_set):
            linenos = [l.lineno for l in self.source_lines]
            min_lineno, max_lineno = min(linenos), max(linenos)
            raise ValueError(f"Error in SingleTrace.check_segment_duplicates(): there are {len(segments_set)} non-duplicate segments out of {len(segments)} total segments for lines {min_lineno}-{max_lineno}.")

    def check_segments_overlap(self):
        """
        Checks for overlapping segments within the trace.
        """
        for s1idx, segment1 in enumerate(self.segments):
            for s2idx, segment2 in enumerate(self.segments):
                if s1idx == s2idx:
                    # don't check for self-segment intersections
                    continue

                if segment1.xy0.almost_equal(segment2.xy0) < 1e-4 or \
                    segment1.xy0.almost_equal(segment2.xy1) < 1e-4 or \
                    segment1.xy1.almost_equal(segment2.xy0) < 1e-4 or \
                    segment1.xy1.almost_equal(segment2.xy1) < 1e-4:
                    # don't check for intersections with segments that share one of the end points
                    continue

                # check if these segments are parallel
                angle_diff = abs(segment1.angle - segment2.angle)
                if angle_diff < geo.ZERO_THRESH or 2*np.pi - angle_diff < geo.ZERO_THRESH:
                    continue

                intersection = segment1.intersection(segment2)
                if intersection is not None:
                    linenos = [l.lineno for l in self.source_lines]
                    raise ValueError(f"Error in SingleTrace.check_segments_overlap(): segments {s1idx} ({segment1}) and {s2idx} ({segment2}) overlap at [{intersection}].\n"
                                     + f"\tLines: {min(linenos)}-{max(linenos)}")

    def get_trace_corner(self, xy_pnt: Point, segment_idx: int, segment: LineSegment) -> TraceCorner | None:
        """
        Get the trace corner for the given xy point.
        Builds the trace corner as necessary.

        Note: returns None for the first and last xy points,
        since there is no corner at the ends of the trace.

        Parameters
        ----------
        xy_pnt : Point
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
        if segment_idx == 0 and xy_pnt == segment.xy0:
            return None
        if segment_idx == len(self.segments)-1 and xy_pnt == segment.xy1:
            return None

        if xy_pnt not in self.xypnt_trace_corners:
            if xy_pnt.almost_equal(segment.xy1):
                segment_a, segment_b = segment, self.segments[segment_idx+1]
            elif xy_pnt.almost_equal(segment.xy0):
                segment_a, segment_b = self.segments[segment_idx-1], segment
            else:
                raise RuntimeError("Error in SingleTrace.get_trace_corner(): " + "xy_pnt != segment.xy1 and xy_pnt != segment.xy2")

            self.xypnt_trace_corners[xy_pnt] = TraceCorner(self, (segment_a, segment_b), self.bend_radius)

        return self.xypnt_trace_corners[xy_pnt]

    def get_xypnt_vtk_verticies(self, xy_pnt: Point, segment_idx: int, segment: LineSegment) -> VtkPointGroup:
        """
        Get the vertices for the given end-point of a trace segment.

        Parameters
        ----------
        xy_pnt : Point
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
                if segment.xy0.almost_equal(xy_pnt):
                    return corner.get_vtk_group(corner.n_points-1)
                elif segment.xy1.almost_equal(xy_pnt):
                    return corner.get_vtk_group(0)
                else:
                    raise RuntimeError("Error in SingleTrace.get_xypnt_vtk_verticies(): " + f"expected the xy_point to be at the beginning or end of the given segment, " + f"but {xy_pnt=} and {segment.xy0=} and {segment.xy1=}!")

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
        xy_locs = { "a": self.segments[0].xy0, "b": self.segments[-1].xy1 }

        # find the component through holes that most closely match this instance
        closest_pins: dict[str, Pin] = { "a": None, "b": None }
        for component in components:
            component = component.get_transformed()
            for pin in component.shape.pins:
                for end in ["a", "b"]:
                    closest_pin, xy_loc = closest_pins[end], xy_locs[end]

                    if pin.location.distance(xy_loc) < Pin.through_hole_diameter():
                        if closest_pin is None:
                            closest_pin = pin
                        else:
                            dist_pin = pin.location.distance(xy_loc)
                            dist_closest = closest_pin.location.distance(xy_loc)
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
        xy_a, xy_b = segment.xy0, segment.xy1

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

    def remove_xypnt(self, xy_pnt: Point) -> tuple[tuple[list[LineSegment], TraceCorner], list[LineSegment]]:
        """ Removes the given xy_pnt.

        Parameters
        ----------
        xy_pnt : Point
            The point to be removed.

        Returns
        -------
        old_segments_and_corner: tuple[tuple[list[LineSegment], TraceCorner]
            The old segments and trace corner that were removed.
        new_segments: list[LineSegment]]
            The new segments that were added.
        """
        # deal with the segments
        old_segments, new_segments = Path.remove_xypnt(self, xy_pnt)

        # remove the old trace corner, if any
        trace_corner = None
        if xy_pnt in self.xypnt_trace_corners:
            trace_corner = self.xypnt_trace_corners[xy_pnt]
            del self.xypnt_trace_corners[xy_pnt]
        
        return (old_segments, trace_corner), new_segments

    def cleanup(self):
        """ Fixes and/or detects problems with the board that would make it difficult to boolean. """
        def cleanup_short_traces():
            while len(self.segments) >= 2:
                found_short_trace = False

                for segment in self.segments:
                    if segment.length < g.TRACE_CORNER_RADIUS:
                        linenos = [fl.lineno for fl in self.source_lines]
                        print(f"Found short trace segment in lines {min(linenos)}-{max(linenos)}")

                        # replace this xy_pair with the intersection of the previous and next lines
                        prev_segments = list(filter(lambda s: s != segment, self.segments_at_xypnt(segment.xy0)))
                        next_segments = list(filter(lambda s: s != segment, self.segments_at_xypnt(segment.xy1)))
                        if len(prev_segments) == 0:
                            # no previous line, replace with the next line
                            (old_segments, old_trace_corner), new_segments = self.remove_xypnt(segment.xy1)
                        elif len(next_segments) == 0:
                            # no next line, replace with the previous line
                            (old_segments, old_trace_corner), new_segments = self.remove_xypnt(segment.xy0)
                        else:
                            prev_segment = prev_segments[0]
                            new_segments: list[LineSegment] = []
                            for i, next_segment in enumerate(next_segments):
                                # get the new point to insert
                                prev_line = Line.from_two_points(prev_segment.xy0, prev_segment.xy1)
                                next_line = Line.from_two_points(next_segment.xy0, next_segment.xy1)
                                new_xy_pnt = prev_line.intersection(next_line)
                                
                                # remove this segment
                                if i == 0:
                                    (old_segments, old_trace_corner), new_segments = self.remove_xypnt(segment.xy0)

                                # TODO fix this
                                # # get the segment to be replaced
                                # new_segment = list(filter(lambda s: s.xy2 == segment.xy2, new_segments))
                                # assert len(new_segment) > 0
                                # new_segment = new_segment[0]

                                # # update the segment with an in-between point
                                # self.insert_xypnt(segment.xy1, new_segment)
                        
                        # force removal of stale trace corners
                        for segment in old_segments:
                            for pnt in segment.xy_points:
                                if pnt in self.xypnt_trace_corners:
                                    del self.xypnt_trace_corners[pnt]
                        for segment in new_segments:
                            for pnt in segment.xy_points:
                                if pnt in self.xypnt_trace_corners:
                                    del self.xypnt_trace_corners[pnt]

                        found_short_trace = True
                        break
                
                if not found_short_trace:
                    break
        
        cleanup_short_traces()

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
        if len(layer) == 1 and layer[0].v.strip().startswith("LINE "):
            # single-value layer
            pass
        elif len(layer) <= 1:
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
            # multiple FLines matching this layer for this route
            assert route[0].v.startswith("ROUTE ")
            pre_layer.insert(0, route[0])
            if layer_helper.end_matches(layer[-1]):
                post_layer.insert(0, layer.pop())

        return pre_routes + pre_route + pre_layer, layer, post_layer + post_route + post_routes
    
    @classmethod
    def _parse_trace_lines(cls, trace: list[FLine]) -> tuple[list[Point], list[_TraceLine], list[Via]]:
        # parse the lines for this route+layer
        segment_lines = list(filter(lambda l: l.v.startswith("LINE "), trace))
        xy_points_orig: list[Point] = []
        edges: list[_TraceLine] = []
        for segment_line in segment_lines:

            x0, y0, x1, y1 = tuple(map(float, segment_line.v.strip()[5:].split(" ")))
            point0, point1 = Point(x0, y0), Point(x1, y1)
            if point0 not in xy_points_orig:
                xy_points_orig.append(point0)
            if point1 not in xy_points_orig:
                xy_points_orig.append(point1)

            point0_idx, point1_idx = xy_points_orig.index(point0), xy_points_orig.index(point1)
            edges.append(_TraceLine(segment_line, point0_idx, point1_idx))
        xy_points: list[Point] = [Point(in2mm(p.x), in2mm(p.y)) for p in xy_points_orig]

        # parse the vias for this route
        vias, l = Via.from_cad_file(trace)
        inner_vias: dict[Point, Via] = {}
        while len(vias) > 0:
            for via in vias:
                inner_vias[via.location] = via
                print(f"Found inner via for single trace at line {via.source_lines[0].lineno}")
            vias, l = Via.from_cad_file(l)

        return xy_points, edges, inner_vias

    @classmethod
    def _sort_edges_in_groups(cls, edge_groups: list[list[_TraceLine]]):
        # order the edges in the edge groups
        for edge_group_idx, edge_group in enumerate(copy.copy(edge_groups)):
            if len(edge_group) == 1:
                assert edge_group[0].is_solo
                continue
            assert len(_TraceLine.only_solos(edge_group)) == 0

            end_a, end_b = edge_group[0], edge_group[-1]
            group_inner_edges = edge_group[1:-1]
            new_edge_group = [end_a]
            while len(group_inner_edges) > 0:
                for edge in copy.copy(group_inner_edges):
                    if (edge.xy0_idx in new_edge_group[-1].xy_idxs) or (edge.xy1_idx in new_edge_group[-1].xy_idxs):
                        new_edge_group.append(edge)
                        group_inner_edges.remove(edge)
                        break
            new_edge_group.append(end_b)
            edge_groups[edge_group_idx] = new_edge_group

    @classmethod
    def _order_xy_points_in_edges(cls, edge_groups: list[list[_TraceLine]]):
        # orient adjacent edges
        for edge_group in edge_groups:
            if len(edge_group) == 1:
                continue

            for edge_idx, edge in enumerate(edge_group[:-1]):
                next_edge = edge_group[edge_idx]
                if edge.xy0_idx in next_edge.xy_idxs:
                    edge.xy0_idx, edge.xy1_idx = edge.xy1_idx, edge.xy0_idx

            # orient the last edge
            end_b = edge_group[-1]
            prev_edge = edge_group[-2]
            if end_b.xy1_idx in prev_edge.xy_idxs:
                end_b.xy0_idx, end_b.xy1_idx = end_b.xy1_idx, end_b.xy0_idx

    @classmethod
    def _group_edges(cls, xy_points: list[Point], edges: list[_TraceLine]) -> list[list[_TraceLine]]:
        """
        Group edges that constitute a single trace (there can be multiple traces per route+layer).
        To do this we look for all edges that join to each other.
        """
        if len(edges) == 1:
            return [edges]
        if not ALLOW_MULTIPLE_TRACES_PER_ROUTE:
            return [edges]
        
        edge_groups: list[list[_TraceLine]] = []

        # find all end edges
        for edge_idx, edge in enumerate(edges):
            n_matches = 0
            matching_edges, m0, m1 = "", [], []
            edge.joined_ends = []
            for other_edge in filter(lambda e: e != edge, edges):
                if (edge.xy0_idx in other_edge.xy_idxs) or (edge.xy1_idx in other_edge.xy_idxs):
                    n_matches += 1
                    matching_edges += str(other_edge) + "\n\t\t"
                    if (edge.xy0_idx in other_edge.xy_idxs):
                        m0.append(other_edge)
                    else:
                        m1.append(other_edge)

            # check for and handle branching traces
            if len(m0) > 1 or len(m1) > 1:
                n_ends, n_inner = 0, 0

                # for each end, count the number of branching traces off of that end
                for end_idx, m, other_m in [(0, m0, m1), (1, m1, m0)]:
                    if len(m) == 0:
                        n_ends += 1
                    elif len(m) == 1:
                        n_inner += 1
                        edge.joined_ends.append(end_idx)
                    elif len(m) == 2:
                        # two branching traces, choose one of them to be the end
                        edge2line = lambda e: Line.from_two_points(xy_points[e.xy0_idx], xy_points[e.xy1_idx])
                        l0, l1, l2 = edge2line(edge), edge2line(m[0]), edge2line(m[1])
                        l01a, l02a, l12a = l0.angle_between(l1), l0.angle_between(l2), l1.angle_between(l2)
                        zero2pi = lambda a: np.pi if a < 1e-6 else a
                        l01a, l02a, l12a = zero2pi(l01a), zero2pi(l02a), zero2pi(l12a)

                        if l01a == max([l01a, l02a, l12a]):
                            # don't count edge 3 as one of the joined edges
                            m.remove(m[1])
                            n_inner += 1
                            edge.joined_ends.append(end_idx)
                        elif l02a == max([l01a, l02a, l12a]):
                            # don't count edge 2 as one of the joined edges
                            m.remove(m[0])
                            n_inner += 1
                            edge.joined_ends.append(end_idx)
                        else:
                            # don't count this edge as one of the joined edges
                            n_ends += 1
                    else: # len(m) > 2:
                        raise RuntimeError(
                            "Branching traces with more than three intersections at a single location aren't currently handled."
                            + f"\n\tSource edge:\n\t\t{edge}\n\tMatching edges:\n\t\t{matching_edges}"
                        )

                # assign the type to the edge
                if n_ends == 2:
                    assert n_inner == 0
                    assert len(edge.joined_ends) == 0
                    edge.is_solo = True
                elif n_ends == 1:
                    assert n_inner == 1
                    assert len(edge.joined_ends) == 1
                    edge.is_end = True
                else:
                    assert n_inner == 2
                    assert len(edge.joined_ends) == 2
                    edge.is_inner_edge = True

            elif n_matches == 0:
                edge.is_solo = True
            elif n_matches == 1:
                edge.is_end = True
                edge.joined_ends.append(0 if len(m0) == 1 else 1)
            elif n_matches == 2:
                assert len(m0) == 1 and len(m1) == 1
                edge.is_inner_edge = True
                edge.joined_ends.append(0)
                edge.joined_ends.append(1)

        assert len(_TraceLine.uncategorized(edges)) == 0
        assert len(_TraceLine.only_ends(edges)) % 2 == 0
        assert len(_TraceLine.only_ends(edges) + _TraceLine.only_solos(edges) + _TraceLine.only_inners(edges)) == len(edges)

        # Reassign edge end points to use points that are _almost_ equal,
        # in case some ends almost match but don't quite.
        end_points: list[Point] = []
        for edge in _TraceLine.only_ends(edges):
            assert len(edge.unjoined_points(xy_points)) == 1
            old_end_point = edge.unjoined_points(xy_points)[0]

            new_end_point = old_end_point
            for end_point in end_points:
                if end_point.almost_equal(old_end_point, delta=1e-3):
                    new_end_point = end_point

            if new_end_point == old_end_point:
                end_points.append(new_end_point)
            else:
                if 0 in edge.unjoined_ends:
                    edge.xy0_idx = xy_points.index(new_end_point)
                else:
                    edge.xy1_idx = xy_points.index(new_end_point)

        # all solo edges are, by definition, an edge group
        for edge in _TraceLine.only_solos(edges):
            edge_groups.append([edge])

        # get the edges for each group
        end_edges = _TraceLine.only_ends(edges)
        inner_edges = _TraceLine.only_inners(edges)
        if len(end_edges) > 0:
            while True:
                end_a = end_edges.pop()
                edge_group = [end_a]
                end_b: _TraceLine = None
                latest_xy_idx = end_a.xy0_idx if 0 == end_a.joined_ends[0] else end_a.xy1_idx

                # find the group inner edges
                group_xy_idxs = end_a.joined_idxs
                found_inner_edge = True
                while found_inner_edge:
                    found_inner_edge = False
                    for edge in copy.copy(inner_edges):
                        if latest_xy_idx in edge.xy_idxs:
                            edge_group.append(edge)
                            latest_xy_idx = edge.xy0_idx if edge.xy1_idx == latest_xy_idx else edge.xy1_idx
                            group_xy_idxs += edge.xy_idxs
                            inner_edges.remove(edge)
                            found_inner_edge = True
                            break

                # find the other end
                for edge in end_edges:
                    if latest_xy_idx in edge.xy_idxs:
                        assert end_b is None
                        end_b = edge
                        group_xy_idxs += edge.xy_idxs
                end_edges.remove(end_b)
                edge_group.append(end_b)

                # add the group
                assert len(set(group_xy_idxs)) == len(edge_group)
                edge_groups.append(edge_group)

                # check if we've found all the edge groups
                if len(sum(edge_groups, start=[])) >= len(edges):
                    break

        # sanity check
        assert len(sum(edge_groups, start=[])) == len(edges)
        assert len(end_edges) == 0
        assert len(inner_edges) == 0

        # # debugging
        # for edge in edges:
        #     print(f"LINE   {xy_points_orig[edge[0]]}   {xy_points_orig[edge[1]]}")

        return edge_groups

    @classmethod
    def _group_and_sort_edges(cls, xy_points: list[Point], edges: list[_TraceLine]) -> list[list[_TraceLine]]:
        """
        Group edges that constitute a single trace (there can be multiple traces per route+layer).
        """
        edge_groups = cls._group_edges(xy_points, edges)
        cls._sort_edges_in_groups(edge_groups)
        cls._order_xy_points_in_edges(edge_groups)
        return edge_groups

    @classmethod
    def from_cad_file(cls, cad_lines: list[FLine], shape: PipeShape=None, bend_radius: float=None) -> tuple[list["SingleTrace"], list[FLine]]:
        """
        Creates zero or more SingleTrace object from CAD file lines.

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
            A tuple containing zero or more SingleTrace instances and the remaining lines from the CAD file.
        """
        # get the lines
        pre_trace, trace, post_trace = cls.get_lines_for_next_trace(cad_lines)
        if len(trace) == 0:
            return [], pre_trace + post_trace

        # get the layer name
        name_lines = list(filter(lambda l: l.v.strip().startswith("LAYER "), trace))
        assert len(name_lines) <= 1
        layer_name = name_lines[0].v.split("LAYER ")[1].strip() if len(name_lines) > 0 else "TOP"

        # parse the lines for this route+layer
        xy_points, edges, inner_vias = cls._parse_trace_lines(trace)

        # group edges into traces
        edge_groups = cls._group_and_sort_edges(xy_points, edges)
        ret: list[SingleTrace] = []

        for edge_group in edge_groups:
            source_lines = [e.fline for e in edge_group]
            instance = cls(source_lines, layer_name, xy_points, edge_group, shape, bend_radius, inner_vias=inner_vias)
            ret.append(instance)

        return ret, pre_trace + post_trace
    
    def _get_segments_with_through_holes(self) -> tuple[LineSegment]:
        """
        Gets the segments of the trace, with the first and last segment
        adjusted to match the pins on either end.

        Returns
        -------
        list[LineSegment]
            Updated segments including adjusted through-holes.
        """
        segments = list(self.segments)
        pins = self._get_pins_ajusted()

        # Recalculate the starting and ending segments depending on
        # if we should be adding through-holes for the traces.
        if pins["a"] is not None:
            segments[0] = LineSegment(pins["a"].location, segments[0].xy1)

        if pins["b"] is not None:
            segments[-1] = LineSegment(segments[-1].xy0, pins["b"].location)

        return tuple(segments)

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

            x, y = pin.location.x, pin.location.y
            y += Pin.through_hole_diameter() / 2 + Pin.via_diameter() / 2
            adjusted_pin = Pin(None, "", "", Point(x, y), pin.layer, pin.is_pad)
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
        xy_points = [s.xy0 for s in segments]
        xy_points.append(segments[-1].xy1)

        # build the segments verts and cells
        for segment_idx, segment in enumerate(segments):
            self.segment_to_vtk(trace_polydata, segment_idx, segment)

        # build the corners
        for xy_pnt in xy_points:
            segment = list(filter(lambda s: xy_pnt in [s.xy0, s.xy1], segments))[0]
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
            ax.arrow(seg.x0, seg.y0, seg.x1-seg.x0, seg.y1-seg.y0, color="teal", head_width=.3)


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

    test_trace = SingleTrace(trace_points, trace_edges, DEFAULT_PIPE_SHAPE(g.DEFAULT_WIRE_DIAMETER), allow_overlap=True)
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