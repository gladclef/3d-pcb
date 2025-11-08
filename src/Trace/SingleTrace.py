import copy
import re

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
from Trace.TraceLine import TraceLine
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
                 route: str,
                 layer: str,
                 segments: list[TraceLine | LineSegment],
                 shape: PipeShape=None,
                 bend_radius: float=1,
                 allow_overlap=False,
                 inner_vias: dict[Point, Via]=None,
                 branch_vias: dict[Point, Via]=None):
        """
        Initializes a SingleTrace object.

        Parameters
        ----------
        source_lines : list[FLine]
            The lines that were parsed in order to create this instance.
        route : str
            The name of the route that this instance came from (eg GND).
        layer : str
            The name of the layer that this instance came from (eg TOP, BOTTOM).
        segments : list[TraceLine | LineSegment]
            List of line segments. Each segment can be defined as either
            a TraceLine or a LineSegment object.
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
        branch_vias : list, optional
            Vias automatically generated at branch points.
        """
        linenos = [l.lineno for l in source_lines]
        print(f"Creating {self.__class__.__name__} instance on layer {layer} from lines ({min(linenos)+1}-{max(linenos)+1})")# +
        #      ":\n\t" + "\n\t".join([l.v.rstrip() for l in source_lines]))
        super().__init__(source_lines, segments, shape)

        # set some defaults
        if bend_radius is None:
            bend_radius = g.TRACE_CORNER_RADIUS

        self.route = route
        """ The name of the route that this instance came from (eg GND). """
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
        self.pins: dict[int, Pin] = { 0: None, 1: None }
        """ The pins that indicate where to put vias at either end of this trace. """
        self.inner_vias: dict[Point, Via] = inner_vias or {}
        """ Vias defined as part of a route. """
        self.branch_vias: dict[Point, Via] = branch_vias or {}
        """ Vias automatically generated at branch points. """

        self.check_segment_duplicates(self.segments)
        if not allow_overlap:
            self.check_segments_overlap()

    def adjust_ends_for_conductive_filament(self):
        """
        For branching vias and inner vias:
            shortens the geometry of the ends of the trace so as to minimally overlap with the vias.
        For end point pins:
            extends the geometry of the ends of the trace so as to fully overlap with the pins.
        """
        # d = Via.via_diameter()

        # get the end points of this trace
        xy0 = self.segments[0].xy0
        xy1 = self.segments[-1].xy1

        # get the list of inner vias
        inner_vias = list(self.inner_vias.values())
        inner_vias += list(self.branch_vias.values())
        # for trace in all_traces:
        #     if trace != self:
        #         if trace.route == self.route:
        #             for inner_via_pnt, inner_via in trace.inner_vias.values():
        #                 if inner_via_pnt.almost_equal(xy0, d) or \
        #                    inner_via_pnt.almost_equal(xy1, d):
        #                     inner_vias[inner_via_pnt] = inner_via
        # if len(inner_vias) == 0:
        #     return
        
        # get the list of end point pins
        end_point_pins: list[Pin] = []
        for pin in self.pins.values():
            if pin.location not in [via.location for via in inner_vias]:
                end_point_pins.append(pin)
        
        # for each inner via, limit the length of the trace to the intersection of the via
        for is_inner, vias in [(True, inner_vias), (False, end_point_pins)]:
            for via in vias:
                if self.segments[0].xy0.distance(via.location) < self.segments[-1].xy1.distance(via.location):
                    is_via_at_xy0 = True
                    segments_out2in = self.segments
                else:
                    is_via_at_xy0 = False
                    segments_out2in = list(reversed(self.segments))

                # Circular via with radius "r", intersecting with trace with radius "o":
                #
                #                  _/|\         +
                #               _/   | \        |
                #       h    _/      |  ]       |
                #         _/         |   }      |
                #      _/            |   [      o
                #   _/               |    |     |
                # /   α              |    |     |
                # -------------------+----|     +
                # 
                # {-----------r-----------}
                # {---------l--------|-a--}
                #
                # given o, r
                # h = r,   l + a = r
                # l^2 + o^2 = h^2
                # 
                # l = sqrt(h^2 - o^2)
                o = self.shape.radius
                if isinstance(via, Via):
                    r = via.via_diameter() / 2
                else:
                    pin: Pin = via
                    r = pin.through_hole_diameter() / 2
                h = r
                l = np.sqrt(h**2 - o**2)

                # Don't adjust anything in the case that there is no segment that is r distance away.
                if sum([s.length for s in self.segments]) < r:
                    continue
                
                # find the segment at distance r from the center of the via
                segs_length = 0
                segments_to_remove: list[LineSegment] = []
                for s in segments_out2in:
                    segs_length += s.length
                    if segs_length >= r:
                        segment = s
                        break
                    else:
                        segments_to_remove.append(s)
                assert segment.xy0 in self.xy_points
                assert segment.xy1 in self.xy_points
                assert segment in self.segments
                
                # remove segments between the center of the via and the found segment
                if len(segments_to_remove) > 0:
                    for s in segments_to_remove:
                        if is_via_at_xy0:
                            self.remove_xypnt(s.xy0)
                        else:
                            self.remove_xypnt(s.xy1)
                assert segment.xy0 in self.xy_points
                assert segment.xy1 in self.xy_points
                assert segment in self.segments

                # order the segment, so that we know which end is on the inner side of the trace
                original_segment = segment
                if is_via_at_xy0:
                    segment_in2out = segment.reversed()
                    seg_inner_pnt = segment.xy1
                    seg_outer_pnt = segment.xy0
                else:
                    segment_in2out = segment
                    seg_inner_pnt = segment.xy0
                    seg_outer_pnt = segment.xy1
                assert seg_inner_pnt.distance(via.location) > seg_outer_pnt.distance(via.location)
                assert segment_in2out.xy0 == seg_inner_pnt
                assert segment_in2out.xy1 == seg_outer_pnt
                
                # adjust the segment
                if is_inner:

                    # change the segment so that it ends at distance l
                    seg_dist_to_via = seg_inner_pnt.distance(via.location)
                    l_pnt = segment_in2out.distance_along_line(seg_dist_to_via - l, seg_inner_pnt)
                    new_pnts = (l_pnt, None) if is_via_at_xy0 else (None, l_pnt)
                    new_segment = self.change_segment_points(segment, new_pnts[0], new_pnts[1])
                else: # end point pin

                    # change the segment so that it overlaps the entire pin
                    pin: Pin = via
                    seg_dist_to_via = seg_inner_pnt.distance(pin.location)
                    r2_pnt = segment_in2out.distance_along_line(seg_dist_to_via + r, seg_inner_pnt)
                    new_pnts = (r2_pnt, None) if is_via_at_xy0 else (None, r2_pnt)
                    new_segment = self.change_segment_points(original_segment, new_pnts[0], new_pnts[1])
                
                # remove any leftover trace corners
                if via.location in self.xypnt_trace_corners:
                    del self.xypnt_trace_corners[via.location]

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
        xy_locs = { 0: self.segments[0].xy0, 1: self.segments[-1].xy1 }

        # find the component through holes that most closely match this instance
        closest_pins: dict[int, Pin] = { 0: None, 1: None }
        for component in components:
            component = component.get_transformed()
            for pin in component.shape.pins:
                for end in [0, 1]:
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
        for end in [0, 1]:
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

    def change_segment_points(self, segment: LineSegment, new_xy0: Point = None, new_xy1: Point = None) -> LineSegment:
        new_segment = super().change_segment_points(segment, new_xy0, new_xy1)

        for old_pnt, new_pnt in [(segment.xy0, new_xy0), (segment.xy1, new_xy1)]:
            if new_pnt is None:
                continue

            if old_pnt in self.inner_vias:
                via = self.inner_vias[old_pnt]
                self.inner_vias[new_pnt] = Via(new_pnt, via.name, via.source_lines)
                del self.inner_vias[old_pnt]

            if old_pnt in self.branch_vias:
                via = self.branch_vias[old_pnt]
                self.branch_vias[new_pnt] = Via(new_pnt, via.name, via.source_lines)
                del self.branch_vias[old_pnt]

            for end in list(self.pins.keys()):
                pin = self.pins[end]
                self.pins[end] = Pin(pin.parent, pin.pin_description, pin.pad_name, new_pnt, pin.layer, pin.is_pad)
            
            if old_pnt in self.xypnt_trace_corners:
                del self.xypnt_trace_corners[old_pnt]
        
        print(f"change_segment_points() from\n\t[{segment.xy0},{segment.xy1}] to\n\t[{new_segment.xy0},{new_segment.xy1}]")
        return new_segment

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
            return cad_lines, [], [], None
        if len(routes) == 1:
            # just the end matcher "$ENDROUTES"
            assert routes[0].v.strip() == "$ENDROUTES"
            return pre_routes, [], post_routes, None
        if len(routes) == 2:
            if routes[0].v.strip() == "$ROUTES" and routes[1].v.strip() == "$ENDROUTES":
                return pre_routes, [], post_routes, None
        assert routes[-1].v.strip() in ["ROUTE", "$ENDROUTES"]

        # break up on route
        pre_route, route, post_route = route_helper.get_next_region(routes)
        if len(route) == 0:
            return cad_lines, [], [], None
        if len(route) == 1:
            # just the end matcher "ROUTE" or "$ENDROUTES"
            assert route[0].v.strip() in ["ROUTE", "$ENDROUTES"]
            # check if there are any more routes
            pre_trace, trace, post_trace, route_name = cls.get_lines_for_next_trace(pre_routes + pre_route + post_route + post_routes)
            if len(trace) != 0:
                return pre_trace, trace, post_trace, route_name
            else:
                return pre_routes + pre_route, [], post_route + post_routes, route_name
        assert route[0].v.strip().startswith("ROUTE")
        assert route[-1].v.strip().startswith("ROUTE") or route[-1].v.strip() == "$ENDROUTES"
        route_name = re.match(r"ROUTE[ \t]*(.*)", route[0].v.strip()).groups()[0]
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
            pre_trace, trace, post_trace, route_name = cls.get_lines_for_next_trace(pre_routes + pre_route + post_route + post_routes)
            if len(trace) != 0:
                return pre_trace, trace, post_trace, route_name
            else:
                return pre_routes + pre_route + pre_layer, [], post_layer + post_route + post_routes, route_name
        else:
            # multiple FLines matching this layer for this route
            assert route[0].v.startswith("ROUTE ")
            pre_layer.insert(0, route[0])
            if layer_helper.end_matches(layer[-1]):
                post_layer.insert(0, layer.pop())

        return pre_routes + pre_route + pre_layer, layer, post_layer + post_route + post_routes, route_name
    
    @classmethod
    def _parse_trace_lines(cls, trace: list[FLine]) -> tuple[list[TraceLine], dict[Point, Via]]:
        # parse the lines for this route+layer
        segment_lines = list(filter(lambda l: l.v.startswith("LINE "), trace))
        xy_points_orig: list[Point] = []
        edges: list[TraceLine] = []
        for segment_line in segment_lines:

            x0, y0, x1, y1 = tuple(map(float, segment_line.v.strip()[5:].split(" ")))
            point0, point1 = Point(x0, y0), Point(x1, y1)
            if point0 not in xy_points_orig:
                xy_points_orig.append(point0)
            if point1 not in xy_points_orig:
                xy_points_orig.append(point1)

            edges.append(TraceLine(segment_line, point0, point1))

        # parse the vias for this route
        vias, l = Via.from_cad_file(trace)
        inner_vias: dict[Point, Via] = {}
        while len(vias) > 0:
            for via in vias:
                inner_vias[via.location] = via
                print(f"Found inner via for single trace at line {via.source_lines[0].lineno}")
            vias, l = Via.from_cad_file(l)

        return edges, inner_vias

    @classmethod
    def _sort_edges_in_groups(cls, edge_groups: list[list[TraceLine]]):
        # order the edges in the edge groups
        for edge_group_idx, edge_group in enumerate(copy.copy(edge_groups)):
            if len(edge_group) == 1:
                assert edge_group[0].is_solo
                continue
            assert len(TraceLine.only_solos(edge_group)) == 0

            end_a, end_b = edge_group[0], edge_group[-1]
            group_inner_edges = edge_group[1:-1]
            new_edge_group = [end_a]
            while len(group_inner_edges) > 0:
                for edge in copy.copy(group_inner_edges):
                    if (edge.xy0 in new_edge_group[-1].xy_points) or (edge.xy1 in new_edge_group[-1].xy_points):
                        new_edge_group.append(edge)
                        group_inner_edges.remove(edge)
                        break
            new_edge_group.append(end_b)
            edge_groups[edge_group_idx] = new_edge_group

    @classmethod
    def _order_xy_points_in_edges(cls, edge_groups: list[list[TraceLine]]):
        # orient adjacent edges
        for edge_group in edge_groups:
            if len(edge_group) == 1:
                continue

            for edge_idx, edge in enumerate(edge_group[:-1]):
                next_edge = edge_group[edge_idx]
                if edge.xy0 in next_edge.xy_points:
                    edge.xy0, edge.xy1 = edge.xy1, edge.xy0

            # orient the last edge
            end_b = edge_group[-1]
            prev_edge = edge_group[-2]
            if end_b.xy1 in prev_edge.xy_points:
                end_b.xy0, end_b.xy1 = end_b.xy1, end_b.xy0

    @classmethod
    def _group_edges(cls, edges: list[TraceLine]) -> list[list[TraceLine]]:
        """
        Group edges that constitute a single trace (there can be multiple traces per route+layer).
        To do this we look for all edges that join to each other.
        """
        if len(edges) == 1:
            edges[0].is_solo = True
            return [edges]
        if not ALLOW_MULTIPLE_TRACES_PER_ROUTE:
            return [edges]
        
        edge_groups: list[list[TraceLine]] = []

        # find all end edges
        for edge_idx, edge in enumerate(edges):
            n_matches = 0
            matching_edges, m0, m1 = "", [], []
            edge.joined_ends = []
            for other_edge in filter(lambda e: e != edge, edges):
                if (edge.xy0 in other_edge.xy_points) or (edge.xy1 in other_edge.xy_points):
                    n_matches += 1
                    matching_edges += str(other_edge) + "\n\t\t"
                    if (edge.xy0 in other_edge.xy_points):
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
                        edge2line = lambda e: Line.from_two_points(e.xy0, e.xy1)
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
                    edge.is_branch = True
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

        assert len(TraceLine.uncategorized(edges)) == 0
        assert len(TraceLine.only_ends(edges)) % 2 == 0
        assert len(TraceLine.only_ends(edges) + TraceLine.only_solos(edges) + TraceLine.only_inners(edges)) == len(edges)

        # Reassign edge end points to use points that are _almost_ equal,
        # in case some ends almost match but don't quite.
        end_points: list[Point] = []
        for edge in TraceLine.only_ends(edges):
            assert len(edge.unjoined_points) == 1
            old_end_point = edge.unjoined_points[0]

            new_end_point = old_end_point
            for end_point in end_points:
                if end_point.almost_equal(old_end_point, delta=1e-3):
                    new_end_point = end_point

            if new_end_point == old_end_point:
                end_points.append(new_end_point)
            else:
                if 0 in edge.unjoined_ends:
                    edge.xy0 = new_end_point
                else:
                    edge.xy1 = new_end_point

        # all solo edges are, by definition, an edge group
        for edge in TraceLine.only_solos(edges):
            edge_groups.append([edge])

        # get the edges for each group
        end_edges = TraceLine.only_ends(edges)
        inner_edges = TraceLine.only_inners(edges)
        if len(end_edges) > 0:
            while True:
                end_a = end_edges.pop()
                edge_group = [end_a]
                end_b: TraceLine = None
                latest_xy = end_a.xy0 if 0 == end_a.joined_ends[0] else end_a.xy1

                # find the group inner edges
                group_xy = end_a.joined_points
                found_inner_edge = True
                while found_inner_edge:
                    found_inner_edge = False
                    for edge in copy.copy(inner_edges):
                        if latest_xy in edge.xy_points:
                            edge_group.append(edge)
                            latest_xy = edge.xy0 if edge.xy1 == latest_xy else edge.xy1
                            group_xy += edge.xy_points
                            inner_edges.remove(edge)
                            found_inner_edge = True
                            break

                # find the other end
                for edge in end_edges:
                    if latest_xy in edge.xy_points:
                        assert end_b is None
                        end_b = edge
                        group_xy += edge.xy_points
                end_edges.remove(end_b)
                edge_group.append(end_b)

                # add the group
                assert len(set(group_xy)) == len(edge_group)
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
    def _group_and_sort_edges(cls, edges: list[TraceLine]) -> list[list[TraceLine]]:
        """
        Group edges that constitute a single trace (there can be multiple traces per route+layer).
        """
        edge_groups = cls._group_edges(edges)
        cls._sort_edges_in_groups(edge_groups)
        cls._order_xy_points_in_edges(edge_groups)
        return edge_groups

    @classmethod
    def _vias_at_branches(cls, edge_groups: list[list[TraceLine]]) -> dict[Point, Via]:
        ret: dict[Point, Via] = {}

        for edge_group in edge_groups:
            for edge in edge_group:
                if edge.is_branch:
                    if 0 in edge.joined_ends:
                        pnt = edge.xy1
                    else:
                        assert 1 in edge.joined_ends
                        pnt = edge.xy0
                    ret[pnt] = Via(pnt, "Branch", [edge.fline])

        return ret

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
        pre_trace, trace, post_trace, route_name = cls.get_lines_for_next_trace(cad_lines)
        if len(trace) == 0:
            return [], pre_trace + post_trace

        # get the layer name
        name_lines = list(filter(lambda l: l.v.strip().startswith("LAYER "), trace))
        assert len(name_lines) <= 1
        layer_name = name_lines[0].v.split("LAYER ")[1].strip() if len(name_lines) > 0 else "TOP"

        # parse the lines for this route+layer
        edges, inner_vias = cls._parse_trace_lines(trace)

        # group edges into traces
        edge_groups = cls._group_and_sort_edges(edges)

        # add vias at branch points
        branch_vias = cls._vias_at_branches(edge_groups)

        ret: list[SingleTrace] = []
        for edge_group in edge_groups:
            source_lines = [e.fline for e in edge_group]
            instance = cls(source_lines, route_name, layer_name, edge_group,
                           shape, bend_radius, inner_vias=inner_vias, branch_vias=branch_vias)
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
        if pins[0] is not None:
            segments[0] = LineSegment(pins[0].location, segments[0].xy1)

        if pins[1] is not None:
            segments[-1] = LineSegment(segments[-1].xy0, pins[1].location)

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
        for end in [0, 1]:
            pin = self.pins[end]
            if pin is None:
                continue

            dist = Pin.through_hole_diameter() / 2 + Via.via_diameter() / 2

            def get_segment_at_dist(segments: list[LineSegment]):
                seg_lengths = 0
                for segment in segments:
                    seg_lengths += segment.length
                    if seg_lengths >= dist:
                        break
                return segment

            if end == 0:
                segment = get_segment_at_dist(self.segments)
                assert pin.location.distance(segment.xy0) < Pin.through_hole_diameter() + Via.via_diameter()
                adjusted_loc = segment.distance_along_line(dist, pin.location)
            else:
                segment = get_segment_at_dist(reversed(self.segments))
                assert pin.location.distance(segment.xy1) < Pin.through_hole_diameter() + Via.via_diameter()
                adjusted_loc = segment.reversed().distance_along_line(dist, pin.location)

            adjusted_pin = Pin(None, "", "", adjusted_loc, pin.layer, pin.is_pad)
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
        
        # add the vias
        vias = list(self.inner_vias.values()) + list(self.branch_vias.values())
        for via in vias:
            via.to_vtk(via_polydata)

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

        for via in self.inner_vias.values():
            via.draw(ax)

        for via in self.branch_vias.values():
            via.draw(ax)


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