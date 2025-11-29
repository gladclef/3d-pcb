from abc import ABC, abstractmethod
import copy
import re
from typing import Literal

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


class SingleTrace(AbstractTrace, ABC):
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
        self._enable_stcheck_structure_is_valid = False
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
        self.pins: dict[Literal[0,1], Pin] = {}
        """ The pins that indicate where to put vias at either end of this trace. """
        self.inner_vias: dict[Point, Via] = inner_vias or {}
        """ Vias defined as part of a route. """
        self.branch_vias: dict[Point, Via] = branch_vias or {}
        """ Vias automatically generated at branch points. """

        self.check_segment_duplicates(self.segments)
        if not allow_overlap:
            self.check_segments_overlap()

        self._enable_stcheck_structure_is_valid = True
        self._stcheck_structure_is_valid()
   
    def _stcheck_structure_is_valid(self):
        if not self._enable_stcheck_structure_is_valid:
            return
        super()._check_structure_is_valid()
        assert (0 not in self.pins) or (self.pins[0].location == self.segments[0].xy0), "Pin0 location does not match segment0 location!"
        assert (1 not in self.pins) or (self.pins[1].location == self.segments[-1].xy1), "Pin1 location does not match segmentN location!"
        assert all([pnt in self.xy_points for pnt in self._xypnt_vtk_verticies.keys()]), "VTK point not in self.xy_points!"
        assert all([pnt in self.xy_points for pnt in self.xypnt_trace_corners.keys()]), "Trace corner not in self.xy_points!"
        assert all([pnt in self.xy_points for pnt in self.inner_vias.keys()]), "Inner via not in self.xy_points!"
        assert all([pnt in self.xy_points for pnt in self.branch_vias.keys()]), "Branch via not in self.xy_points!"
        assert all([pnt == self.inner_vias[pnt].location for pnt in self.inner_vias.keys()]), "Inner via has a misconfigured location!"
        assert all([pnt == self.branch_vias[pnt].location for pnt in self.branch_vias.keys()]), "Branch via has a misconfigured location!"
        assert all([isinstance(s, TraceLine) for s in self.segments])

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
        # sanity check
        assert len(self.pins) == 0

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
        
        # adjust the segments so that they end at the new component pad locations
        self._enable_stcheck_structure_is_valid = False
        if 0 in self.pins:
            (s0, s1), (ps, ns) = self.change_segment_points(self.segments[0], self.pins[0].location, None, force=True)
            assert s0 and s1
        if 1 in self.pins:
            (s0, s1), (ps, ns) = self.change_segment_points(self.segments[-1], None, self.pins[1].location, force=True)
            assert s0 and s1
        self._enable_stcheck_structure_is_valid = True
        
        self._stcheck_structure_is_valid()

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

    def _unset_xypnt(self, xy_pnt: Point):
        if xy_pnt in self.xypnt_trace_corners:
            del self.xypnt_trace_corners[xy_pnt]

    def remove_xypnt(self, xy_pnt: Point) -> tuple[bool, tuple[list[LineSegment], list[LineSegment]]]:
        # Don't remove this point if it's the same as a pin location
        if (0 in self.pins) and (self.pins[0].location == xy_pnt):
            return False, ([], [])

        # Remove trace corners
        self._unset_xypnt(xy_pnt)
        
        # Deal with moving inner vias and branch vias
        adjoining_pnt = None
        def get_adjoining_point():
            segment = list(filter(lambda s: s is not None, self.segments_at_xypnt(xy_pnt)))[0]
            if segment.xy0 == xy_pnt:
                return segment.xy1
            else:
                assert segment.xy1 == xy_pnt
                return segment.xy0
        for vias in self.inner_vias, self.branch_vias:
            if xy_pnt in vias:
                via = vias[xy_pnt]
                del vias[xy_pnt]

                adjoining_pnt = adjoining_pnt or get_adjoining_point()
                if adjoining_pnt not in vias:
                    # move the via
                    vias[adjoining_pnt] = via
                    via.location = adjoining_pnt
                else:
                    # just discard the via
                    pass

        old_segments, new_segments = super().remove_xypnt(xy_pnt)
        self._stcheck_structure_is_valid()
        return True, (old_segments, new_segments)

    def remove_segment(
            self,
            segment: TraceLine,
            modify_adjoining_segments = True,
            max_distance_from_center = 2.0
        ) -> tuple[bool, tuple[TraceLine, TraceLine]]:
        # Don't remove this segment if it joins one of the pins
        if segment == self.segments[0]:
            if 0 in self.pins:
                if self.pins[0].location == segment.xy0:
                    return False, None, None
        if segment == self.segments[-1]:
            if 1 in self.pins:
                if self.pins[1].location == segment.xy1:
                    return False, None, None

        # Remove trace corners
        self._unset_xypnt(segment.xy0)
        self._unset_xypnt(segment.xy1)
        
        prev_segment, next_segment = super().remove_segment(segment, modify_adjoining_segments, max_distance_from_center)

        # Deal with moving inner vias and branch vias
        for vias in self.inner_vias, self.branch_vias:
            for xy_pnt in segment.xy_points:
                if xy_pnt in vias:
                    via = vias[xy_pnt]
                    del vias[xy_pnt]

                    new_pnt = prev_segment.xy1
                    if new_pnt not in vias:
                        # move the via
                        vias[new_pnt] = via
                        via.location = new_pnt
                    else:
                        # just discard the via
                        pass

        self._stcheck_structure_is_valid()

        return True, (prev_segment, next_segment)

    def change_segment_points(
            self,
            segment: LineSegment,
            new_xy0: Point = None,
            new_xy1: Point = None,
            force = False
        ) -> tuple[tuple[bool, bool], tuple[LineSegment|None, LineSegment|None]]:
        """ Moves the points for the given segment and related components to
        the new locations.

        Parameters
        ----------
        segment : LineSegment
            The segment to be moved.
        new_xy0 : Point, optional
            The new point to move the start of the segment to, or None. By default None
        new_xy1 : Point, optional
            The new point to move the end of the segment to, or None. By default None
        force : bool, optional
            True to force the new point to be applied, regardless of if it is
            the end pin point. If True then it is assumed that the calling code
            will manage the end pins. By default False.

        Returns
        -------
        xy0_success, xy1_success: tuple[bool, bool]
            True if the xy0 or xy1 points were changed, respectively.
        prev_segment, next_segment: tuple[LineSegment|None, LineSegment|None]
            The previous and next segments that are next to the given segment.
        """
        original_pnts = segment.xy_points

        # If the pins are set then don't change these locations
        success_val_0, success_val_1 = True, True
        if not force:
            if new_xy0 is not None:
                if segment == self.segments[0]:
                    if 0 in self.pins:
                        if self.pins[0].location != new_xy0:
                            new_xy0 = None
                            success_val_0 = False
            if new_xy1 is not None:
                if segment == self.segments[-1]:
                    if 1 in self.pins:
                        if self.pins[1].location != new_xy0:
                            new_xy1 = None
                            success_val_1 = False

        # update the locations
        prev_segment, next_segment = super().change_segment_points(segment, new_xy0, new_xy1)

        # Remove trace corners
        for pnt in original_pnts:
            if pnt in self.xypnt_trace_corners:
                del self.xypnt_trace_corners[pnt]

        # Deal with moving inner vias and branch vias
        for old_pnt, new_pnt in [(original_pnts[0], new_xy0), (original_pnts[1], new_xy1)]:
            if new_pnt is None:
                continue

            if old_pnt in self.inner_vias:
                if new_pnt not in self.inner_vias:
                    via = self.inner_vias[old_pnt]
                    self.inner_vias[new_pnt] = Via(new_pnt, via.name, via.source_lines)
                del self.inner_vias[old_pnt]

            if old_pnt in self.branch_vias:
                if new_pnt not in self.branch_vias:
                    via = self.branch_vias[old_pnt]
                    self.branch_vias[new_pnt] = Via(new_pnt, via.name, via.source_lines)
                del self.branch_vias[old_pnt]
            
            if old_pnt in self.xypnt_trace_corners:
                del self.xypnt_trace_corners[old_pnt]
        
        # print(f"change_segment_points() from\n\t{original_pnts} to\n\t[{segment.xy0},{segment.xy1}]")
        self._stcheck_structure_is_valid()
        return (success_val_0, success_val_1), (prev_segment, next_segment)

    def cleanup(self):
        """ Fixes and/or detects problems with the board that would make it difficult to boolean. """
        def cleanup_short_traces():
            while len(self.segments) >= 2:
                found_short_trace = False

                for segment in self.segments:
                    if segment.length < g.TRACE_CORNER_RADIUS:
                        linenos = [fl.lineno for fl in self.source_lines]
                        print(f"Found short trace segment in lines {min(linenos)}-{max(linenos)}")

                        # replace this segment with the intersection of the previous and next segments
                        prev_segment = self.get_previous_segment(segment)
                        next_segment = self.get_next_segment(segment)
                        if prev_segment is None and next_segment is None:
                            # only one segment to this trace, do nothing
                            break
                        if prev_segment is None:
                            # no previous line, replace with the next line
                            success, (old_segments, new_segments) = self.remove_xypnt(segment.xy1)
                            assert success
                        elif next_segment is None:
                            # no next line, replace with the previous line
                            success, (old_segments, new_segments) = self.remove_xypnt(segment.xy0)
                            assert success
                        else:
                            # There is a previous and next segment, remove this segment
                            # and update the adjoining segments with their intersection point.
                            success, (previous_segment, next_segment) = self.remove_segment(segment)
                            assert success

                        self._stcheck_structure_is_valid()

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
        
        # add all via lines
        for line in layer:
            if line.v.startswith("VIA "):
                if line not in layer:
                    layer.append(line)

        pre_lines = pre_routes + pre_route + pre_layer
        target_lines = layer
        post_lines = post_layer + post_route + post_routes
        return pre_lines, target_lines, post_lines, route_name
    
    @classmethod
    def _parse_trace_lines(cls, trace: list[FLine]) -> tuple[list[TraceLine], dict[Point, Via]]:
        # parse the lines for this route+layer
        segment_lines = list(filter(lambda l: l.v.startswith("LINE "), trace))
        xy_points: set[Point] = set()
        edges: list[TraceLine] = []
        for segment_line in segment_lines:

            x0, y0, x1, y1 = tuple(map(float, segment_line.v.strip()[5:].split(" ")))
            x0, y0, x1, y1 = in2mm(x0), in2mm(y0), in2mm(x1), in2mm(y1)
            point0, point1 = Point(x0, y0), Point(x1, y1)
            xy_points.add(point0)
            xy_points.add(point1)

            edges.append(TraceLine(point0, point1, segment_line))

        # parse the vias for this route
        vias, l = Via.from_cad_file(trace)
        inner_vias: dict[Point, Via] = {}
        while len(vias) > 0:
            for via in vias:
                if via.location in xy_points:
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

            line_a, line_b = edge_group[0], edge_group[-1]
            group_inner_lines = edge_group[1:-1]
            new_edge_group = [line_a]
            while len(group_inner_lines) > 0:
                for line in copy.copy(group_inner_lines):
                    if (line.xy0 in new_edge_group[-1].xy_points) or (line.xy1 in new_edge_group[-1].xy_points):
                        new_edge_group.append(line)
                        group_inner_lines.remove(line)
                        break
            new_edge_group.append(line_b)
            edge_groups[edge_group_idx] = new_edge_group

    @classmethod
    def _order_xy_points_in_edges(cls, edge_groups: list[list[TraceLine]]):
        # orient adjacent edges
        for edge_group in edge_groups:
            if len(edge_group) == 1:
                continue

            for line_idx, line in enumerate(edge_group[:-1]):
                next_line = edge_group[line_idx+1]
                if line.xy0 in next_line.xy_points:
                    assert line.xy1 not in next_line.xy_points, f"Two lines share the same points:\n\t{line.fline}\n\t{next_line.fline}"
                    edge_group[line_idx] = line.reversed()

            # orient the last edge
            end_b = edge_group[-1]
            prev_line = edge_group[-2]
            if end_b.xy1 in prev_line.xy_points:
                assert end_b.xy0 not in prev_line.xy_points
                edge_group[-1] = end_b.reversed()

            assert all([edge_group[i].xy1 == edge_group[i+1].xy0 for i in range(len(edge_group)-1)])

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
    def _vias_at_branches(cls, edge_groups: list[list[TraceLine]]) -> list[dict[Point, Via]]:
        ret: list[dict[Point, Via]] = []

        for edge_group in edge_groups:
            group_vias = {}
            ret.append(group_vias)

            for edge in edge_group:
                if edge.is_branch:
                    if 0 in edge.joined_ends:
                        pnt = edge.xy1
                    elif 1 in edge.joined_ends:
                        pnt = edge.xy0
                    else:
                        raise RuntimeError(f"Unexpected value in joined_ends ({edge.joined_ends=})")
                    group_vias[pnt] = Via(pnt, "Branch", [edge.fline])

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
        for edge_group, group_branch_vias in zip(edge_groups, branch_vias):
            source_lines = [e.fline for e in edge_group]
            instance = cls(source_lines, route_name, layer_name, edge_group,
                           shape, bend_radius, inner_vias=inner_vias, branch_vias=group_branch_vias)
            ret.append(instance)

        return ret, pre_trace + post_trace
    
    @staticmethod
    def remove_short_traces(traces: list["SingleTrace"]):
        for trace in copy.copy(traces):
            trace_tot_length = sum([s.length for s in trace.segments])
            if trace_tot_length < Via.via_diameter() / 2 + Pin.through_hole_diameter() / 2:
                traces.remove(trace)

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
        if 0 in self.pins:
            segments[0] = LineSegment(pins[0].location, segments[0].xy1)

        if 1 in self.pins:
            segments[-1] = LineSegment(segments[-1].xy0, pins[1].location)

        return tuple(segments)

    def _get_pins_ajusted(self) -> dict[Literal[0, 1], Pin|None]:
        return self.pins

    @abstractmethod
    def to_vtk(self, trace_polydata: vtk.vtkPolyData, via_polydata: vtk.vtkPolyData):
        pass

    @abstractmethod
    def draw(self, ax: maxis.Axis):
        pass
