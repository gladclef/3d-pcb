import copy
from typing import Literal

import matplotlib.axis as maxis
import numpy as np
import vtk

from Component.Pin import Pin
from Component.Via import Via
from Geometry.Point import Point
from Geometry.LineSegment import LineSegment
from Trace.SingleTrace import SingleTrace
from tool.units import *


class SingleTraceWire(SingleTrace):
    def _get_pins_ajusted(self) -> dict[Literal[0, 1], Pin|None]:
        """
        Adjusts end point pin locations to be offset from the component pins.

        The distance of the offset is based on through-hole diameter and via diameter.

        Returns
        -------
        dict[int, Union[Pin, None]]
            Adjusted pins for the trace ends.
        """
        ret = copy.copy(self.pins)

        # Recalculate the pin locations to be offset from the component through holes
        for end, pin in self.pins.items():
            dist = Pin.through_hole_diameter() / 2 + Via.via_diameter() / 2

            def get_segment_at_dist(start: Point, segments: list[LineSegment], end):
                for segment in segments:
                    seg_end = segment.xy0 if end == 0 else segment.xy1
                    if start.distance(seg_end) >= dist:
                        break
                return segment

            if end == 0:
                segment = get_segment_at_dist(pin.location, self.segments, 1)
                assert pin.location.distance(segment.xy0) < Pin.through_hole_diameter() + Via.via_diameter()
                adjusted_loc = segment.distance_along_line(dist, pin.location)
            else:
                segment = get_segment_at_dist(pin.location, reversed(self.segments), 0)
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