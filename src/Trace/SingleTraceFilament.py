import matplotlib.axis as maxis
import numpy as np
import vtk

from Component.Pin import Pin
from Component.Via import Via
from Geometry.LineSegment import LineSegment
from Trace.SingleTrace import SingleTrace
from tool.units import *


class SingleTraceFilament(SingleTrace):
    def _adjust_ends_for_conductive_filament(self):
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
                
                while True:
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
                                success, (ps, ns) = self.remove_xypnt(s.xy1)
                                assert success
                            else:
                                success, (ps, ns) = self.remove_xypnt(s.xy0)
                                assert success

                    if len(segments_to_remove) == 0:
                        break
                assert segment.xy0 in self.xy_points
                assert segment.xy1 in self.xy_points

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

        # Adjust the end points to be at the correct length from the pins
        self._adjust_ends_for_conductive_filament()

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
        for pin in self.pins.values():
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
        # Adjust the end points to be at the correct length from the pins
        self._adjust_ends_for_conductive_filament()

        # Recalculate the starting and ending segments depending on
        # if we should be adding through-holes for the traces.
        segments = self._get_segments_with_through_holes()

        for pin in self.pins.values():
            if pin is not None:
                pin.draw(ax)

        for seg in segments:
            ax.arrow(seg.x0, seg.y0, seg.x1-seg.x0, seg.y1-seg.y0, color="teal", head_width=.3)

        for via in self.inner_vias.values():
            via.draw(ax)

        for via in self.branch_vias.values():
            via.draw(ax)
