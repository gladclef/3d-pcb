import math

import numpy as np
import pyvista
from scipy.spatial.transform import Rotation
import vtk
from vtk.util import numpy_support # type: ignore
from vtkbool import vtkBool

from Trace.AbstractTrace import AbstractTrace
from Trace.AbstractVtkPointTracker import AbstractVtkPointTracker as PntInc
import Geometry.geometry_tools as geo
from Trace.PipeShape import PipeShape
from Geometry.LineSegment import LineSegment
from Trace.TraceCorner import TraceCorner
from Trace.VtkPointGroup import VtkPointGroup
import tool.vtk_tools as vt


class SingleTrace(AbstractTrace):
    """
    Represents a singular wire trace on a PCB.

    This is limited to a single continuous line. If there are
    junctions they should be built up of multiple SingleTrace objects
    using the more general Trace class.
    """
    def __init__(self, xy_points: list[tuple[float, float]], segments: list[tuple[int, int]] | list[LineSegment], shape: PipeShape, bend_radius: float=1):
        super().__init__(xy_points, segments, shape)

        self._xypnt_vtk_verticies: dict[int, VtkPointGroup] = {}
        """ Dictionary from xy point index to vtk points. """
        self.xypnt_trace_corners: dict[int, TraceCorner] = {}
        """ Dictionary from xy point index to trace corners. """
        self.bend_radius = bend_radius

        self.check_segment_duplicates(self.segments)
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
        vt.join_with_quads(polydata, va.vtk_idx_0, vb.vtk_idx_0, len(va))
    
    def inc_vtk_indicies(self, start: int, cnt: int):
        for segment_points in self._xypnt_vtk_verticies.values():
            segment_points.inc_vtk_indicies(start, cnt)

    def to_vtk(self):
        # Create vtkPolyData and set points and cells
        vtk_points = vtk.vtkPoints()
        polydata = vtk.vtkPolyData()
        vtk_cells = vtk.vtkCellArray()
        polydata.SetPoints(vtk_points)
        polydata.SetPolys(vtk_cells)

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
        
        # Assign point normals
        vt.calculate_point_normals(polydata)
        
        return polydata

if __name__ == "__main__":
    from Trace.PipeShape import PipeBasicBox
    from tool.units import *

    trace_points = [
        (0, 0),
        (1, 0),
        (2, 1),
    ]
    trace_edges = [
        (0, 1),
        (1, 2)
    ]

    test_trace = SingleTrace(trace_points, trace_edges, PipeBasicBox(awg2mm(26)))
    vtk_test_trace = test_trace.to_vtk()
    print(f"{vtk_test_trace.GetNumberOfPoints()=}")
    vt.save_to_vtk(vtk_test_trace, "test_trace.vtk")

    test_trace_mesh = pyvista.PolyData(vtk_test_trace)
    pyvista.global_theme.allow_empty_mesh = True
    test_trace_mesh.plot(show_edges=True, opacity=1, show_vertices=True)
    # pyvista.PolyDataFilters.plot_normals(test_trace_mesh, mag=0.5, flip=False, faces=False, show_edges=True, opacity=0.95, show_verticies=True)