import math

import numpy as np
import pyvista
from scipy.spatial.transform import Rotation
import vtk
from vtk.util import numpy_support # type: ignore
from vtkbool import vtkBool

from AbstractTrace import AbstractTrace
from AbstractVtkPointTracker import AbstractVtkPointTracker as PntInc
import geometry as geo
from PipeShape import PipeShape
from Segment import Segment
import vtk_tools as vt

class SegmentPoints(PntInc):
    def __init__(self, xyz_points: np.ndarray, vtk_indices: list[int] = None):
        if vtk_indices is None:
            vtk_indices = [None for i in range(xyz_points.shape[0])]

        self.xyz_points = xyz_points
        self.vtk_indices = vtk_indices

    def inc_vtk_indicies(self, start: int, cnt: int):
        for i in range(len(self.vtk_indices)):
            if self.vtk_indices[i] >= start:
                self.vtk_indices[i] += cnt
    
    def __len__(self):
        return self.xyz_points.shape[0]

class SingleTrace(AbstractTrace):
    """
    Represents a singular wire trace on a PCB.

    This is limited to a single continuous line. If there are
    junctions they should be built up of multiple SingleTrace objects
    using the more general Trace class.
    """
    def __init__(self, xy_points: list[tuple[float, float]], segments: list[tuple[int, int]] | list[Segment], shape: PipeShape):
        super().__init__(xy_points, shape)
        self.segments = [(s if isinstance(s, Segment) else Segment(self, s)) for s in segments]
        print(repr(self.segments))
        self.segment_vertices: dict[int, SegmentPoints] = {}
        """ Dictionary from segment index to points. """

        self.check_segment_duplicates(self.segments)
        self.check_segments_overlap()
    
    def check_segment_duplicates(self, segments: list[Segment]):
        segment_tuples = [s.xy_point_indicies for s in segments]
        segments_set = set(segment_tuples)
        if len(segment_tuples) != len(segments_set):
            raise ValueError(f"Error in SingleTrace.check_segment_duplicates(): there are {len(segments_set)} non-duplicate segments out of {len(segments)} total segments.")

    def check_segments_overlap(self):
        for e1idx, segment1 in enumerate(self.segments):
            for e2idx, segment2 in enumerate(self.segments):
                if e1idx == e2idx:
                    # don't check for self-segment intersections
                    continue

                if segment1.xy0 == segment2.xy0 or \
                    segment1.xy0 == segment2.xy1 or \
                    segment1.xy1 == segment2.xy0 or \
                    segment1.xy1 == segment2.xy1:
                    # don't check for intersections with segments that share one of the end points
                    continue

                # check if these segments are parallel
                angle_diff = abs(segment1.angle - segment2.angle)
                if angle_diff < 0.001 or 2*np.pi - angle_diff < 0.001:
                    continue

                p1a, p1b = self.xy_points[segment1[0]], self.xy_points[segment1[1]]
                p2a, p2b = self.xy_points[segment2[0]], self.xy_points[segment2[1]]
                intersection = geo.line_segments_intersection((p1a, p1b), (p2a, p2b))
                if intersection is not None:
                    raise ValueError(f"Error in SingleTrace.check_segments_overlap(): segments {segment1} (a={p1a}, b={p1b}) and {segment2} (a={p2a}, b={p2b}) overlap at [{intersection}].")

    def get_segment_vertices(self, xy_point: tuple[float, float], angle: float, top_or_bottom='top') -> np.ndarray:
        """ Get the vertices for the end-cap of a given trace segment, at the given point in the xy plane.

        Parameters
        ----------
        xy_point : tuple[float, float]
            Either end of an trace segment, in the xy plane.
        angle : float
            The angle of the trace segment in radians, as determined by its two end points.
        top_or_bottom : str, optional
            If the vertices are part of a top or bottom trace.
        """

        # get the cross section from the shape
        xz_points = np.array(self.shape.normalize_points())
        xyz_points = np.array([[x, 0, z] for x, z in xz_points])
        
        # apply rotation and translation
        r = Rotation.from_euler('z', angle+np.pi/2)
        xyz_rotated: np.ndarray = r.apply(xyz_points)
        xyz_translated: np.ndarray = xyz_rotated + np.array([xy_point[0], xy_point[1], 0])
        
        return xyz_translated
    
    def segment_to_vtk(self, polydata: vtk.vtkPolyData, segment: tuple[int, int]):
        vtk_points: vtk.vtkPoints = polydata.GetPoints()
        vtk_cells: vtk.vtkCellArray = polydata.GetPolys()

        # get core segment values
        xy_idx_a, xy_idx_b = segment[0], segment[1]
        xy_a, xy_b = self.xy_points[xy_idx_a], self.xy_points[xy_idx_b]
        angle = self.get_setment_angle(segment)

        # get the vertices at either end of the trace segment
        if xy_idx_a not in self.segment_vertices:
            self.segment_vertices[xy_idx_a] = SegmentPoints(self.get_segment_vertices(xy_a, angle))
        if xy_idx_b not in self.segment_vertices:
            self.segment_vertices[xy_idx_b] = SegmentPoints(self.get_segment_vertices(xy_b, angle))
        va = self.segment_vertices[xy_idx_a]
        vb = self.segment_vertices[xy_idx_b]

        # Add to vtk points
        for verticies in [va, vb]:
            for i in range(len(verticies)):
                if verticies.vtk_indices[i] is None:
                    xyz_points: np.ndarray = verticies.xyz_points[i]
                    vtk_points.InsertNextPoint(xyz_points.tolist())
                    verticies.vtk_indices[i] = vtk_points.GetNumberOfPoints()-1
        
        # build the sides along the length of the trace
        pa0 = va.vtk_indices[0]
        pb0 = vb.vtk_indices[0]
        for i in range(len(va)):
            j = 0 if i == len(va) - 1 else i+1

            quad = vtk.vtkQuad()
            quad.GetPointIds().SetId(0, pa0+i)
            quad.GetPointIds().SetId(1, pa0+j)
            quad.GetPointIds().SetId(2, pb0+j)
            quad.GetPointIds().SetId(3, pb0+i)
            vtk_cells.InsertNextCell(quad)
    
    def inc_vtk_indicies(self, start: int, cnt: int):
        for segment_points in self.segment_vertices.values():
            segment_points.inc_vtk_indicies(start, cnt)

    def to_vtk(self):
        # Create vtkPolyData and set points and cells
        vtk_points = vtk.vtkPoints()
        polydata = vtk.vtkPolyData()
        vtk_cells = vtk.vtkCellArray()
        polydata.SetPoints(vtk_points)
        polydata.SetPolys(vtk_cells)

        # build the segment pipes
        for seg_idx, segment in enumerate(self.segments):
            # build the segment verts and cells
            self.segment_to_vtk(polydata, segment)
            
        # Close the ends
        xy_idx_a = 0
        xy_idx_b = max(list(self.segment_vertices.keys()))
        for xy_idx in [xy_idx_a, xy_idx_b]:
            vertices = self.segment_vertices[xy_idx]
            point_id_0 = vertices.vtk_indices[0]
            new_vtk_points = self.shape.add_vtk_cells(polydata, point_id_0)
            self.inc_vtk_indicies(point_id_0 + len(vertices), len(new_vtk_points))
            for new_vtk_id in new_vtk_points:
                xyz_point = list(vtk_points.GetPoint(new_vtk_id))
                vertices.xyz_points = np.concat((vertices.xyz_points, np.array([xyz_point])), axis=0)
                vertices.vtk_indices.append(new_vtk_id)
        
        # Assign point normals
        vt.calculate_point_normals(polydata)
        
        return polydata

if __name__ == "__main__":
    from PipeShape import PipeBasicBox
    from units import *

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
    # pyvista.PolyDataFilters.plot_normals(test_trace_mesh, mag=0.5, flip=False, faces=False, show_edges=True, opacity=0.95, show_vertices=True)