import math

import numpy as np
import pyvista
from scipy.spatial.transform import Rotation
import vtk
from vtk.util import numpy_support # type: ignore
from vtkbool import vtkBool

from PipeShape import PipeShape
import vtk_tools as vt

class XYZ:
    def __init__(self, x: float, y: float, z: float):
        self._val_tuple = x, y, z

    @property
    def val(self) -> np.ndarray:
        return np.array(self._val_tuple)

    @property
    def x(self) -> float:
        return self.val[0]
    
    @property
    def y(self) -> float:
        return self.val[1]
    
    @property
    def z(self) -> float:
        return self.val[2]
    
    @staticmethod
    def as_np_array(xyz_list: list["XYZ"]) -> np.ndarray:
        ret = []
        for v in xyz_list:
            ret.append([*v._val_tuple])
        return np.array(ret)
    
    def __hash__(self):
        return hash(self._val_tuple)

    def __eq__(self, value):
        if not isinstance(value, XYZ):
            return False
        return self._val_tuple == value._val_tuple
    
    def __repr__(self):
        return f"<{self.x:0.02f},{self.y:0.02f},{self.z:0.02f}>"

class Trace:
    """ Represents a wire trace on a PCB. """
    def __init__(self, xy_points: list[tuple[float, float]], segments: list[tuple[int, int]], shape: PipeShape):
        self.xy_points = xy_points
        self.segments = segments
        self.shape = shape

        self.check_segments_overlap()
        self.check_segments_overlap()
    
    def check_segments_overlap(self):
        segments_set = set(self.segments)
        if len(self.segments) != len(segments_set):
            raise ValueError(f"Error in Trace.check_segments_overlap(): there are {len(segments_set)} non-duplicate segments out of {len(self.segments)} total segments.")

    def get_setment_angle(self, segment: tuple[int, int]) -> float:
        pa, pb = self.xy_points[segment[0]], self.xy_points[segment[1]]
        ang = math.atan2(pb[1] - pa[1], pb[0] - pa[0])
        return ang

    def check_segments_overlap(self):
        for e1idx, segment1 in enumerate(self.segments):
            for e2idx, segment2 in enumerate(self.segments):
                if e1idx == e2idx:
                    # don't check for self-segment intersections
                    continue
                if segment1[0] == segment2[0] or segment1[0] == segment2[1] or segment1[1] == segment2[0] or segment1[1] == segment2[1]:
                    # don't check for intersections with segments that share one of the end points
                    continue

                # check if these segments are parallel
                ang1, ang2 = self.get_setment_angle(segment1), self.get_setment_angle(segment2)
                if abs(ang1 - ang2) < 0.001:
                    continue

                # check if these segments overlap
                def line_intersection(p0, p1, p2, p3):
                    s1_x = p1[0] - p0[0]
                    s1_y = p1[1] - p0[1]
                    s2_x = p3[0] - p2[0]
                    s2_y = p3[1] - p2[1]

                    s = (-s1_y * (p0[0] - p2[0]) + s1_x * (p0[1] - p2[1])) / (-s2_x * s1_y + s2_y * s1_x)
                    t = ( s2_x * (p0[1] - p2[1]) - s2_y * (p0[0] - p2[0])) / (-s2_x * s1_y + s2_y * s1_x)

                    if s >= 0 and s <= 1 and t >= 0 and t <= 1:
                        # Collision detected
                        i_x = p0[0] + (t * s1_x)
                        i_y = p0[1] + (t * s1_y)
                        return True, i_x, i_y

                    return False, None, None

                p1a, p1b = self.xy_points[segment1[0]], self.xy_points[segment1[1]]
                p2a, p2b = self.xy_points[segment2[0]], self.xy_points[segment2[1]]
                intersection, i_x, i_y = line_intersection(p1a, p1b, p2a, p2b)
                if intersection:
                    raise ValueError(f"Error in Trace.check_segments_overlap(): segments {segment1} (a={p1a}, b={p1b}) and {segment2} (a={p2a}, b={p2b}) overlap at [{i_x}, {i_y}].")

    def get_segment_vertices(self, xy_point: tuple[float, float], angle: float, top_or_bottom='top'):
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
        xyz_rotated = r.apply(xyz_points)
        xyz_translated = xyz_rotated + np.array([xy_point[0], xy_point[1], 0])

        # convert to XYZ instances
        xyz_objs: list[XYZ] = []
        for i in range(len(xyz_points)):
            xyz_vals = xyz_translated[i]
            xyz_objs.append(XYZ(*xyz_vals))
        
        return xyz_objs

    def to_vtk(self):
        # build the segment pipes
        trace_segments: list[vtk.vtkPolyData] = []
        for seg_idx, segment in enumerate(self.segments):
            pa, pb = self.xy_points[segment[0]], self.xy_points[segment[1]]
            angle = self.get_setment_angle(segment)

            # get the vertices at either end of the trace segment
            va = self.get_segment_vertices(pa, angle)
            vb = self.get_segment_vertices(pb, angle)

            # Create vtkPoints
            xyz_array = XYZ.as_np_array(va + vb)
            xyz_array[:,2] += seg_idx*0.1
            vtk_points = vtk.vtkPoints()
            vtk_points.SetData(numpy_support.numpy_to_vtk(xyz_array, deep=True))

            # Create vtkPolyData and set points
            polydata = vtk.vtkPolyData()
            polydata.SetPoints(vtk_points)
            vtk_cells = vtk.vtkCellArray()
            polydata.SetPolys(vtk_cells)

            # Close the ends
            pa0, pb0 = 0, len(va)
            pnts_tot = vtk_points.GetNumberOfPoints()
            assert pnts_tot == len(va) + len(vb)
            self.shape.add_vtk_cells(polydata, pa0, add_edges=False, fill_with_triangles=True)
            pb0 += vtk_points.GetNumberOfPoints() - pnts_tot
            self.shape.add_vtk_cells(polydata, pb0, add_edges=False, fill_with_triangles=True)
            
            # build the sides along the length of the trace
            for i in range(len(va)):
                j = 0 if i == len(va) - 1 else i+1

                quad = vtk.vtkQuad()
                quad.GetPointIds().SetId(0, i)
                quad.GetPointIds().SetId(1, j)
                quad.GetPointIds().SetId(2, pb0+j)
                quad.GetPointIds().SetId(3, pb0+i)
                vtk_cells.InsertNextCell(quad)
            
            # Assign point normals
            vt.calculate_point_normals(polydata)

            if True:
                # Register the trace segment
                trace_segments.append(polydata)
            else:
                # print("PipeShape")
                # print(polydata)
                # print("")

                polydata = pyvista.Sphere(radius=2, center=(*pa, 0))
                # print("Sphere")
                # print(polydata)
                # print(polydata.GetCellData())
                # print(polydata.GetPointData())
                # print("")
                trace_segments.append(polydata)
                polydata = pyvista.Sphere(radius=2, center=(*pb, 0))
                trace_segments.append(polydata)
        
        # build up a union of all the segment pipes
        ret = None
        for trace_segment in trace_segments:
            vt.triangulate_quads(trace_segment)
            pv_trace_segment = pyvista.PolyData(trace_segment)
            # print(pv_trace_segment)
            if ret is None:
                ret = pv_trace_segment
            else:
                # boolean_op = vtk.vtkBooleanOperationPolyDataFilter()
                # boolean_op.SetOperationToUnion()
                # boolean_op.SetInputData(0, ret)
                # boolean_op.SetInputData(1, pv_trace_segment)
                # ret = ret.boolean_intersection(pv_trace_segment, tolerance=1e-1000)
                # ret = ret.merge(pv_trace_segment, tolerance=1e-100)

                bf = vtkBool.vtkPolyDataBooleanFilter()
                bf.SetDebug(True)
                bf.SetInputData(0, ret)
                bf.SetInputData(1, pv_trace_segment)
                bf.SetOperModeToUnion()
                bf.Update()
                ret = bf.GetOutput()
        
        return ret

if __name__ == "__main__":
    from PipeShape import PipeBasicBox
    from units import *

    trace_points = [
        (0, 0),
        (4, 0)
    ]
    trace_segments = [
        (0, 1)
    ]
    test_trace = Trace(trace_points, trace_segments, PipeBasicBox(awg2mm(26)))

    vt.save_to_vtk(test_trace.to_vtk(), "test_trace.vtk")

    test_trace_mesh: pyvista.PolyData = pyvista.read("test_trace.vtk")
    pyvista.global_theme.allow_empty_mesh = True
    test_trace_mesh.plot(show_edges=True, opacity=0.9, show_vertices=True)