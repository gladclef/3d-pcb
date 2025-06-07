import copy
import math

import numpy as np
import vtk
from vtk.util import numpy_support # type: ignore

from units import *
import vtk_tools as vt

BREADBOARD_SPACING=in2mm(0.1)
NOZZLE_DIAMETER=0.4
LAYER_HEIGHT=0.1

class PipeShape:
    def __init__(self, xz_pointz: list[list[float]], symetric_about_x=False, symetric_about_z=True):
        global BREADBOARD_SPACING
        global NOZZLE_DIAMETER
        global LAYER_HEIGHT
        self.BREADBOARD_SPACING = BREADBOARD_SPACING
        self.NOZZLE_DIAMETER = NOZZLE_DIAMETER
        self.LAYER_HEIGHT = LAYER_HEIGHT
        
        self.symetric_about_x = symetric_about_x
        self.symetric_about_z = symetric_about_z
        self._xz_pointz = xz_pointz
        self.xz_pointz = self.normalize_points()
    
    @property
    def diameter(self) -> float:
        xz_points = self.normalize_points()
        xs = [p[0] for p in xz_points]
        xmin, xmax = np.min(xs), np.max(xs)
        return xmax - xmin

    def normalize_points(self) -> list[list[float]]:
        """ Returns the actual xz points of the shape, with properties applied. """
        ret = copy.deepcopy(self._xz_pointz)

        # add new points, reflected about the x axis
        if self.symetric_about_x:
            new_points: list[list[float]] = []
            for pnt in reversed(ret):
                pnt = [pnt[0], -pnt[1]]
                new_points.append(pnt)
            ret += new_points
            
        # add new points, reflected about the z axis
        if self.symetric_about_z:
            new_points: list[list[float]] = []
            for pnt in reversed(ret):
                pnt = [-pnt[0], pnt[1]]
                new_points.append(pnt)
            ret += new_points

        # Adjust the points so that the highest point is at z=0
        z_max = np.max([z for x, z in ret])
        ret = [[x, z-z_max] for x, z in ret]
        
        z_max = np.max([z for x, z in ret])
        assert z_max == 0

        # Adjust the points so that the shape is centered around x=0
        x_min = np.min([x for x, y in ret])
        x_max = np.max([x for x, y in ret])
        x_offset = x_min + ((x_max - x_min) / 2)
        ret = [[x-x_offset, y] for x, y in ret]
        
        x_min = np.min([x for x, y in ret])
        x_max = np.max([x for x, y in ret])
        assert -x_min == x_max

        # # Reorder the points so that they are in polar angle order,
        # # starting with the point closest to angle=0.
        # points_with_angles = [(math.atan2(p[1], p[0]), p) for p in ret]
        # sorted_points_with_angles = sorted(points_with_angles)
        # ret = [point for _, point in sorted_points_with_angles]
        
        return ret
    
    def to_vtk(self, mode: str = None) -> vtk.vtkPolyData:
        # mirror and normalize
        xz_pointz = self.normalize_points()

        # Add a y axis
        xyz_points = [[x, 0, z] for x, z in xz_pointz]

        # Create vtkPoints
        vtk_points = vtk.vtkPoints()
        vtk_points.SetData(numpy_support.numpy_to_vtk(xyz_points, deep=True))

        # Create vtkPolyData and set points
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(vtk_points)
        vtk_cells = vtk.vtkCellArray()
        polydata.SetPolys(vtk_cells)

        self.add_vtk_cells(polydata, 0, mode)

        return polydata

    def add_vtk_cells(self, polydata: vtk.vtkPolyData, starting_point_id: int, mode="triangles") -> list[int]:
        p0_id = starting_point_id
        npoints = len(self.normalize_points())
        vtk_points: vtk.vtkPoints = polydata.GetPoints()
        vtk_cells: vtk.vtkCellArray = polydata.GetPolys()

        # set defaults
        if mode is None:
            mode = "triangles"

        # get the set of points that this shape belongs to
        xyz_points: list[tuple[float, float, float]] = []
        for i in range(npoints):
            j = i + p0_id
            xyz_points.append(vtk_points.GetPoint(j))

        if mode == "triangles":
            # Add a point in the middle for all triangles to adjoin to
            xs, ys, zs = [x for x, y, z in xyz_points], [y for x, y, z in xyz_points], [z for x, y, z in xyz_points]
            z_min = np.min(zs)
            x_avg, y_avg = np.average(xs), np.average(ys)
            mid_point = tuple([x_avg, y_avg, z_min/2])
            xyz_points.append(mid_point)
            mid_point_id = p0_id + npoints
            vt.insert_points(polydata, mid_point_id, mid_point)

            # Create triangles
            for i in range(npoints):
                p1 = p0_id+i
                p2 = p0_id if i == npoints-1 else p1+1
                triangle = vtk.vtkTriangle()
                triangle.GetPointIds().SetId(0, p1)
                triangle.GetPointIds().SetId(1, p2)
                triangle.GetPointIds().SetId(2, mid_point_id)
                vtk_cells.InsertNextCell(triangle)
            
            return [mid_point_id]

        elif mode == "edges":
            # Create lines
            for i in range(npoints):
                p1 = p0_id+i
                p2 = p0_id if i == npoints-1 else p1+1
                line = vtk.vtkLine()
                line.GetPointIds().SetId(0, p1)
                line.GetPointIds().SetId(1, p2)
                vtk_cells.InsertNextCell(line)
            
            return []

        else:
            raise RuntimeError(f"In PipeShape.add_vtk_cells(): unknown mode, expected either \"edges\" or \"triangles\" but got \"{mode}\".")
    
class PipeBasicBox(PipeShape):
    """
    Box shape that is slightly wider than the wire diameter,
    with an opening at the top that is only as wide as the wire.
    """
    def __init__(self, wire_diameter: float, lip_height: float=None, max_box_width: float=None):
        """
        Parameters
        ----------
        wire_diameter : float
            The size of the wire in millimeters.
        lip_height : float, optional
            The height of the lip (aka opening) from the surface to the box, by default LAYER_HEIGHT*2.
        max_box_width : float, optional
            The maximum width of the internal box structure, by default BREADBOARD_SPACING-NOZZLE_DIAMETER.
        """
        global BREADBOARD_SPACING
        global NOZZLE_DIAMETER
        global LAYER_HEIGHT

        wire_radius = wire_diameter / 2
        xz_pointz: list[list[float]] = []

        # assign default values
        lip_height = lip_height if lip_height is not None else LAYER_HEIGHT*2
        max_box_width = max_box_width if max_box_width is not None else BREADBOARD_SPACING-NOZZLE_DIAMETER

        # opening
        xz_pointz.append([wire_radius, 0])
        xz_pointz.append([wire_radius, -lip_height])

        # box
        box_radius = min(max(wire_diameter*1.5, wire_diameter+NOZZLE_DIAMETER), max_box_width)
        box_height = wire_diameter+LAYER_HEIGHT
        xz_pointz.append([box_radius, -lip_height])
        xz_pointz.append([box_radius, -lip_height-box_height])

        super().__init__(xz_pointz, symetric_about_x=False, symetric_about_z=True)
    
class PipeDebug(PipeShape):
    """
    Box with 4 sides. Generally used for debugging.
    """
    def __init__(self, side_length: float):
        """
        Parameters
        ----------
        side_length : float
            The size of the box in millimeters.
        """
        xz_pointz: list[list[float]] = []

        # opening
        xz_pointz.append([side_length/2, 0])
        xz_pointz.append([side_length/2, side_length])

        super().__init__(xz_pointz, symetric_about_x=False, symetric_about_z=True)

