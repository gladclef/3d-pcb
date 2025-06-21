import numpy as np
from scipy.spatial.transform import Rotation
import vtk
from vtk.util import numpy_support # type: ignore

def save_to_vtk(polydata: vtk.vtkPolyData, file_path_name_ext: str):
    # Write the polydata to a .vtk file
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(file_path_name_ext)
    writer.SetInputData(polydata)
    writer.Write()

def insert_points(polydata: vtk.vtkPolyData, insert_start_idx: int, values: tuple[float, float, float] | list[tuple] | vtk.vtkPoints):
    vtk_points: vtk.vtkPoints = polydata.GetPoints()
    vtk_cells: vtk.vtkCellArray = polydata.GetPolys()
    npoints_orig = vtk_points.GetNumberOfPoints()

    # normalize the inputs
    if isinstance(values, vtk.vtkPoints):
        vtk_vals = values
    else:
        original_values = values
        if isinstance(values, tuple):
            values = []
            values.append(original_values)
        vtk_vals = vtk.vtkPoints()
        vtk_vals.SetData(numpy_support.numpy_to_vtk(values, deep=True))
    ninserts = vtk_vals.GetNumberOfPoints()

    # increase the length of the vtk_points
    for i in range(ninserts):
        vtk_points.InsertNextPoint((0, 0, 0))
    
    # move the existing points down
    dstStart = insert_start_idx + ninserts
    n = npoints_orig - insert_start_idx
    srcStart = insert_start_idx
    vtk_points.InsertPoints(dstStart, n, srcStart, vtk_points)

    # add in the new points
    vtk_points.InsertPoints(insert_start_idx, ninserts, 0, vtk_vals)

    # adjust the cell point ids
    for cell_idx in range(vtk_cells.GetNumberOfCells()):
        pnt_id_list = vtk.vtkIdList()
        vtk_cells.GetCellAtId(cell_idx, pnt_id_list)
        for pnt_idx in range(vtk_cells.GetCellSize(cell_idx)):
            pnt_id = pnt_id_list.GetId(pnt_idx)
            if pnt_id >= insert_start_idx:
                vtk_cells.ReplaceCellPointAtId(cell_idx, pnt_idx, pnt_id+ninserts)

def adjoin_with_quads(polydata: vtk.vtkPolyData, a_points_idx: int, b_points_idx: int, num_points: int):
    """
    Join two groups A and B with equal numbers of verticies with quads.
    This function is most useful for extruding a shape.

    Parameters
    ----------
    polydata : vtk.vtkPolyData
        The object containing the verticies.
    a_points_idx : int
        The starting vertex index for the A group.
    b_points_idx : int
        The starting vertex index for the B group.
    num_points : int
        How many verticies are in the A group. It is assumed that
        the B group has the same number of verticies.
    """
    # vtk_points: vtk.vtkPoints = polydata.GetPoints()
    vtk_cells: vtk.vtkCellArray = polydata.GetPolys()
        
    # build the sides along the length of the trace
    for i in range(num_points):
        j = 0 if i == num_points - 1 else i+1

        quad = vtk.vtkQuad()
        quad.GetPointIds().SetId(0, a_points_idx+i)
        quad.GetPointIds().SetId(1, a_points_idx+j)
        quad.GetPointIds().SetId(2, b_points_idx+j)
        quad.GetPointIds().SetId(3, b_points_idx+i)
        vtk_cells.InsertNextCell(quad)

def triangulate_quads(polydata: vtk.vtkPolyData):
    """ Find all cells with 4 point ids and replace them with two cells with 3 ids each. """
    polys = polydata.GetPolys()
    new_polys = vtk.vtkCellArray()

    # Iterate through the existing cells
    for i in range(polys.GetNumberOfCells()):
        cell = vtk.vtkIdList()
        polys.GetCellAtId(i, cell)

        if cell.GetNumberOfIds() == 4:
            # Quad found, create two triangles from it
            triangle1 = vtk.vtkIdList()
            triangle2 = vtk.vtkIdList()

            # First triangle: points 0, 1, 2
            triangle1.InsertNextId(cell.GetId(0))
            triangle1.InsertNextId(cell.GetId(1))
            triangle1.InsertNextId(cell.GetId(2))

            # Second triangle: points 0, 2, 3
            triangle2.InsertNextId(cell.GetId(0))
            triangle2.InsertNextId(cell.GetId(2))
            triangle2.InsertNextId(cell.GetId(3))

            # Add triangles to new_polys
            new_polys.InsertNextCell(triangle1)
            new_polys.InsertNextCell(triangle2)
        else:
            # Not a quad, just copy the original cell
            new_polys.InsertNextCell(cell)

    # Replace the old polys with new_polys in polydata
    polydata.SetPolys(new_polys)

def set_point_normals(polydata: vtk.vtkPolyData, new_normals: vtk.vtkFloatArray, start_idx=-1, end_idx=-1):
    point_data = polydata.GetPointData()
    point_normals: vtk.vtkDoubleArray = point_data.GetNormals()
    if point_normals is None:
        point_normals = vtk.vtkDoubleArray()
        point_normals.SetNumberOfComponents(3)
    if point_normals.GetNumberOfTuples() < polydata.GetNumberOfPoints():
        point_normals.SetNumberOfTuples(polydata.GetNumberOfPoints())

    if start_idx < 0:
        start_idx = 0
    if end_idx < 0:
        end_idx = min(start_idx + new_normals.GetNumberOfTuples(), polydata.GetNumberOfPoints())

    for i in range(end_idx - start_idx):
        point_normals.SetTuple3(i + start_idx, *new_normals.GetTuple3(i))
    
    point_data.SetNormals(point_normals)

def calculate_point_normals(polydata, start_idx=-1, end_idx=-1, from_point: tuple[float, float, float]=None):
    """
    Calculate point normals that point away from the centroid of the polydata.

    Args:
        polydata (vtk.vtkPolyData): The input polydata with points.
        start_idx (int): The range of vtk points to calculate point normals for, inclusive.
        end_idx (int): The range of vtk points to calculate point normals for, exclusive.

    Returns:
        None. Modifies the input polydata to include normals.
    """

    if not isinstance(polydata, vtk.vtkPolyData):
        raise TypeError("Input must be a vtkPolyData object")

    # Get the number of points
    num_points = polydata.GetNumberOfPoints()
    if num_points == 0:
        return
    
    # Get the range
    start_idx = start_idx if start_idx > -1 else 0
    end_idx = end_idx if end_idx > -1 else num_points
    to_calculate = end_idx - start_idx

    # Create an array to store normals
    normals = vtk.vtkDoubleArray()
    normals.SetName("Normals")
    normals.SetNumberOfComponents(3)
    normals.SetNumberOfTuples(to_calculate)

    # Get the points as a numpy array for easier manipulation
    points = np.array([
        polydata.GetPoint(i) for i in range(start_idx, end_idx)
    ])

    # Choose the ray casting point
    # (calculating normals is hard, just use rays from the centroid)
    if from_point is None:
        centroid = np.mean(points, axis=0)
        from_point = centroid
    from_point = np.array(from_point)

    # Calculate normals pointing away from centroid
    for i, point in enumerate(points):
        normal = point - from_point
        normal /= np.linalg.norm(normal)  # Normalize the vector to unit length
        normals.SetTuple3(i, *normal)

    # Assign the normals array to the polydata's point data
    set_point_normals(polydata, normals, start_idx, end_idx)

def new_polydata() -> vtk.vtkPolyData:
    # Create vtkPolyData and set points and cells
    vtk_points = vtk.vtkPoints()
    polydata = vtk.vtkPolyData()
    vtk_cells = vtk.vtkCellArray()
    polydata.SetPoints(vtk_points)
    polydata.SetPolys(vtk_cells)
    return polydata

def verts_and_cells(polydata: vtk.vtkPolyData) -> tuple[vtk.vtkPoints, vtk.vtkCellArray]:
    vtk_points: vtk.vtkPoints = polydata.GetPoints()
    vtk_cells: vtk.vtkCellArray = polydata.GetPolys()
    return vtk_points, vtk_cells

def rotate(polydata: vtk.vtkPolyData, rot: Rotation):
    verts, cells = verts_and_cells(polydata)

    for vi in range(verts.GetNumberOfPoints()):
        x, y, z = verts.GetPoint(vi)
        rotated = rot.apply(np.array([x, y, z]))
        rx, ry, rz = rotated[0], rotated[1], rotated[2]
        verts.SetPoint(vi, (rx, ry, rz))

def translate(polydata: vtk.vtkPolyData, translation: tuple[float, float, float]):
    verts, cells = verts_and_cells(polydata)

    for vi in range(verts.GetNumberOfPoints()):
        x, y, z = verts.GetPoint(vi)
        tx, ty, tz = x+translation[0], y+translation[1], z+translation[2]
        verts.SetPoint(vi, (tx, ty, tz))

def join(polydata1: vtk.vtkPolyData, polydata2: vtk.vtkPolyData, *other_polydatas: vtk.vtkPolyData):
    """ Merges polydata2 into polydata1 by adding new verticies and associated cells. """
    v1, c1 = verts_and_cells(polydata1)
    v2, c2 = verts_and_cells(polydata2)

    # insert the points
    start_idx = v1.GetNumberOfPoints()
    for vi in range(v2.GetNumberOfPoints()):
        v1.InsertNextPoint(v2.GetPoint(vi))
    set_point_normals(polydata1, polydata2.GetPointData().GetNormals(), start_idx)

    # insert cells
    for ci in range(c2.GetNumberOfCells()):
        cell = vtk.vtkIdList()
        c2.GetCellAtId(ci, cell)

        if not hasattr(cell, 'GetPointId'):
            old_poly = cell
            new_poly = vtk.vtkPolygon()
            new_poly.GetPointIds().SetNumberOfIds(cell.GetNumberOfIds())
            for i in range(cell.GetNumberOfIds()):
                id = old_poly.GetId(i)+start_idx
                new_poly.GetPointIds().SetId(i, id)
            c1.InsertNextCell(new_poly)

        elif cell.GetNumberOfIds() == 3:
            # insert triangles
            tri_old: vtk.vtkTriangle = cell
            tri_new = vtk.vtkTriangle()
            tri_new.GetPointIds().SetId(0, tri_old.GetPointId(0)+start_idx)
            tri_new.GetPointIds().SetId(1, tri_old.GetPointId(1)+start_idx)
            tri_new.GetPointIds().SetId(2, tri_old.GetPointId(2)+start_idx)
            c1.InsertNextCell(tri_new)

        elif cell.GetNumberOfIds() == 4:
            # insert quads
            quad_old: vtk.vtkQuad = cell
            quad_new = vtk.vtkQuad()
            quad_new.GetPointIds().SetId(0, quad_old.GetPointId(0)+start_idx)
            quad_new.GetPointIds().SetId(1, quad_old.GetPointId(1)+start_idx)
            quad_new.GetPointIds().SetId(2, quad_old.GetPointId(2)+start_idx)
            quad_new.GetPointIds().SetId(3, quad_old.GetPointId(3)+start_idx)
            c1.InsertNextCell(quad_new)
        
        else:
            raise RuntimeError(f"In vtk_tools.join(): unknown cell type with {cell.GetNumberOfIds()} ids")

    # join the other polydatas
    for other_polydata in other_polydatas:
        join(polydata1, other_polydata)