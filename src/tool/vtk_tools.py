import numpy as np
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

def join_with_quads(polydata: vtk.vtkPolyData, a_points_idx: int, b_points_idx: int, num_points: int):
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

def calculate_point_normals(polydata):
    """
    Calculate point normals that point away from the centroid of the polydata.

    Args:
        polydata (vtk.vtkPolyData): The input polydata with points.

    Returns:
        None. Modifies the input polydata to include normals.
    """

    if not isinstance(polydata, vtk.vtkPolyData):
        raise TypeError("Input must be a vtkPolyData object")

    # Get the number of points
    num_points = polydata.GetNumberOfPoints()
    if num_points == 0:
        return

    # Create an array to store normals
    normals = vtk.vtkDoubleArray()
    normals.SetName("Normals")
    normals.SetNumberOfComponents(3)
    normals.SetNumberOfTuples(num_points)

    # Get the points as a numpy array for easier manipulation
    points = np.array([
        polydata.GetPoint(i) for i in range(num_points)
    ])

    # Calculate centroid (average point location)
    centroid = np.mean(points, axis=0)

    # Calculate normals pointing away from centroid
    for i, point in enumerate(points):
        normal = point - centroid
        normal /= np.linalg.norm(normal)  # Normalize the vector to unit length
        normals.SetTuple3(i, *normal)

    # Assign the normals array to the polydata's point data
    polydata.GetPointData().SetNormals(normals)
