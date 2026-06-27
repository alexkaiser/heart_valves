
import os, vtk
import math
import numpy as np
import json 
from vtkmodules.util.numpy_support import vtk_to_numpy

from scipy.io import loadmat, savemat




def read_vtp(filename):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()

def write_vtk_polydata(poly, fn):
    """
    This function writes a vtk polydata to disk
    """
    # print('Writing vtp with name:', fn)
    if (fn == ""):
        return 0
    _ , extension = os.path.splitext(fn)
    if extension == '.vtk':
        writer = vtk.vtkPolyDataWriter()
    elif extension == '.stl':
        writer = vtk.vtkSTLWriter()
    elif extension == '.vtp':
        writer = vtk.vtkXMLPolyDataWriter()
    elif extension == '.obj':
        writer = vtk.vtkOBJWriter()
    else:
        raise ValueError("Incorrect extension" + extension)
    writer.SetInputData(poly)
    writer.SetFileName(fn)
    writer.Update()
    writer.Write()  
    return


def pca_2d_projection(points3d):
    X = np.array(points3d)
    X_centered = X - np.mean(X, axis=0)
    cov = np.cov(X_centered.T)
    eigvals, eigvecs = np.linalg.eigh(cov)
    idx = np.argsort(eigvals)[::-1]
    principal_components = eigvecs[:, idx[:2]]
    projected = X_centered @ principal_components
    return projected


def _segments_intersect_2d(a1, a2, b1, b2):
    def orient(p, q, r):
        return (q[0]-p[0])*(r[1]-p[1]) - (q[1]-p[1])*(r[0]-p[0])
    o1 = orient(a1,a2,b1)
    o2 = orient(a1,a2,b2)
    o3 = orient(b1,b2,a1)
    o4 = orient(b1,b2,a2)
    return (o1*o2 < 0) and (o3*o4 < 0)


def order_points(points):
    """
    Orders a closed loop of 3D points by:
      1) Starting at the point farthest from the centroid,
      2) Greedily walking to nearest neighbors,
      3) Skipping any choice that would cause a crossing in the best-fit 2D projection.
    """
    pts = [tuple(p) for p in points] 
    n = len(pts)
    if n < 3:
        return pts + pts[:1]  

    # 1) Find best-fit plane via PCA
    X = np.array(pts)
    uv = pca_2d_projection(pts)


    # 2) Pick starting point = farthest from centroid
    centroid = np.mean(X, axis=0)
    dists = np.linalg.norm(X - centroid, axis=1)
    start_idx = int(np.argmax(dists))
    ordered_3d = [pts.pop(start_idx)]
    ordered_2d = [tuple(uv[start_idx])]
    uv = np.delete(uv, start_idx, axis=0)

    # 3) Walk
    while pts:
        last3 = np.array(ordered_3d[-1])
        last2 = ordered_2d[-1]

        candidates = sorted(
            [(math.dist(last3, pts[i]), i) for i in range(len(pts))],
            key=lambda x: x[0]
        )

        for _, ci in candidates:
            cand3 = pts[ci]
            cand2 = tuple(uv[ci])

            crosses = False
            for j in range(len(ordered_2d)-1):
                if _segments_intersect_2d(
                    ordered_2d[j], ordered_2d[j+1], 
                    last2, cand2
                ):
                    crosses = True
                    break
            if not crosses:
                ordered_3d.append(pts.pop(ci))
                ordered_2d.append(cand2)
                uv = np.delete(uv, ci, axis=0)
                break
        else:
            _, ci = candidates[0]
            ordered_3d.append(pts.pop(ci))
            ordered_2d.append(tuple(uv[ci]))
            uv = np.delete(uv, ci, axis=0)

    # 4) Close loop
    # ordered_3d.append(ordered_3d[0])

    return ordered_3d


def convertPointsToVTKPolyData(points):
    vtkPointsObj = vtk.vtkPoints()

    if isinstance(points, vtk.vtkPoints):
        num_points = points.GetNumberOfPoints()
        for i in range(num_points):
            vtkPointsObj.InsertNextPoint(points.GetPoint(i))
    elif isinstance(points, list):
        for p in points:
            vtkPointsObj.InsertNextPoint(p)
    else:
        raise TypeError("Input must be a list of (x, y, z) tuples or a vtk.vtkPoints object.")

    num_points = vtkPointsObj.GetNumberOfPoints()
    polyLine = vtk.vtkPolyLine()
    polyLine.GetPointIds().SetNumberOfIds(num_points + 1)

    for i in range(num_points):
        polyLine.GetPointIds().SetId(i, i)
    # polyLine.GetPointIds().SetId(num_points, 0)  # Close the loop

    cells = vtk.vtkCellArray()
    cells.InsertNextCell(polyLine)

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(vtkPointsObj)
    polyData.SetLines(cells)

    return polyData


def resamplePointsOnMesh(initial_curve_points, num_landmark_points, mesh_polydata, extended_baffle_bool = None):
    # 1) Order and check curve points
    orderedCurvePoints = order_points(initial_curve_points)
    # 2)Convert to VTK PolyData
    vtkCurvePolyData = convertPointsToVTKPolyData(orderedCurvePoints)

    curvePoints = vtkCurvePolyData.GetPoints()
    numberOfCurvePoints = curvePoints.GetNumberOfPoints() if curvePoints else 0 

    # 3) Compute total curve length
    curveLength = 0 
    for i in range(1, numberOfCurvePoints):
        p0 = curvePoints.GetPoint(i-1)
        p1 = curvePoints.GetPoint(i)
        curveLength += math.sqrt((p1[0] - p0[0]) ** 2 + (p1[1] - p0[1]) ** 2 + (p1[2] - p0[2]) ** 2)

    samplingDistance = curveLength / (num_landmark_points-0.1)
    
    # 4) Create and sample the spline
    spline = vtk.vtkParametricSpline()
    spline.SetPoints(vtkCurvePolyData.GetPoints())
    
    splineFunctionSource = vtk.vtkParametricFunctionSource()
    splineFunctionSource.SetParametricFunction(spline)
    splineFunctionSource.SetUResolution(int(curveLength / samplingDistance))
    splineFunctionSource.Update()

    resampledPolyData = splineFunctionSource.GetOutput()
    resampledPoints = resampledPolyData.GetPoints()

   # 5) Project each resampled point onto the mesh surface
    resampledPoints = resampledPolyData.GetPoints()
    projectedPoints = vtk.vtkPoints()

    implicit_distance = vtk.vtkImplicitPolyDataDistance()
    implicit_distance.SetInput(mesh_polydata)

    for i in range(resampledPoints.GetNumberOfPoints()):
        p = resampledPoints.GetPoint(i)
        closest = [0.0, 0.0, 0.0]
        implicit_distance.EvaluateFunctionAndGetClosestPoint(p, closest)

        projectedPoints.InsertNextPoint(closest)

    return projectedPoints


def get_vtu_points(filepath):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(filepath)
    reader.Update()
    points = reader.GetOutput().GetPoints()
    return points


def extract_positions(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)

    positions = []

    def find_positions(obj):
        if isinstance(obj, dict):
            for key, value in obj.items():
                if key == "position":
                    positions.append(value)
                else:
                    find_positions(value)
        elif isinstance(obj, list):
            for item in obj:
                find_positions(item)

    find_positions(data)

    return np.array(positions)


def vtk_points_to_numpy(vtk_points):
    vtk_array = vtk_points.GetData()   # vtkDataArray
    np_array = vtk_to_numpy(vtk_array) # convert to NumPy
    return np_array


def example():
    mesh_outputs_folder = os.path.join(os.path.dirname(__file__), "outputs", "mesh_outputs")
    os.makedirs(mesh_outputs_folder, exist_ok=True)
    
    inputs_folder = os.path.join(os.path.dirname(__file__), "inputs")
    curve_points_path = os.path.join(inputs_folder, "contour_points.vtu")
    print(curve_points_path)
    
   
    initialCurvePoints = np.array(get_vtu_points(curve_points_path).GetData()).tolist()
    numberOfCurveLandmarkPoints = 600

    lumen_mesh_path = os.path.join("lumen.vtp")
    lumen_pd = read_vtp(lumen_mesh_path)

    resampled_bdy_pts = resamplePointsOnMesh(initialCurvePoints, numberOfCurveLandmarkPoints, lumen_pd)
    resampled_bdy_pts_pd = convertPointsToVTKPolyData(resampled_bdy_pts)

    write_vtk_polydata(resampled_bdy_pts_pd, os.path.join(mesh_outputs_folder, "target_curve.vtp"))
   

if __name__ == "__main__":

    historical_3 = True 
    if historical_3: 

        mesh_names = ['1_aorta_native_original_cropped_cap_pt25.vtp', 
                      '6_aorta_1pt98_original_cropped_pt25.vtp',
                      '8_aorta_1pt7_original_cropped_pt25.vtp']

        mesh_suffixes = ['_native', 
                      '_1pt98',
                      '_1pt7']

        data = loadmat("ring_points_hist_3.mat")

        names = ['ring_pts_LR_cusp_sim_coords_mm', 'ring_pts_Non_cusp_sim_coords_mm']

        arrays_processed = {}


        for lumen_mesh, mesh_suffix in zip(mesh_names, mesh_suffixes):

            lumen_mesh_path = os.path.join(lumen_mesh)
            lumen_pd = read_vtp(lumen_mesh_path)

            for basename in names:

                initialCurvePoints = data[basename]

                initialCurvePoints_tuples = [tuple(row) for row in initialCurvePoints]

                print("initialCurvePoints = ", initialCurvePoints)
                print("initialCurvePoints_tuples = ", initialCurvePoints_tuples)

                initialCurvePoints_pd = convertPointsToVTKPolyData(initialCurvePoints_tuples)
                write_vtk_polydata(initialCurvePoints_pd, os.path.join(basename + "_initial.vtp"))
       

                numberOfCurveLandmarkPoints = 50

                resampled_bdy_pts = resamplePointsOnMesh(initialCurvePoints, numberOfCurveLandmarkPoints, lumen_pd)
                resampled_bdy_pts_pd = convertPointsToVTKPolyData(resampled_bdy_pts)

                # print(resampled_bdy_pts)

                write_vtk_polydata(resampled_bdy_pts_pd, os.path.join(basename + "_to_mesh.vtp"))

       
                print(f"{basename} = 0.1 * [")

                resampled_bdy_pts_numpy = vtk_points_to_numpy(resampled_bdy_pts)

                for row in vtk_points_to_numpy(resampled_bdy_pts):
                    print(",".join(map(str, row)))

                print("]';\n")

                arrays_processed[basename + mesh_suffix] = resampled_bdy_pts_numpy

        # save all 
        for name in arrays_processed:
            arrays_processed[name] = arrays_processed[name].astype(np.float64)

        savemat("ring_points_hist_3_processed.mat", arrays_processed)


    ross_1 = False 
    if ross_1: 
        lumen_mesh_path = os.path.join("meshes_stl/lv_aorta_ross_1.vtp")
        lumen_pd = read_vtp(lumen_mesh_path)

        for filename in os.listdir("."):
            if filename.startswith("aortic_annulus") and filename.endswith(".json"):

                basename = filename.split('.')[0]

                file_path = os.path.join(".", filename)

                initialCurvePoints = extract_positions(file_path) 
                # print(initialCurvePoints)

                initialCurvePoints_tuples = [tuple(row) for row in initialCurvePoints]
                initialCurvePoints_pd = convertPointsToVTKPolyData(initialCurvePoints_tuples)
                write_vtk_polydata(initialCurvePoints_pd, os.path.join(basename + "_initial.vtp"))
       

                numberOfCurveLandmarkPoints = 50

                resampled_bdy_pts = resamplePointsOnMesh(initialCurvePoints, numberOfCurveLandmarkPoints, lumen_pd)
                resampled_bdy_pts_pd = convertPointsToVTKPolyData(resampled_bdy_pts)

                # print(resampled_bdy_pts)

                write_vtk_polydata(resampled_bdy_pts_pd, os.path.join(basename + "_to_mesh.vtp"))

       
                print(f"{basename} = 0.1 * [")

                resampled_bdy_pts_numpy = vtk_points_to_numpy(resampled_bdy_pts)

                for row in vtk_points_to_numpy(resampled_bdy_pts):
                    print(",".join(map(str, row)))

                print("]';\n")



