import vtk
from vtk import *

 #write normals
points = vtk.vtkPoints()
verts = vtk.vtkCellArray()
polydata_pts = vtk.vtkPolyData()
    
pointNormalsArray = vtk.vtkDoubleArray()
pointNormalsArray.SetNumberOfComponents(3)

for i in range(clipSeedIds.GetNumberOfIds()):

    seedId = clipSeedIds.GetId(i)

    locator = vtk.vtkPointLocator()
    locator.SetDataSet(clippedSurface)
    locator.BuildLocator()

    seedPoint = self.Surface.GetPoint(seedId)
    seedPointId = locator.FindClosestPoint(seedPoint)

    planeEstimator = vtkvmtk.vtkvmtkPolyDataNormalPlaneEstimator()
    planeEstimator.SetInputData(clippedSurface)
    planeEstimator.SetOriginPointId(seedPointId)
    planeEstimator.Update()

    plane = vtk.vtkPlane()
    plane.SetOrigin(planeEstimator.GetOrigin())
    plane.SetNormal(planeEstimator.GetNormal())

    #testing the normal vector
    
    origin_pt = planeEstimator.GetOrigin()
    normal_pt = planeEstimator.GetNormal()
    
    id = points.InsertNextPoint(origin_pt)
    verts.InsertNextCell(1) 
    verts.InsertCellPoint(id)
    pointNormalsArray.InsertNextTuple(normal_pt)
    


#Add the points to the polydata container
polydata_pts.SetPoints(points)
polydata_pts.SetVerts(verts)

# Add the normals to the points in the polydata
#polydata_pts.GetCellData().SetNormals(pointNormalsArray)
polydata_pts.GetPointData().SetNormals(pointNormalsArray)
polydata_pts.Modified()


writer = vtk.vtkXMLPolyDataWriter()
writer.SetDataModeToAscii()
writer.SetFileName("UnstructuredPoints.vtp")
writer.SetInputData(polydata_pts)
writer.Write()
