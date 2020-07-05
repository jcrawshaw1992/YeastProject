import vtk
import pdb

# Generate two sets of points for this example.
points0 = vtk.vtkPoints()
points1 = vtk.vtkPoints()

points0.InsertNextPoint(1., 0., 0.)
points0.InsertNextPoint(1., 1., 0.)

points1.InsertNextPoint(1., 0., 1.)
points1.InsertNextPoint(1., 1., 1.)

polydata10 = vtk.vtkPolyData()
polydata11 = vtk.vtkPolyData()
polydata10.SetPoints(points0)
polydata11.SetPoints(points1)
#-------------------------------------

# Create a set of points that joints the two sets of points from polydata10 and polydata11.

pointsJoined = vtk.vtkPoints()

for i in range(polydata10.GetNumberOfPoints()):
    pointsJoined.InsertNextPoint(polydata10.GetPoint(i))

for i in range(polydata11.GetNumberOfPoints()):
    pointsJoined.InsertNextPoint(polydata11.GetPoint(i))

# Initialize a polydata object and set the joint point set.
polydata = vtk.vtkPolyData()
polydata.SetPoints(pointsJoined)

# Create verticies so the points can be visualized
verts = vtk.vtkCellArray()
for i in range(polydata.GetNumberOfPoints()):
    verts.InsertNextCell(1) # Create a 1 dimensional cell
    verts.InsertCellPoint(i) # Append point i to the verts array
polydata.SetVerts(verts)

polydata.Modified()
# polydata.GetPointData().SetScalars(Radii)
writer = vtk.vtkXMLPolyDataWriter()
# writer.SetDataModeToBinary()
writer.SetDataModeToAscii()
writer.SetFileName("JointedPoints.vtp")
writer.SetInputData(polydata)

writer.Write()

