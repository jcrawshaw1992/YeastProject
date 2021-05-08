
# ReadVtu.py
import meshio



# mesh = meshio.read(file)
# # mesh.points, mesh.cells, mesh.point_data, ...


import numpy
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as VN


file = "/Volumes/Hardrive/Projects/FSI/UpperBranchCollapse/CollectedResults/surface-traction_5.5.vtu"


reader = vtkUnstructuredGridReader()
reader.SetFileName(file)
reader.ReadAllVectorsOn()
reader.ReadAllScalarsOn()
reader.Update()

# data = reader.GetOutput()
# potential = data.GetPointData().GetScalars("potential")

# print(type(potential))
