import subprocess
import vtk
import shutil
import os
from xml.etree import ElementTree
import glob
import numpy as np
import time
import pdb
import string
import math
import sys


# working_directory = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory/SmallCylinders/Radius1/'

hemelb_setup_exe = 'env PYTHONPATH=/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools:/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool:/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters:$PYTHONPATH /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool/scripts/hemelb-setup-nogui'
hemelb_setup_exe1 = 'env PYTHONPATH=/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool:$PYTHONPATH /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool/scripts/hemelb-setup-nogui'

# _______________________
#  Convert vtu to the stl i need 
#  _______________________

vtu_reader = vtk.vtkXMLUnstructuredGridReader()
vtu_reader.SetFileName( '~/Plexus.vtu')

extract_surface_filter = vtk.vtkDataSetSurfaceFilter()
extract_surface_filter.AddInputConnection(vtu_reader.GetOutputPort())

# Write out the data in unstructured grid format
stl_writer = vtk.vtkSTLWriter()
stl_writer.SetInputConnection(extract_surface_filter.GetOutputPort())
stl_writer.SetFileName( '~/config.stl')
stl_writer.Write()
print "  Doneies    "

# inputFile = "//Users/jcrawshaw/Documents/ChasteWorkingDirectory/ShrunkPlexusWithLongInlets/SetUpData/PlexusMesh.vtu"
# outputFile = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/ShrunkPlexus/SetUpData/config.stl"

# working_directory = '/Users/jcrawshaw/docker-polnet-master/NewMesh/'
# mesh = meshio.read(working_directory + 'NewReMeshedPlexus.')
# print "  Convert stl to vtu    "
#     # Read the VTU file from disk
# stl_reader = vtk.vtkSTLUnstructuredGridReader()
# stl_reader.SetFileName(working_directory + 'NewReMeshedPlexus.stl')

# extract_surface_filter = vtk.vtkDataSetSurfaceFilter()
# extract_surface_filter.AddInputConnection(vtu_reader.GetOutputPort())

# # Write out the data in unstructured grid format
# stl_writer = vtk.vtkSTLWriter()
# stl_writer.SetInputConnection(extract_surface_filter.GetOutputPort())
# stl_writer.SetFileName(working_directory + 'ScalledMesh.stl')
# stl_writer.Write()
# print "  Doneies    "



# working_directory = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory/LargestCylinder/'

# print "  Convert vtu to stl    "
#     # Read the VTU file from disk
# vtu_reader = vtk.vtkXMLUnstructuredGridReader()
# vtu_reader.SetFileName(working_directory + 'config.vtu')

# extract_surface_filter = vtk.vtkDataSetSurfaceFilter()
# extract_surface_filter.AddInputConnection(vtu_reader.GetOutputPort())

# # Write out the data in unstructured grid format
# stl_writer = vtk.vtkSTLWriter()
# stl_writer.SetInputConnection(extract_surface_filter.GetOutputPort())
# stl_writer.SetFileName(working_directory + 'config.stl')
# stl_writer.Write()


# import numpy
# from vtk import vtkStructuredPointsReader
# from vtk.util import numpy_support as VN

# reader = vtkStructuredPointsReader()

# # working_directory = '/Users/jcrawshaw/docker-polnet-master/NewMesh/'
# # mesh = meshio.read(working_directory + 'NewReMeshedPlexus.')x
# reader.ReadAllVectorsOn()
# reader.ReadAllScalarsOn()
# reader.Update()

# data = reader.GetOutput()

# dim = data.GetDimensions()
# vec = list(dim)
# vec = [i-1 for i in dim]
# vec.append(3)

# u = VN.vtk_to_numpy(data.GetCellData().GetArray('velocity'))
# b = VN.vtk_to_numpy(data.GetCellData().GetArray('cell_centered_B'))

# u = u.reshape(vec,order='F')
# b = b.reshape(vec,order='F')

# x = zeros(data.GetNumberOfPoints())
# y = zeros(data.GetNumberOfPoints())
# z = zeros(data.GetNumberOfPoints())

# for i in range(data.GetNumberOfPoints()):
#         x[i],y[i],z[i] = data.GetPoint(i)

# x = x.reshape(dim,order='F')
# y = y.reshape(dim,order='F')
# z = z.reshape(dim,order='F')