import subprocess
import vtk
import shutil
import os
from xml.etree import ElementTree
import glob
from argparse import ArgumentParser
import numpy as np
import time
import matplotlib.pyplot as plt
import csv
import pdb
import string
import math


# working_directory = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory/SmallCylinders/Radius1/'

hemelb_setup_exe = 'env PYTHONPATH=/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools:/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool:/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters:$PYTHONPATH /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool/scripts/hemelb-setup-nogui'
hemelb_setup_exe1 = 'env PYTHONPATH=/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool:$PYTHONPATH /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool/scripts/hemelb-setup-nogui'

# _______________________
#  COnvert vtu to the stl i need 
#  _______________________



working_directory = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory/BifurcationScalled/SetUpData/'
print "  Convert vtu to stl    "
    # Read the VTU file from disk
vtu_reader = vtk.vtkXMLUnstructuredGridReader()
vtu_reader.SetFileName(working_directory + 'configChaste.vtu')

extract_surface_filter = vtk.vtkDataSetSurfaceFilter()
extract_surface_filter.AddInputConnection(vtu_reader.GetOutputPort())

# Write out the data in unstructured grid format
stl_writer = vtk.vtkSTLWriter()
stl_writer.SetInputConnection(extract_surface_filter.GetOutputPort())
stl_writer.SetFileName(working_directory + 'configChaste.stl')
stl_writer.Write()
print "  Doneies    "



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