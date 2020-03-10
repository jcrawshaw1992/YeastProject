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
sys.path.insert(1, '/Applications/ParaView-5.4.1-822-g597adef982.app/Contents/Python/')
from paraview.simple import *



hemelb_setup_exe = 'env PYTHONPATH=/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools:/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool:/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters:$PYTHONPATH /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool/scripts/hemelb-setup-nogui'
hemelb_setup_exe1 = 'env PYTHONPATH=/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool:$PYTHONPATH /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool/scripts/hemelb-setup-nogui'

# _______________________
#  Convert vtu to the stl i need 
#  _______________________


inputFile = "//Users/jcrawshaw/Documents/ChasteWorkingDirectory/ShrunkPlexus/SetUpData/PlexusMesh.vtu"
outputFile = "//Users/jcrawshaw/Documents/ChasteWorkingDirectory/ShrunkPlexus/SetUpData/config.stl"
print ': Converting ', inputFile, '  ->  ', outputFile
reader = LegacyVTKReader( FileNames= inputFile )
writer = XMLUnstructuredGridWriter()
writer.FileName = outputFile
writer.UpdatePipeline()
Delete(reader)
Delete(writer)

pdb.set_trace()

working_directory = '//Users/jcrawshaw/Documents/ChasteWorkingDirectory/ShrunkPlexus/SetUpData/'
print "  Convert vtu to stl    "
    # Read the VTU file from disk
vtu_reader = vtk.vtkXMLUnstructuredGridReader()
vtu_reader.SetFileName(working_directory + 'ScalledMesh.vtu')

extract_surface_filter = vtk.vtkDataSetSurfaceFilter()
extract_surface_filter.AddInputConnection(vtu_reader.GetOutputPort())

# Write out the data in unstructured grid format
stl_writer = vtk.vtkSTLWriter()
stl_writer.SetInputConnection(extract_surface_filter.GetOutputPort())
stl_writer.SetFileName(working_directory + 'ScalledMesh.stl')
stl_writer.Write()
print "  Doneies    "


