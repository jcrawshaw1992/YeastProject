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
import MIRTK
# import vtkCommonCorePython
# # export PYTHONPATH="$PYTHONPATH:usr/local/lib/python2.7/site-packages/vtk/"
# sys.path.insert(1,"/usr/local/Cellar/vtk/8.2.0_11/lib/python3.8/site-packages/")
# sys.path.insert(1,"/usr/local/Cellar/vtk/8.2.0_11/lib/python3.8/site-packages/vtkmodules/vtkCommonCorePython.so")
# sys.path.insert(1, '/Applications/ParaView-5.4.1-822-g597adef982.app/Contents/Python/')
# from paraview.simple import *


#PATH=
#  29 export PATH ="/Library/Frameworks/Python.framework/Versions/3.7/bin:${PATH}"


# # Directory = "/Users/jcrawshaw/Downloads/Nonsymmetric/"
# Directory = "/home/vascrem/MeshCollection/IdealisedNetwork/IdealMeshWIthCentralCollapse/Nonsymmetric/"
# vtuFile=  Directory+"Mesh7.vtk"
# Clipped_Mesh = Directory+"MeshClipped7.vtk"
     
# inputFile = Clipped_Mesh
# outputFile = vtuFile
# print ' Converting ', inputFile, '  ->  ', outputFile
# reader = LegacyVTKReader( FileNames= inputFile )
# writer = XMLUnstructuredGridWriter()
# writer.FileName = outputFile
# writer.UpdatePipeline()
# Delete(reader)
# Delete(writer)

# print "Done :) "


