#!/usr/bin/env python
#
import subprocess
import vtk
# import shutil
import os
# from xml.etree import ElementTree
# import glob
from argparse import ArgumentParser
# import time
# import matplotlib.pyplot as plt
# import csv
import pdb
import string
import math
import vmtk
from vmtk import pypes

from vtk import vtkXMLUnstructuredGridWriter as Writer, VTK_VERSION
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
import numpy
import meshio


from FileConverter import vtuTostl


    # ---------------------------------
  
def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")


if __name__=="__main__":
   

    OriginalGeometry = "/Users/jcrawshaw/Documents/testoutput/CylinderCollection/BunchOfCylinder/40/Deformed/results_from_time_0/mesh_0.vtu"# "/Users/jcrawshaw/Documents/testoutput/CylinderCollection/TestRemeshDeformingCylinder/InitialMesh/results_from_time_0/mesh_0.vtu"
    RemeshedGeometry =  "/Users/jcrawshaw/Documents/testoutput/CylinderCollection/BunchOfCylinder/40/Deformed/Remeshed.stl"# "/Users/jcrawshaw/Documents/testoutput/CylinderCollection/TestRemeshDeformingCylinder/InitialMesh/Remeshed.stl"
    OriginalInSTLFormat = "/Users/jcrawshaw/Documents/testoutput/CylinderCollection/BunchOfCylinder/40/Deformed/Deformed.stl"# "/Users/jcrawshaw/Documents/testoutput/CylinderCollection/TestRemeshDeformingCylinder/InitialMesh/Deformed.stl"
   
    RemeshedGeometry_0 =  "/Users/jcrawshaw/Documents/testoutput/CylinderCollection/TestRemeshDeformingCylinder/InitialMesh/Remeshed_0.5.stl"
    RemeshedGeometry_1 =  "/Users/jcrawshaw/Documents/testoutput/CylinderCollection/TestRemeshDeformingCylinder/InitialMesh/Remeshed_0.05.stl"
    RemeshedGeometry_2 =  "/Users/jcrawshaw/Documents/testoutput/CylinderCollection/TestRemeshDeformingCylinder/InitialMesh/Remeshed_0.005.stl"
    RemeshedGeometry_3 =  "/Users/jcrawshaw/Documents/testoutput/CylinderCollection/TestRemeshDeformingCylinder/InitialMesh/Remeshed_0.1.stl"
    
    ifile = vtuTostl(OriginalGeometry, OriginalInSTLFormat)

    command = 'vmtksurfaceremeshing -ifile '+OriginalInSTLFormat +' -iterations 50 -area 0.5 -ofile ' +RemeshedGeometry
    subprocess.call(command, shell=True)


    print "-------------------------------------"
    print "------ Remeshing Complete  ----------"
    print "-------------------------------------"
    print "-------------  Finito  --------------"
    print "-------------------------------------"






  
  
