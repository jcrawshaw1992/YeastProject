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
   
    MeshSize =  [20,30,40,50,55,60,65,70,75,80,85]
    # Areas = [0.4, 0.1,0.05, 0.04,0.03,0.028,0.024, 0.02,0.015,0.013, 0.0115, 0.01, 0.009,0.0085, 0.008 , 0.0075,0.007,0.006]
    Areas = [0.0015] #[0.04,0.03, 0.02,0.015,0.013, 0.0115, 0.01, 0.009,0.0085, 0.008 , 0.007,0.006,0.005,0.004,0.003,0.001 , 0.0027,0.0024,0.002]


# % Need Refinment to go to 5000

    for i in MeshSize:
        for j in Areas :

            InitalMesh =  "/Users/jcrawshaw/Documents/testoutput/CylinderCollection/BunchOfCylinder/"+ str(i)+ "/Inital/results_from_time_0/mesh_0.vtu"
            DeformedMesh = "/Users/jcrawshaw/Documents/testoutput/CylinderCollection/BunchOfCylinder/" +str(i)+ "/Deformed/results_from_time_0/mesh_0.vtu"


            InitalMesh_stl =  "/Users/jcrawshaw/Documents/testoutput/CylinderCollection/BunchOfCylinder/"+ str(i)+ "/Inital.stl"
            DeformedMesh_stl = "/Users/jcrawshaw/Documents/testoutput/CylinderCollection/BunchOfCylinder/" +str(i)+ "/Deformed.stl"

            RemeshedMesh = "/Users/jcrawshaw/Documents/testoutput/CylinderCollection/BunchOfCylinder/" +str(i)+ "/Remeshed"+str(j) + ".stl"
            vtuTostl(InitalMesh, InitalMesh_stl)
            vtuTostl(DeformedMesh, DeformedMesh_stl)
            command = "vmtksurfaceremeshing -ifile "+DeformedMesh_stl+" -iterations 30 -area " + str(j)+  " -ofile " + RemeshedMesh
            subprocess.call(command, shell=True)


    print "-------------------------------------"
    print "------ Remeshing Complete  ----------"
    print "-------------------------------------"
    print "-------------  Finito  --------------"
    print "-------------------------------------"




  
  
