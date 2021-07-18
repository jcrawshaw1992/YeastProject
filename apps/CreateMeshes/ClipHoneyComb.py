import vtk
from vtk.util import numpy_support
import vtk
import numpy as np
from subprocess import call


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
# import vmtk
# from vmtk import pypes

from vtk import vtkXMLUnstructuredGridWriter as Writer, VTK_VERSION
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
import numpy
import meshio

def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")


def clip_surface_with_plane(vtkFile, plane_origin, plane_normal, save_name):

    reader = vtk.vtkPolyDataReader() 

    reader.SetFileName(vtkFile)
    

    plane = vtk.vtkPlane()
    
    plane.SetOrigin(*plane_origin)
    plane.SetNormal(*plane_normal)


    clipper = vtk.vtkClipPolyData()
    
    clipper.SetInputConnection(reader.GetOutputPort())

    clipper.SetClipFunction(plane)
    clipper.SetValue(0.1)
    clipper.GenerateClipScalarsOn()
    clipper.GenerateClippedOutputOn()
    clipper.SetValue(0)
    

    writer = vtk.vtkPolyDataWriter()
    
    writer.SetFileName(save_name)
    
    writer.SetInputConnection(clipper.GetOutputPort())
    
    writer.Write()



def getBoundedCylinderFunction(center, axis, radius, centerTop,centerBottom):
   # Create cut cylinder. The vtkCylinder is just an infinite cylinder. To
   # add two ends, combine an implicit cylinder (vtkCylinder) with two implicit
   # plane functions (vtkPlane) through vtkImplicitBoolean.
   
   cylinderImplicit = vtk.vtkCylinder()
   cylinderImplicit.SetCenter( center)
   cylinderImplicit.SetAxis(axis)
   cylinderImplicit.SetRadius(radius)

   print center, axis,radius
#    A = (1,0,0)*axis
#    print A
#    centerTop = center + height/2.*axis
#    centerBottom = center + height/2.*axis
   print centerBottom, centerTop

   plane1 = vtk.vtkPlane()
   plane1.SetOrigin(centerTop)
   plane1.SetNormal(np.array(cylinderImplicit.GetAxis()))
   plane2 = vtk.vtkPlane()
   plane2.SetOrigin(centerBottom)
   plane2.SetNormal(-np.array(cylinderImplicit.GetAxis()))

   implicitBoolean = vtk.vtkImplicitBoolean()
   implicitBoolean.AddFunction(cylinderImplicit)
   implicitBoolean.AddFunction(plane1)
   implicitBoolean.AddFunction(plane2)
   implicitBoolean.SetOperationTypeToIntersection()

   return implicitBoolean


def getBoundedPlaneFunction(plane_origin, plane_normal, plane_origin2nd, plane_normal2nd):
   # Create cut cylinder. The vtkCylinder is just an infinite cylinder. To
   # add two ends, combine an implicit cylinder (vtkCylinder) with two implicit
   # plane functions (vtkPlane) through vtkImplicitBoolean.
   
    plane = vtk.vtkPlane()
    plane.SetOrigin(*plane_origin)
    plane.SetNormal(*plane_normal)

    plane1 = vtk.vtkPlane()
    plane1.SetOrigin(*plane_origin2nd)
    plane1.SetNormal(*plane_normal2nd)
  
    implicitBoolean = vtk.vtkImplicitBoolean()
    implicitBoolean.AddFunction(plane)
    implicitBoolean.AddFunction(plane1)
    implicitBoolean.SetOperationTypeToIntersection()

    return implicitBoolean




def doClipCylinderSimple(vtkFile, cylinderImplicit, Save_name):
    # Apply cut using vtkClipPolyData.

    reader = vtk.vtkPolyDataReader() 

    reader.SetFileName(vtkFile)
    clipper = vtk.vtkClipPolyData()
    
    clipper.SetInputConnection(reader.GetOutputPort())

    # clipper.SetClipFunction(plane)
    clipper.SetClipFunction(cylinderImplicit)
    # clipper.SetValue(0.1)
    clipper.GenerateClipScalarsOn()
    clipper.GenerateClippedOutputOn()
    clipper.SetValue(0)
    

    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(Save_name)
    
    writer.SetInputConnection(clipper.GetOutputPort())
    
    writer.Write()


if __name__=="__main__":
    print "B"
    file_basename = '/Users/jcrawshaw/Downloads/InitialHoneycombNetworkForFSI/Mesh_Scaled.vtk'
    Save_name ='/Users/jcrawshaw/Downloads/InitialHoneycombNetworkForFSI/MeshClipped.vtk'



    plane_origin = (0.13303499668450677 , 0,0)
    plane_normal = (1,0,0)
    clip_surface_with_plane(file_basename, plane_origin, plane_normal, Save_name)

    plane_origin = (0.24804051175590053 , 0.05852047082094186,0)
    plane_normal = (-0.7020425109198389,-0.7020425109198389,0)
    clip_surface_with_plane(Save_name, plane_origin, plane_normal, Save_name)

    plane_origin = (0.24804051175590053 , -0.05852047082094186,0)
    plane_normal = (-0.7020425109198389,0.7020425109198389,0)
    clip_surface_with_plane(Save_name, plane_origin, plane_normal, Save_name)


 




#  cd  ~/Documents/CGAL-5.0.2/examples/Polygon_mesh_processing/;

# ./isotropic_remeshing_ForChaste -input /Users/jcrawshaw/Downloads/InitialHoneycombNetworkForFSI/Bifucation.off -output /Users/jcrawshaw/Downloads/InitialHoneycombNetworkForFSI/BifucationRemeshed.off -target_edge_length 0.002 -iterations 5