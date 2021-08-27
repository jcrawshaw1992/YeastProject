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

    file_basename = "/Volumes/Hardrive/Projects/MeshCollection/VascularNetworks/Network/VascularNetwork.vtk"
    Save_name = "/Volumes/Hardrive/Projects/MeshCollection/VascularNetworks/Network/ClippedVascularNetwork.vtk" 


    
    Origin = (10.920486962892474, 116.3874796679927, -0.6276400308686206 )
    Normal = (0.9477460065412346, -0.31902584258272776,0.000137293563571154 )
    Origin2nd = (21.813720423837108, 130.53001283165923, -2.525749995223146)
    Normal2nd = (0,1,0)

    Plane = getBoundedPlaneFunction(Origin, Normal, Origin2nd, Normal2nd)
    doClipCylinderSimple(file_basename, Plane, Save_name)


    Origin = (149.6697601199692,22.71236044521207, -0.7836004395217067  )
    Normal = (0.08436010716094508, 0.9964010238344646, -0.008268737595575627 )

    Origin2nd = (134.34712340504694,10.19651459001711, -15.157072663314448 )
    Normal2nd = (-1,0,0)

    Plane = getBoundedPlaneFunction(Origin, Normal, Origin2nd, Normal2nd)
    doClipCylinderSimple(Save_name, Plane, Save_name)

   clip_surface_with_plane(Save_name, (37.076811980226985, 214.0089896119276, -0.6576238345166555 ), (0.7917293866548619,  -0.610797056164505, 0.009567365775952187), Save_name)




# #  cd  ~/Documents/CGAL-5.0.2/examples/Polygon_mesh_processing/;

# # ./isotropic_remeshing_ForChaste -input /Volumes/Hardrive/Projects/MeshCollection/VascularNetworks/Network/ClippedNetwork.off -output /Volumes/Hardrive/Projects/MeshCollection/VascularNetworks/Network/Remeshed.off -target_edge_length 1000 -iterations 5

# ./isotropic_remeshing_ForChaste -input /Volumes/Hardrive/Projects/MeshCollection/Plexus/Mesh.off -output /Volumes/Hardrive/Projects/MeshCollection/Plexus/Remesh.off -target_edge_length 5 -iterations 5