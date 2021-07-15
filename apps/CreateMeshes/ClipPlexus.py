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
    file_basename = '/Users/jcrawshaw/Downloads/Plexus.vtk'
    Save_name ='/Users/jcrawshaw/Downloads/PlexusClipped.vtk'

    Origin = (0.04445307579525465/0.00006684491*1.29, 0.05162259742566157/0.00006684491*1.29, -0.001990588305433403/0.00006684491*1.29)
    Normal = (-0.8105631234548122,-0.5853191641674113,0.019720521127546863)

    Origin2nd = (0.05/0.00006684491*1.29,  0.04/0.00006684491*1.29, -0.0006/0.00006684491*1.29)
    Normal2nd = (0.6382554028840824, -0.7675489817154649, 0.059148975957581465)
    

    plane_origin = (0.04127266890226009/0.00006684491*1.29,00.01885961520968956/0.00006684491*1.29, 0.0007232045314343617/0.00006684491*1.29)
    plane_normal = (0.057879903604992705,0.997388116052124,0.04320720676775996)

 

    Plane = getBoundedPlaneFunction(Origin, Normal, Origin2nd, Normal2nd)
    doClipCylinderSimple(file_basename, Plane, Save_name)
    clip_surface_with_plane(Save_name, plane_origin, plane_normal, Save_name)

    # clip_surface_with_plane(Save_name, (0.045098203892165734,0.05064189958380424,0.0011589943289266515), (0.8081276713886205,-0.5856195847394546,-0.06308223763940606), Save_name)
    # clip_surface_with_plane(Save_name, (0.03656403371703994,0.049415287040463565,0.0005315046730622214), (0.6955018850657018,-0.7184806801140912,0.007914554494637781), Save_name)
    clip_surface_with_plane(Save_name, (0.028133553209430767,0.04062601993172238,0.0010395808490631854), (0.983017698694772,-0.1828596366438072,-0.015445301529249871), Save_name)

    clip_surface_with_plane(Save_name, (0.03817188881578499, 0.019092013302893857, -0.0003001335411293224), (0.03619867499024044,0.9989466830122159,-0.02819894373), Save_name)
    # clip_surface_with_plane(Save_name, (1162.8189537477685, 573.6140054014987, -1.5705496953176452), (-0.6612833830623196,0.7436053470763407, 0.09876930234201606), Save_name)
    clip_surface_with_plane(Save_name, (524.0685139327954, 701.4029880692041, 2.251895091415726), (0.9976171514538998, -0.05816030877723643, 0.03711330769337815), Save_name)
    



    # clip_surface_with_plane(Save_name, (0.047632004415947796, 0.018907967775871823, 0.0015750103541767799), (-0.6183178411933676,0.7845743897013229,0.04610937309019025), Save_name)
    # clip_surface_with_plane(Save_name, (0.06539797148540963, 0.025939351832405394, 0.0038007231009339957), (-0.8188719509487834,-0.5731501839108796,0.030783025066716225), Save_name)
    clip_surface_with_plane(Save_name, (0.027592706932313197, 0.03520549399450725, 0.0007548304121448202), (0.9930115746535606,-0.1034408880150438, -0.05681544940168478), Save_name)
    # 
    # print "B"

    # # convertFile(Save_name, ScalledMeshVTU)




#  cd  ~/Documents/CGAL-5.0.2/examples/Polygon_mesh_processing/;

# ./isotropic_remeshing_ForChaste -input /Users/jcrawshaw/Downloads/PlexusClipped.off -output /Users/jcrawshaw/Downloads/PlexusRemeshed.off -target_edge_length 10 -iterations 5