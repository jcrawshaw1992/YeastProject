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
import vmtk
from vmtk import pypes

from vtk import vtkXMLUnstructuredGridWriter as Writer, VTK_VERSION
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
import numpy
import meshio

def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")



def convertFile(vtkFile, VTUfile):
    if os.path.isfile(vtkFile):
        basename = os.path.basename(vtkFile)
        print "converting VTK to VTU"

        reader = vtk.vtkPolyDataReader() 
   
        reader.SetFileName(vtkFile)
        reader.Update()
        output = reader.GetOutput()

        Nodes = []
        ElementList = []
        # Add the first Node so this isnt an emtpy list
        Nodes.append(output.GetCell(0).GetPoints().GetPoint(0))
        
        # loop over all of the elements in the vtk polydata file 
        # For each element, the polydata has recorded the location of each node,
        # but not the node index. As such here we loop over the elements 
        # and select out the node, which are put in the Node list. 
        # As I loop over the elements I will also record which node indices are in each element
        for i in range(output.GetNumberOfCells()):
            Element = []
            pts = output.GetCell(i).GetPoints()  
            np_pts = numpy.array([pts.GetPoint(i) for i in range(pts.GetNumberOfPoints())]) 
            
            # Saving each of the nodes in this element
            for j in [0,1,2]:
                AddNodeToList = 1 
                # want to add the Node to the Nodes list, But need to check if it is in the list yet
                # Do this by looping over list and seeing if it is in there
                for k in Nodes: 
                    if(k == pts.GetPoint(j)):
                        AddNodeToList =0
                if (AddNodeToList ==1):
                    # print 'Adding Node'
                        # Filling in the array storing Node locations
                    Nodes.append(pts.GetPoint(j))
                NodeIndex = Nodes.index(pts.GetPoint(j))
                Element.append(NodeIndex)
                # This works to fill in Nodes
            ElementList.append(Element)
        # pdb.set_trace()

        # Save the nodes as points and elements for the meshio writer
        points = numpy.array(Nodes)
        elements = {
        "triangle": numpy.array(ElementList
        )
        }    

        meshio.write_points_cells(
        VTUfile,
        points,
        elements,
        # Optionally provide extra data on points, cells, etc.
        # point_data=point_data,
        # cell_data=cell_data,
        # field_data=field_data
        )
        print 'Finished'



def clip_surface_with_plane(vtkFile, plane_origin, plane_normal, save_name):

    reader = vtk.vtkPolyDataReader() 

    reader.SetFileName(vtkFile)

    plane = vtk.vtkPlane()
    
    plane.SetOrigin(*plane_origin)
    plane.SetNormal(*plane_normal)


    clipper = vtk.vtkClipPolyData()
    
    clipper.SetInputConnection(reader.GetOutputPort())
    
    clipper.SetClipFunction(plane)
    clipper.GenerateClipScalarsOn()
    clipper.GenerateClippedOutputOn()
    clipper.SetValue(0)
    

    writer = vtk.vtkPolyDataWriter()
    
    writer.SetFileName(save_name)
    
    writer.SetInputConnection(clipper.GetOutputPort())
    
    writer.Write()


file_basename = '/Users/jcrawshaw/docker-polnet-master/PlexusNotScalled/PlexusMesh.vtk'
Save_name = "/Users/jcrawshaw/docker-polnet-master/PlexusNotScalled/PlexusMesh_Clipped.vtk"
ScalledMeshVTU = "/Users/jcrawshaw/docker-polnet-master/PlexusNotScalled/Plexus.vtu"

#clip_surface_with_plane(file_basename, (500, 400, 10), (1, 0, 0))
#clip_surface_with_plane(file_basename, (2000, 1200, 10), (-1, 0, 0))
#clip_surface_with_plane(file_basename, (472, 420, 0), (0, -11, 0))
clip_surface_with_plane(file_basename, (584.0682462724309, 1220.205,-46.025), (0.671,-0.7325,0.1142), Save_name)


clip_surface_with_plane(Save_name, (188.38,740.09,23.53), (0.999,-0.0158,0.015), Save_name)


clip_surface_with_plane(Save_name, (490.7664298534873, 350.7211701339496, 6.666), ( 0.7067859059916219, 0.7065452213525341, -0.0353204370796983), Save_name)
clip_surface_with_plane(Save_name, (739.05,181.48,3.609), ( 0.259,0.96,-0.065), Save_name)
clip_surface_with_plane(Save_name, (1234.1761995532954, 231.67,11.24), ( -0.762,0.646,0.015), Save_name)
clip_surface_with_plane(Save_name, (1261.6921084158262, 875.832, 23.36), (-0.74,-0.67,-0.0219 ), Save_name)

clip_surface_with_plane(Save_name, (1246.6,1177.32, -2.20), (-0.86, -0.505, -0.0042 ), Save_name)

print "B"

convertFile(Save_name, ScalledMeshVTU)



        
    



