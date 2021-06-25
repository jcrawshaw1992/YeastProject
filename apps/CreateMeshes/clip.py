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
    clipper.SetValue(0.001)
    clipper.GenerateClipScalarsOn()
    clipper.GenerateClippedOutputOn()
    clipper.SetValue(0)
    

    writer = vtk.vtkPolyDataWriter()
    
    writer.SetFileName(save_name)
    
    writer.SetInputConnection(clipper.GetOutputPort())
    
    writer.Write()





if __name__=="__main__":
    print "B"
    file_basename = '/Users/jcrawshaw/docker-polnet-master/Plexus/Plexus_.vtk'
    Save_name ='/Users/jcrawshaw/docker-polnet-master/Plexus/Plexus_Clipped.vtk'
    ScalledMeshVTU = '/Users/jcrawshaw/docker-polnet-master/Plexus/Plexus.vtu'


    clip_surface_with_plane(file_basename, (584.0682462724309, 1220.205,-46.025), (0.671,-0.7325,0.1142), Save_name)
    clip_surface_with_plane(Save_name, (724.0,1033, 7.07), ( 0.6,-0.77,0.021024), Save_name)
    clip_surface_with_plane(Save_name, (515,725.22,-3.9552), (1,0,0), Save_name) 
    clip_surface_with_plane(Save_name, (790.36,335.33,-9.04), (0.1633,0.9849,-0.0558), Save_name)
    clip_surface_with_plane(Save_name, (1074.24,387.48,7.1176), (-0.718,0.6957,-0.0060), Save_name)

    clip_surface_with_plane(Save_name, (566.662,419.26,-6.4678), (0.749,0.6611,0.034), Save_name)
    clip_surface_with_plane(Save_name, (923.2201,999.48,64.698), (-0.8,-0.59,0.066), Save_name) # This is the double inlet clip 


    clip_surface_with_plane(Save_name, (589,449,-10.92), (0.74,0.66,0.07), Save_name)


    print "B"

    convertFile(Save_name, ScalledMeshVTU)


