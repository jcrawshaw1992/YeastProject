import subprocess
import vtk
import shutil
import os
from xml.etree import ElementTree
import glob
import time
import pdb
import string
import math
import sys


import pdb
# import vmtk
# from vmtk import pypes

from vtk import vtkXMLUnstructuredGridWriter as Writer, VTK_VERSION
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
import numpy
import meshio



def vtuTostl(vtuFile, stlfile):
# _______________________
#  Convert vtu to the stl i need 
#  _______________________
    print "  Convert vtu to stl    "
        # Read the VTU file from disk
    vtu_reader = vtk.vtkXMLUnstructuredGridReader()
    vtu_reader.SetFileName(vtuFile)

    extract_surface_filter = vtk.vtkDataSetSurfaceFilter()
    extract_surface_filter.AddInputConnection(vtu_reader.GetOutputPort())

    # Write out the data in unstructured grid format
    stl_writer = vtk.vtkSTLWriter()
    stl_writer.SetInputConnection(extract_surface_filter.GetOutputPort())
    stl_writer.SetFileName(stlfile)
    stl_writer.Write()
    print "  Doneies    "
    return None



def vtkTovtu(vtkFile, VTUfile):
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


def stlTovtu(stlFile, VTUfile):
    


    mesh = meshio.read(
        stlFile,  # string, os.PathLike, or a buffer/open file
        file_format="stl"  # optional if filename is a path; inferred from extension
    )
    # mesh.points, mesh.cells, mesh.cells_dict, ...
    # now figure out how to write, shouldnt be too hard
    reader = vtk.vtkSTLReader()
    reader.SetFileName(stlFile)
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


    #     print 'Finished'

if __name__=="__main__":
    
    Directory = "/Users/jcrawshaw/Downloads/CollapseWithAngleVariation/"
    # # # # # Set up the points for the centerlines and write into a file to be read in cpp 

    # Collapse = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4]#, 0.3, 0.2, 0.1]
    # for i in Collapse: 
    #     vtuFile=  Directory+"Mesh"+str(int(10*i))+".vtk"
    #     Clipped_Mesh = Directory+"MeshClipped"+str(int(10*i))+".vtk"

    #     vtkTovtu(Clipped_Mesh, vtuFile)

    # print "Done :) "



    # stlTovtu(Directory+"mesh.stl", Directory+"converted.vtu")
    vtkTovtu(Directory+"MeshClipped_A5_C0.55.vtk", Directory+"converted2.vtu")



