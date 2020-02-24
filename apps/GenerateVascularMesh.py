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


    # ---------------------------------
  
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



if __name__=="__main__":
    parser = ArgumentParser(description='Get the radius file ')
    
    parser.add_argument('--ifile', default='/Users/jcrawshaw/docker-polnet-master/Plexus_Course_2/', type=str, help='Need to supply a input folder')
    args = parser.parse_args()    
    # ' ------- Setting up args ------- '

    Directory  = args.ifile #+ 'results_from_time_0/'
    CenterLines_filename = args.ifile + 'PlexusCenterlines.vtp'
    SmootherCenterlinesFile = args.ifile + 'PlexusCenterlines_Smoothed.vtp'
    ScalledCenterLines = args.ifile + 'PlexusCenterlines_Smoothed.vtp' # args.ifile + 'PlexusCenterlines_Scalled.vtp'

    ScalledMeshVTK = args.ifile + 'InialMesh.vtk'
    ScalledMeshVTK_2 = args.ifile + 'PlexusMesh.vtk'
    ScalledMeshVTU = args.ifile + 'Plexus_2.vtu'
    
    
    
    AdaptiveRadius = 1
    ResampleCenterlinesFile =0
    Scalling = 10

    CenterlinesGenerationNeeded =0
    

    if CenterlinesGenerationNeeded:
        pdb.set_trace()
        print ' You should  not need this if you have a centerlines file, go get the centerlines file properly '
        # Generating a centerlines file form thse original mesh (the one from James and Miguel)
        GetCenterlinesCommand = 'vmtk vmtknetworkextraction -ifile ' + Directory + 'config.stl' + ' -ofile ' + CenterLines_filename+"NewlyGenerated"
        subprocess.call(GetCenterlinesCommand, shell=True)

    # Resamlpe the centerlines smooth centerlines with a moving average filter -->    vmtkcenterlineresampling -ifile PlexusInversed.vtp -length 20 -ofile resampledCenterlines.vtp  
    if ResampleCenterlinesFile:
        SmoothCenterlinesCommond = 'vmtkcenterlineresampling -ifile '+ CenterLines_filename + ' -ofile '+ SmootherCenterlinesFile
        subprocess.call(SmoothCenterlinesCommond, shell=True)

        SmoothCenterlinesCommond = 'vmtkcenterlineresampling -ifile '+ SmootherCenterlinesFile + ' -ofile '+ SmootherCenterlinesFile
        subprocess.call(SmoothCenterlinesCommond, shell=True)

  #     # Read in the centerlines file and edit the radius -> this is saved in a new file, which will be read in to generate the new mesh
    if AdaptiveRadius:
       
        reader = vtk.vtkXMLPolyDataReader()

        #Read in the centerlines file
        reader.SetFileName(CenterLines_filename)
        reader.Update()

        # Set up the vtk writer for the edited centerlines data 
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(ScalledCenterLines)
    
        point_data = reader.GetOutput().GetPointData()

        assert point_data.GetArrayName(0) == 'Radius', "VTP file doesn't contain array Radius"
        radii = point_data.GetArray(0)
        
        # loop over the radius and save scale
        for i in range(radii.GetSize()):
            radius = radii.GetValue(i)/2
            reader.GetOutput().GetPointData().GetArray(0).SetValue(i,radius)

        # Write the edited vtk data into the new centerlines file 
        writer.SetInputData(reader.GetOutput())
        writer.Write()
       

    print "Developing Mesh"

    # With the scalled radii generate a new mesh from with adapted centerlines file  -- here the discretisation dimension (i.e nunber of nodes in each axis) is currently 200 200 200, but might need changing 
    command = 'vmtk vmtkcenterlinemodeller -ifile ' + ScalledCenterLines +' -radiusarray Radius -dimensions 150 150 150 --pipe vmtkmarchingcubes -ofile '+ ScalledMeshVTK # -handle Self'# -dimensions '+3   /Users/jcrawshaw/docker-polnet-master/GeneratingShrunkMesh/Original2.vtk
    subprocess.call(command, shell=True)

    
    # # pause()
    # The Mesh is currently dense and messy, remesh to get a nicer mesh, can control the target size of each element
    command = 'vmtksurfaceremeshing -ifile '+ScalledMeshVTK +' -iterations 6  -area 125 -ofile ' +ScalledMeshVTK_2
    subprocess.call(command, shell=True)

    # # Need to convert the mesh from a vtk to a vtu, this step is critical 
    convertFile(ScalledMeshVTK_2, ScalledMeshVTU)
    # ScalledMeshVTUN = args.ifile + 'PlexusNew.vtu'


# #     # # Remove the ends from the mesh 
#     RemoveInletCaps = 'vmtkmeshclipper -ifile '+ScalledMeshVTU +' -ofile '+ScalledMeshVTUN
#     subprocess.call(RemoveInletCaps, shell=True)



    print "-------------------------------------"
    print "- Scalled Mesh has been generated  -"
    print "-------------------------------------"
    print "-------------  Finito  --------------"
    print "-------------------------------------"






