#!/usr/bin/env python
#
import subprocess
import vtk
import os
import pdb
import string
import math as m
import vmtk
from vmtk import pypes
import clip
import array as A
import numpy as np
import meshio
import vtk


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
        print output.GetNumberOfCells()
        
        # loop over all of the elements in the vtk polydata file 
        # For each element, the polydata has recorded the location of each node, but not the node index. As such here we loop over the elements   and select out the node, which are put in the Node list.  As I loop over the elements I will also record which node indices are in each element
        for i in range(output.GetNumberOfCells()):
            Element = []
            pts = output.GetCell(i).GetPoints()  
            np_pts = np.array([pts.GetPoint(i) for i in range(pts.GetNumberOfPoints())]) 
            
            # Saving each of the nodes in this element
            for j in [0,1,2]:
                AddNodeToList = 1 
                # want to add the Node to the Nodes list, But need to check if it is in the list yet,  Do this by looping over list and seeing if it is in there
                for k in Nodes: 
                    #if Nodes.count(pts.GetPoint(j)) > 0:
                    if(k == pts.GetPoint(j)):
                        AddNodeToList =0
                        break
                if (AddNodeToList ==1):
                    # print 'Adding Node'
                        # Filling in the array storing Node locations
                    Nodes.append(pts.GetPoint(j))
                NodeIndex = Nodes.index(pts.GetPoint(j))
                Element.append(NodeIndex)
                # This works to fill in Nodes
            ElementList.append(Element)
    
        # Save the nodes as points and elements for the meshio writer
        points = np.array(Nodes)
        elements = {
        "triangle": np.array(ElementList
        )
        }    

        meshio.write_points_cells(
        VTUfile,
        points,
        elements,
        )

        print 'Finished'

def UpdateNodesAndEdges(NewNodes,Nodes,Edges):
    # Index = [[]]
    Index = [ 0 ,0 ]
    for i in [0,1]:
        if NewNodes[i] in Nodes :
            # THis node is already in the list
            Index[i] = Nodes.index(NewNodes[i])# Get the correct index
        else:
            Nodes.append(NewNodes[i]) # Add to the list 
            Index[i] = len(Nodes)-1# Last index is the needed one
    Edges.append(Index)
    return Nodes,Edges

def CreateIdealSkeleton(Directory, Generations,GenerationsX,Height,HorizonatalEdgeLength, theta, Translation):

    # Symetric about the yaxis, Standard length of Horizontal edges, y is the vertical hight
    L = HorizonatalEdgeLength
    y = Height

    Nodes = [[]]
    Edges= []
    X_0 = 0

    Range = range(- (Generations -1)/2, (Generations -1)/2+1)

    for j in range(GenerationsX):
        for i in Range:
            # Do a whole hex in one shot
            NewNodes = [[X_0,i*y], [X_0+L,i*y]] # Horizontal compotent on left
            Nodes,Edges = UpdateNodesAndEdges(NewNodes,Nodes,Edges)
    
            # LastNode = [L,i*y]
            A = y/(2*m.tan(theta))+(X_0+L)
            NewNodes1 = [[X_0+L,i*y], [A, y/2+i*y ]] # Left diag
            NewNodes2 = [[X_0+L,i*y], [A, -y/2+i*y ]] # Left diag
            Nodes,Edges = UpdateNodesAndEdges(NewNodes1,Nodes,Edges)
            Nodes,Edges = UpdateNodesAndEdges(NewNodes2,Nodes,Edges)

            NewNodes1 = [ [A, y/2+i*y ], [A+L, y/2+i*y ]] # horixontal 
            NewNodes2 = [[A, -y/2+i*y ], [A+L, -y/2+i*y ]] # horixontal 
            Nodes,Edges = UpdateNodesAndEdges(NewNodes1,Nodes,Edges)
            Nodes,Edges = UpdateNodesAndEdges(NewNodes2,Nodes,Edges)

            A2 = A + y/(2*m.tan(theta))
            NewNodes1 = [[A+L, y/2+i*y ], [ A2+L,i*y ]] # Left diag
            NewNodes2 = [[A+L, -y/2+i*y ], [ A2+L,i*y ]] # Left diag
            Nodes,Edges = UpdateNodesAndEdges(NewNodes1,Nodes,Edges)
            Nodes,Edges = UpdateNodesAndEdges(NewNodes2,Nodes,Edges)

        maxInColumns = np.amax(Nodes, axis=0)
        X_0 = maxInColumns[0]

    # THe outlet region
    for i in Range:
        NewNodes = [[X_0,i*y], [X_0+L,i*y]] # Horizontal compotent on left
        Nodes,Edges = UpdateNodesAndEdges(NewNodes,Nodes,Edges)
    MaxX = X_0+L
  

    with open(Directory+ 'CenterlinePoints.txt', 'w') as filehandle:
        filehandle.writelines("%s\n" % points for points in Nodes)

    # Set up the edges and write into a file 
    with open(Directory+ 'CenterlineEdges.txt', 'w') as filehandle:
        filehandle.writelines("%s\n" % Edges for Edges in Edges)
    return MaxX
    


if __name__=="__main__":
  
    server = 1
    if server: 
        CPP_Centerlines_vtp_writer = "/home/vascrem/Chaste/projects/VascularRemodelling/build/optimised/GenerateIdealVascularMesh/Test_VTP_writerRunner"
        Directory = "/home/vascrem/MeshCollection/IdealisedNetwork/DelfatedNetwork/"
    else:
        CPP_Centerlines_vtp_writer = "/Users/jcrawshaw/Chaste/projects/VascularRemodelling/build/optimised/GenerateIdealVascularMesh/Test_VTP_writerRunner"
        Directory = "/Users/jcrawshaw/Documents/Projects/MeshCollection/IdealiseNetworks/DelfatedNetwork/"

    VTK_Mesh = Directory+"Meshinital.vtk"
    VTK_MeshRefined = Directory+"MeshRefined.vtk"
    Clipped_Mesh = Directory+"MeshClipped.vtk"

    GenerationsHeigh = 2
    Height = 28
    HorizonatalEdgeLength = 20
    GenerationsLong = 2
    alpha = m.pi/4

    MaxX = CreateIdealSkeleton(Directory, GenerationsHeigh, GenerationsLong, Height, HorizonatalEdgeLength, alpha, 0)

    CenterLines_filename = Directory + "Centerlines.vtp"
    print CenterLines_filename

    # # # # # # # # read in the centerlines points into cpp to generate the centerlines.vtp file
    command = CPP_Centerlines_vtp_writer + ' -ofile ' + CenterLines_filename + ' -CenterlinePoints ' +Directory+ 'CenterlinePoints.txt -CenterlineEdges ' + Directory +'CenterlineEdges.txt -Radius 0.1' 
    subprocess.call(command, shell=True)

    # # #  # interpolate the points in the centerlines file, this will reduce the refinment needed in the centerline modeller
    SmoothCenterlinesCommond = 'vmtkcenterlineresampling -ifile '+ CenterLines_filename + ' -length 0.7 -ofile '+ CenterLines_filename
    subprocess.call(SmoothCenterlinesCommond, shell=True)
    
    
    print "Developing Mesh"
    # With the scalled radii generate a new mesh from with adapted centerlines file  -- here the discretisation dimension (i.e nunber of nodes in each axis) is currently 200 200 200, but might need changing 
    command = 'vmtk vmtkcenterlinemodeller -ifile ' + CenterLines_filename +' -radiusarray Radius -dimensions 150 150 150 --pipe vmtkmarchingcubes -ofile '+ VTK_Mesh
    subprocess.call(command, shell=True)
    
    VTK_MeshRefined1 = Directory+"FirstRefined.vtk"
    # The Mesh is currently dense and messy, remesh to get a nicer mesh, can control the target size of each element
    command = 'vmtksurfaceremeshing -ifile '+VTK_Mesh +' -iterations 5 -edgelength 0.001 -elementsizemode "edgelength" -ofile ' + VTK_MeshRefined1
    subprocess.call(command, shell=True)

    VTK_MeshRefined2 = Directory+"SecondRefined2.vtk"
    # The Mesh is currently dense and messy, remesh to get a nicer mesh, can control the target size of each element
    command = 'vmtksurfaceremeshing -ifile '+VTK_MeshRefined1 +' -iterations 10 -edgelength 0.02 -elementsizemode "edgelength" -ofile ' + VTK_MeshRefined2
    subprocess.call(command, shell=True)
    

    # VTK_MeshRefined3 = Directory+"ThirdRefined.vtk"
    # # # The Mesh is currently dense and messy, remesh to get a nicer mesh, can control the target size of each element
    # command = 'vmtksurfaceremeshing -ifile '+VTK_MeshRefined2 +' -iterations 5 -edgelength 0.05 -elementsizemode "edgelength" -ofile ' + VTK_MeshRefined3
    # subprocess.call(command, shell=True)

    
   
    Clipped_Mesh = Directory+"Clipped.vtk"
    clip.clip_surface_with_plane(VTK_MeshRefined2,(11,0,0), (1,0,0), Clipped_Mesh)
    clip.clip_surface_with_plane(Clipped_Mesh,(170,0,0), (-1,0,0), Clipped_Mesh)


    print "-------------------------------------"
    print "-------------  Finito  --------------"
    print "-------------------------------------"




    print "-------------------------------------"
    print "-------------  Finito  --------------"
    print "-------------------------------------"

