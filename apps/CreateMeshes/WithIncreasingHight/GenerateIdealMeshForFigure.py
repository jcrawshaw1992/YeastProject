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
import array as A
# from vtk import vtkXMLUnstructuredGridWriter as Writer, VTK_VERSION
# from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
import numpy as np
import meshio
from vtk.util import numpy_support
import vtk
import clip
import math
import ConvertVTKtoSTL
from os import path

# from vtk.numpy_interface import dataset_adapter as dsa



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

def CreateIdealSkeleton(Directory, Generations,GenerationsX,Height,HorizonatalEdgeLength, theta, Translation, CollapseBranch):

    # Symetric about the yaxis, Standard length of Horizontal edges, y is the vertical hight
    L = HorizonatalEdgeLength
    y = Height

    Nodes = [[]]
    Edges= []
    X_0 = 0
    print -Generations+2
    Range = range(-Generations+2,2)
    print Range

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

    LeftCollpaseBoundary = 0
    RightCollpaseBoundary = 1

    # for i in range(1,len(Nodes)):
    #     Nodes[i] += [0,7.0]
    LeftBound = 3.0
    RightBound = 5.0
    Branch = []
    NumberOfNodes = len(Nodes)
    # print Edges
    counter = 0

    Branch = []
    NumberOfNodes = len(Nodes)
    # print Edges
    counter = 0
    if CollapseBranch:
        for i in Nodes:
            
            if len(i)>1:
                counter =  counter +1
                if (i[1] == 1.4 and  i[0]>LeftBound and i[0]< RightBound):
 
                    Branch.append(counter)
        Edges.remove(Branch)
        

    

    Boundaryies  = [LeftCollpaseBoundary, RightCollpaseBoundary]
    with open(Directory+ 'CenterlinePoints.txt', 'w') as filehandle:
        filehandle.writelines("%s\n" % points for points in Nodes)

    # Set up the edges and write into a file 
    with open(Directory+ 'CenterlineEdges.txt', 'w') as filehandle:
        filehandle.writelines("%s\n" % Edges for Edges in Edges)
    return Boundaryies
    

if __name__=="__main__":
  
    print " The mesh generated by this code is for a HemeLB mesh, not okay for Chaste"
    CPP_Centerlines_vtp_writer = "/Users/jcrawshaw/Chaste/projects/VascularRemodelling/build/optimised/GenerateIdealVascularMesh/Test_VTP_writer_WithCollapseRunner"
    CPP_Centerlines_vtp_writer = "/home/vascrem/Chaste/projects/VascularRemodelling/build/optimised/GenerateIdealVascularMesh/Test_VTP_writer_WithCollapseRunner"

    Directory = "/Users/jcrawshaw/docker-polnet-master/"
    Directory = "/data/vascrem/MeshCollection/IdealisedNetwork/ForFigure/"
    if path.isdir(Directory)==0:
        os.mkdir(Directory)
        
    # # # # Set up the points for the centerlines and write into a file to be read in cpp 
    L = 7.8 # Total length of network


    GenerationsLong = 1
    GenerationsHeigh = 1
    Height =1.4

    Levels = [1]#,2,3,4]
    scaleY = [0.75, 0.6, 0.51, 0.45] 
    scaleX = [0.7, 0.535, 0.43, 0.35] 
    RightBound = [ 10.7,14.1, 17.5,21]

    
    counter = -1
    alpha = m.pi/4

    counter = counter +1
    LevelsDirectory = Directory
    if path.isdir(LevelsDirectory)==0:
        os.mkdir(LevelsDirectory)

    HorizonatalEdgeLength =1


    i = 1.0

    CenterLines_filename = LevelsDirectory + "Centerlines_"+str(i)+".vtp"
            #  CreateIdealSkeleton(Directory, Generations,GenerationsX,Height,HorizonatalEdgeLength, theta, Translation, CollapseBranch):
    GHigh = GenerationsHeigh
    GLong = GenerationsLong
    HorizonatalEdgeLength = 1#(7+j*2) *0.14285714285
    Height = 1.4
    # Boundaries = Fun(LevelsDirectory, GHigh, GLong, Height, HorizonatalEdgeLength, alpha, 0, 0)
    Boundaries = CreateIdealSkeleton(LevelsDirectory, GHigh, GLong, Height, HorizonatalEdgeLength, alpha, 0, 0)

    # ---- Read points into the vtp writer to generate a centerlines.vtp file -------------# 
    command = CPP_Centerlines_vtp_writer + ' -ofile ' + CenterLines_filename + ' -CenterlinePoints ' +LevelsDirectory+ 'CenterlinePoints.txt -CenterlineEdges ' + LevelsDirectory +'CenterlineEdges.txt -Radius 0.2 -LeftBound '+str(Boundaries[0]+0.05)+' -RightBound '+str(Boundaries[1]-0.05)+' -Collapse ' + str(i) + ' -Branch ' + str(1)
    subprocess.call(command, shell=True)

    VTK_Mesh = LevelsDirectory+"Meshinital_"+str(i)+".vtp"
    ScalledVTK_Mesh = LevelsDirectory+"Scalled_"+str(i)+".vtp"

    VTK_Meshremeshed = LevelsDirectory+"Mesh_"+str(i)+".vtk"
    Clipped_Mesh = LevelsDirectory+"MeshClipped_"+str(i)+".vtk"
    Scalled_Mesh = LevelsDirectory+"MeshScalled_"+str(i)+".vtk"

    Outputstl = LevelsDirectory+"Mesh_"+str(i)+".stl" 



    # ----  Generate a mesh from the centerlines file -------------# 

    # command = 'vmtk vmtkcenterlinemodeller -ifile ' + CenterLines_filename +' -radiusarray Radius -dimensions 40 40 40 --pipe vmtkmarchingcubes -ofile '+ VTK_Mesh
    # subprocess.call(command, shell=True)

    # The Mesh is currently dense and messy, remesh to get a nicer mesh, can control the target size of each element

    # command = 'vmtksurfaceremeshing -ifile '+VTK_Mesh +' -iterations 5 -edgelength 0.09 -elementsizemode "edgelength" -ofile ' + VTK_Meshremeshed
    # subprocess.call(command, shell=True)



    # ----  Clip Edges -------------# 
    print "Jess is good"
    print Clipped_Mesh
    clip.clip_surface_with_plane(VTK_Meshremeshed,(0.4,0,0), (1,0,0), Clipped_Mesh)
    clip.clip_surface_with_plane(Clipped_Mesh,(4,0,0), (-1,0,0), Clipped_Mesh)
    print 'about to convert'

    # ----  Convert files vtk to stl :)   -------------#    

    Outputstl = LevelsDirectory+"ScalledMesh"+str(i)+".stl"
    ConvertVTKtoSTL.convert(Clipped_Mesh, Outputstl)
