#!/usr/bin/env python
#
import subprocess
import vtk
import os
import pdb
import string
import math as m
import numpy as np
import meshio
import array as A
# ##
import vmtk
from vmtk import pypes
# from vtk import vtkXMLUnstructuredGridWriter as Writer, VTK_VERSION
# from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
from vtk.util import numpy_support
import vtk
import clip

from vtk.numpy_interface import dataset_adapter as dsa



# ---------------------------------
  
def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

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
  
    # for i in range(1,len(Nodes)):
    #     Nodes[i] += [0,7.0]
    LeftBound = 3.0
    RightBound = 5.0
    Branch = []
    NumberOfNodes = len(Nodes)
    # print Edges
    counter = 0
    if CollapseBranch:
        for i in Nodes:
            
            if len(i)>1:
                counter =  counter +1
                if (i[1] == 0 and  i[0]>LeftBound and i[0]< RightBound):
 
                    Branch.append(counter)
    Edges.remove(Branch)


    # If the connected branches here
    
    for i in [1,0]:
        EdgesWithEndA = [ E for E in Edges if E[1-i]==Branch[i]]
        BiA = EdgesWithEndA[0]
        BiA = [ E for E in BiA if E!=Branch[i]][0]
        BiB = EdgesWithEndA[1]
        BiB = [ E for E in BiB if E!=Branch[i]][0]

        Edges.remove(EdgesWithEndA[0]), Edges.remove(EdgesWithEndA[1])
        NewNodes1 = np.array(Nodes[Branch[i]]) - 1.0/6* (np.array(Nodes[Branch[i]]) - np.array(Nodes[BiA]) )
        NewNodes2 = np.array(Nodes[Branch[i]]) - 1.0/6* (np.array(Nodes[Branch[i]]) - np.array(Nodes[BiB]) )
        Nodes.append(NewNodes1)
        Nodes.append(NewNodes2)

        Edges.append([EdgesWithEndA[0][0],NumberOfNodes]), Edges.append([EdgesWithEndA[0][1],NumberOfNodes])
        Edges.append([EdgesWithEndA[1][0],NumberOfNodes+1]), Edges.append([EdgesWithEndA[1][1],NumberOfNodes+1])

        NumberOfNodes = len(Nodes)
    # NumberOfNodes = len(Nodes)
    NumbNewNodes = 3

    NewNodes1 = np.array(Nodes[Branch[0]]) + 1.0*(np.array(Nodes[Branch[1]]) - np.array(Nodes[Branch[0]]))/6.0
    LeftCollpaseBoundary = NewNodes1[0]
    Nodes.append(NewNodes1)
    N1 = NumberOfNodes
    NewNodes2 =  np.array(Nodes[Branch[0]]) + 2.0*(np.array(Nodes[Branch[1]]) - np.array(Nodes[Branch[0]]))/6.0
    Nodes.append(NewNodes2)
    N2 = NumberOfNodes+1
    NewNodes3 =  np.array(Nodes[Branch[0]]) + 3.0*(np.array(Nodes[Branch[1]]) - np.array(Nodes[Branch[0]]))/6.0
    Nodes.append(NewNodes3)
    N3 = NumberOfNodes+2
    NewNodes4 =  np.array(Nodes[Branch[0]]) + 4.0*(np.array(Nodes[Branch[1]]) - np.array(Nodes[Branch[0]]))/6.0
    Nodes.append(NewNodes4)
    N4 = NumberOfNodes+3
    NewNodes5 =  np.array(Nodes[Branch[0]]) + 5.0*(np.array(Nodes[Branch[1]]) - np.array(Nodes[Branch[0]]))/6.0
    RightCollpaseBoundary = NewNodes5[0]
    Nodes.append(NewNodes5)
    N5 = NumberOfNodes+4
    NewEdges = [[Branch[0],N1], [N1,N2],  [N2,N3], [N3,N4], [N4,N5], [N5,Branch[1] ]]

    for j in NewEdges:
        Edges.append(j)



 # Set the collapsing branch 
    Boundaryies  = [LeftCollpaseBoundary, RightCollpaseBoundary]

    LeftBound = 3.0
    RightBound = 5.0
    Branch = []
    NumberOfNodes = len(Nodes)
    # print Edges
    counter = 0
    if CollapseBranch:
        for i in Nodes:
            if len(i)>1:
                counter =  counter +1
                if (i[1] < -0.5 and  i[0]>LeftBound and i[0]< RightBound):
                    Branch.append(counter)
    print Branch
    Edges.remove(Branch)

    NumbNewNodes = 3

    NewNodes1 = np.array(Nodes[Branch[0]]) + 1.0*(np.array(Nodes[Branch[1]]) - np.array(Nodes[Branch[0]]))/6.0
    LeftCollpaseBoundary = NewNodes1[0]
    Nodes.append(NewNodes1)
    N1 = NumberOfNodes
    NewNodes2 =  np.array(Nodes[Branch[0]]) + 2.0*(np.array(Nodes[Branch[1]]) - np.array(Nodes[Branch[0]]))/6.0
    Nodes.append(NewNodes2)
    N2 = NumberOfNodes+1
    NewNodes3 =  np.array(Nodes[Branch[0]]) + 3.0*(np.array(Nodes[Branch[1]]) - np.array(Nodes[Branch[0]]))/6.0
    Nodes.append(NewNodes3)
    N3 = NumberOfNodes+2
    NewNodes4 =  np.array(Nodes[Branch[0]]) + 4.0*(np.array(Nodes[Branch[1]]) - np.array(Nodes[Branch[0]]))/6.0
    Nodes.append(NewNodes4)
    N4 = NumberOfNodes+3
    NewNodes5 =  np.array(Nodes[Branch[0]]) + 5.0*(np.array(Nodes[Branch[1]]) - np.array(Nodes[Branch[0]]))/6.0
    RightCollpaseBoundary = NewNodes5[0]
    Nodes.append(NewNodes5)
    N5 = NumberOfNodes+4
    NewEdges = [[Branch[0],N1], [N1,N2],  [N2,N3], [N3,N4], [N4,N5], [N5,Branch[1] ]]

    for j in NewEdges:
        Edges.append(j)

    
    with open(Directory+ 'CenterlinePoints.txt', 'w') as filehandle:
        filehandle.writelines("%s\n" % points for points in Nodes)

    # Set up the edges and write into a file 
    with open(Directory+ 'CenterlineEdges.txt', 'w') as filehandle:
        filehandle.writelines("%s\n" % Edges for Edges in Edges)
    return Boundaryies
    


if __name__=="__main__":
    print " The mesh generated by this code is for a HemeLB mesh, not okay for Chaste"
    CPP_Centerlines_vtp_writer = "/Users/jcrawshaw/Chaste/projects/VascularRemodelling/build/optimised/GenerateIdealVascularMesh/Test_VTP_writer_WithCollapseAndGrowthRunner"
    CPP_Centerlines_vtp_writer = "/home/vascrem/Chaste/projects/VascularRemodelling/build/optimised/GenerateIdealVascularMesh/Test_VTP_writer_WithCollapseAndGrowthRunner"

    Directory = "/Users/jcrawshaw/docker-polnet-master/"
    # /Users/jcrawshaw/docker-polnet-master/IdealiseNetworks/IdealMeshWIthCentralCollapse/Nonsymmetric/"
    Directory = "/home/vascrem/MeshCollection/IdealisedNetwork/NonSymmetricCollapseAndGrowth/"

    try:
        os.mkdir(Directory)
    except OSError:
        print ("Creation of the directory %s failed" % Directory)
    else:
        print ("Successfully created the directory %s " % Directory)
    
    # # # # Set up the points for the centerlines and write into a file to be read in cpp 

    GenerationsHeigh = 2
    Height =1.4
    HorizonatalEdgeLength =1
    GenerationsLong = 2
    alpha = m.pi/4
    CollapseBranch =1
    ExpandBranch = 1
     
    Boundaries = CreateIdealSkeleton(Directory, GenerationsHeigh, GenerationsLong, Height, HorizonatalEdgeLength, alpha, 0, CollapseBranch)
    print Boundaries
    Collapse =[0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0]
    for i in Collapse:
        CenterLines_filename = Directory + "Centerlines"+str(int(10*i))+".vtp"
        print CenterLines_filename
        ExpandingFactor = 1-i+1
        # ---- Read points into the vtp writer to generate a centerlines.vtp file -------------# 
        command = CPP_Centerlines_vtp_writer + ' -ofile ' + CenterLines_filename + ' -CenterlinePoints ' +Directory+ 'CenterlinePoints.txt -CenterlineEdges ' + Directory +'CenterlineEdges.txt -Radius 0.2 -LeftBound '+str(Boundaries[0]-0.2)+' -RightBound '+str(Boundaries[1]+0.2)+' -Collapse ' + str(i) + ' -Expand '+ str(ExpandingFactor) 
        subprocess.call(command, shell=True)

        # print "Jess is great"
        # # # ----  Generate a mesh from the centerlines file -------------# 
        VTK_Mesh = Directory+"Meshinital"+str(int(10*i))+".vtk"
        print VTK_Mesh
        if i == 0.1 or i ==  0.2 or i ==  0:
            command = 'vmtk vmtkcenterlinemodeller -ifile ' + CenterLines_filename +' -radiusarray Radius -dimensions 210 200 200 --pipe vmtkmarchingcubes -ofile '+ VTK_Mesh
            subprocess.call(command, shell=True)
        else:
            command = 'vmtk vmtkcenterlinemodeller -ifile ' + CenterLines_filename +' -radiusarray Radius -dimensions 110 110 110 --pipe vmtkmarchingcubes -ofile '+ VTK_Mesh
            subprocess.call(command, shell=True)
        
        VTK_Meshremeshed = Directory+"Mesh"+str(int(10*i))+".vtk"
       
        if i == 0.1:
            # VTK_Mesh= Directory+"MeshClipped1Cource.vtk"
            command = 'vmtksurfaceremeshing -ifile '+VTK_Mesh +' -iterations 3 -edgelength 0.02 -elementsizemode "edgelength" -ofile ' + VTK_Meshremeshed
            subprocess.call(command, shell=True)
        elif i < 0.35:
            command = 'vmtksurfaceremeshing -ifile '+VTK_Mesh +' -iterations 3 -edgelength 0.04 -elementsizemode "edgelength" -ofile ' + VTK_Meshremeshed
            subprocess.call(command, shell=True)
        else: #if i == 0.6 or i ==  0.5 or i ==  0.4 or i ==  0.3:
            command = 'vmtksurfaceremeshing -ifile '+VTK_Mesh +' -iterations 3 -edgelength 0.05 -elementsizemode "edgelength" -ofile ' + VTK_Meshremeshed
            subprocess.call(command, shell=True)
    
        # ----  Clip Edges -------------# 
        Clipped_Mesh = Directory+"MeshClipped"+str(int(10*i))+".vtk"
        clip.clip_surface_with_plane(VTK_Meshremeshed,(0.4,0,0), (1,0,0), Clipped_Mesh)
        clip.clip_surface_with_plane(Clipped_Mesh,(7.4,0,0), (-1,0,0), Clipped_Mesh)