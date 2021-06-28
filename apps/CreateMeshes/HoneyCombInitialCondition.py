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
from argparse import ArgumentParser
# from vtk.numpy_interface import dataset_adapter as dsa



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

def CreateIdealSkeleton(Directory, Generations,GenerationsX,Height,HorizonatalEdgeLength, theta, Translation, CollapseBranch, Bounds):

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
    LeftBound = Bounds[0]
    RightBound = Bounds[1]
    Branch = []
    NumberOfNodes = len(Nodes)
    counter = 0
    if CollapseBranch:
        for i in Nodes:
            if len(i)>1:
                counter =  counter +1
                if (i[1] ==1.4 and  i[0]>LeftBound-0.01 and i[0]< RightBound+0.01):
                    Branch.append(counter)
    # pdb.set_trace()
    Edges.remove(Branch)
    
    RightCollpaseBoundary =  Nodes[Branch[1]][0]
    LeftCollpaseBoundary = Nodes[Branch[0]][0]
    NumberOfIntervals = 100
    for k in range(1, NumberOfIntervals):
        # NewNodes =  np.array(Nodes[Branch[0]]) + k*(np.array(Nodes[Branch[1]]) - np.array(Nodes[Branch[0]]))/NumberOfIntervals
        NewNodes =  np.array(Nodes[Branch[0]]) + k*(np.array(Nodes[Branch[1]]) - np.array(Nodes[Branch[0]]))/NumberOfIntervals

        Nodes.append(NewNodes)

        NodeNumber = NumberOfNodes+k
        if k == 1:
            LeftCollpaseBoundary = NewNodes[0]
            NewEdge = [Branch[0],NodeNumber]
            Edges.append(NewEdge)
            
        elif k == NumberOfIntervals-1:
            NewEdge = [NodeNumber-1, Branch[1]]
            Edges.append(NewEdge)
            RightCollpaseBoundary = NewNodes[0]
        else: 
            NewEdge = [NodeNumber-1, NodeNumber]
            Edges.append(NewEdge)


    RightCollpaseBoundary =  Nodes[Branch[1]][0]
    LeftCollpaseBoundary = Nodes[Branch[0]][0]
        
        

    Boundaryies  = [LeftCollpaseBoundary, RightCollpaseBoundary]
    with open(Directory+ 'CenterlinePoints.txt', 'w') as filehandle:
        filehandle.writelines("%s\n" % points for points in Nodes)

    # Set up the edges and write into a file 
    with open(Directory+ 'CenterlineEdges.txt', 'w') as filehandle:
        filehandle.writelines("%s\n" % Edges for Edges in Edges)
    return Boundaryies
    


def CreateNewFolder(Directory):
    if path.isdir(Directory)==0:
        os.mkdir(Directory)


if __name__=="__main__":
    print 'Hi'
    CPP_Centerlines_vtp_writer = "/home/vascrem/Chaste/projects/VascularRemodelling/build/optimised/GenerateIdealVascularMesh/Test_VTP_writer_WithCollapseRunner"

    # # Define arguments to be parsed
    # parser = ArgumentParser(description='Create a honeycomb mesh')
    # parser.add_argument('-Directory', dest='Directory', type=str)
    # parser.add_argument('-Collapse',  dest='Collapse',  type=str)
    # parser.add_argument('-Length',   dest='Length',     type=str)
    # parser.add_argument('-Angle',    dest='Angle',      type=str)
    # parser.add_argument('-WaitFileGeneration', dest='WaitFileGeneration',type=str)
       
    # # Parse arguments (this will create args.flow_to_vtu etc. variables)
    # args = parser.parse_args()
    i = "1"
    Directory = "/data/vascrem/MeshCollection/InitialHoneycombNetworkForFSI/"
    CreateNewFolder(Directory) 

    HorizonatalEdgeLength =1 
    Alpha = m.pi/4# float(args.Angle)
    Radius = '0.02'

    # --------------------------
    Height =1.4
    y = 0.7
    xia = y/m.tan(Alpha)

    Bound1  = 2* HorizonatalEdgeLength+2*xia 
    Bound2  = 3* HorizonatalEdgeLength+2*xia 
    Bounds = [Bound1, Bound2]

    # --------------------------
  
    End = 5*HorizonatalEdgeLength + 4* xia
    ClippingRegion = 0.2 *HorizonatalEdgeLength
    ends = [ClippingRegion,End-ClippingRegion]

    # --------------------------

    CenterLines_filename = Directory + "Centerlines_"+i+".vtp"
    Boundaries = CreateIdealSkeleton(Directory, 3, 2, Height, HorizonatalEdgeLength, Alpha, 0, 1, Bounds)
    
    # # ---- Read points into the vtp writer to generate a centerlines.vtp file -------------# 
    command = CPP_Centerlines_vtp_writer + ' -ofile ' + CenterLines_filename + ' -CenterlinePoints ' +Directory+ 'CenterlinePoints.txt -CenterlineEdges ' + Directory +'CenterlineEdges.txt -Radius 0.02 -LeftBound '+str(Boundaries[0]+0.001)+' -RightBound '+str(Boundaries[1]-0.01)+' -Collapse ' + i + ' -Branch 1'
    
    subprocess.call(command, shell=True)

    # --------------------------
    VTK_Mesh = Directory+"Meshinital_"+i+".vtp"
    VTK_Meshremeshed = Directory+"Mesh_"+i+".vtk"
    Clipped_Mesh = Directory+"MeshClipped_"+i+".vtk"
    
    Outputstl = Directory+"Mesh_"+i+".stl" 

    # ---- Interpolate the points in the centerlines file, this will reduce the refinment needed in the centerline modeller -------------# 
    SmoothCenterlinesCommond = 'vmtkcenterlineresampling -ifile '+ CenterLines_filename + ' -length 0.01 -ofile '+ CenterLines_filename
    subprocess.call(SmoothCenterlinesCommond, shell=True)

    # # ----  Generate a mesh from the centerlines file -------------# 

    command = 'vmtk vmtkcenterlinemodeller -ifile ' + CenterLines_filename +' -radiusarray Radius -dimensions 210 200 200 --pipe vmtkmarchingcubes -ofile '+ VTK_Mesh
    subprocess.call(command, shell=True)

    # The Mesh is currently dense and messy, remesh to get a nicer mesh, can control the target size of each element

    command = 'vmtksurfaceremeshing -ifile '+VTK_Mesh +' -iterations 5 -edgelength 0.01 -elementsizemode "edgelength" -ofile ' + VTK_Meshremeshed
    subprocess.call(command, shell=True)

    # ----  Clip Edges -------------# 

    clip.clip_surface_with_plane(VTK_Meshremeshed,(ends[0],0,0), (1,0,0), Clipped_Mesh)
    clip.clip_surface_with_plane(Clipped_Mesh,(ends[1],0,0), (-1,0,0), Clipped_Mesh)
    
    # ----  Convert files vtk to stl :)   -------------#    

    Outputstl = Directory+"Mesh_PreScale"+i+".stl"
    ConvertVTKtoSTL.convert(Clipped_Mesh, Outputstl)

    Outputvtu = Directory+"Mesh_PreScale"+i+".vtu"
    command = 'meshio-convert  ' + Outputstl + ' ' +Outputvtu
    subprocess.call(command, shell=True)

    ScaledMesh = Directory+"Mesh_Scaled"+i+".vtu"
    # # ---- Interpolate the points in the centerlines file, this will reduce the refinment needed in the centerline modeller -------------# 
    Scale = 'vmtkmeshscaling -ifile '+ Outputvtu + ' -scale 0.20  --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ ScaledMesh
    subprocess.call(Scale, shell=True)
    print "Meshio covert 1 "

    ScaledMeshstl = Directory+"ScaledMesh."+i+".stl"
    convert = 'meshio-convert '+ ScaledMesh +'  '+ ScaledMeshstl
    subprocess.call(convert, shell=True)

    print "Meshio covert 2 "
    # ----  Scale files  -------------# 
    Scale = 'vmtkmeshscaling -ifile '+ Outputvtu + ' -scale 0.20  --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ ScaledMesh
    subprocess.call(Scale, shell=True)

    print "Meshio covert 2 "
    # ----  Convert files  -------------# 
    convert = 'meshio-convert '+ ScaledMesh +'  '+ ScaledMeshstl
    subprocess.call(convert, shell=True)

    print "----  Initial mesh generated  ------------- "
    


