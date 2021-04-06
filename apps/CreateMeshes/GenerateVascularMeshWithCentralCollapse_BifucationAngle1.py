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
    


if __name__=="__main__":
  
    print " The mesh generated by this code is for a HemeLB mesh, not okay for Chaste"
    CPP_Centerlines_vtp_writer = "/Users/jcrawshaw/Chaste/projects/VascularRemodelling/build/optimised/GenerateIdealVascularMesh/Test_VTP_writer_WithCollapseRunner"
    CPP_Centerlines_vtp_writer = "/home/vascrem/Chaste/projects/VascularRemodelling/build/optimised/GenerateIdealVascularMesh/Test_VTP_writer_WithCollapseRunner"

    Directory = "/Users/jcrawshaw/docker-polnet-master/"
    Directory = "/data/vascrem/MeshCollection/IdealisedNetwork/AngleVariation_3X3NetworkExtended/"
    if path.isdir(Directory)==0:
        os.mkdir(Directory)
        
    
    # # # # Set up the points for the centerlines and write into a file to be read in cpp 
    L = 7.8 # Total length of network
    GenerationsHeigh = 3
    Height =1.4
    OldHorizonatalEdgeLength =1
    x1 = 1
    x2 = 0.7000000000000002
    y = 0.7
    GenerationsLong = 2

    print "I love Jessie"
    #    --- 
    Alpha = [m.pi/6, m.pi/3, m.pi/2.2] 
    Alpha = [m.pi/6] 
    FileLabels = ['/PI_6/', '/PI_3/', '/PI_2.2/'] 
    FileLabels = ['/PI_6/'] 
    counter = -1
    for alpha in Alpha:
        counter = counter+1
        AngleDirectory = Directory+FileLabels[counter]
        if path.isdir(AngleDirectory)==0:
            os.mkdir(AngleDirectory)

        HorizonatalEdgeLength = 0.2*(L - 4*y/m.tan(alpha))
        xia = y/m.tan(alpha)

        Bound1  = 2* HorizonatalEdgeLength+2*xia 
        Bound2  = 3* HorizonatalEdgeLength+2*xia 
        Bounds = [Bound1, Bound2]
        print Bounds
        ends = [0.4,7.4]
    # scp ~/Chaste/projects/VascularRemodelling/apps/GenerateVascularMeshWithCentralCollapse_BifucationAngle.py vascrem@6200rs-osborne-l.science.unimelb.edu.au:///home/vascrem//Chaste/projects/VascularRemodelling/apps/GenerateVascularMeshWithCentralCollapse_BifucationAngle.py
    
        # // need to get the right collapse points here matlabd 
        # Collapse = [0, 0.1227, 0.2248, 0.3178, 0.4170]
        Collapse = [ 0.5124, 0.6119, 0.7080, 0.8059, 0.9032, 1.0000]
        for i in Collapse:
            CenterLines_filename = AngleDirectory + "Centerlines_"+str(i)+".vtp"

            Boundaries = CreateIdealSkeleton(AngleDirectory, GenerationsHeigh, GenerationsLong, Height, HorizonatalEdgeLength, alpha, 0, 1, Bounds)
            
            # ---- Read points into the vtp writer to generate a centerlines.vtp file -------------# 
            command = CPP_Centerlines_vtp_writer + ' -ofile ' + CenterLines_filename + ' -CenterlinePoints ' +AngleDirectory+ 'CenterlinePoints.txt -CenterlineEdges ' + AngleDirectory +'CenterlineEdges.txt -Radius 0.2 -LeftBound '+str(Boundaries[0])+' -RightBound '+str(Boundaries[1])+' -Collapse ' + str(i) + ' -Branch ' + str(1)
            subprocess.call(command, shell=True)


            VTK_Mesh = AngleDirectory+"Meshinital_"+str(i)+".vtp"
            VTK_Meshremeshed = AngleDirectory+"Mesh_"+str(i)+".vtk"
            Clipped_Mesh = AngleDirectory+"MeshClipped_"+str(i)+".vtk"
            
            Outputstl = AngleDirectory+"Mesh_"+str(i)+".stl" 

            # ---- Interpolate the points in the centerlines file, this will reduce the refinment needed in the centerline modeller -------------# 
            SmoothCenterlinesCommond = 'vmtkcenterlineresampling -ifile '+ CenterLines_filename + ' -length 0.01 -ofile '+ CenterLines_filename
            subprocess.call(SmoothCenterlinesCommond, shell=True)
        
            # # ----  Generate a mesh from the centerlines file -------------# 
        
            if i == 0.1 or i ==  0.2 or i ==  0.3:
                command = 'vmtk vmtkcenterlinemodeller -ifile ' + CenterLines_filename +' -radiusarray Radius -dimensions 210 200 200 --pipe vmtkmarchingcubes -ofile '+ VTK_Mesh
                subprocess.call(command, shell=True)
            else:
                command = 'vmtk vmtkcenterlinemodeller -ifile ' + CenterLines_filename +' -radiusarray Radius -dimensions 130 130 130 --pipe vmtkmarchingcubes -ofile '+ VTK_Mesh
                subprocess.call(command, shell=True)

            # The Mesh is currently dense and messy, remesh to get a nicer mesh, can control the target size of each element
            if i == 0.1:
                # VTK_Mesh= Directory+"MeshClipped1Cource.vtk"
                command = 'vmtksurfaceremeshing -ifile '+VTK_Mesh +' -iterations 5 -edgelength 0.03 -elementsizemode "edgelength" -ofile ' + VTK_Meshremeshed
                subprocess.call(command, shell=True)
            elif i < 0.35:
                command = 'vmtksurfaceremeshing -ifile '+VTK_Mesh +' -iterations 5 -edgelength 0.01 -elementsizemode "edgelength" -ofile ' + VTK_Meshremeshed
                subprocess.call(command, shell=True)
            else:
                command = 'vmtksurfaceremeshing -ifile '+VTK_Mesh +' -iterations 5 -edgelength 0.03 -elementsizemode "edgelength" -ofile ' + VTK_Meshremeshed
                subprocess.call(command, shell=True)

            # ----  Clip Edges -------------# 
            print "Jess is good"
            print Clipped_Mesh
            clip.clip_surface_with_plane(VTK_Meshremeshed,(ends[0],0,0), (1,0,0), Clipped_Mesh)
            clip.clip_surface_with_plane(Clipped_Mesh,(ends[1],0,0), (-1,0,0), Clipped_Mesh)
            print 'about to convert'
            # ----  Convert files vtk to stl :)   -------------#    

            Outputstl = AngleDirectory+"Mesh_PreScale"+str(i)+".stl"
            ConvertVTKtoSTL.convert(Clipped_Mesh, Outputstl)


            Outputvtu = AngleDirectory+"Mesh_PreScale"+str(i)+".vtu"
            command = 'meshio-convert  ' + Outputstl + ' ' +Outputvtu
            subprocess.call(command, shell=True)

            ScaledMesh = AngleDirectory+"Mesh_Scaled"+str(i)+".vtu"
            # # ---- Interpolate the points in the centerlines file, this will reduce the refinment needed in the centerline modeller -------------# 
            Scale = 'vmtkmeshscaling -ifile '+ Outputvtu + ' -scale 0.20  --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ ScaledMesh
            subprocess.call(Scale, shell=True)

            ScaledMeshstl = AngleDirectory+"ScaledMesh."+str(i)+".stl"
            convert = 'meshio-convert '+ ScaledMesh +'  '+ ScaledMeshstl
            subprocess.call(convert, shell=True)
            print "Done one"
        


            print 'Scale'
            # ----  Scale files  -------------# 
            Scale = 'vmtkmeshscaling -ifile '+ Outputvtu + ' -scale 0.20  --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ ScaledMesh
            subprocess.call(Scale, shell=True)

            # ----  Convert files  -------------# 
            convert = 'meshio-convert '+ ScaledMesh +'  '+ ScaledMeshstl
            subprocess.call(convert, shell=True)

    print " Finished Creating files "
        
        