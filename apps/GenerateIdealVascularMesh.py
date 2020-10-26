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
# from vtk.util import numpy_support
import vtk
import clip

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

def CreateIdealSkeleton(Directory, Generations,GenerationsX,Height,HorizonatalEdgeLength, theta):

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
  
    CPP_Centerlines_vtp_writer = "/home/vascrem/Chaste/projects/VascularRemodelling/build/optimised/GenerateIdealVascularMesh/Test_VTP_writerRunner"

    Directory = "//home/vascrem/MeshCollection/IdealisedNetwork/"
    CenterLines_filename = Directory + "Centerlines.vtp"
   
    VTK_Mesh = Directory+"Meshinital.vtk"
    Clipped_Mesh = Directory+"MeshClippedThird.vtk"
    VTK_MeshRefined = Directory+"MeshRefined_Second.vtk"
    VTK_MeshRefinedSec = Directory+"MeshRefined_Third.vtk"
    VTU_Mesh = Directory+"Mesh3.vtu"

    # # # # Set up the points for the centerlines and write into a file to be read in cpp 

    GenerationsHeigh = 3
    Height =2
    HorizonatalEdgeLength =1.3
    GenerationsLong = 2
    alpha = m.pi/3
    MaxX = CreateIdealSkeleton(Directory, GenerationsHeigh, GenerationsLong, Height, HorizonatalEdgeLength, alpha)
    # # # # # read in the centerlines points into cpp to generate the centerlines.vtp file
    command = CPP_Centerlines_vtp_writer + ' -ofile ' + CenterLines_filename + ' -CenterlinePoints ' +Directory+ 'CenterlinePoints.txt -CenterlineEdges ' + Directory +'CenterlineEdges.txt -Radius 0.1' 
    subprocess.call(command, shell=True)


    # #  # interpolate the points in the centerlines file, this will reduce the refinment needed in the centerline modeller
    SmoothCenterlinesCommond = 'vmtkcenterlineresampling -ifile '+ CenterLines_filename + ' -length 0.01 -ofile '+ CenterLines_filename
    subprocess.call(SmoothCenterlinesCommond, shell=True)

    print "Developing Mesh"
    With the scalled radii generate a new mesh from with adapted centerlines file  -- here the discretisation dimension (i.e nunber of nodes in each axis) is currently 200 200 200, but might need changing 
    command = 'vmtk vmtkcenterlinemodeller -ifile ' + CenterLines_filename +' -radiusarray Radius -dimensions 180 180 150 --pipe vmtkmarchingcubes -ofile '+ VTK_Mesh
    subprocess.call(command, shell=True)
    
    # # # pause()
    # #The Mesh is currently dense and messy, remesh to get a nicer mesh, can control the target size of each element
    # command = 'vmtksurfaceremeshing -ifile '+VTK_MeshRefined +' -iterations 5 -edgelength 0.075 -elementsizemode "edgelength" -ofile ' +VTK_MeshRefinedSec
    # subprocess.call(command, shell=True)


    # # clip.clip_surface_with_plane(VTK_MeshRefined,(0,-5.5,0), (0,1,0), Clipped_Mesh)
    # # clip.clip_surface_with_plane(Clipped_Mesh,(0,5.5,0), (0,-1,0), Clipped_Mesh)

    # # clip.clip_surface_with_plane(VTK_MeshRefined,(0.39,0,0), (1,0,0), Clipped_Mesh)
    # # clip.clip_surface_with_plane(Clipped_Mesh,(7.5,0,0), (-1,0,0), Clipped_Mesh)

    # clip.clip_surface_with_plane(VTK_MeshRefinedSec,(0.25,0,0), (1,0,0), Clipped_Mesh)
    # clip.clip_surface_with_plane(Clipped_Mesh,(MaxX-0.25,0,0), (-1,0,0), Clipped_Mesh)

    # convertFile(Clipped_Mesh , VTU_Mesh)

    # print "-------------------------------------"
    # print "-------------  Finito  --------------"
    # print "-------------------------------------"









# def CreateIdealSkeleton_Orig(Directory, Generations,Height, Inlets,Outlets, HorizonatalEdgeLength, alpha):

#     # Standard length of Horizontal edges
#     L = HorizonatalEdgeLength

#     # Symetric about the yaxis
#     Nodes = np.array([[0,0]])
#     # Inlet Node
#     y = 2*L*m.sin(alpha)
#     x = L*m.cos(alpha)
 
#     for i in range(Inlets):
#         i = i-Inlets/2
#         # Do a whole hex in one shot
#         NewNodes = np.array([[0,i*y], [L,i*y]])
#         Nodes =  np.append(Nodes, NewNodes,axis=0)

#         LastNode = [L,i*y]

#         NewNodes1 = np.array([[x+LastNode[0], y/2+LastNode[1]], [x+LastNode[0], -y/2+LastNode[1]]])
#         NewNodes2 = np.array([[1+x+LastNode[0], y/2+LastNode[1]], [1+x+LastNode[0], -y/2+LastNode[1]]])
#         NewNodes3 = np.array([[1+2*x+LastNode[0], 0+LastNode[1]]])

#         Nodes =  np.append(Nodes, NewNodes1,axis=0)
#         Nodes =  np.append(Nodes, NewNodes2,axis=0)
#         Nodes =  np.append(Nodes, NewNodes3,axis=0)

    
#     Nodes = np.delete(Nodes, 0, 0)
#     maxInColumns = np.amax(Nodes, axis=0)
    
    
#     MaxX = maxInColumns[0]
#     MaxY = maxInColumns[1]
    
#     for i in range(Generations):
#         NewNode =  np.array([[MaxX+L  ,0]])
#         Nodes =  np.append(Nodes, NewNode,axis=0)
#         NewNode =  np.array([[MaxX  ,0]])
#         Nodes =  np.append(Nodes, NewNode,axis=0)

#         if(Inlets % 2) != 0:  # Odd inlets =
#             for j in range(Height+1):
#                 j=(j+1)/2
#                 # On the left
#                 NewNode =  np.array([[MaxX, y*j ]])
#                 Nodes =  np.append(Nodes, NewNode,axis=0)

#                 NewNode =  np.array([[MaxX - x, y*j-y/2 ]])
#                 Nodes =  np.append(Nodes, NewNode,axis=0)

#                 NewNode =  np.array([[MaxX, -y*j ]])
#                 Nodes =  np.append(Nodes, NewNode,axis=0)

#                 NewNode =  np.array([[MaxX - x, -y*j+y/2 ]])
#                 Nodes =  np.append(Nodes, NewNode,axis=0)

#                 # On the right
#                 NewNode =  np.array([[MaxX+L, y*j ]])
#                 Nodes =  np.append(Nodes, NewNode,axis=0)

#                 NewNode =  np.array([[MaxX + L+x, y*j-y/2 ]])
#                 Nodes =  np.append(Nodes, NewNode,axis=0)

#                 NewNode =  np.array([[MaxX+L, -y*j ]])
#                 Nodes =  np.append(Nodes, NewNode,axis=0)


#                 NewNode =  np.array([[MaxX + L+x, -y*j+y/2 ]])
#                 Nodes =  np.append(Nodes, NewNode,axis=0)

#             maxInColumns = np.amax(Nodes, axis=0)
#             MaxX = maxInColumns[0] +x+L
#             MaxY = maxInColumns[1]
#         else:
#             # Get the minimum values of each column i.e. along axis 0
#             MinY = np.amin(Nodes, axis=0)[1]
#             YMedian = (MaxY +MinY)/2
#             MaxX = MaxX
            

#             NewNode =  np.array([[MaxX+x+L, YMedian] ,[MaxX+L, YMedian-0.5*y]])
#             Nodes =  np.append(Nodes, NewNode,axis=0)

#             for j in range(Height):
#                 j=j+1
#                 # On the left
#                 NewNode =  np.array([[MaxX, 0.5*y+y*j+YMedian], [MaxX, YMedian-(0.5*y+y*j)]])
#                 Nodes =  np.append(Nodes, NewNode,axis=0)

#                 NewNode =  np.array([[MaxX-x, YMedian + y*j], [MaxX-x, YMedian-(y*j)]])
#                 Nodes =  np.append(Nodes, NewNode,axis=0)

#                  # On the right
#                 NewNode =  np.array([[MaxX+x+L, YMedian-y*j] , [MaxX+x+L, YMedian+j*y]   ])
#                 Nodes =  np.append(Nodes, NewNode,axis=0)

#                 NewNode =  np.array([[MaxX+L, YMedian-y*j-0.5*y] , [MaxX+L, YMedian+j*y+0.5*y]   ])
#                 Nodes =  np.append(Nodes, NewNode,axis=0)

#             maxInColumns = np.amax(Nodes, axis=0)
#             MaxX = maxInColumns[0] +x+L
#             MaxY = maxInColumns[1]
    
#     if  (Inlets % 2) != 0: 
#         MaxX = maxInColumns[0]-x
#         Midline = 0

#         for k in range (1):
#             for i in range(Outlets-2):
#                 i=i-1

#                 NewNode = np.array([[MaxX,y*i ], [MaxX +x,y*i +y/2]])
#                 Nodes =  np.append(Nodes, NewNode,axis=0)

#                 NewNode = np.array([[MaxX +x,y*i-y/2], [MaxX+x+L,y*i+y/2]])
#                 Nodes =  np.append(Nodes, NewNode,axis=0)
#                 # # # print "ABC"
#                 NewNode = np.array([[MaxX+x+L,-y*i-y/2]])
#                 Nodes =  np.append(Nodes, NewNode,axis=0)

#                 # NewNode = np.array([[MaxX+2*x+2*L,-y*i]])
#                 # Nodes =  np.append(Nodes, NewNode,axis=0)

#                 Vein = MaxX+L+x
    
#             # maxInColumns = np.amax(Nodes, axis=0)
#             # MaxX = maxInColumns[0]

#     else:
#             MaxX = maxInColumns[0]
#             NewNode = np.array([[MaxX+L,YMedian]])
#             N  = [MaxX+L,YMedian] 
#             Nodes =  np.append(Nodes, NewNode,axis=0)
#             for i in range(Outlets-1):
#                 i=i+1
#                 NewNode1 = np.array([[N[0]+x,N[1]+y/2*i], [N[0]+x,N[1]-y/2*i]])
#                 Nodes =  np.append(Nodes, NewNode1,axis=0)

#                 NewNode1 = np.array([[N[0],N[1]+y*i], [N[0],N[1]-y*i]])
#                 Nodes =  np.append(Nodes, NewNode1,axis=0)

#                 NewNode1 = np.array([[N[0]+x+L,N[1]+y/2*i], [N[0]+x+L,N[1]-y/2*i]])
#                 Nodes =  np.append(Nodes, NewNode1,axis=0)
            

#     Nodes = np.unique(Nodes,axis=0)
#     # Now try get the edges -- Going to loop through the Nodes and try find what is nearest 
#     # Junk edge to set up array, is removed after Edges are filled
#     Edges= np.array([[-100,-100]])
#     Nodelist = []
#     for i in Nodes:
        
#         CurrentIndex = np.where((Nodes ==i).all(axis=1))[0][0]
#                             # Nodes along the horizontal                    # Nodes digagonally left                                                                                 # Nodes digagonally right
#         NearestNodes = (   ((Nodes >= i-[L, 0]) & (Nodes <= i+[L, 0])) | ( (Nodes >= i+[L*m.cos(alpha), L*m.sin(alpha)]) & (Nodes <= i+[L*m.cos(alpha), -L*m.sin(alpha)])  ) |    ( (Nodes >= i+[L*m.cos(alpha), -L*m.sin(alpha)]) & (Nodes <= i+[L*m.cos(alpha), L*m.sin(alpha)])  )     ).all(axis=1)
#         NearestNodes = np.where(NearestNodes)[0] # np.concatenate((NearestNodes1,NearestNodes2, NearestNodes3))
        
#         for j in NearestNodes:
#             if CurrentIndex != j:
#                  Edges =  np.append(Edges, [[CurrentIndex,j]],axis=0)
     
#     Edges = np.unique(Edges,axis=0)
#     Edges = np.delete(Edges, 0, 0)
#     # Edges = np.unique(Edges,axis=0)



#     # try add the inlets and outlets here 

#     # EdgeLine1= np.array([[11, -1.73205081], [15, 1.73205081]])
#     EdgeLine1= np.array([[0, 6], [0,-6]])
#     Nodes =  np.append(Nodes , EdgeLine1,axis=0)

#     Edges =  np.append(Edges, [[len(Nodes)-2,len(Nodes)-1]],axis=0)

#     EdgeLine1= np.array([[Vein, 6], [Vein,-6]])
#     Nodes =  np.append(Nodes , EdgeLine1,axis=0)

#     Edges =  np.append(Edges, [[len(Nodes)-2,len(Nodes)-1]],axis=0)

#     with open(Directory+ 'CenterlinePoints.txt', 'w') as filehandle:
#         filehandle.writelines("%s\n" % points for points in Nodes)

#     # Set up the edges and write into a file 
#     with open(Directory+ 'CenterlineEdges.txt', 'w') as filehandle:
#         filehandle.writelines("%s\n" % Edges for Edges in Edges)
#     # pause()

