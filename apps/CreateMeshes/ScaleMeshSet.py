#!/usr/bin/env python
#
import subprocess
import vtk
import os
import pdb
import string
import math as m
import vmtk
# from vmtk import pypes
import array as A
from vtk import vtkXMLUnstructuredGridWriter as Writer, VTK_VERSION
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
import numpy as np
import meshio
from vtk.util import numpy_support
import vtk
import clip
import math

from vtk.numpy_interface import dataset_adapter as dsa

if __name__=="__main__":
  
    print " Scale set of meshes "
    MeshDirectory = "/data/vascrem/MeshCollection/IdealisedNetwork/CollapseOf3By3Network/UpperBranchAround6_3/"

    Collapse = [0]
    # Label = [0.585,0.595,0.615,0.605]
    Label = [0.5961, 0.5971,  0.6098, 0.6117, 0.6119, 0.6127, 0.6137, 0.6146]
    Label =[ 0.5924]
    # Collapse = [0,1,2,3,4,5,6,7]

    for i in Collapse:

        OriginalMeshFile = MeshDirectory+"mesh_"+str(i)+".vtu"  

        ScaledMesh = MeshDirectory+"ScaledMesh."+str(Label[i])+".vtu"
        # # ---- Interpolate the points in the centerlines file, this will reduce the refinment needed in the centerline modeller -------------# 
        Scale = 'vmtkmeshscaling -ifile '+ OriginalMeshFile + ' -scale 0.20  --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ ScaledMesh
        subprocess.call(Scale, shell=True)

        
        ScaledMeshstl = MeshDirectory+"ScaledMesh."+str(Label[i])+".stl"
        convert = 'meshio-convert '+ ScaledMesh +'  '+ ScaledMeshstl
        subprocess.call(convert, shell=True)
        print "Done one"
        


            
    print " Finished "
    

        # ----  Turn turn it to a stl  -------------# 
        # STL_Mesh = Directory+"Mesh"+str(int(10*i))+".stl"
        # command = "meshio-convert " + Clipped_Mesh +" " +STL_Mesh
        # subprocess.call(command, shell=True)
        #  # ----  Turn turn it to a stl  -------------# 
        # VTU_Mesh = Directory+"Mesh"+str(int(10*i))+".vtu"
        # command = "meshio-convert " + STL_Mesh +" " +VTU_Mesh
        # subprocess.call(command, shell=True)
        