#!/usr/bin/env python
#
import subprocess
import vtk
import os
import pdb
import string
import math as m
# import vmtk
# from vmtk import pypes
import array as A
# from vtk import vtkXMLUnstructuredGridWriter as Writer, VTK_VERSION
# from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
# import numpy as np
import meshio
# from vtk.util import numpy_support
# import vtk
# import clip
import math

# from vtk.numpy_interface import dataset_adapter as dsa

if __name__=="__main__":
  
    print " Convert stl to vtu"
    print " Jess is good"
    Directory = "/Users/jcrawshaw/Downloads/meshes/"

    # Collapse = ['0.6','0.55','0.56','0.57','0.58','0.59','0.61','0.62','0.63','0.64','0.65']

    Collapse = [0,1,2,3,4]


    for i in Collapse:

        ScaledMeshstl = Directory+"mesh."+str(i)+".stl"
        OriginalMeshFile = Directory+"mesh_"+str(i)+".vtu"
        convert = 'meshio-convert '+ ScaledMeshstl +'  '+ OriginalMeshFile
        subprocess.call(convert, shell=True)

    print " Finished "
    

