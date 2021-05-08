# #!/usr/bin/env python
# # vasc_res_sim.py
# import subprocess
# import shutil
# import os
# import glob
# import numpy as np
# import time
# import pdb
# import string
# import math
# import sys
# from os import path
# import os
# import shutil
# from xml.etree import ElementTree
# import tempfile 


import numpy
from stl import mesh

def GetTheDetailsOfTheMesh(MeshFile):


        # Using an existing stl file:
        your_mesh = mesh.Mesh.from_file('/Users/jcrawshaw/Documents/testoutput/VariableEdgeLength/Scalling2/mesh_0.stl')

        # Or creating a new mesh (make sure not to overwrite the `mesh` import by
        # naming it `mesh`):
        MaxPoint = 0
        MinPoint =100

        for i in range(0,len(your_mesh.points)):
            vector = your_mesh.points[i]
            for i in [0,3,6]:
                # print vector[i]
                if vector[i]<MinPoint:
                    MinPoint = vector[i]
                elif vector[i]>MaxPoint:
                    MaxPoint = vector[i]

        Length = MaxPoint - MinPoint


        # So now the boundaries are

        LeftBoundary =  MinPoint + Length/100
        RightBoundary =  MaxPoint - Length/100
        Seed = your_mesh.points[200][0:3]
        print Length
        Result = [LeftBoundary,RightBoundary, Seed[0],Seed[1],Seed[2]]
        print Result

        return Result        
         
         
if __name__=="__main__":
    MeshFile ='/Users/jcrawshaw/Documents/testoutput/VariableEdgeLength/Scalling2/mesh_0.stl'
    Result = GetTheDetailsOfTheMesh(MeshFile)
    print Result