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

if __name__=="__main__":
  
    MeshDirectory = "/data/vascrem/MeshCollection/IdealisedNetwork/CollapseOf3By3Network/UpperBranch/"
       
    # Collapse = ['1','2','3','4','5','6','7','8','9','10','0','5.9','6','6.1','5.5','5.6','5.7','5.8','6.2','6.3','6.4','6.5']
    Collapse = ['0','0.1','0.2','0.3','0.4','0.5','0.55','0.56','0.57','0.58','0.59','0.6','0.61','0.62','0.63','0.64','0.65','0.7','0.8','0.9','1',]
    
    for i in Collapse:
        
        vtuFile = MeshDirectory+"ScaledMesh."+i+".vtu"
        Input = MeshDirectory+"TempScaledMesh."+i+".vtk"

        # command = 'meshio-convert ' +vtuFile + ' '+Input
        # subprocess.call(command, shell=True)

        VTK_Meshremeshed = MeshDirectory+"MeshClipped"+str(i)+".vtk"
        Output = MeshDirectory+"TempClippedMesh."+i+".vtk"
        # clip.clip_surface_with_plane(VTK_Meshremeshed,(1.0704000057295497/0.2,0,-1/0.2), (-1,0,0), Output)
        # clip.clip_surface_with_plane(Output,(0.7436938976296382/0.2,0.00235322834737794/0.2,-0.00377883645399601/0.2), (1,0,0), Output)
        # clip.clip_surface_with_plane(Output,(1/0.2,-0.06/0.2,-0.01/0.2), (0,1,0), Output)

        stlFile = MeshDirectory+"ClippedMesh."+i+".vtu"
        command = 'meshio-convert ' +Output + ' '+stlFile
        subprocess.call(command, shell=True)

    print "-------------------------------------"
    print "-------------  Finito  --------------"
    print "-------------------------------------"