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
  
    server = 1
    if server: 
        CPP_Centerlines_vtp_writer = "/home/vascrem/Chaste/projects/VascularRemodelling/build/optimised/GenerateIdealVascularMesh/Test_VTP_writerRunner"
        Directory = "/home/vascrem/MeshCollection/IdealisedNetwork/DeflatedNetwork_2/"
    else:
        CPP_Centerlines_vtp_writer = "/Users/jcrawshaw/Chaste/projects/VascularRemodelling/build/optimised/GenerateIdealVascularMesh/Test_VTP_writerRunner"
        Directory = "/Users/jcrawshaw/Documents/Projects/MeshCollection/IdealiseNetworks/DelfatedNetwork/"

    
    Input = Directory+"Clipped.vtk"
    Output = Directory+"SmallBifucation.vtk"
    clip.clip_surface_with_plane(Input,(80.5,-14,0), (1,0,0), Output)
    clip.clip_surface_with_plane(Output,(96.2276112912127,-18,0), (0,-1,0), Output)
    clip.clip_surface_with_plane(Output,(95,-20,-0.05), (-0.7,-0.7,0), Output)
    clip.clip_surface_with_plane(Output,(95,-34,-0.05), (-0.7,0.7,0), Output)

    print "-------------------------------------"
    print "-------------  Finito  --------------"
    print "-------------------------------------"