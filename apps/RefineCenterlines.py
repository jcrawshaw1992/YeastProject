#!/usr/bin/env python
#
import subprocess
import vtk
import os
from argparse import ArgumentParser
import pdb
import string
import math
import vmtk
from vmtk import pypes

from vtk import vtkXMLUnstructuredGridWriter as Writer, VTK_VERSION
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
import numpy
import meshio


    # ---------------------------------
  
def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

if __name__=="__main__":
    parser = ArgumentParser(description='Get the radius file ')
    
    parser.add_argument('--ifile', default='/Users/jcrawshaw/docker-polnet-master/PlexusNotScalled/', type=str, help='Need to supply a input folder')
    args = parser.parse_args()    
    # ' ------- Setting up args ------- '

    CenterLines_filename = args.ifile + 'PlexusCenterlines.vtp'
    SmootherCenterlinesFile = args.ifile + 'PlexusCenterlinesRefined.vtp'
       # Resamlpe the centerlines smooth centerlines with a moving average filter -->    vmtkcenterlineresampling -ifile PlexusInversed.vtp -length 20 -ofile resampledCenterlines.vtp  
    SmoothCenterlinesCommond = 'vmtkcenterlineresampling -ifile '+ CenterLines_filename + ' -length 0.05 -ofile '+ SmootherCenterlinesFile
    subprocess.call(SmoothCenterlinesCommond, shell=True)

    # SmoothCenterlinesCommond = 'vmtkcenterlineresampling -ifile '+ SmootherCenterlinesFile + ' -length 0.001 -ofile '+ SmootherCenterlinesFile
    # subprocess.call(SmoothCenterlinesCommond, shell=True)
    # subprocess.call(SmoothCenterlinesCommond, shell=True)
    # subprocess.call(SmoothCenterlinesCommond, shell=True)
    # subprocess.call(SmoothCenterlinesCommond, shell=True)






