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
import numpy as np
import meshio
# from vtk.util import numpy_support
import vtk
import clip
from argparse import ArgumentParser


# ---------------------------------
  
def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")        

if __name__=="__main__":

    parser = ArgumentParser(description='Get the radius file ')
    
    parser.add_argument('--ifile', default='/Users/jcrawshaw/docker-polnet-master/IdealiseNetworks/3_by_4/MeshClipped4.vtk', type=str, help='Need to supply a input folder')
    parser.add_argument('--ofile', default='/Users/jcrawshaw/docker-polnet-master/IdealiseNetworks/3_by_4/MeshCoarse.vtk', type=str, help='Need to supply a input folder')
    parser.add_argument('--edgelength', default=0.5, type= float, help='Need edge lenght')
    parser.add_argument('--iters', default=10, type=int, help='Need to supply a input folder')
    args = parser.parse_args()    
    # ' ------- Setting up args ------- '

    Directory  = args.ifile #+ 'results_from_time_0/'
    InputFile = args.ifile
    OutputFile =  args.ofile
    EdgeLength = str(args.edgelength)
    Iterations = str(args.iters)

    #The Mesh is currently dense and messy, remesh to get a nicer mesh, can control the target size of each element
    command = 'vmtksurfaceremeshing -ifile '+InputFile +' -iterations '+Iterations+' -edgelength '+EdgeLength+' -elementsizemode "edgelength" -ofile ' +OutputFile
    subprocess.call(command, shell=True)

    print "-------------------------------------"
    print "-------------  Finito  --------------"
    print "-------------------------------------"






