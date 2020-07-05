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

    parser = ArgumentParser(description='Remesh current geometry ')
    parser.add_argument('--ifile', default='/Users/jcrawshaw/Documents/testoutput/RemeshingComparison/WithRemeshing4th/', type=str, help='Need to supply a input folder')                                       
    args = parser.parse_args()    
#     # ' ------- Setting up args ------- '

    Directory  = args.ifile
    CenterLines_filename = args.ifile + 'PlexusCenterlines.vtp'
    SmootherCenterlinesFile = args.ifile + 'PlexusCenterlines_Smoothed.vtp'
    ScalledCenterLines = args.ifile + 'PlexusCenterlines_Smoothed.vtp' # args.ifile + 'PlexusCenterlines_Scalled.vtp'

    ScalledMeshVTK = args.ifile + 'NewMesh.vtk'
    ScalledMeshVTK_2 = args.ifile + 'SmoothedMesh.vtk'
    ScalledMeshVTU = args.ifile + 'Plexus_2.vtu'



    # COnvert geom to stl 
    vtuFile = Directory + 'results_from_time_0/mesh_1300.vtu'
    print " Convert vtu to stl "
    print vtuFile
    # # Read the VTU file from disk
    # vtu_reader = vtk.vtkXMLUnstructuredGridReader()
    # vtu_reader.SetFileName(vtuFile)
    # extract_surface_filter = vtk.vtkDataSetSurfaceFilter()
    # extract_surface_filter.AddInputConnection(vtu_reader.GetOutputPort())
    # # Write out the data in unstructured grid format
    # stl_writer = vtk.vtkSTLWriter()
    # stl_writer.SetInputConnection(extract_surface_filter.GetOutputPort())
    # stl_writer.SetFileName(Directory+'New.stl')
    # stl_writer.Write()


    
    # # # Generating a centerlines file from mesh
    GetCenterlinesCommand = 'vmtk vmtknetworkextraction -ifile ' + Directory + 'New.stl' + ' -ofile '+ CenterLines_filename
    subprocess.call(GetCenterlinesCommand, shell=True)

    # # Resamlpe the centerlines smooth centerlines with a moving average filter -->    vmtkcenterlineresampling -ifile PlexusInversed.vtp -length 20 -ofile resampledCenterlines.vtp  

    # SmoothCenterlinesCommond = 'vmtkcenterlineresampling -ifile '+ CenterLines_filename + ' -ofile '+ SmootherCenterlinesFile
    SmoothCenterlinesCommond = 'vmtkcenterlineresampling -ifile '+ CenterLines_filename + ' -length 0.0001 -ofile '+ CenterLines_filename
    subprocess.call(SmoothCenterlinesCommond, shell=True)


    # print "Developing Mesh"

    # With the scalled radii generate a new mesh from with adapted centerlines file  -- here the discretisation dimension (i.e nunber of nodes in each axis) is currently 200 200 200, but might need changing 
    command = 'vmtk vmtkcenterlinemodeller -ifile ' + CenterLines_filename +' -radiusarray Radius -dimensions 150 150 150 --pipe vmtkmarchingcubes -ofile '+ ScalledMeshVTK # -handle Self'# -dimensions '+3   /Users/jcrawshaw/docker-polnet-master/GeneratingShrunkMesh/Original2.vtk
    subprocess.call(command, shell=True)


         # # pause()
    # The Mesh is currently dense and messy, remesh to get a nicer mesh, can control the target size of each element
    command = 'vmtksurfaceremeshing -ifile '+ScalledMeshVTK +' -iterations 30  -area 100 -ofile ' +ScalledMeshVTK_2
    subprocess.call(command, shell=True)


    print "-------------  Finito  --------------"
    print "-------------------------------------"






