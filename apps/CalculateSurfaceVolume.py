#!/usr/bin/env python
# vasc_res_sim.py
import subprocess
import vtk
import shutil
import os
from xml.etree import ElementTree
import glob
from argparse import ArgumentParser
import numpy as np
import time
import matplotlib.pyplot as plt
import csv
import pdb
import string
import math
from vmtk import pypes
import matplotlib.pyplot as plt



Iterations=1



chaste_setup_exe = '/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/build/optimised/TestSetupFlowInInitallyCollapsedPipeRunner'
chaste_run_exe = '/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/build/optimised/TestRunFlowInInitallyCollapsedPipeRunner'

hemelb_setup_exe = 'env PYTHONPATH=/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool:$PYTHONPATH /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool/scripts/hemelb-setup-nogui'

radii_over_time = [] # empty array



# def pause():
#     programPause = raw_input("Press the <ENTER> key to continue...")

    # ---------------------------------
  
def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

# Creats a file callsed centerlines iter_number  in the file where the config.stl is saved. needs to be pointed at the stl folder 
def GetVolume(filename): 
    # print "**************   Compute Volume  **************  "
    
    # vmtk vmtksurfacemassproperties -ifile 'config.stl'

    command = 'vmtk vmtksurfacemassproperties -ifile ' + filename
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=None, shell=True)

    #Launch the shell command:
    output = process.communicate()
    
    SAIndex = output[0].find('SurfaceArea')

    VIndex = output[0].find('Volume')
    SIndex = output[0].find('ShapeIndex')

    
    SurfaceArea = float(output[0][SAIndex+14:VIndex-5])
    Volume =  float(output[0][VIndex+9:SIndex-5])
    return [Volume, SurfaceArea]
    # pdb.set_trace()

    # print "    vmtk_compute_stl_radii DONE   "


def vtu2stl(file,WorkingDirectory , OutPutDirectory, iter):

    # print "  Convert vtu to stl    "
    # Read the VTU file from disk
    vtu_reader = vtk.vtkXMLUnstructuredGridReader()
    vtu_reader.SetFileName(WorkingDirectory +file)

    extract_surface_filter = vtk.vtkDataSetSurfaceFilter()
    extract_surface_filter.AddInputConnection(vtu_reader.GetOutputPort())
    # Write out the data in unstructured grid format
    stl_writer = vtk.vtkSTLWriter()
    stl_writer.SetInputConnection(extract_surface_filter.GetOutputPort())
    stl_writer.SetFileName(OutPutDirectory + 'mesh_' +str(iter)+'.stl')
    stl_writer.Write()


    # print "  stl created     "

   # ---------------------------------
 
# '~/Documents/testoutput/PlexusDeformation/Homogeneous_HighPressure/results_from_time_0/'
  

if __name__=="__main__":
    parser = ArgumentParser(description='Get the radius file ')
    parser.add_argument('--ifile', default='/Users/jcrawshaw/Documents/testoutput/PlexusDeformation/Homogeneous_HighPressure/', type=str, help='Need to supply a input folder')
    args = parser.parse_args()
    
    # print ' ------- Setting up args ------- '
    # print args.ifile
    # print args.ofile
    WorkingDirectory  = args.ifile + 'results_from_time_0/'
    OutPutDirectory  = args.ifile
    # Need to turn the vtus into stls 
    iter = 1
    for file in os.listdir(WorkingDirectory):
        if file.startswith("mesh"):
            vtu2stl(file, WorkingDirectory , OutPutDirectory,iter  )
            iter = iter +1
    
    # # Now take these stls and find the volume

    Data = []
    for j in range(1, iter):
        filename = OutPutDirectory + 'mesh_'+str(j) +'.stl' 

        Data.append(GetVolume(filename))


    # Now I need to write the csv file to save this data 

    print "makeing csv"
    # in this csv file the left array will be volume, and the right is surface area
    csvfile = open(OutPutDirectory + 'VolumeAndSurfaceAreaOfCapillary.csv', 'wb')
    writer = csv.writer(csvfile)
    for j in range(0, iter-1): 
        writer.writerow(Data[j])
    csvfile.close()

    print "---Finito---"







    
   