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




Iterations=1



chaste_setup_exe = '/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/build/optimised/TestSetupFlowInInitallyCollapsedPipeRunner'
chaste_run_exe = '/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/build/optimised/TestRunFlowInInitallyCollapsedPipeRunner'

hemelb_setup_exe = 'env PYTHONPATH=/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool:$PYTHONPATH /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool/scripts/hemelb-setup-nogui'

radii_over_time = [] # empty array



# def pause():
#     programPause = raw_input("Press the <ENTER> key to continue...")

    # ---------------------------------
  

# Creats a file callsed centerlines iter_number  in the file where the config.stl is saved. needs to be pointed at the stl folder 
def GetRadius(CenterLines_filename, NewOutput): 
    print "**************   vmtk_compute_stl_radii  **************  "
    
    # CenterLines_filename = working_directory + 'centerlines' +  str(iter_number) + '.vtp'
    # NewOutput = working_directory + 'NewFile' +  str(iter_number) + '.vtp'

    # Get the center line file 
    # # command = 'vmtk vmtknetworkextraction -ifile ' + working_directory + 'config.stl' + ' -ofile ' + CenterLines_filename 
    # subprocess.call(command, shell=True)

    # pdb.set_trace()
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(CenterLines_filename)
    reader.Update()
    point_data = reader.GetOutput().GetPointData()
    assert point_data.GetArrayName(0) == 'Radius', "VTP file doesn't contain array Radius"
    radii = point_data.GetArray(0)
    # radii_over_time.append(np.array([radii.GetValue(i) for i in range(radii.GetSize())]))
    #  pause():
    print "    vmtk_compute_stl_radii DONE   "

  

if __name__=="__main__":
    parser = ArgumentParser(description='Get the radius file ')
    parser.add_argument('--Geometryfile', type=str, help='Need to supply a config.stl file (directory)')
    parser.add_argument('--ifile', type=str, help='Need to supply a centerlines.vtp file (directory)')
    parser.add_argument('--ofile', default='~/Documents/ChasteWorkingDirectory/Radius.vtp', type=str, help='Need to supply an output directory for the new radius file(directory)')
    args = parser.parse_args()
    # print ' ------- Setting up args ------- '
    # print args.ifile
    # print args.ofile
    CenterLines_filename = args.ifile
    reader = vtk.vtkXMLPolyDataReader()
    writer = vtk.vtkXMLPolyDataWriter()
    reader.SetFileName(CenterLines_filename)
    writer.SetFileName(CenterLines_filename)

    reader.Update()
    point_data = reader.GetOutput().GetPointData()
    assert point_data.GetArrayName(0) == 'Radius', "VTP file doesn't contain array Radius"
    radii = point_data.GetArray(0)
    # print radii
    Scalling =10
    radii_over_time = []
    radii_over_time.append([radii.GetValue(i)/Scalling for i in range(radii.GetSize())])
    # print radii_over_time[0]
    print len(radii_over_time[0])
    writer.SetRadius(radii_over_time[0])


    

    # point_data.SetupOutputData = radii_over_time[0]

    # reader.Update()
    # # assert point_data.GetArrayName(1) == 'Radius', "VTP file doesn't contain array Radius"
    # radii = point_data.GetArray(0)

    # NewArray = []
    # NewArray.append([radii.GetValue(i) for i in range(radii.GetSize())])
    # print NewArray


    


    # point_data.GetArra(0)  = radii_over_time[0]

    # vmtk vmtkcenterlinemodeller -ifile ~/Documents/ChasteWorkingDirectory/VesselPlexus/JessNewCenterLinesFile.vtp -radiusarray 1
   
    # arr = np.array(radii_over_time[0]) 
    # print arr
    # # np.savetxt('array_hf.csv', [arr], delimiter=',', fmt='%d' , header='A Sample 2D Numpy Array :: Header', footer='This is footer')
    # Save Numpy array to csv
    # np.savetxt('array.csv', radii_over_time[0], delimiter=',', fmt='%d')

    # # # # # try generate the surface or the mesh 
    # # # # theRadius =  str(radii_over_time[0])
    # command = 'vmtk vmtkcenterlinemodeller -ifile ' + CenterLines_filename + ' -radiusarray Radius --pipe vmtkmarchingcubes -ofile ' + args.ofile
    # # command = 'vmtk vmtkcenterlinemodeller -ifile ' + CenterLines_filename + ' -radiusarray ' +str(radii_over_time[0])+' -ofile ' + args.ofile
    # subprocess.call(command, shell=True)


    print "---Finito---"





    
   