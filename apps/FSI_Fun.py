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


hemelb_setup_exe = 'env PYTHONPATH=/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool:$PYTHONPATH /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool/scripts/hemelb-setup-nogui'


def SetWorkingDirectory(WorkingDirectory):
    global working_directory
    working_directory = WorkingDirectory

def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

    # ---------------------------------
  

# Creats a file callsed centerlines iter_number  in the file where the config.stl is saved. needs to be pointed at the stl folder 
def vmtk_compute_stl_radii(radii_over_time, iter_number): 
    print "**************   vmtk_compute_stl_radii  **************  "
    
    CenterLines_filename = working_directory + 'centerlines' +  str(iter_number) + '.vtp'
    NewOutput = working_directory + 'NewFile' +  str(iter_number) + '.vtp'

    # Get the center line file 
    command = 'vmtk vmtknetworkextraction -ifile ' + working_directory + 'config.stl' + ' -ofile ' + CenterLines_filename 
    subprocess.call(command, shell=True)

    # pdb.set_trace()
    reader = vtk.vtkXMLPolyDataReader()
    #Read in the centerlines file
    reader.SetFileName(CenterLines_filename)
    reader.Update()
    point_data = reader.GetOutput().GetPointData()
    assert point_data.GetArrayName(0) == 'Radius', "VTP file doesn't contain array Radius"
    radii = point_data.GetArray(0)

    radii_over_time.append(np.array([radii.GetValue(i) for i in range(radii.GetSize())]))
    RadArr = radii_over_time[len(radii_over_time)-1]
    RadArr = RadArr[RadArr != 0]
    # MaxRadi = np.amax(radii)
    MinRadi = np.amin(RadArr)
    # np.mean(RadArr)
    
    # pdb.set_trace()
    return MinRadi

    

def generate_flow_vtus(working_directory):
    print "  generate_flow_vtus   "

    # Make the vtu 
    shutil.copyfile(working_directory + 'config.xml', working_directory + 'config2.xml')
    command = 'python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/GmyUnstructuredGridReader.py' + ' ' + working_directory + 'config2.xml'
    subprocess.call(command, shell=True)

    # command = 'python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py' + ' ' +  working_directory +'config2.vtu' +' ' +  working_directory+ 'results/Extracted/surface-pressure.xtr' +' ' +  working_directory+ 'results/Extracted/wholegeometry-velocity.xtr'+' ' +  working_directory+ 'results/Extracted/surface-traction.xtr'+' ' +  working_directory+ 'results/Extracted/surface-tangentialprojectiontraction.xtr'
    command = 'python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py' + ' ' +  working_directory +'config2.vtu' +' ' +  working_directory+ 'results/Extracted/surface-pressure.xtr' +' ' +  working_directory+ 'results/Extracted/wholegeometry-velocity.xtr'+' ' +  working_directory+ 'results/Extracted/surface-traction.xtr'
   
    subprocess.call(command, shell=True)

    print "generate_flow_vtus DONE    "
   # ---------------------------------
 
def run_hemelb_setup():
    print "HemeLB setting up"
 
   # TODO: It seems that when you use the --stl option below, the voxel size specified in .pr2 is ignored and the setup tool guesses one for you
  
    voxel_size =  0.017262106761336327 # 0.15e-3
    # VoxelSize = 0.48330092430114746
 
    heme_profile_filename = working_directory + 'config.pr2' 
    ConfigDirectory = working_directory + 'config.stl'
   
    command = hemelb_setup_exe + ' ' + heme_profile_filename #+ ' --voxel ' + str(voxel_size) # + ' --geometry ' + working_directory + 'config.gmy' + ' --xml ' + working_directory + 'config.xml'

    subprocess.call(command, shell=True)
    print "HemeLB set up is DONE"

    # ---------------------------------
 

def update_xml_file(period):

    print "update_xml_file"
    # Load automatically generated XML file
    filename = working_directory + 'config.xml'
    
    tree = ElementTree.parse(filename)
    root = tree.getroot()

    # Add monitoring of incompressibility and convergence
    monitoring = ElementTree.SubElement(root, 'monitoring')
    ElementTree.SubElement(monitoring, 'incompressibility') 
    convergence = ElementTree.SubElement(monitoring, 'steady_flow_convergence', {'tolerance': '1e-3', 'terminate': 'false'})
    ElementTree.SubElement(convergence, 'criterion', {'type': 'velocity', 'value': '0.01', 'units': 'm/s'})
    
    # Add definition of properties to be extracted
    extr = ElementTree.SubElement(root, 'properties')

    # surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-tangentialprojectiontraction.xtr'})
    # ElementTree.SubElement(surface, 'geometry', type='surface')
    # ElementTree.SubElement(surface, 'field', type='tangentialprojectiontraction')

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(period), 'file': 'surface-traction.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='traction')

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(period), 'file': 'surface-tractions.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='traction')
    ElementTree.SubElement(surface, 'field', type='tangentialprojectiontraction')    

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(period), 'file': 'surface-pressure.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='pressure')

    wholegeometry = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(period), 'file': 'wholegeometry-velocity.xtr'})
    ElementTree.SubElement(wholegeometry, 'geometry', type='whole')
    ElementTree.SubElement(wholegeometry, 'field', type='velocity')

    # Save XML file to disk
    tree.write(filename)
    print "update_xml_file DONE"

    # ---------------------------------

def update_pr2_file(Radius_New):
    filename = working_directory + 'config.pr2'
    # Find the location of node 200, which is my seed, must make this more generailable so I can select my seeding node
    Nodesfilename = working_directory + 'ChasteMeshes/results.viznodes'
    MutationCells = working_directory + 'ChasteMeshes/MutationStates'
    f = open(Nodesfilename, "r")
    fMut = open(MutationCells, "r")
    NodeLocationsList = f.readlines() # Open the file
    NodeMutations = fMut.readlines()

    NodeInFluidAndStructureSimulations = (NodeMutations[0][3:len(NodeMutations[0])-2].index('4'))+1 # Mutation state 2 detects nodes outside the the fluid simulation/ outside the boundaries, so 
    # a node in both is has a mutation state 4. Also cant be node 0 (the first node) because the x coord is attached to the time mark, annoying and will throw problems

    # The first two items are the time and '\t', so remove them from the count
    XCoordIndex = 3*(NodeInFluidAndStructureSimulations)
    print XCoordIndex

    # Last Line of the NodeLocation list has the Node locations at the last timestep  x1 y1 z1 x2 y2 z2 ext  
    NodeLocationAtLastTime = NodeLocationsList[-1] 
    NodeLocationAtLastTime= NodeLocationAtLastTime.split(" ") # Formate correctly

    # Get the Node Location of the 200th node
    # Node200Locafindtion = NodeLocationAtLastTime[597:600] 
    
    NodeLocation = NodeLocationAtLastTime[XCoordIndex:XCoordIndex+3]
    NodeLocation[0] = str(float(NodeLocation[0])*1e3)
    NodeLocation[1] = str(float(NodeLocation[1])*1e3)
    NodeLocation[2] = str(float(NodeLocation[2])*1e3)
    
        
    Seed_New = 'SeedPoint: {x: '+ NodeLocation[0] +', y: '+ NodeLocation[1]  + ', z: '+ NodeLocation[2] +'} \n'
    
    # pdb.set_trace()
    # Find the lines that the radius and the seed are defined 
    file = open(filename).readlines()
    Radius_Old=""

    for line in file:
        if 'Radius:' in line:                                                                                         
            Radius_Old = line 
        if 'SeedPoint:' in line:
            Seed_Old = line 
        if 'VoxelSize:' in line:
            Voxel_Old = line 
        if 'TimeStepSeconds:' in line:
            Time_Old = line 

    NewRad = '  Radius: '+ str(Radius_New+4.5) +'\n'
   
    # print  'New Seed: '+  Seed_New 
    dx = 2*Radius_New/15
    print Radius_New
    dt = 0.1/4*dx*dx
    print dx
    print dt
    VoxelSize = 'VoxelSize: ' + str(dx)
    Time_new = 'TimeStepSeconds: '+  str(dt) +'\n'
     
    s = open(filename).read()
    s = s.replace(Seed_Old, Seed_New)
    s = s.replace(Time_Old, Time_new)
    s = s.replace(Voxel_Old, VoxelSize)
    for line in file:
        if 'Radius:' in line:                                                                                         
            Radius_Old = line 
            s = s.replace(line, NewRad)

    
    # s = s.replace(Radius_Old, newRadius)
    f = open(filename, 'w')
    f.write(s)
    f.close()
  
    print "Updated Seed and Radius in .pr2 file"
    return dt
    
    # ------------------------------------------

def run_hemelb():

    
    print "**************   hemelb step  **************   "
    shutil.rmtree(working_directory + 'results/', ignore_errors=True)
    xml_config_file = working_directory + 'config.xml'
    command = 'mpirun -np 2 hemelb -in ' + xml_config_file
    subprocess.call(command, shell=True)
    KillCommand = 'Stop-Process -Name "hemelb" -Force'

    # ------------------------------------------

def WriteReadme(ElasticShearModulus, AreaDilationModulus, membrane_constant, Area_constant , Time_step, Sampling, duration, working_directory ):
    print "**************   WriteReadme  **************   "

    file1 = open(working_directory +"readme.txt","w") 
    file1.write(" --- NOTE - PIPE NOT INITIALLY COLLAPSED IN THIS SIMULATION, NEED TO RENAME FOLDER --- \n\n")
    file1.write(" --- Membrane Constants --- \n\n")
    file1.write( ElasticShearModulus+ "\n" )
    file1.write( AreaDilationModulus +"\n" )
    file1.write( membrane_constant +"\n" )
    file1.write( Area_constant +"\n\n" )
    file1.write( " -- Chaste parameters -- \n\n" )
    file1.write( Time_step +"\n"  )
    file1.write( Sampling +"\n" )
    file1.write( duration +"\n" )
    file1.close() 

    # ------------------------------------------




# # Create the directory 
#     if not os.path.exists(working_directory):
#         os.makedirs(working_directory)