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




Iterations=1



chaste_setup_exe = '/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/build/optimised/TestSetupFlowInInitallyCollapsedPipeRunner'
chaste_run_exe = '/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/build/optimised/TestRunFlowInInitallyCollapsedPipeRunner'



hemelb_setup_exe = 'env PYTHONPATH=/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool:$PYTHONPATH /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool/scripts/hemelb-setup-nogui'

radii_over_time = [] # empty array



def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

    # ---------------------------------
  

# Creats a file callsed centerlines iter_number  in the file where the config.stl is saved. needs to be pointed at the stl folder 

def generate_flow_vtus(timestep):
    print "  generate_flow_vtus   "

    # Make the vtu 
    shutil.copyfile(working_directory + 'config.xml', working_directory + 'config2.xml')
    command = 'python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/GmyUnstructuredGridReader.py' + ' ' + working_directory + 'config.xml'
    subprocess.call(command, shell=True)


    # command = 'python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py' + ' ' +  working_directory +'config2.vtu' +' ' +  working_directory+ 'results/Extracted/surface-pressure.xtr' +' ' +  working_directory+ 'results/Extracted/wholegeometry-velocity.xtr'+' ' +  working_directory+ 'results/Extracted/surface-traction.xtr'+' ' +  working_directory+ 'results/Extracted/surface-tangentialprojectiontraction.xtr'
    command = 'python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py' + ' ' +  working_directory +'config.vtu' +' ' +  working_directory+ 'results/Extracted/surface-pressure.xtr' +' ' +  working_directory+ 'results/Extracted/wholegeometry-velocity.xtr'  
    #+' ' +  working_directory+ 'results/Extracted/surface-traction.xtr'
   
    subprocess.call(command, shell=True)

    print "generate_flow_vtus DONE    "


    # ---------------------------------
  
    

def vtu2stl(iter):

    print "  Convert vtu to stl    "
    # Read the VTU file from disk
    vtu_reader = vtk.vtkXMLUnstructuredGridReader()
    vtu_reader.SetFileName(working_directory + 'config.vtu')

    extract_surface_filter = vtk.vtkDataSetSurfaceFilter()
    extract_surface_filter.AddInputConnection(vtu_reader.GetOutputPort())

    # Write out the data in unstructured grid format
    stl_writer = vtk.vtkSTLWriter()
    stl_writer.SetInputConnection(extract_surface_filter.GetOutputPort())
    stl_writer.SetFileName(working_directory + 'config.stl')
    stl_writer.Write()



    print "  stl created     "

   # ---------------------------------
 
def run_hemelb_setup():
    print "HemeLB setting up"
 
   # TODO: It seems that when you use the --stl option below, the voxel size specified in .pr2 is ignored and the setup tool guesses one for you
  
    # voxel_size =  0.017262106761336327 # 0.15e-3
    # voxel_size = 1.222e-6 
    # voxel_size= 0.48330092430114746
 
    # heme_profile_filename = working_directory + 'config.pr2' 
    # ConfigDirectory = working_directory + 'config.stl'
   
    # command = hemelb_setup_exe + ' ' + heme_profile_filename + ' --voxel ' + str(voxel_size)  + ' --geometry ' + working_directory + 'config.gmy' + ' --xml ' + working_directory + 'config.xml' + ' --stl ' + ConfigDirectory 

    # subprocess.call(command, shell=True)
    # print "HemeLB set up is DONE"

    voxel_size = 0.15e-5

    heme_profile_filename = working_directory +  'config.pr2' 
    command = hemelb_setup_exe + ' ' + heme_profile_filename + ' --stl ' + working_directory + 'config.stl' + ' --voxel ' + str(voxel_size) + ' --geometry ' + working_directory + 'config.gmy' + ' --xml ' + working_directory + 'config.xml'
    subprocess.call(command, shell=True)

    

    # ---------------------------------
 

def update_xml_file(iter_num, num_iters):

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

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(500), 'file': 'surface-traction.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='traction')

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(500), 'file': 'surface-tractions.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='traction')
    ElementTree.SubElement(surface, 'field', type='tangentialprojectiontraction')    

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(500), 'file': 'surface-pressure.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='pressure')

    wholegeometry = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(500), 'file': 'wholegeometry-velocity.xtr'})
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
    f = open(Nodesfilename, "r")
    NodeLocationsList = f.readlines() # Open the file

    # Last Line of the NodeLocation list has the Node locations at the last timestep  x1 y1 z1 x2 y2 z2 ext  
    NodeLocationAtLastTime = NodeLocationsList[-1] 
    NodeLocationAtLastTime= NodeLocationAtLastTime.split(" ") # Formate correctly

    # Get the Node Location of the 200th node
    Node200Location = NodeLocationAtLastTime[597:600] 
 

    Seed_New = 'SeedPoint: {x: '+ Node200Location[0] +', y: '+ Node200Location[1] + ', z: '+ Node200Location[2] +'} \n'

    # pdb.set_trace()
    # Find the lines that the radius and the seed are defined 
    file = open(filename).readlines()
    Radius_Old=""

    for line in file:
        if 'Radius:' in line:                                                                                         
            Radius_Old = line 
        if 'SeedPoint:' in line:
            Seed_Old = line 

    newRadius = '  Radius: '+ str(Radius_New+5.5) +'\n'
   
    print  'New Seed: '+  Seed_New 

    s = open(filename).read()
    s = s.replace(Seed_Old, Seed_New)
    # s = s.replace(Radius_Old, newRadius)
    f = open(filename, 'w')
    f.write(s)
    f.close()

    print "Updated Seed and Radius in .pr2 file"

    
    # ---------------------------------

def run_hemelb():
    print "**************   hemelb step  **************   "
    shutil.rmtree(working_directory + 'results/', ignore_errors=True)
    xml_config_file = working_directory + 'config.xml'
    command = 'mpirun -np 1 hemelb -in ' + xml_config_file
    subprocess.call(command, shell=True)

    # ------------------------------------------

def run_hemelb2():

    
    print "**************   hemelb step  **************   "
    shutil.rmtree(working_directory + 'results/', ignore_errors=True)
    xml_config_file = working_directory + 'config.xml'
    command = 'mpirun -np 2 hemelb -in ' + xml_config_file
    subprocess.call(command, shell=True)

    # ------------------------------------------

    # ------------------------------------------


if __name__=="__main__":
    working_directory = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory/VesselPlexus/' 
    data_path =working_directory + 'SetUpData/'
    
    
    shutil.copy(data_path +'config.stl' , working_directory +'config.stl' )
    # shutil.copy(data_path +'config.xml' , working_directory +'config.xml' )
    # shutil.copy(data_path +'config.vtu' , working_directory +'config.vtu' )
    shutil.copy(data_path +'config.pr2' , working_directory +'config.pr2' )
    

# Create the directory 
    if not os.path.exists(working_directory):
        os.makedirs(working_directory)

# **** Set up what you need for the collected chaste mesh folder


    # Define arguments to be parsed
    parser = ArgumentParser(description='Run a vascular remodelling simulation')
    parser.add_argument('--num_iterations', default=Iterations, type=int, help='Number of Hemelb/Chaste iterations to be run (optional, default is 5).')
    parser.add_argument('--compute_radii', dest='compute_radii', action='store_true', help='Use VMTK to compute axis radii.')
    parser.add_argument('--output_postfix', dest='output_postfix', default='', help='This string will be added to ChasteWorkingDirectory to get the output folder.')
    parser.add_argument('--div_threshold', dest='div_threshold', default=1e10, help='This specifies the length that edges will divide at. (Defaults to no division, but 6e-4 is good for mm meshes')
    parser.add_argument('--mesh_scale', dest='mesh_scale', default=1e-3, help='This specifies what to scale the mesh by so that all distances are in meters (defaults to mm).')
    args = parser.parse_args()
 
    
    start_time = 0
    # duration = 0.3
    for iter in range(args.num_iterations):
    
        print "Fluid in Plexus"
        
        #
        # Step 2: Run HemeLB setup
        #
        
        print "Run HemeLB setup "

        run_hemelb_setup()

        print "Update xml file"

        # Step 2a: Update xml

        update_xml_file(iter, args.num_iterations)
    
        print "Run HemeLB"
        # Step 3: HemeLB simulation
        run_hemelb2() 
        print'\n HemeLB simulation complete \n'

        print'\n Generate flow vtus \n'

        generate_flow_vtus(3)
        
  
        # # pdb.set_trace()
        # vtu2stl(iter)


    print '\n ********* ------ Completed ------ ********* \n'

        # ---------------------------------