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





Iterations=20

working_directory = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory/DeformingCylindricalMembraneWithFluidFlow/' 
data_path =working_directory + 'SetUpData/'

chaste_setup_exe = '/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/build/optimised/TestSetupFlowInHetroPipeRunner'
chaste_run_exe = '/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/build/optimised/TestRunFlowInHetroPipeRunner'

hemelb_setup_exe = 'env PYTHONPATH=/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool:$PYTHONPATH /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool/scripts/hemelb-setup-nogui'

radii_over_time = [] # empty array



def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

    # ---------------------------------
  

# Creats a file callsed centerlines iter_number  in the file where the config.stl is saved. needs to be pointed at the stl folder 
def vmtk_compute_stl_radii(iter_number): 
    print "**************   vmtk_compute_stl_radii  **************  "
    output_filename = working_directory + 'centerlines' +  str(iter_number) + '.vtp'

    command = 'vmtk vmtknetworkextraction -ifile ' + working_directory + 'config.stl' + ' -ofile ' + output_filename
 
    subprocess.call(command, shell=True)
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(output_filename)
    reader.Update()
    point_data = reader.GetOutput().GetPointData()
    assert point_data.GetArrayName(0) == 'Radius', "VTP file doesn't contain array Radius"
    radii = point_data.GetArray(0)
    radii_over_time.append(np.array([radii.GetValue(i) for i in range(radii.GetSize())]))
    print "    vmtk_compute_stl_radii DONE   "

    # ---------------------------------

    

def generate_flow_vtus(timestep):
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
  
    voxel_size =  0.017262106761336327 # 0.15e-3
    # VoxelSize = 0.48330092430114746
 
    heme_profile_filename = working_directory + 'config.pr2' 
    ConfigDirectory = working_directory + 'config.stl'
   
    command = hemelb_setup_exe + ' ' + heme_profile_filename #+ ' --voxel ' + str(voxel_size) # + ' --geometry ' + working_directory + 'config.gmy' + ' --xml ' + working_directory + 'config.xml'

    subprocess.call(command, shell=True)
    print "HemeLB set up is DONE"

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

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-traction.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='traction')

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-tractions.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='traction')
    ElementTree.SubElement(surface, 'field', type='tangentialprojectiontraction')    

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-pressure.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='pressure')

    wholegeometry = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'wholegeometry-velocity.xtr'})
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
    for line in f.readlines(): # get the last entry in the data, this will be the cell locations at the last time step 
        A = line
    f.close()
    a = A.split(" ")[:-1]
    Node200Location = a[597:600]

    Seed_New = 'SeedPoint: {x: '+Node200Location[0]+', y: '+Node200Location[1]+ ', z: '+Node200Location[2]+'} \n'


    # Find the lines that the radius and the seed are defined 
    file = open(filename).readlines()
    Radius_Old=""

    for line in file:
        if 'Radius:' in line:                                                                                         
            Radius_Old = line 
        if 'SeedPoint:' in line:
            Seed_Old = line 

    newRadius = '  Radius: '+ str(Radius_New+5.5) +'\n'
    # Seed_New = 'SeedPoint: {x: 0.0, y: ' + str(Radius_New)+ ', z: 0.5} \n'


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


if __name__=="__main__":
    shutil.copy(data_path +'config.stl' , working_directory +'config.stl' )
    shutil.copy(data_path +'config.xml' , working_directory +'config.xml' )
    shutil.copy(data_path +'config.vtu' , working_directory +'config.vtu' )
    shutil.copy(data_path +'config.pr2' , working_directory +'config.pr2' )
    

# Create the directory 
    if not os.path.exists(working_directory):
        os.makedirs(working_directory)

# **** Set up what you need for the collected chaste mesh folder

    newpath = working_directory +'ChasteMeshes/' 
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    else:
        shutil.rmtree(newpath)
        os.makedirs(newpath)

    Hemepath = working_directory +'HemeLBFluid/' 
    if not os.path.exists(Hemepath):
        os.makedirs(Hemepath)
    else:
        shutil.rmtree(Hemepath)
        os.makedirs(Hemepath)

    # pr2_scripts.Format_pr2_files(working_directory)
    # pause()
   
    duration =0.2
    dt= 0.002
    SamplingTimestepMultiple = 100
    TimeStep = dt * SamplingTimestepMultiple
    EndMeshFileNumber = duration /dt
    # pdb.set_trace()
    range(0, 600,100)
    list = range(0, int(EndMeshFileNumber + 1), SamplingTimestepMultiple)
    DurationSeconds = 2
    Hemelist = range(5000, DurationSeconds*10000 + 1, 5000)
    i=0
    k=0
    

    # Define arguments to be parsed
    parser = ArgumentParser(description='Run a vascular remodelling simulation')
    parser.add_argument('-o', dest='overwrite_results', action='store_false', help='Allow overwriting of existing results folder.')
    parser.add_argument('--hemelb_vtu', dest='flow_to_vtu', action='store_true', help='Create VTU files with flow results.')
    parser.add_argument('--num_iterations', default=Iterations, type=int, help='Number of Hemelb/Chaste iterations to be run (optional, default is 5).')
    parser.add_argument('--compute_radii', dest='compute_radii', action='store_true', help='Use VMTK to compute axis radii.')
    parser.add_argument('--output_postfix', dest='output_postfix', default='', help='This string will be added to ChasteWorkingDirectory to get the output folder.')
    parser.add_argument('--div_threshold', dest='div_threshold', default=1e10, help='This specifies the length that edges will divide at. (Defaults to no division, but 6e-4 is good for mm meshes')
    parser.add_argument('--mesh_scale', dest='mesh_scale', default=1e-3, help='This specifies what to scale the mesh by so that all distances are in meters (defaults to mm).')
    args = parser.parse_args()
    
    # Parse arguments (this will create args.flow_to_vtu etc. variables)
    
    os.environ['CHASTE_TEST_OUTPUT'] = working_directory
    
    #
    # Step 1: Run preliminary Chaste setup
    #
    mesh_filename = working_directory + 'config.vtu' # Use proper path concatentation
    xml_filename = working_directory +'config.xml' # Use proper path concatentation
    # ElasticShearModulus = 0
    print ' ------- Setting up Chaste simulation ------- '

     # Hyperelastic energy 
    ElasticShearModulus = ' -ElasticShearModulus ' +  str(0.5 * 1e-4) #  1e-4 Dont change this, stops werid stuff happening 
    AreaDilationModulus = ' -AreaDilationModulus ' + str(1e-9) # str(5e-9) this one controls the dilation of the cylinder, so can be played with 
    # Bending energy 
    membrane_constant = ' -membrane_constant ' + str(1* 1e-15)
    # Area energy 
    Area_constant = ' -Area_constant ' + str(1 * 1e-6) # str(5 * 1e-6)
    chaste_setup_call = chaste_setup_exe +' -mesh ' + mesh_filename + ' -xml ' + xml_filename + ' -div_threshold ' +  str(args.div_threshold) + ' -mesh_scale ' +  str(args.mesh_scale) + ' -wss_threshold ' + str(0.025)  + ElasticShearModulus + AreaDilationModulus + membrane_constant + Area_constant  
    subprocess.call(chaste_setup_call, shell=True)
  
    print ' ------- Chaste is set up ------- '
   
    
   # Step 1.a: Prepare mesh for HemeLB setup tool and compute radii along the axes of the domain (if requested)
    
    vtu2stl(0)

    # 

    start_time = 0
    duration = 0.5
    for iter in range(args.num_iterations):
        # pdb.set_trace()
        # print "About to update pr2 file"
        # # 
        # print 'Iteration ' + str(iter)
        # #
        # # Step 2: Run HemeLB setup
        # #
        
        # print "About to run set up "

        # run_hemelb_setup()
        # print "About to update xml file"
  
        # # Step 2a: Update xml
 
        # update_xml_file(iter, args.num_iterations)
        # print "About to  run HemeLB"


        # # Step 3: HemeLB simulation
        # run_hemelb() 
        # print'\n HemeLB simulation complete \n'
  
        # Step 4: Run Chaste with the tractions previously computed
        traction_filename = working_directory + 'results/Extracted/surface-tractions.xtr'
        command = chaste_run_exe + ' -start_time ' + str(start_time) + ' -duration ' + str(duration) + ' -traction_file ' +  traction_filename + ' -mesh_scale ' +  str(args.mesh_scale)
        print 'Chaste step'
        subprocess.call(command, shell=True)
        
        # start_time += duration 
    
    
        # Step 5: Convert Chaste output (vtu file) into the input of the next iteration (stl file required by hemelb setup tool)
        vtu2stl(iter)
       
        # Step 5a: compute radii
        vmtk_compute_stl_radii(iter)

        # Update the radius of the tube*, this can then be used to update the pr2 file to enable the flow to be correctly calculated 
        LengthOfArray = len(radii_over_time[iter-1])        
        
        # Step 6: Collect Chaste Mesh outputs 
        print start_time
        if (float(start_time).is_integer()):
            start_time = int(start_time)
            print 'yes'

        print start_time

        # need an iff here looking for an int 
        
        OldFilePath = working_directory +'results_from_time_' + str(start_time) + '/'
        shutil.move(OldFilePath +'results.viznodes' , newpath+'results.viznodes')

        for j in list:
            shutil.move(OldFilePath +'mesh_'+ str(j) +'.vtu' , newpath+'mesh_'+ str(i) +'.vtu'  )
            shutil.move(OldFilePath +'results_'+ str(j) +'.vtu' , newpath+'results_'+ str(i) +'.vtu'  )   
            i+=1
        # shutil.rmtree(OldFilePath) 

        start_time += duration
       
        # Radius_New = radii_over_time[iter-1][LengthOfArray-1] 
        # print Radius_New
        # update_pr2_file(Radius_New)
       

        # generate_flow_vtus(1)
        # OldHemePath = working_directory +'results/Extracted/'
        # for j in Hemelist:
        #     shutil.move(OldHemePath +'surface-pressure_'+ str(j) +'.vtu' , Hemepath +'surface-pressure_'+ str(k) +'.vtu' )
        #     shutil.move(OldHemePath +'surface-traction_'+ str(j) +'.vtu' , Hemepath +'surface-traction_'+ str(k) +'.vtu' )
        #     shutil.move(OldHemePath +'wholegeometry-velocity_'+ str(j) +'.vtu' , Hemepath +'wholegeometry-velocity_'+ str(k) +'.vtu' )
        #     # shutil.move(OldHemePath +'surface-tangentialprojectiontraction_'+ str(j) +'.vtu' , Hemepath +'surface-tangentialprojectiontraction_'+ str(k) +'.vtu' )
        #     k+=1


    

    
    csvfile = open(working_directory + 'radii_over_time.csv', 'wb')
    writer = csv.writer(csvfile)
    for iter_num, radii in enumerate(radii_over_time):
        time = iter_num*duration # The end of the n-th iteration
        plt.plot(radii, label='t={}'.format(time))
        line = [time]
        line.extend(radii)
        writer.writerow(line)
    plt.legend()
    plt.savefig(working_directory + 'radii_over_time.png')
    csvfile.close()
    

    print '\n ********* ------ Completed ------ ********* \n'





#   Old stuff
  # output_folders = glob.glob(working_directory + 'results_from_time_*')

    # for folder in output_folders:
    #     folder_timestep = folder.split('results_from_time_')[-1]
    #     if float(folder_timestep) == timestep:
    #         shutil.copyfile(working_directory +  'results/Extracted/surface-pressure_2000.vtu', folder + '/surface-pressure.vtu')
    #         shutil.copyfile(working_directory +  'results/Extracted/surface-traction_2000.vtu', folder + '/surface-traction.vtu')
    #         shutil.copyfile(working_directory + 'results/Extracted/wholegeometry-velocity_2000.vtu', folder + '/wholegeometry-velocity.vtu')
    #         shutil.copyfile(working_directory + 'results/Extracted/surface-tangentialprojectiontraction_2000.vtu', folder + '/surface-tangentialprojectiontraction.vtu')
  
