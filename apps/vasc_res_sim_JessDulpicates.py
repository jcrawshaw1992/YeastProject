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
def vmtk_compute_stl_radii(iter_number): 
    print "**************   vmtk_compute_stl_radii  **************  "
    
    CenterLines_filename = working_directory + 'centerlines' +  str(iter_number) + '.vtp'
    NewOutput = working_directory + 'NewFile' +  str(iter_number) + '.vtp'

    # Get the center line file 
    command = 'vmtk vmtknetworkextraction -ifile ' + working_directory + 'config.stl' + ' -ofile ' + CenterLines_filename 
    subprocess.call(command, shell=True)

    # pdb.set_trace()
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(CenterLines_filename)
    reader.Update()
    point_data = reader.GetOutput().GetPointData()
    assert point_data.GetArrayName(0) == 'Radius', "VTP file doesn't contain array Radius"
    radii = point_data.GetArray(0)
    radii_over_time.append(np.array([radii.GetValue(i) for i in range(radii.GetSize())]))
    print "    vmtk_compute_stl_radii DONE   "

    

def generate_flow_vtus(timestep):
    print "  generate_flow_vtus   "

    # Make the vtu 
    shutil.copyfile(working_directory + 'config.xml', working_directory + 'config2.xml')
    command = 'python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/GmyUnstructuredGridReader.py' + ' ' + working_directory + 'config2.xml'
    subprocess.call(command, shell=True)

    # command = 'python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py' + ' ' +  working_directory +'config2.vtu' +' ' +  working_directory+ 'results/Extracted/surface-pressure.xtr' +' ' +  working_directory+ 'results/Extracted/wholegeometry-velocity.xtr'+' ' +  working_directory+ 'results/Extracted/surface-traction.xtr'+' ' +  working_directory+ 'results/Extracted/surface-tangentialprojectiontraction.xtr'
    command = 'python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py' + ' ' +  working_directory +'config2.vtu' +' ' +  working_directory+ 'results/Extracted/surface-pressure.xtr' +' ' +  working_directory+ 'results/Extracted/wholegeometry-velocity.xtr'+' ' +  working_directory+ 'results/Extracted/surface-traction.xtr'
   
    subprocess.call(command, shell=True)

    print "generate_flow_vtus DONE"


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

    # surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(5000), 'file': 'surface-tangentialprojectiontraction.xtr'})
    # ElementTree.SubElement(surface, 'geometry', type='surface')
    # ElementTree.SubElement(surface, 'field', type='tangentialprojectiontraction')

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(900), 'file': 'surface-traction.xtr'}) # ElementTree.SubElement(extr, 'propertyoutput', {'period': str(10000), 'file': 'surface-traction.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='traction')

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(900), 'file': 'surface-tractions.xtr'}) # surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(10000), 'file': 'surface-tractions.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='traction')
    ElementTree.SubElement(surface, 'field', type='tangentialprojectiontraction')    

    surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(900), 'file': 'surface-pressure.xtr'}) # surface = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(10000), 'file': 'surface-pressure.xtr'})
    ElementTree.SubElement(surface, 'geometry', type='surface')
    ElementTree.SubElement(surface, 'field', type='pressure')

    wholegeometry = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(90), 'file': 'wholegeometry-velocity.xtr'}) #wholegeometry = ElementTree.SubElement(extr, 'propertyoutput', {'period': str(10000), 'file': 'wholegeometry-velocity.xtr'})
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

def WriteReadme(ElasticShearModulus, AreaDilationModulus, membrane_constant, Area_constant , Time_step, Sampling, duration, Scalling ):
    print "**************   WriteReadme  **************   "

    file1 = open(working_directory +"readme.txt","w") 
    file1.write(" --- NOTE - PIPE NOT INITIALLY COLLAPSED IN THIS SIMULATION, NEED TO RENAME FOLDER --- \n\n")
    file1.write( Scalling +"\n" )
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


if __name__=="__main__":
    working_directory = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory/CollapsedPlexues/' 
    data_path =working_directory + 'SetUpData/'
    
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


    duration =  0.01# changed from 0.001
    dt = 0.001 # changed from 0.0001
    SamplingTimestepMultiple = 5
    TimeStep = dt * SamplingTimestepMultiple
    EndMeshFileNumber = duration /dt
    list = range(0, int(EndMeshFileNumber), SamplingTimestepMultiple)
 
    ResultsCounter=0
    MeshCounter=0
    PressureCounter =0
    

    # Define arguments to be parsed
    parser = ArgumentParser(description='Run a vascular remodelling simulation')
    parser.add_argument('--num_iterations', default=Iterations, type=int, help='Number of Hemelb/Chaste iterations to be run (optional, default is 5).')
    parser.add_argument('--compute_radii', dest='compute_radii', action='store_true', help='Use VMTK to compute axis radii.')
    parser.add_argument('--output_postfix', dest='output_postfix', default='', help='This string will be added to ChasteWorkingDirectory to get the output folder.')
    parser.add_argument('--div_threshold', dest='div_threshold', default=1e10, help='This specifies the length that edges will divide at. (Defaults to no division, but 6e-4 is good for mm meshes')
    parser.add_argument('--mesh_scale', dest='mesh_scale', default=1e-3, help='This specifies what to scale the mesh by so that all distances are in meters (defaults to mm).')
    args = parser.parse_args()
    
    # Parse arguments (this will create args.flow_to_vtu etc. variables)   
    os.environ['CHASTE_TEST_OUTPUT'] = working_directory  
    #
    # Step 3: Run preliminary Chaste setup
    #
    mesh_filename = working_directory + 'config.vtu' # Use proper path concatentation
    xml_filename = working_directory +'config.xml' # Use proper path concatentation



    print ' ------- Setting up Chaste simulation ------- '

    ElasticShearModulus = ' -ElasticShearModulus '+   str(5.0119e-05) #  +   str(4.4e-07) # str(5e-8) #str(6e-8) #  1e-4 Dont change this, stops werid stuff happening 
    AreaDilationModulus = ' -AreaDilationModulus ' + str(1e-05) # str(0.9e-11) # str(5e-9) this one controls the dilation of the cylinder, so can be played with 
    membrane_constant = ' -membrane_constant ' + str( 0) # str(100)
    Area_constant = ' -Area_constant ' +  str(1e-05) # str(0.9e-15) # str(8 * 1e-17) # str(5 * 1e-6)
    Time_step = ' -dt ' + str(dt)
    Sampling = ' -SamplingTimestepMultiple ' + str(SamplingTimestepMultiple)
    duration_of_iteration = " -Duration " + str(duration)

    Scalling =1 
    ScallingFactor = " -Scalling " + str(Scalling)

    WriteReadme(ElasticShearModulus, AreaDilationModulus, membrane_constant, Area_constant , Time_step, Sampling, duration_of_iteration , ScallingFactor)

    chaste_setup_call = chaste_setup_exe +' -mesh ' + mesh_filename + ' -xml ' + xml_filename + ' -div_threshold ' +  str(args.div_threshold) + ' -mesh_scale ' +  str(args.mesh_scale) + ' -wss_threshold ' + str(0.025)  + ElasticShearModulus + AreaDilationModulus + membrane_constant + Area_constant + Time_step  +Sampling +ScallingFactor
    subprocess.call(chaste_setup_call, shell=True)
    print ' ------- Chaste is set up ------- '
   
    
   # Step 3.a: Prepare mesh for HemeLB setup tool and compute radii along the axes of the domain (if requested)
    
    vtu2stl(0)

    # 

    start_time = 0
    FluidSimulation = 1
    # duration = 0.3
    for iter in range(args.num_iterations):
        print 'Iteration ' + str(iter)
        if FluidSimulation ==1:
            print "Fluid Pressure with membrane Forces"
            
            start_time = round(start_time, 4)
            print start_time
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
            generate_flow_vtus(3)
  
  
        # # Step 4: Run Chaste with the tractions previously computed
        traction_filename = working_directory + 'results/Extracted/surface-tractions.xtr'
        # traction_filename = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/DeformingCylindricalBasicHetroWallTesting/results/Extracted/surface-tractions.xtr"
        
        command = chaste_run_exe + ' -start_time ' + str(start_time) + ' -duration ' + str(duration) + ' -traction_file ' +  traction_filename + ' -mesh_scale ' +  str(args.mesh_scale)
        print 'Chaste step'
        subprocess.call(command, shell=True)

        
        # # Step 6: Collect Chaste Mesh outputs 

        if start_time == int(start_time):
            start_time = int(start_time)
            print start_time
          
        pdb.set_trace()
        OldFilePath = working_directory +'results_from_time_' + str(start_time) + '/'
        shutil.move(OldFilePath +'results.viznodes' , newpath+'results.viznodes')

        for j in list:
            shutil.move(OldFilePath  + 'mesh_'+ str(j) +'.vtu' ,newpath + 'mesh_'+ str(MeshCounter) +'.vtu')
            MeshCounter+=1
            shutil.move(OldFilePath  + 'results_'+ str(j) +'.vtu' ,newpath + 'results_'+ str(ResultsCounter) +'.vtu')
            ResultsCounter+=1
        start_time += duration 

        # Update the radius of the tube*, this can then be used to update the pr2 file to enable the flow to be correctly calculated 
        # Step 5: Convert Chaste output (vtu file) into the input of the next iteration (stl file required by hemelb setup tool)
        if FluidSimulation ==1:
            # vtu2stl(iter)

        # Step 5a: compute radii
            # vmtk_compute_stl_radii(iter)
            # Radius_New = radii_over_time[iter][-1]  # LengthOfArray = len(radii_over_time[iter])    
        # pdb.set_trace()
        # print Radius_New
            # update_pr2_file(Radius_New) 
            generate_flow_vtus(3)
            pause()
        # shutil.rmtree(OldFilePath)   
        #  Get the two fluid outputs you need  
            directory = working_directory +'results/Extracted/'
            print "Jess is great"
            for file in os.listdir(directory):
                if file.endswith(".vtu"):
                    if file.startswith("surface-pressure_"):
                        PressureFile = directory + file 
                    if file.startswith("wholegeometry-velocity_"):
                        VelocityFile = directory + file  
            # Make sure the number of fluid outputs is equal to the number of solid outputs 
            for j in list:
                shutil.copy(VelocityFile , Hemepath  +'wholegeometry-velocity_'+  str(PressureCounter) +'.vtu')
                shutil.copy(PressureFile  , Hemepath + 'surface-pressure_'+  str(PressureCounter) +'.vtu')
                PressureCounter+=1

  

    
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
    

    file1 = open(working_directory +"readme.txt","a") 
    file1.write("\n\n ---- Sucessful ----")
    file1.close() 

    print '\n ********* ------ Completed ------ ********* \n'

  



    # --------------

    # Need to go through this stuff and take what I need
     # there is problems, its asking for an input directory 

    #  CenterLines_filename = working_directory + 'centerlines' +  str(2) + '.vtp'
    # command = 'vmtkcenterlines -ifile '  + working_directory + 'config.stl' + ' -ofile ' + CenterLines_filename
    # subprocess.call(command, shell=True)
    #     # vmtkcenterlinestonumpy


    # ---------------
    #     command =  'vmtk vmtkcenterlinestonumpy -ifile ' + CenterLines_filename + ' -celltopoint ' + str(1) + ' -ofile ' + NewOutput 
    # subprocess.call(command, shell=True)
    #     pdb.set_trace()
        

    #     # Find the points along the center line
    #     CenterlineAttributes_FileName =  working_directory + 'CenterLineAttributes' +  str(iter_number) + '.vtp'
    #     command = 'vmtkcenterlineattributes -ifile ' + CenterLines_filename + ' -ofile '+ CenterlineAttributes_FileName
    #     subprocess.call(command, shell=True)
    #     pdb.set_trace()


    # NewOutput = working_directory + 'NewFile' +  str(iter_number) + '.vtp'
    # #  vmtkcenterlines -ifile foo.vtp -ofile foo_centerlines.vtp
    # command = 'vmtkcenterlines -ifile ' + CenterLines_filename + ' -ofile ' +NewOutput 
    # subprocess.call(command, shell=True)
    # pdb.set_trace()

    # ---------------------------------