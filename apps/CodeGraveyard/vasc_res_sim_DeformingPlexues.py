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

from vasc_res_sim import *

Iterations=100
working_directory = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory/BifurcationSecond/' 
chaste_setup_exe = '/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/build/optimised/TestSetupFlowInPlexusRunner'

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


duration =  0.001
dt = 0.00001 # changed from 0.0001
SamplingTimestepMultiple = 5*20
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

ElasticShearModulus = ' -ElasticShearModulus ' +    str(4.4e-9) #+ +   str(4.4e-07) #  +    #  1e-4 Dont change this, stops werid stuff happening 
AreaDilationModulus = ' -AreaDilationModulus '+   str(0.9e-11) # +   str(0.9e-11) 
membrane_constant = ' -membrane_constant ' +    str(5) #+    str( 100 ) # str(3 * 1e-14)
Area_constant = ' -Area_constant ' +    str(0.9e-15) #+    str(0.9e-15) # str(8e-11) # str(8 * 1e-17) # str(5 * 1e-6)
Time_step = ' -dt ' + str(dt)
Sampling = ' -SamplingTimestepMultiple ' + str(SamplingTimestepMultiple)
duration_of_iteration = " -Duration " + str(duration)

Scalling =1
ScallingFactor = " -Scalling " + str(Scalling)

WriteReadme(working_directory, ElasticShearModulus, AreaDilationModulus, membrane_constant, Area_constant , Time_step, Sampling, duration_of_iteration , ScallingFactor)

chaste_setup_call = chaste_setup_exe +' -mesh ' + mesh_filename + ' -xml ' + xml_filename + ' -div_threshold ' +  str(args.div_threshold) + ' -mesh_scale ' +  str(args.mesh_scale) + ' -wss_threshold ' + str(0.025)  + ElasticShearModulus + AreaDilationModulus + membrane_constant + Area_constant + Time_step  +Sampling +ScallingFactor
subprocess.call(chaste_setup_call, shell=True)
# pause()
print ' ------- Chaste is set up ------- '


# Step 3.a: Prepare mesh for HemeLB setup tool and compute radii along the axes of the domain (if requested)

vtu2stl(working_directory,0)

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

        run_hemelb_setup(working_directory, 1)

        print "Update xml file"

        # Step 2a: Update xml

        update_xml_file(working_directory, iter, args.num_iterations)
    
        print "Run HemeLB"
        # Step 3: HemeLB simulation
        run_hemelb2(working_directory) 
        print'\n HemeLB simulation complete \n'
    

    # # Step 4: Run Chaste with the tractions previously computed
    traction_filename = working_directory + 'results/Extracted/surface-tractions.xtr'
    # traction_filename = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/DeformingCylindricalBasicHetroWallTesting/results/Extracted/surface-tractions.xtr"
    
    command = chaste_run_exe + ' -start_time ' + str(start_time) + ' -duration ' + str(duration) + ' -traction_file ' +  traction_filename + ' -mesh_scale ' +  str(args.mesh_scale)
    print 'Chaste step'
    subprocess.call(command, shell=True)

    
    # Step 6: Collect Chaste Mesh outputs 

    if start_time == int(start_time):
        start_time = int(start_time)
        print start_time
        
    # pdb.set_trace()
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
    # FluidSimulation = 1
    if FluidSimulation ==1:
        vtu2stl(working_directory, iter)

    # Step 5a: compute radii
        # vmtk_compute_stl_radii(working_directory , iter)
        # Radius_New = radii_over_time[iter][-1]  # LengthOfArray = len(radii_over_time[iter])    
    # pdb.set_trace()
    # print Radius_New
        # update_pr2_file(working_directory, Radius_New) 
        generate_flow_vtus(working_directory, 3)

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


    # pause()

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