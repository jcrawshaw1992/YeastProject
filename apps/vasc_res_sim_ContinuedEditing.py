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
from FileConverter import vtuTostl
import FSI_Fun

Iterations=30
# SimulationsList = ['/Users/jcrawshaw/Documents/ChasteWorkingDirectory/IncreasedFluidPressure/' ,'/Users/jcrawshaw/Documents/ChasteWorkingDirectory/IncreasedFluidPressure/' ]
chaste_setup_exe = '/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/build/optimised/HemeLBChasteTests/TestSetupFlowInVesselRunner'
chaste_run_exe = '/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/build/optimised/HemeLBChasteTests/TestRunFlowInVesselRunner'
hemelb_setup_exe = 'env PYTHONPATH=/Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool:$PYTHONPATH /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/setuptool/scripts/hemelb-setup-nogui'
radii_over_time = [] # empty array


if __name__=="__main__":

    working_directory = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory/ShrunkPlexusWithLongInlets/' 
    FSI_Fun.SetWorkingDirectory(working_directory)  
    data_path =working_directory + 'SetUpData/'

    files = [data_path +'config.stl', data_path +'config.xml', data_path +'config.vtu', data_path +'config.pr2']
    for f in files:
        shutil.copy(f, working_directory)

# *** Set up Chaste Folder to collect the deformation data
    ChastePath = working_directory +'ChasteMeshes/' 
    if not os.path.exists(ChastePath):
        os.makedirs(ChastePath)
    else:
        shutil.rmtree(ChastePath)
        os.makedirs(ChastePath)


    HemeLBPath = working_directory +'HemeLBFluid/' 
    if not os.path.exists(HemeLBPath):
        os.makedirs(HemeLBPath)
    else:
        shutil.rmtree(HemeLBPath)
        os.makedirs(HemeLBPath)
 
    duration = 0.25
    dt = 0.02 # changed from 0.0001
    SamplingTimestepMultiple = 5 # changed from 10
    # EndMeshFileNumber = duration /dt
    list = range(0, int((duration/dt))+1, SamplingTimestepMultiple)
    ResultsCounter =  0
  
    # Step 3: Run preliminary Chaste setup
    mesh_scale = 1e-3
    os.environ['CHASTE_TEST_OUTPUT'] = working_directory  
    chaste_setup_call = chaste_setup_exe +' -mesh ' + working_directory + 'config.vtu -xml ' + working_directory +'config.xml -mesh_scale ' +  str(mesh_scale) + ' -dt ' + str(dt)  +' -SamplingTimestepMultiple ' + str(SamplingTimestepMultiple)
    subprocess.call(chaste_setup_call, shell=True)
    # print ' ------- Chaste is set up ------- '
   
    
    # vtuTostl(working_directory + 'config.vtu', working_directory + 'config.stl')

    start_time = 0.0
    for iter in range(Iterations):
        print 'Iteration ' + str(iter)

        start_time = round(start_time, 4)

        # Step 2: Run HemeLB setup
        FSI_Fun.run_hemelb_setup()

        # Step 2a: Update xml
        FSI_Fun.update_xml_file(700)

        # Step 3: HemeLB simulation
        FSI_Fun.run_hemelb() 
        # pdb.set_trace()
        # Step 4: Run Chaste with the tractions previously computed
        command = chaste_run_exe + ' -start_time ' + str(start_time) + ' -duration ' + str(duration) + ' -traction_file ' +  working_directory + 'results/Extracted/surface-tractions.xtr' + ' -mesh_scale ' +  str(mesh_scale)
        subprocess.call(command, shell=True)
        
        # Update the radius of the tube*, this can then be used to update the pr2 file to enable the flow to be correctly calculated 
        # Step 5: Convert Chaste output (vtu file) into the input of the next iteration (stl file required by hemelb setup tool)

        vtuTostl(working_directory + 'config.vtu', working_directory + 'config.stl')
        FSI_Fun.generate_flow_vtus(working_directory)

        # Step 5: Collect Chaste Mesh outputs 
        if start_time == int(start_time):
            start_time = int(start_time)
        
        OldFilePath = working_directory +'results_from_time_' + str(start_time) + '/'
        for j in list:
            shutil.copy(OldFilePath  + 'mesh_'+ str(j) +'.vtu' ,ChastePath + 'mesh_'+ str(ResultsCounter) +'.vtu')
            shutil.copy(OldFilePath  + 'results_'+ str(j) +'.vtu' ,ChastePath + 'results_'+ str(ResultsCounter) +'.vtu')
            ResultsCounter+=1
        shutil.copy(OldFilePath +'results.viznodes' , ChastePath+'results.viznodes')
        shutil.copy(OldFilePath +'results.vizmutationstates' , ChastePath+'MutationStates')
        # shutil.rmtree(OldFilePath)  

        Numb_FluidFiles =int(math.floor((duration/dt)/SamplingTimestepMultiple)+1)

        # Collect the fluid vtus
        directory = working_directory +'results/Extracted/'
        for file in os.listdir(directory):
            if file.startswith("surface-pressure_") and file.endswith(".vtu"):
                for j in range(1, Numb_FluidFiles+1 , 1):
                    FluidOutputNumber =  j  +3*(iter)
                    shutil.copy(directory + file , HemeLBPath + 'surface-pressure_'+  str(FluidOutputNumber) +'.vtu')

            if file.startswith("wholegeometry-velocity_") and file.endswith(".vtu"):
                for j in range(1, Numb_FluidFiles+1 , 1):
                    FluidOutputNumber =  j  +3*(iter)
                    shutil.copy(directory + file , HemeLBPath + 'wholegeometry-velocity_'+  str(FluidOutputNumber) +'.vtu')

        # Step 5a: compute radii
        FSI_Fun.vmtk_compute_stl_radii(radii_over_time, iter)
        Radius_New = radii_over_time[iter][-1]  # LengthOfArray = len(radii_over_time[iter]) 
        FSI_Fun.update_pr2_file(Radius_New) 
        start_time += duration     
      


    
  

    
    csvfile = open(working_directory + 'radii_over_time.csv', 'a')
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




        # #  Get the two fluid outputs you need  
        # directory = working_directory +'results/Extracted/'
        # # Move the fluids files into one place for ordered and safe keeping -- No, have the one fluids save for however many times you need it 
        # # FluidOutputInterval = 250
        # for file in os.listdir(directory):
        #     if file.endswith(".vtu"):
        #         if file.startswith("surface-pressure_"):
        #             Index1 = file.find('_',1,len(file)) +1
        #             Index2 = file.find('.',1,len(file))
        #             FluidOutputNumber =  int(file[Index1:Index2])/FluidOutputInterval  +3*(iter)

        #             PressureFile = directory + file 
        #             shutil.copy(PressureFile, HemeLBPath + 'surface-pressure_'+  str(FluidOutputNumber) +'.vtu')
        #         if file.startswith("wholegeometry-velocity_"):

        #             Index1 = file.find('_',1,len(file)) +1
        #             Index2 = file.find('.',1,len(file))
        #             FluidOutputNumber =  int(file[Index1:Index2])/FluidOutputInterval  +3*(iter)

        #             VelocityFile = directory + file  
        #             shutil.copy(VelocityFile , HemeLBPath  +'wholegeometry-velocity_'+ str(FluidOutputNumber)  +'.vtu')
