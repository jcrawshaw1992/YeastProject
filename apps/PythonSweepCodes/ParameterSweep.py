#!/usr/bin/env python
# vasc_res_sim.py
import subprocess
import shutil
import os
import glob
from argparse import ArgumentParser
import numpy as np
import time
import pdb
import string
import math
import sys

Iterations=30
# SimulationsList = ['/Users/jcrawshaw/Documents/ChasteWorkingDirectory/IncreasedFluidPressure/' ,'/Users/jcrawshaw/Documents/ChasteWorkingDirectory/IncreasedFluidPressure/' ]
chaste_run_exe = '/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/build/optimised/HemeLBChasteTests/TestRunFlowInVesselRunner'



if __name__=="__main__":
    GenerateRunner =1
    if GenerateRunner ==1:
        command = " /home/vascrem/Chaste scons b=GccOpt projects/VascularRemodelling/test/ParameterSweep/TestMembraneParameters.hpp"
#         subprocess.call(command, shell=True)
        # print"Need to generate runner"

#     CompletedAreaParameter = [6, 6.5, 7,7.5, 8]
#     CompletedDilationParameter =[6, 6.5, 7,7.5, 8]
#     CompletedDeformationParamter = [6, 6.5, 7,7.5, 8]

#     AreaParameter = [4, 4.5,5, 5.5,6, 6.5, 7,7.5, 8, 8,8.5,9]
#     DilationParameter =[4, 4.5,5,5.5, 6, 6.5, 7,7.5, 8,8.5,9]
#     DeformationParamter = [6, 6.5,5.5, 7,7.5, 8, 8,8.5,9]

#     # for i in AreaParameter:
#     #     for j in DilationParameter:
#     #         for k in DeformationParamter:
#     #             if ((i in CompletedAreaParameter) & (j in CompletedDilationParameter) & (k in CompletedDeformationParamter)):
#     #                 print "skip"
#     #             else:
#                       command = chaste_run_exe + ' -start_time ' + str(start_time) + ' -duration ' + str(duration) + ' -traction_file ' +  working_directory + 'results/Extracted/surface-tractions.xtr' + ' -mesh_scale ' +  str(mesh_scale)
# #         subprocess.call(command, shell=True)
#                      subprocess.Popen([sys.executable, "/home/vascrem/Chaste/projects/VascularRemodelling/apps/CollectSweepResults.py"])
#                      print"somethingelse"
#     #                 print" open new screen with naming convention"
#     #                 print(" Set Chaste to run in this screen ")
#     #                 print("Leave the screen  ")

#     print '\n ********* ------ Completed ------ ********* \n'







#     working_directory = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory/ShrunkPlexusWithLongInlets/' 
#     FSI_Fun.SetWorkingDirectory(working_directory)  
#     data_path =working_directory + 'SetUpData/'

#     files = [data_path +'config.stl', data_path +'config.xml', data_path +'config.vtu', data_path +'config.pr2']
#     for f in files:
#         shutil.copy(f, working_directory)

# # *** Set up Chaste Folder to collect the deformation data
#     ChastePath = working_directory +'ChasteMeshes/' 
#     if not os.path.exists(ChastePath):
#         os.makedirs(ChastePath)
#     else:
#         shutil.rmtree(ChastePath)
#         os.makedirs(ChastePath)

#     duration = 0.2
#     dt = 0.002 # changed from 0.0001
#     SamplingTimestepMultiple = 10 # changed from 10
#     # EndMeshFileNumber = duration /dt
#     list = range(0, int((duration/dt))+1, SamplingTimestepMultiple)
#     ResultsCounter =  0
 
 #         # Step 4: Run Chaste with the tractions previously computed
#         command = chaste_run_exe + ' -start_time ' + str(start_time) + ' -duration ' + str(duration) + ' -traction_file ' +  working_directory + 'results/Extracted/surface-tractions.xtr' + ' -mesh_scale ' +  str(mesh_scale)
#         subprocess.call(command, shell=True)
    

