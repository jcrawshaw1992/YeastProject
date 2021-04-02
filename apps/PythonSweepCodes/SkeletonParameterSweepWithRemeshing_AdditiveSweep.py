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
import os.path
from os import path

def RunSimulations(i,j, chaste_run_exe, TerminalOutputFolder):

    DeformationParamter = [  6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]

    dt = ' -dt 50'
    EndTime = 50
    SamplingTimestepMultiple = ' -SamplingTimestepMultiple 500'

    K =5.5
    # Input1 = chaste_run_exe +' -AreaParameter '+str(i)+' -DilationParameter '+str(j)+' -DeformationParamter ' +str(K)+   dt + SamplingTimestepMultiple +' -EndTime ' + str(EndTime) 
    # Input2 = TerminalOutputFolder+'AreaParameter_'+str(i)+'_DilationParameter'+str(j)+'_DeformationParamter_' +str(K)+'.txt'
    # Input3 = TerminalOutputFolder+'WaitFile.txt'
    # subprocess.call(['./RunChaste', Input1,Input2,Input3] )

    ParameterSet = "Param_" + str(i) + "_DilationParam_" +  str(j) +  "_DeformationParam_" +  str(K)
    ArchieveFile = "ParameterSweepWithRemeshing/Cylinder/" + ParameterSet
    StartTime = EndTime
    # run the first parameter sweep, then all the others are just addition on to the first :) 
    for k in DeformationParamter:
        AdjustEndTime = 20
        Input1 = chaste_run_exe +' -AreaParameter '+str(i)+' -DilationParameter '+str(j)+' -DeformationParamter ' +str(k)+ ' -ArchieveFile '+ArchieveFile + ' -StartTime '+str(StartTime) +' -EndTime ' + str(AdjustEndTime) 
        Input2 = TerminalOutputFolder+'AreaParameter_'+str(i)+'_DilationParameter_'+str(j)+'_DeformationParamter_' +str(k)+'.txt'
        Input3 = TerminalOutputFolder+'WaitFile.txt'
        subprocess.call(['./RunChaste', Input1,Input2,Input3] )
        # subprocess.Popen(['./RunChaste', Input1,Input2,Input3 ])
        StartTime = StartTime + AdjustEndTime
        ParameterSet = "Param_" + str(i) + "_DilationParam_" +  str(j) +  "_DeformationParam_" +  str(k)
        ArchieveFile = "ParameterSweepWithRemeshing/Cylinder/" + ParameterSet

if __name__=="__main__":

    # Define arguments to be parsed
    parser = ArgumentParser(description='Run a vascular remodelling simulation')
    parser.add_argument('-i', dest='I',type=float, default='-10', help='Need Area parameter ')
    parser.add_argument('-j', dest='J',type=float, default='-10', help='Need Dilation parameter  ')
    parser.add_argument('-chaste_run_exe', dest='Chaste_run_exe',type=str, default='/home/vascrem/Chaste/projects/VascularRemodelling/build/optimised/ParameterSweep/TestMembraneParametersWithRemeshingRunner ', help='Set the archieve files for chaste to continue from ')
    parser.add_argument('-TerminalOutputFolder', dest='TerminalOutput',type=str, default="/data/vascrem/testoutput/ParameterSweep/Cylinder/SweepTerminalOutputs/", help='Set the archieve files for chaste to continue from ')
    # Parse arguments (this will create args.flow_to_vtu etc. variables)
    args = parser.parse_args()


    chaste_run_exe = args.Chaste_run_exe
    TerminalOutputFolder = args.TerminalOutput

    t0 = time.time()
    GenerateRunner =1
    if GenerateRunner ==1:
        command = "cd ~/Chaste && scons b=GccOpt projects/VascularRemodelling/test/ParameterSweep/TestMembraneParametersWithRemeshing.hpp"
        subprocess.call(command, shell=True)
    Server = 1
    if Server ==1:
       chaste_run_exe = '/home/vascrem/Chaste/projects/VascularRemodelling/build/optimised/ParameterSweep/TestMembraneParametersWithRemeshingRunner '
        TerminalOutputFolder = "/data/vascrem/testoutput/ParameterSweepWithRemeshing/Cylinder/SweepTerminalOutputs/"
    else:
        chaste_run_exe =  '/Users/jcrawshaw/Chaste/projects/VascularRemodelling/build/optimised/ParameterSweep/TestMembraneParametersWithRemeshingRunner '
        TerminalOutputFolder = "/Users/jcrawshaw/Documents/testoutput/ParameterSweepWithRemeshing/Cylinder/SweepTerminalOutputs/"

    if args.I<0: # Set up the 30 simulation here
        t0 = time.time()
        
        DilationParameter = [ 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
        AreaParameter = [6, 8, 10]
        for j in DilationParameter:
            for i in AreaParameter:
                p = os.popen('python SkeletonParameterSweepWithRemeshing_AdditiveSweep.py -i '+str(i) +' -j '+str(j) ) # this takes to line 89 all at once; 30 at a time
                print(p.read())
    else: # And here is where everything is actually ran in parallel
        RunSimulations(args.I,args.J,chaste_run_exe, TerminalOutputFolder)


    

    
    