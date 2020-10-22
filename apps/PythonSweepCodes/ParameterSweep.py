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

if __name__=="__main__":
    GenerateRunner =1
    if GenerateRunner ==1:
        command = "cd ~/Chaste && scons b=GccOpt projects/VascularRemodelling/test/ParameterSweep/TestMembraneParameters.hpp"
        subprocess.call(command, shell=True)
    Server = 1
    if Server ==1:
        chaste_run_exe = '~/Chaste/projects/VascularRemodelling/build/optimised/ParameterSweep/TestMembraneParametersRunner'
        TerminalOutputFolder = "/data/vascrem/testoutput/ParameterSweep/Cylinder/SweepTerminalOutputs/"
    else:
        chaste_run_exe =  '/Users/jcrawshaw/Chaste/projects/VascularRemodelling/build/optimised/ParameterSweep/TestMembraneParametersRunner '
        TerminalOutputFolder = "/Users/jcrawshaw/Documents/testoutput/ParameterSweep/Cylinder/SweepTerminalOutputs/"

    
    # TerminalOutputFolder = "/data/vascrem/testoutput/ParameterSweep/Cylinder/SweepTerminalOutputs/"
    if path.isdir(TerminalOutputFolder)==0:
        os.mkdir(TerminalOutputFolder)

    # Set up the bash scrip that runs Chaste
    command = "chmod 700 RunChaste"
    subprocess.call(command, shell=True)

    CompletedAreaParameter = [6, 6.5, 7,7.5, 8]
    CompletedDilationParameter =[6, 6.5, 7,7.5, 8]
    CompletedDeformationParamter = [6, 6.5, 7,7.5, 8]

    AreaParameter = [5, 5.5,6, 6.5, 7,7.5, 8, 8,8.5,9, 8.5,10]
    DilationParameter =[5,5.5, 6, 6.5, 7,7.5, 8,8.5,9,8.5, 10]
    DeformationParamter = [5,5.5, 6, 6.5, 7,7.5, 8, 8,8.5,9,10]

    AreaParameter = [5]
    DilationParameter =[5]
    DeformationParamter = [6]

    dt = ' -dt 0.0001'
    SamplingTimestepMultiple = ' -SamplingTimestepMultiple 10000'

    for i in AreaParameter:
        for j in DilationParameter:
            for k in DeformationParamter:
                if ((i in CompletedAreaParameter) & (j in CompletedDilationParameter) & (k in CompletedDeformationParamter)):
                    print "skip"
                else:
                    print i 
                    print j
                    print k 
                    if (i<6 or j<6 or k<6) :
                        Input1 = chaste_run_exe +' -AreaParameter '+str(i)+' -DilationParameter '+str(j)+' -DeformationParamter ' +str(k)+   dt + SamplingTimestepMultiple
                        Input2 = TerminalOutputFolder+'AreaParameter_'+str(i)+'_DilationParameter+'+str(j)+'_DeformationParamter_' +str(k)+'.txt'
                        subprocess.Popen(['./RunChaste', Input1,Input2])
                    else:
                        Input1 = chaste_run_exe +' -AreaParameter '+str(i)+' -DilationParameter '+str(j)+' -DeformationParamter ' +str(k)  
                        Input2 = TerminalOutputFolder+'AreaParameter_'+str(i)+'_DilationParameter'+str(j)+'_DeformationParamter_' +str(k)+'.txt'
                        subprocess.Popen(['./RunChaste', Input1,Input2])

    print '\n ********* ------ Completed ------ ********* \n'



    # #     ArchivedDirectory = ' -ArchivedDirectory /data/vascrem/testoutput/ParameterSweep/Cylinder/Parameteres2/Param_6_DilationParam_6_DeformationParam_6/'
    # # StartingParameterForSlowIncrease = " -StartingParameterForSlowIncrease 1e-6"

    #     NewEndTime = ' -NewEndTime 45'
    # EndTime = ' -EndTime 30'
