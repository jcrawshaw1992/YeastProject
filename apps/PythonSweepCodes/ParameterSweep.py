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
        chaste_run_exe = '/home/vascrem/Chaste/projects/VascularRemodelling/build/optimised/ParameterSweep/TestMembraneParametersRunner '
        TerminalOutputFolder = "/data/vascrem/testoutput/ParameterSweep/Cylinder/SweepTerminalOutputs/"
    else:
        chaste_run_exe =  '/Users/jcrawshaw/Chaste/projects/VascularRemodelling/build/optimised/ParameterSweep/TestMembraneParametersRunner '
        TerminalOutputFolder = "/Users/jcrawshaw/Documents/testoutput/ParameterSweep/Cylinder/SweepTerminalOutputs/"

    if path.isdir(TerminalOutputFolder)==0:
        os.mkdir(TerminalOutputFolder)

    subprocess.call("chmod 700 RunChaste", shell=True)

    Completedarameters = [5.5, 6, 6.5, 7,7.5, 8]
    AreaParameter = [5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
    DilationParameter = [5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
    DeformationParamter = [5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]

    dt = ' -dt 0.0001'
    EndTime = ' -NewEndTime 15' 
    SamplingTimestepMultiple = ' -SamplingTimestepMultiple 10000'
    Parallel = 16
    SleepyTime =50
    AvaliablePaths = range(Parallel)
    print AvaliablePaths
    for i in AreaParameter:
        for j in DilationParameter:
            for k in DeformationParamter:
                if ((i in Completedarameters) & (j in Completedarameters) & (k in Completedarameters)):
                    continue
                Core = AvaliablePaths[0]
                if (i<6 or j<6 or k<6) :
                    Input1 = chaste_run_exe +' -AreaParameter '+str(i)+' -DilationParameter '+str(j)+' -DeformationParamter ' +str(k)+   dt + SamplingTimestepMultiple +EndTime
                    Input2 = TerminalOutputFolder+'AreaParameter_'+str(i)+'_DilationParameter+'+str(j)+'_DeformationParamter_' +str(k)+'.txt'
                    Input3 = TerminalOutputFolder+'WaitFile'+str(Core)+'.txt'
                    subprocess.Popen(['./RunChaste', Input1,Input2,Input3 ])
                else:
                    Input1 = chaste_run_exe +' -AreaParameter '+str(i)+' -DilationParameter '+str(j)+' -DeformationParamter ' +str(k)  + EndTime
                    Input2 = TerminalOutputFolder+'AreaParameter'+str(i)+'_DilationParameter'+str(j)+'_DeformationParamter' +str(k)+'.txt'
                    Input3 = TerminalOutputFolder+'WaitFile'+str(Core)+'.txt'
                    subprocess.Popen(['./RunChaste', Input1,Input2,Input3 ])
                AvaliablePaths.remove(Core) 
                # # Check if all positions are taken
                while len(AvaliablePaths) ==0:
                    time.sleep(SleepyTime)
                    # print "Awake and checking for spare cores" 
                    for P in range(Parallel):
                        OutputFile = TerminalOutputFolder+'WaitFile'+str(P)+'.txt'
                        if path.exists(OutputFile):
                            AvaliablePaths.append(P)
                            os.remove(OutputFile)
                    # if len(AvaliablePaths) >0:
                        # print AvaliablePaths
                        # print "Have found a spare core or two :-) " 


    print '\n ********* ------ Completed ------ ********* \n'