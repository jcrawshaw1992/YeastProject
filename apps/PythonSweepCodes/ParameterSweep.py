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
    t0 = time.time()
    GenerateRunner =0
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

    AreaParameter = [5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
    DilationParameter = [ 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
    DeformationParamter = [ 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
    Completed = [ [5.5, 5.5]  , [5.5, 6], [5.5, 6.5, 5.5] , [5.5, 6.5, 6], [5.5, 6.5, 6.5], [5.5, 6.5, 7] ] 


    AreaParameter = [8, 8.5, 9, 9.5, 10]
    DilationParameter = [ 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
    DeformationParamter = [ 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]


    
    dt = ' -dt 0.0005'
    EndTime = ' -NewEndTime 15' 
    SamplingTimestepMultiple = ' -SamplingTimestepMultiple 2000'
    Parallel = 14
    SleepyTime = 120
    AvaliablePaths = range(Parallel)
    print AvaliablePaths
    for i in AreaParameter:
        for j in DilationParameter:
            for k in DeformationParamter:
                Core = AvaliablePaths[0]
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
                    if len(AvaliablePaths) >0:
                        print AvaliablePaths, "Have found a spare core or two :-) " 
                        print time.time() - t0, "seconds time"


    # print '\n ********* ------ Completed Sweep ------ ********* \n'
    # print '\n ********* ------ Collect results ------ ********* \n'

    # NewFolder = "/data/vascrem/testoutput/ParameterSweep/Cylinder/CollectedResults/"
    # OldFolder = "/data/vascrem/testoutput/ParameterSweep/Cylinder/Parameteres3/"
    # os.mkdir(NewFolder)

    # Parameter1 = [5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
    # Parameter2 = [ 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
    # Parameter3 = [ 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]

    # for i in Parameter1:
    #     for j in Parameter2:
    #         for k in Parameter3:
    #             Oldfile = OldFolder + "Param_" + str(i) + "_DilationParam_" + str(j) + "_DeformationParam_"+ str(k) + "/results_from_time_30/results.viznodes" 
    #             if path.exists(Oldfile):
    #                 NewFile = NewFolder + "Area_" + str(i) + "_Dil_" + str(j) + "_Def_" + str(k)+".viznodes" 
    #                 shutil.copy(Oldfile, NewFile)


    # OldFolder = "/data/vascrem/testoutput/ParameterSweep/Cylinder/Parameteres2/"

    # Parameter1 = [5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
    # Parameter2 = [ 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
    # Parameter3 = [ 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]

    # for i in Parameter1:
    #     for j in Parameter2:
    #         for k in Parameter3:
    #             Oldfile = OldFolder + "Param_" + str(i) + "_DilationParam_" + str(j) + "_DeformationParam_"+ str(k) + "/results_from_time_30/results.viznodes" 
    #             if path.exists(Oldfile):
    #                 NewFile = NewFolder + "Area_" + str(i) + "_Dil_" + str(j) + "_Def_" + str(k)+".viznodes" 
    #                 shutil.copy(Oldfile, NewFile)

   
    print '\n ********* ------ Finished ------ ********* \n'
   

