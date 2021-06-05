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
    GenerateRunner =1
    if GenerateRunner ==1:
        command = "cd ~/Chaste && scons b=GccOpt projects/VascularRemodelling/test/ParameterSweep/TestMembraneParameterSweep.hpp"
        subprocess.call(command, shell=True)

    chaste_run_exe = '/home/vascrem/Chaste/projects/VascularRemodelling/build/optimised/ParameterSweep/TestMembraneParameterSweepRunner '
    TerminalOutputFolder = "/data/vascrem/testoutput/MembraneParameterSweep/Cylinder/NewTerminalOutputTesting/"

    if path.isdir(TerminalOutputFolder)==0:
        os.mkdir(TerminalOutputFolder)

    Parameters = [[7,  9,  9.5] , [7,  9.5, 9.5],  [7, 10,  9.5],  [8,  9,  9.5],  [8,  9.5, 9.5],  [8, 10,  9.5],  [9,  9,  9.5],  [9,  9.5, 9.5],  [9, 10,  9.5], [10,  9,  9.5], [10,  9.5, 5.5], [10,  9.5, 6.5], [10,  9.5, 7], [10,  9.5, 7.5], [10,  9.5, 8], [10,  9.5, 8.5], [10,  9.5, 9], [10,  9.5, 9.5], [10,  9.5,10], [10, 10,  5.5], [10, 10,  6.5], [10, 10,  7], [10, 10,  7.5], [10, 10,  8], [10, 10,  8.5], [10, 10,  9], [10, 10,  9.5], [10, 10, 10]]


    Parameters = [[7,    8.5,   9.5],    [7,    9.0,   8.5],    [7,    9.5,   8.5],    [7,   10.0,   8.5],    [7,   10.0,  10.0],    [8,    8.5,   9.0],    [8,    8.5,   9.5],    [8,    8.5,  10.0],    [8,    9.0,   8.5],    [8,    9.0,  10.0],    [8,    9.5,   8.5],    [8,    9.5,   9.0],    [8,    9.5,  10.0],    [8,   10.0,   8.5],    [8,   10.0,   9.0],    [8,   10.0,  10.0],    [9,  8.5,   9.0],    [9,  8.5,   9.5],    [9,  8.5,  10.0],    [9,  9.0,   8.5],    [9,  9.0,  10.0],    [9,  9.5,   8.5],    [9,  9.5,   9.0],    [9,  9.5,  10.0],    [9, 10.0,   8.5],    [9, 10.0,   9.0],    [9, 10.0,  10.0],    [10,   5.0,   5.0],    [10,   5.0,   5.5],    [10,   5.0,   6.0],    [10,   5.0,   6.5],    [10,   5.0,   9.5],    [10,   5.0,  10.0],    [10,   8.5,   9.0],    [10,   8.5,   9.5],    [10,   8.5,  10.0],    [10,   9.0,   8.5],    [10,   9.0,  10.0],    [10,   9.5,   5.0],    [10,   9.5,   6.0],    [10,   9.5,  10.0],    [10,  10.0,   5.0],    [10,  10.0,   6.0],    [10,  10.0,  10.0]]


    ParameteresUnder50 =[[7.0,   8.5,   9.5],[ 7.0,   9.0,   8.5],[ 7.0,   9.5,   8.5],[ 7.0,  10.0,   8.5],[ 8.0,   8.5,   9.0],[ 8.0,   8.5,   9.5],[ 8.0,   8.5,  10.0],[ 8.0,   9.0,   8.5],[ 8.0,   9.5,   8.5],[ 8.0,  10.0,   8.5],[ 9.0,   8.5,   9.0],[ 9.0,   8.5,   9.5],[ 9.0,   8.5,  10.0],[ 9.0,   9.0,   8.5],[ 9.0,   9.5,   8.5],[ 9.0,  10.0,   8.5],[10.0,   5.0,   5.0],[10.0,   5.0,   5.5],[10.0,   5.0,   6.0],[10.0,   5.0,   6.5],[10.0,   5.0,   9.5],[10.0,   5.0,  10.0],[10.0,   8.5,   9.0],[10.0,   8.5,   9.5],[10.0,   8.5,  10.0],[10.0,   9.0,   8.5],[10.0,   9.5,   5.0],[10.0,   9.5,   6.0],[10.0,  10.0,   5.0],[10.0,  10.0,   6.0000]]
    ParameteresUnder100 = [[8.0,   9.0,  10.0],[ 8.0,   9.5,   9.0],[ 8.0,  10.0,   9.0],[ 9.0,   9.0,  10.0],[ 9.0,   9.5,   9.0],[ 9.0,  10.0,   9.0],[10.0,   9.0,  10.0000]]
    ParameteresUnder200 =[[ 7.0,  10.0,  10.0],[ 8.0,   9.5,  10.0],[9.0,   9.5,  10.0],[10.0,   9.5,  10.0],[10.0,  10.0,  10.0000]]
    ParameteresUnder300 =[[8, 10, 10], [9, 10, 10]]






    print len(Parameters)
    for I in Parameters:
        i = I[0]
        j = I[1]
        k = I[2]
        dt = ' 0.0001'
        # if i <5.5 or j<5.5 or k<5.5:
        #     SimulationTime = ' 60'
        # elif i <7.3 or j<7.3 or k<7.3:
        #     SimulationTime = ' 70'
        # elif i <7.3 or j<7.3 or k<7.3:
        #     SimulationTime = ' 70'
        # elif i <7.3 or j<7.3 or k<7.3:
        #     dt = ' 0.0001'
        # elif i <7.3 or j<7.3 or k<7.3:
        #     dt = ' 0.0001'
        #     SimulationTime = ' 150'

        SimulationTime = ' 100'

        Input1 = chaste_run_exe +' -AreaParameter '+str(i)+' -DilationParameter '+str(j)+' -DeformationParamter ' +str(k) +' -dt'+dt+' -NewEndTime'+SimulationTime
        Input2 = TerminalOutputFolder+'AreaParameter'+str(i)+'_DilationParameter'+str(j)+'_DeformationParamter' +str(k)+'.txt'
        Input3 = TerminalOutputFolder+'WaitFile'+str(1)+'.txt'
        subprocess.Popen(['./RunChaste', Input1,Input2,Input3 ])


    


    # OldFolder = "/data/vascrem/testoutput/MembraneParameterSweep/Cylinder/Parameters/"
    # NewFolder = "/data/vascrem/testoutput/MembraneParameterSweep/Cylinder/CollectedResults5th/"
    # os.mkdir(NewFolder)

    
   
    # for I in Parameters:
    #     i = I[0]
    #     j = I[1]
    #     k = I[2]
    #     Oldfile = OldFolder + "Param_" + str(i) + "_DilationParam_" + str(j) + "_DeformationParam_"+ str(k) + "/results_from_time_10/results.viznodes" 
    #     if path.exists(Oldfile):
    #         NewFile = NewFolder + "Area_" + str(i) + "_Dil_" + str(j) + "_Def_" + str(k)+".viznodes" 
    #         shutil.copy(Oldfile, NewFile)

   