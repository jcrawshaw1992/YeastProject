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
    TerminalOutputFolder = "/data/vascrem/testoutput/MembraneParameterSweep/Cylinder/SweepTerminalOutputsNew/"

    if path.isdir(TerminalOutputFolder)==0:
        os.mkdir(TerminalOutputFolder)

    Parameters = [[10, 9,  9],   [10, 9, 10],   [9, 9, 9],   [9, 9, 10],   [9, 9.5, 9],   [9, 9.5, 10],   [9, 10, 9],   [9, 10, 10],  [8, 9,  9], [8, 9, 10], [8, 9.5, 9], [8, 9.5, 10], [8, 10, 9], [8, 10, 10], [7, 8.5, 9], [7, 8.5, 10], [7, 9, 9], [7, 9, 10], [7, 9.5, 9], [7, 9.5, 10], [7, 10, 9], [7, 10, 10]]
    print len(Parameters)
    print Parameters
    for I in Parameters:
        i = I[0]
        j = I[1]
        k = I[2]
        dt = ' 0.0005'
        if i <7.3 or j<7.3 or k<7.3:
            SimulationTime = ' 60'
            dt = ' 0.0001'
        else:
            SimulationTime = ' 60'

        Input1 = chaste_run_exe +' -AreaParameter '+str(i)+' -DilationParameter '+str(j)+' -DeformationParamter ' +str(k) +' -dt'+dt+' -NewEndTime'+SimulationTime
        Input2 = TerminalOutputFolder+'AreaParameter'+str(i)+'_DilationParameter'+str(j)+'_DeformationParamter' +str(k)+'.txt'
        Input3 = TerminalOutputFolder+'WaitFile'+str(1)+'.txt'
        subprocess.Popen(['./RunChaste', Input1,Input2,Input3 ])