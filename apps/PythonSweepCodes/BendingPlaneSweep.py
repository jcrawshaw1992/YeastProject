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
        command = "cd ~/Chaste && scons b=GccOpt projects/VascularRemodelling/test/ParameterSweep/TestBendingOnPlaneSweep.hpp"
        subprocess.call(command, shell=True)

    chaste_run_exe = '/home/vascrem/Chaste/projects/VascularRemodelling/build/optimised/ParameterSweep/TestBendingOnPlaneSweepRunner '
    TerminalOutputFolder = "/data/vascrem/testoutput/BendingForceOnBentRectanlge/SweepTerminalOutputs/"

    if path.isdir("/data/vascrem/testoutput/BendingForceOnBentRectanlge")==0:
        os.mkdir("/data/vascrem/testoutput/BendingForceOnBentRectanlge")

    if path.isdir(TerminalOutputFolder)==0:
        os.mkdir(TerminalOutputFolder)

    AspectRatio = [ 0.75,1.5,3] 
    Refinment = [ 7,15,25]
       
    for i in AspectRatio:
        for j in Refinment:
    
            Input1 = chaste_run_exe +' -Nc '+str(j)+' -AspectRatio '+str(i)
            Input2 = TerminalOutputFolder+'AspectRatio_'+str(i)+'_Nc_'+str(j)+'.txt'
            Input3 = TerminalOutputFolder+'WaitFile'+str(1)+'.txt'
            subprocess.Popen(['./RunChaste', Input1,Input2,Input3 ])
   
    print '\n ********* ------ Finished ------ ********* \n'
   

