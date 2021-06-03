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

def RunTheCode(a,b,p):

    Input1 = chaste_run_exe +' -a '+str(a)+' -b '+str(b)+' -p ' +str(p) 
    Input2 = TerminalOutputFolder+'a_'+str(a)+'_b_'+str(b)+'p_' +str(p)+'.txt'
    Input3 = TerminalOutputFolder+'WaitFile'+str(a+b+p)+'.txt'
    subprocess.Popen(['./RunChaste', Input1,Input2,Input3 ])


if __name__=="__main__":
    
    t0 = time.time()
    time.sleep(60*30)
    GenerateRunner =0
    if GenerateRunner ==1:
        command = "cd ~/Chaste && scons b=GccOpt projects/VascularRemodelling/test/ParameterSweep/TestHetroCylinderExampleParameterTesting.hpp"
        subprocess.call(command, shell=True)


    chaste_run_exe =  '/home/vascrem/Chaste/projects//VascularRemodelling/build/optimised/ParameterSweep/TestHetroCylinderExampleParameterTestingRunner '
    TerminalOutputFolder = "/data/vascrem/testoutput/HetroCylinderParameters/SweepTerminalOutputs/"

    if path.isdir("/data/vascrem/testoutput/HetroCylinderParameters")==0:
        os.mkdir("/data/vascrem/testoutput/HetroCylinderParameters")

    if path.isdir(TerminalOutputFolder)==0:
        os.mkdir(TerminalOutputFolder)

    # subprocess.call("chmod 700 RunChaste", shell=True)

    # Altering B
    RunTheCode(8,2.2,12)
    RunTheCode(8,2,12)
    RunTheCode(8,0.5,12)
    RunTheCode(8,4,12)


    # Altering A
    RunTheCode(8,2,10)
    RunTheCode(12,2,12)
    RunTheCode(6,2,14)
    

    # Altering P
    RunTheCode(8,2,10)
    RunTheCode(8,0.5,12)
    RunTheCode(8,4,14)

    
   
    print '\n ********* ------ Finished ------ ********* \n'
   

