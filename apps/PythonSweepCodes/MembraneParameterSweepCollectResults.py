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


    print '\n ********* ------ Completed Sweep ------ ********* \n'
    print '\n ********* ------ Collect results ------ ********* \n'

    OldFolder = "/data/vascrem/testoutput/MembraneParameterSweep/Cylinder/Parameters/"
    NewFolder = "/data/vascrem/testoutput/MembraneParameterSweep/Cylinder/CollectedResults/"
    os.mkdir(NewFolder)

    
    AreaParameter = [ 10, 9, 8,7,6,  5] 
    DilationParameter = [ 9.5, 10, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9]
    DeformationParamter = [ 9.5, 10 , 5 , 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9]

    for i in AreaParameter:
        for j in DilationParameter:
            for k in DeformationParamter:
                Oldfile = OldFolder + "Param_" + str(i) + "_DilationParam_" + str(j) + "_DeformationParam_"+ str(k) + "/results_from_time_10/results.viznodes" 
                if path.exists(Oldfile):
                    NewFile = NewFolder + "Area_" + str(i) + "_Dil_" + str(j) + "_Def_" + str(k)+".viznodes" 
                    shutil.copy(Oldfile, NewFile)

   
    print '\n ********* ------ Finished ------ ********* \n'
   

