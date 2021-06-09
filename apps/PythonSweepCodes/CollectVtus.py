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

    OldFolder = "/Users/jcrawshaw/Documents/testoutput/StepChangeHetroCylinder/Collapse/results_from_time_"
    NewFolder = "/Users/jcrawshaw/Documents/testoutput/StepChangeHetroCylinder/Collapse/CollectedResults/"
    # os.mkdir(NewFolder)
    NumberOfSimulations =  range(9,114)
    NumberOfSimulations = NumberOfSimulations[0::4]
    counter =0
    for i in NumberOfSimulations:
    

        MeshFiles = []
        ResultsFiles = []
        CurrentFolder =OldFolder+str(i)+"/"
        for file in os.listdir(CurrentFolder):
            if file.startswith("mesh_") and file.endswith(".vtu"):
                MeshFiles.append(int(file[5:-4]))
                # print file
            # if file.startswith("results_"):
        MeshFiles.sort()
        for j in MeshFiles:
            if j ==0:
                continue
            Oldfile = CurrentFolder+"mesh_"+str(j)+".vtu"
            NewFile = NewFolder + "mesh_" + str(counter) +".vtu" 
            counter = counter+1
            shutil.copy(Oldfile, NewFile)


        # print MeshFiles

            # if file.startswith("wholegeometry-velocity_"):



    # K = [9,13,17,21,25,29,33]
    # J = [100,200,300,400,500,600,700,800]

    # # K = [37,41,45,49,53]
    # # J = [500,1000,1500,2000,2500,3000,3500,4000]
    # counter = 0
    # OldFolder = "/Users/jcrawshaw/Documents/testoutput/StepChangeHetroCylinder/Collapse/results_from_time_"
    # NewFolder = "/Users/jcrawshaw/Documents/testoutput/StepChangeHetroCylinder/Collapse/CollectedResults/"
    # os.mkdir(NewFolder)
    # Oldfile = OldFolder  +"9/results_0.vtu"    
    # NewFile = NewFolder + "results_" + str(0) +".vtu" 
    # shutil.copy(Oldfile, NewFile)
    # Oldfile = OldFolder  +"9/mesh_0.vtu"    
    # NewFile = NewFolder + "mesh_" + str(0) +".vtu" 
    # shutil.copy(Oldfile, NewFile)
    # counter = counter+1

    # for k in K:
    #     for j in J:
    #         Oldfile = OldFolder  + str(k)+"/mesh_"+str(j)+".vtu"    
    #         NewFile = NewFolder + "mesh_" + str(counter) +".vtu" 
    #         shutil.copy(Oldfile, NewFile)

    #         Oldfile = OldFolder  + str(k)+"/results_"+str(j)+".vtu"    
    #         NewFile = NewFolder + "results_" + str(counter) +".vtu" 
    #         shutil.copy(Oldfile, NewFile)
    #         counter = counter+1
    #         print counter

    # K = [37,41,45,49,53]
    # J = [500,1000,1500,2000,2500,3000,3500,4000]
    
    # for k in K:
    #     for j in J:
    #         Oldfile = OldFolder  + str(k)+"/mesh_"+str(j)+".vtu"    
    #         NewFile = NewFolder + "mesh_" + str(counter) +".vtu" 
    #         shutil.copy(Oldfile, NewFile)

    #         Oldfile = OldFolder  + str(k)+"/results_"+str(j)+".vtu"    
    #         NewFile = NewFolder + "results_" + str(counter) +".vtu" 
    #         shutil.copy(Oldfile, NewFile)
    #         counter = counter+1
    #         print counter

   
    print '\n ********* ------ Finished ------ ********* \n'
   

