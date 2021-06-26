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
 
    OldFolder = "/Users/jcrawshaw/Downloads/PlexusExampleSmall/results_from_time_"
    NewFolder = "/Users/jcrawshaw/Downloads/PlexusExampleSmall/CollectedResults/"
    # os.mkdir(NewFolder)
    # NumberOfSimulations =  range(9,114)
    # NumberOfSimulations = NumberOfSimulations[0::4]
    NumberOfSimulations = [0, 0.1, 0.2,0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,2,2.1,2.2,2.4,2.6,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5,1.6]
    counter =0
    # Meshes
    # for i in NumberOfSimulations:
    #     MeshFiles = []
    #     ResultsFiles = []
    #     CurrentFolder =OldFolder+str(i)+"/"
    #     for file in os.listdir(CurrentFolder):
    #         if file.startswith("results_") and file.endswith(".vtu"):
    #             MeshFiles.append(int(file[5:-4]))
    #             # print file
    #         # if file.startswith("results_"):
    #     MeshFiles.sort()
    #     for j in MeshFiles:
    #         if j ==0:
    #             continue
    #         Oldfile = CurrentFolder+"results_"+str(j)+".vtu"
    #         NewFile = NewFolder + "results_" + str(counter) +".vtu" 
    #         counter = counter+1
    #         shutil.copy(Oldfile, NewFile)


    # Results
    counter = 0
    for i in NumberOfSimulations:
        MeshFiles = []
        ResultsFiles = []
        CurrentFolder =OldFolder+str(i)+"/"
        for file in os.listdir(CurrentFolder):
            if file.startswith("results_") and file.endswith(".vtu"):
                MeshFiles.append(int(file[8:-4]))
                # print file
            # if file.startswith("results_"):
        MeshFiles.sort()
        for j in MeshFiles:
            if j ==0:
                continue
            Oldfile = CurrentFolder+"results_"+str(j)+".vtu"
            NewFile = NewFolder + "results_" + str(counter) +".vtu" 
            counter = counter+1
            shutil.copy(Oldfile, NewFile)

