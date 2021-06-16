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

    OldFolder = "/Users/jcrawshaw/Downloads/PlexusExample/results_from_time_"
    NewFolder = "/Users/jcrawshaw/Downloads/PlexusExample/CollectedResults/"
    os.mkdir(NewFolder)
    # NumberOfSimulations =  range(9,114)
    # NumberOfSimulations = NumberOfSimulations[0::4]
    NumberOfSimulations = [0, 0.2,0.4,0.6,0.8,1]#,1.2,1.7,2.2,2.7,3.2    ]
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

