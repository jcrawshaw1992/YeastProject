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
 
    OldFolder = "/data/vascrem/testoutput/FSICylinder/Medium/Hetro9/"
    NewFolder = OldFolder+"CollectedResults/"
    if path.exists(NewFolder) ==0:
        os.mkdir(NewFolder)
    NewMeshFolder = OldFolder+"CollectedMeshes/"
    if path.exists(NewMeshFolder) ==0:
        os.mkdir(NewMeshFolder)
  

    Simulations = np.arange(20, 181, 10)
    # results.viznodes

    # Results
    counter = 0
    # for i in Simulations:
    #     File = OldFolder+"results_from_time_"+str(i)+"/results.viznodes"
    #     NewFile  = NewFolder + "/results_"+str(i)+".viznodes"
    #     shutil.copy(File, NewFile)

    for i in Simulations:
        MeshFiles = []
        CurrentFolder = OldFolder+"results_from_time_"+str(i)+"/"
        for file in os.listdir(CurrentFolder):
            if file.startswith("mesh_") and file.endswith(".vtu"):
                MeshFiles.append(int(file[5:-4]))
        MeshFiles.sort()
        for j in MeshFiles:
            if j ==0:
                continue
            Oldfile = CurrentFolder+"mesh_"+str(j)+".vtu"
            NewFile = NewMeshFolder + "mesh_" + str(counter) +".vtu" 
            counter = counter+1
            shutil.copy(Oldfile, NewFile)

