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
 
    OldFolder = "/data/vascrem/testoutput/FSISimulations/VascularNetworkFine/GrowingToEqui/ConstantForceArchiving/"
    NewFolder = OldFolder+"CollectedResults/"
    if path.exists(NewFolder) ==0:
        os.mkdir(NewFolder)
    NewMeshFolder = NewFolder+"CollectedMeshes/"
    if path.exists(NewMeshFolder) ==0:
        os.mkdir(NewMeshFolder)
  

    Simulations = np.arange(0, 66,5)
    # results.viznodes

    # Results
    counter = 0
    for i in Simulations:
        MeshFiles = []
        # if i ==10 or i ==11 or i ==12 or i== 10.999999999999996:
        #     CurrentFolder = OldFolder+"results_from_time_"+str(int(i))+"/"
        # else:
        #     print round(i,2)
        #     CurrentFolder = OldFolder+"results_from_time_"+str(round(i,3))+"/"
        CurrentFolder = OldFolder+"results_from_time_"+str(i)+"/"
        for file in os.listdir(CurrentFolder):
            if file.startswith("mesh_") and file.endswith(".vtu"):
                MeshFiles.append(int(file[5:-4]))
        MeshFiles.sort()
       
        for j in MeshFiles:
            Oldfile = CurrentFolder+"mesh_"+str(j)+".vtu"
            NewFile = NewMeshFolder + "mesh_" + str(counter) +".vtu" 
            counter = counter+1
            shutil.copy(Oldfile, NewFile)
