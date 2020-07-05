#!/usr/bin/env python
# vasc_res_sim.py
import subprocess
import vtk
import shutil
import os
from xml.etree import ElementTree
import glob
from argparse import ArgumentParser
import numpy as np
import time
import matplotlib.pyplot as plt
import csv
import pdb
import string
import math
# from FileConverter import vtuTostl
# import FSI_Fun


if __name__=="__main__":
    parser = ArgumentParser(description='Sort HemeLB vtu files ')
    parser.add_argument('-Directory', default = '/Users/jcrawshaw/Documents/ChasteWorkingDirectory/FlowFolder/', type=str, help='Need to supply a input folder')
    parser.add_argument('-CurrentNumberOfFiles', default=10, type=float, help='Need to know which how many files I an adding')
    args = parser.parse_args()    
# # ' ------- Setting up args ------- '
    print "/Users/jcrawshaw/Documents/testoutput/" + args.Directory+ "HemeLB_results_from_time_0/Extracted/"
    Directory  = "/Users/jcrawshaw/Documents/testoutput/" + args.Directory+ "HemeLB_results_from_time_0/Extracted/"
    HemeLBDirectory = "/Users/jcrawshaw/Documents/testoutput/" + args.Directory  + "HemeLB_results_from_time_0/"
    ChasteDirectory = "/Users/jcrawshaw/Documents/testoutput/" +  args.Directory +'results_from_time_0/'
    PriorFinalOutput = args.CurrentNumberOfFiles
    Directory  = "/Users/jcrawshaw/Documents/testoutput/" + args.Directory+ "HemeLBFluid/results_PriorTimeStep/Extracted/"

    # Get the vtu with the largest file number 
    VtuNumbering = []
    CopyTraction = 0
    CopyVelocity = 0
    # Collect the fluid vtus
    for file in os.listdir(Directory):
        if file.startswith("surface-pressure_") and file.endswith(".vtu"):
            VtuNumbering.append(int(file[17:-4]))
        if file.startswith("surface-traction_"):
            CopyTraction =1
        if file.startswith("wholegeometry-velocity_"):
            CopyVelocity =1

    FinalFile =  str(max(VtuNumbering)) 
 
    ChasteNumbering = []
    for file in os.listdir(ChasteDirectory):
        if file.startswith("mesh_") and file.endswith(".vtu"):
            FileNumber = int(file[5:-4])
            if FileNumber > PriorFinalOutput:
                ChasteNumbering.append(FileNumber)


    for j in ChasteNumbering:
        # print j
        shutil.copy(Directory + "surface-pressure_"+ FinalFile +".vtu" , HemeLBDirectory + 'surface-pressure_'+  str(j) +'.vtu')
        if CopyTraction == 1:
            shutil.copy(Directory + "surface-traction_"+ FinalFile +".vtu", HemeLBDirectory + 'surface-traction_'+  str(j) +'.vtu')
        if CopyVelocity == 1:
            shutil.copy(Directory + "wholegeometry-velocity_"+ FinalFile +".vtu", HemeLBDirectory + 'wholegeometry-velocity_'+  str(j) +'.vtu')

    MaxNumber = str(max(ChasteNumbering))
    outF = open(HemeLBDirectory+"CurrentLastFluidOutput.txt", "w")
    outF.write(MaxNumber)
    outF.close()
    # print '\n ********* ------ vtus from last time step sorted ------ ********* \n'