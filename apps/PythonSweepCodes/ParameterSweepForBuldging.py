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
    GenerateRunner =0
    if GenerateRunner ==1:
        command = "cd ~/Chaste && scons b=GccOpt projects/VascularRemodelling/test/ParameterSweep/TestSweepingForBifucationDeformation.hpp"
        subprocess.call(command, shell=True)
    Server = 1
    if Server ==1:
        chaste_run_exe = '/home/vascrem/Chaste/projects/VascularRemodelling/build/optimised/ParameterSweep/TestSweepingForBifucationDeformationRunner '
        TerminalOutputFolder = "/data/vascrem/testoutput/TestForBuldgingInBifucation/SweepTerminalOutputs/"
        mesh_file = "/home/vascrem/MeshCollection/IdealisedNetwork/SimpleBifucation/mesh.vtu"

        if path.isdir("/data/vascrem/testoutput/TestForBuldgingInBifucation/")==0:
            os.mkdir("/data/vascrem/testoutput/TestForBuldgingInBifucation/")
    else:
        chaste_run_exe =  '/Users/jcrawshaw/Chaste/projects/VascularRemodelling/build/optimised/ParameterSweep/TestSweepingForBifucationDeformationRunner '
        TerminalOutputFolder = "/Users/jcrawshaw/Documents/testoutput/TestForBuldgingInBifucation/SweepTerminalOutputs/"
        mesh_file = "/Users/jcrawshaw/Documents/Projects/MeshCollection/SimpleBifucation/mesh.vtu"


    if path.isdir(TerminalOutputFolder)==0:
        os.mkdir(TerminalOutputFolder)

    # subprocess.call("chmod 700 RunChaste", shell=True)



    ParameterSets =[1,2,3,4,5,6,7,8]
    BendingParameters =[1,2,3,4,5,6,7,8,9]
    
    RunSweep = 0
    if RunSweep ==1:
        Parallel = 10
        SleepyTime = 300
        AvaliablePaths = range(Parallel)
        print AvaliablePaths
        for i in ParameterSets:
            for j in BendingParameters:
                Core = AvaliablePaths[0]
                Input1 = chaste_run_exe +' -ParameterSet '+str(i)+' -BendingParameter '+str(j)+ ' -MeshFile ' + mesh_file
                Input2 = TerminalOutputFolder+'ParameterSet'+str(i)+'BendingParameter'+str(j)+'.txt'
                Input3 = TerminalOutputFolder+'WaitFile'+str(Core)+'.txt'
                subprocess.Popen(['./RunChaste', Input1,Input2,Input3 ])
                AvaliablePaths.remove(Core) 
                    # # Check if all positions are taken
                while len(AvaliablePaths) ==0:
                    time.sleep(SleepyTime)
                    # print "Awake and checking for spare cores" 
                    for P in range(Parallel):
                        OutputFile = TerminalOutputFolder+'WaitFile'+str(P)+'.txt'
                        if path.exists(OutputFile):
                            AvaliablePaths.append(P)
                            os.remove(OutputFile)
                    if len(AvaliablePaths) >0:
                        print AvaliablePaths, "Have found a spare core or two :-) " 
                        print time.time() - t0, "seconds time"


    print '\n ********* ------ Completed Sweep ------ ********* \n'
    print '\n ********* ------ Collect results ------ ********* \n'

    CollectResults =1
    if CollectResults ==1:
        NewFolder = "/data/vascrem/testoutput/TestForBuldgingInBifucation/CollectedResults2/"
        OldFolder = "/data/vascrem/testoutput/TestForBuldgingInBifucation/"


        if path.exists(NewFolder)==0:
            os.mkdir(NewFolder)

        for i in ParameterSets:
            for j in BendingParameters:
                Oldfile = OldFolder + "PSet" + str(i) + "_Bending_" + str(j) + "/results_from_time_0/mesh_154000.vtu" 

                print Oldfile
                print path.exists(Oldfile)
                # Oldfile = OldFolder + "PSet" + str(i) + "_Bending_" + str(j) + "/results_from_time_30/results.viznodes" 
                if path.exists(Oldfile):
                    NewFile = NewFolder + "PSet" + str(i) + "_Bending_" + str(j) + ".mesh.vtu"
                    shutil.copy(Oldfile, NewFile)

        print '\n ********* ------ Finished ------ ********* \n'
   

