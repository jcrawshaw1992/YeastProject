#!/usr/bin/env python
# vasc_res_sim.py
# import vtk
import shutil
import os
import string
import math
import os.path
from os import path



if __name__=="__main__":

    # os.mkdir("/data/vascrem/testoutput/InitialInfection/")
    NewFolder = "/Volumes/Hardrive/VASCREM_server/SweepOnSphere/EndPoints/mesh_"
    OldFolder = "/Volumes/Hardrive/VASCREM_server/SweepOnSphere/_Bend_"
    # /Volumes/Hardrive/VASCREM_server/SweepOnSphere/_Bend_8/results_from_time_0/mesh_60000.vtu
    # os.mkdir(NewFolder)
    OneVariable=1
    TwoVariables=0
    ThreeVariabes =0

    if OneVariable ==1:
         Parameter = [ 8,9,10,11,12,13,14,15]
         for i in Parameter:
            Oldfile = OldFolder+ str(i) + "/results_from_time_0/mesh_60000.vtu"
            if path.exists(Oldfile):
                NewFile = NewFolder +str(i) + ".vtu"
                shutil.copy(Oldfile, NewFile)



    if TwoVariables ==1:
        Parameter1 = [ 0.01, 0.05, 0.1 ]
        Parameter2 = [1,2,3,4,5,6,7,8,9, 10]

        for i in Parameter1:
            for j in Parameter2:
                Oldfile = OldFolder + "FractionInfected_" + str(i) + "/" + str(j)+ "/results_from_time_0/covid19celltypes.dat" 
                if path.exists(Oldfile):
                    NewFile = NewFolder + "FractionInfected_" + str(i) + "_" + str(j) + "_viruscelltypes.dat" 
                    shutil.copy(Oldfile, NewFile)


    NewFolder = "/data/vascrem/testoutput/InitialInfection/Random/"
    OldFolder = "/data/vascrem/testoutput/VaryingInitialCovidInfection/Random/"
    # os.mkdir(NewFolder)
    TwoVariables=1
    ThreeVariabes =0

    if TwoVariables==1:
        Parameter1 = [ 0.01, 0.05, 0.1 ]
        Parameter2 = [1,2,3,4,5,6,7,8,9, 10]

        for i in Parameter1:
            for j in Parameter2:
                Oldfile = OldFolder + "FractionInfected_" + str(i) + "/" + str(j)+ "/results_from_time_0/covid19celltypes.dat" 
                if path.exists(Oldfile):
                    NewFile = NewFolder + "FractionInfected_" + str(i) + "_" + str(j) + "_viruscelltypes.dat" 
                    shutil.copy(Oldfile, NewFile)


    NewFolder = "/data/vascrem/testoutput/InitialInfection/Clummped/"
    OldFolder = "/data/vascrem/testoutput/VaryingInitialCovidInfection/Clummped/"
    # os.mkdir(NewFolder)

    ThreeVariabes =1
    if ThreeVariabes==1:
        Parameter1 = [ 0.01, 0.05, 0.1 ]
        Parameter2 = [1,2,3]
        Parameter3 = [1,2,3,4,5,6,7,8,9, 10]

        for i in Parameter1:
            for j in Parameter2:
                for k in Parameter3:
                    Oldfile = OldFolder + "FractionInfected_" + str(i) + "/ClusterDesnity_" + str(j) + "/" + str(k) + "/results_from_time_0/covid19celltypes.dat" 
                    if path.exists(Oldfile):
                        NewFile = NewFolder + "FractionInfected_" + str(i) + "_ClusterDesnity_" + str(j) + "_" + str(k)+"_viruscelltypes.dat" 
                        shutil.copy(Oldfile, NewFile)
    
    
    print '\n ********* ------ Finished ------ ********* \n'
   