#!/usr/bin/env python
# vasc_res_sim.py
import subprocess
import shutil
import os
import glob




if __name__=="__main__":


    # chmod 700 RunHemeLBSweepBash
    # subprocess.call("chmod 700 RunHemeLBSweepBash", shell=True)
    TerminalOutputFolder = '/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3Collapse/AngleVariation_3X3Network/'
    MeshDirectory = "/data/vascrem/MeshCollection/IdealisedNetwork/AngleVariation_3X3Network/"

    print MeshDirectory

    Collapse = ['0.6','0.55','0.56','0.57','0.58','0.59','0.61','0.62','0.63','0.64','0.65']
   
    # Angle == 3: Seed = [2.19658, 1.97595, -0.155648]
    # Angle == 6:  Seed = [2.05519, 1.98863, 0.165396]
    # Angle == 2.2: Seed = [2.11408, 2.11277, 0.199536]

    counter = 0
    BifucationAngles = [ 'PI_2.2/','PI_3/','PI_6/' ] 
    Seeds = [[2.11408, 2.11277, 0.199536] , [2.19658, 1.97595, -0.155648],[2.05519, 1.98863, 0.165396] ]
    for angle in BifucationAngles:
        counter = counter+1
        for i in Collapse:
            # print Seeds[counter]
            print counter
            # mHemeLBDirectory = TerminalOutputFolder+angle+str(float(i)*10)+'/'
            # print mHemeLBDirectory
     
     

            # MeshFile = MeshDirectory+angle+"Mesh_"+i+".stl"
            # print MeshFile
            # dX = 0.08/15*(float(i))/5
            # print dX

           
    


            
            