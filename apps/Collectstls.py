#!/usr/bin/env python
# vasc_res_sim.py
# import vtk
import shutil
import os
import string
import math
import os.path
from os import path
import subprocess



if __name__=="__main__":

    # os.mkdir("/data/vascrem/testoutput/InitialInfection/")
    OldFolder = "/data/vascrem/testoutput/HemeLBSweep/FlowThrough3X3Collapse/UpperBranchFolder/"
    NewFolder = OldFolder +"/Collectedstls/"
    # /Volumes/Hardrive/VASCREM_server/SweepOnSphere/_Bend_8/results_from_time_0/mesh_60000.vtu
    if path.isdir(NewFolder)==0:
        os.mkdir(NewFolder)

    Parameter = [0,1,2,3,4,5,6,7,8,9,10, 5.5,5.6,5.7,5.8,5.9,6.1,6.2,6.3,6.4,6.5]
    for i in Parameter:
        Stls = OldFolder+ str(i) + "/config.stl"
        Newstl = NewFolder +"mesh_"+ str(i) +".stl"
        vtu = NewFolder +"mesh_"+ str(i) +".vtu"
        shutil.copy(Stls , Newstl )

                # ----  Convert files  -------------# 
        # convert = 'meshio-convert '+ Stls +'  '+ vtu
        # subprocess.call(convert, shell=True)

    
    print '\n ********* ------ Finished ------ ********* \n'
   