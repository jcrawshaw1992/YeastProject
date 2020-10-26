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
    
    NewFolder = "/data/vascrem/testoutput/Tissue2dVirusInfection_InfectiousPercent_S/CollectedResults/"
    OldFolder = "/data/vascrem/testoutput/Tissue2dVirusInfection_InfectiousPercent_S/"
    # "InitialInfec_0.05_Diffusion_0.5/6/results_from_time_0/pde_results_ext_virus_4000.vtu
    os.mkdir(NewFolder)
    Parameter1 = [0.05,0.1,0.01]
    Parameter2 =[1,0.1,0.5,0.01,0.05,0.005]
    Parameter3 = [1,2,3,4,5,6,7,8,9, 10]

    for i in Parameter1:
        for j in Parameter2:
            for k in Parameter3:
                Oldfile = OldFolder + "InitialInfec_" + str(i) + "_Diffusion_" + str(j) + "/" + str(k) + "/results_from_time_0/viruscelltypes.dat" 
                if path.exists(Oldfile):
                    NewFile = NewFolder + "InitialInfec_" + str(i) + "_Diffusion_" + str(j) + "_" + str(k)+"_viruscelltypes.dat" 
                    shutil.copy(Oldfile, NewFile)
   
    print '\n ********* ------ Finished ------ ********* \n'
   