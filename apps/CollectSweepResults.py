#!/usr/bin/env python
# vasc_res_sim.py
# import vtk
import shutil
import os
import string
import math



if __name__=="__main__":
    
    NewFolder = "/data/vascrem/testoutput/ParameterSweep/Cylinder/Parameteres2/CollectedResults/"
    os.mkdir(NewFolder)
    AreaParameter = [6, 6.5, 7,7.5, 8]
    DilationParameter =[6, 6.5, 7,7.5, 8]
    DeformationParamter = [6, 6.5, 7,7.5, 8]

    for i in AreaParameter:
        for j in DilationParameter:
            for k in DeformationParamter:
                Oldfile = "/data/vascrem/testoutput/ParameterSweep/Cylinder/Parameteres2/Param_"+  str(i) +"_DilationParam_"+str(j) + "_DeformationParam_" +str(k) +"/results_from_time_30/results.viznodes" 

                NewFile = NewFolder+"Area_"+  str(i) +"_Dil_"+str(j) + "_Def_" +str(k)+".viznodes"
                # print"Hit"
                shutil.copy(Oldfile, NewFile)
   
    print '\n ********* ------ Finished ------ ********* \n'
   