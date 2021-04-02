#!/usr/bin/env python
import subprocess
import os
import pdb
from os import path
import time


def LoopOverSet(AreaParameter,DilationParameter,DeformationParamter,BendingParamter):
    t0 = time.time()
    SleepyTime =1*60
    counter =0

    Folder = "/data/vascrem/testoutput/ParameterSweepOnBifucation1/"
    NewDirectory = "/data/vascrem/testoutput/BifucationCharacterisitcs6/"

    if path.isdir(NewDirectory)==0:
        os.mkdir(NewDirectory)

    for i in AreaParameter:
        for j in DilationParameter:
            for k in DeformationParamter:
                for l in BendingParamter:
                    OutputDirectory = NewDirectory + "Area_"+  str(i)+"_Dil_"+  str(j)+"_Shear_"+  str(k)+"_Bend_"+  str(l)+"/"
                    if ~path.isdir(OutputDirectory):
                        WaitFileGeneration = NewDirectory+'TerminalOutput'+str(k)+'.txt'
                        subprocess.Popen(['./BifucationCharacteristicsSingleParameterSetBash', str(i), str(j), str(k),  str(l), Folder, NewDirectory,WaitFileGeneration])
                        counter = counter+1
                        if counter ==20:
                            print "SleepyTime"
                            counter=0
                            time.sleep(SleepyTime/4)
                            print "1/4th SleepyTime"
                            time.sleep(SleepyTime/4)
                            print "1/2 SleepyTime - If you get to here and nothing is going on, you should stop the code and reduce SleepyTime -- currently 4 minutes"
                            time.sleep(SleepyTime/4+SleepyTime/8)
                            print "3/4th SleepyTime - Need to be done by now "
                            time.sleep(SleepyTime/8)
                
                    


if __name__=="__main__":
    subprocess.call("chmod 700 BifucationCharacteristicsSingleParameterSetBash", shell=True)



    AreaParameter = [6,7,8,9,10]
    DilationParameter =  [5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
    DeformationParamter = [ 5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
    BendingParamter = [7,8,9,10,11,12,13,14,15,16]

    LoopOverSet(AreaParameter,DilationParameter,DeformationParamter,BendingParamter)

    
    