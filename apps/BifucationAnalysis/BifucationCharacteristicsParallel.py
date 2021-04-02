#!/usr/bin/env python
import subprocess
import os
import pdb
from os import path

if __name__=="__main__":

    DeformationParamter = [5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
    DeformationParamter = [6.7,]
    
    

    Folder = "/data/vascrem/testoutput/ParameterSweepOnBifucation1/"
    NewDirectory = "/data/vascrem/testoutput/BifucationCharacterisitcs3/"

    if path.isdir(NewDirectory)==0:
        os.mkdir(NewDirectory)


    for k in DeformationParamter:
        WaitFileGeneration = NewDirectory+'TerminalOutput'+str(k)+'.txt'
        subprocess.Popen(['./BifucationCharacteristicsBash', str(k), Folder, NewDirectory,WaitFileGeneration])
                    
        

    print "Collection is done "
      