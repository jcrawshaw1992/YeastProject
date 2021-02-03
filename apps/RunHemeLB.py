#!/usr/bin/env python
# vasc_res_sim.py
import subprocess
import shutil
import os
from xml.etree import ElementTree
import glob
from argparse import ArgumentParser
import numpy as np
import string

if __name__=="__main__":


    Collapse = ['3']
    for i in Collapse:
        mHemeLBDirectory = "/data/vascrem/testoutput/TestFlowThroughCollapse/"+i+"/HemeLBFluid/"
        mHemeLBPath ="/home/vascrem/hemelb-dev/"
        print "GmyUnstructuredGridReader for ", i
        GmyUnstructuredGridReader = "python " +mHemeLBPath+ "Tools/hemeTools/converters/GmyUnstructuredGridReader.py " + mHemeLBDirectory + "config.xml >nul"
        # subprocess.call(GmyUnstructuredGridReader, shell=True)
        print "GenerateFlowVtus for ", i
        GenerateFlowVtus = "python " +mHemeLBPath+ "Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory + "results/Extracted/surface-pressure.xtr " + mHemeLBDirectory + "results/Extracted/wholegeometry-velocity.xtr " + mHemeLBDirectory + "results/Extracted/surface-traction.xtr >nul"
        subprocess.call(GenerateFlowVtus, shell=True)
    print "Done"




    # # # Currently this code only collects the vtus at the end
    # Collapse = ['2']
    # for i in Collapse:
    #     mHemeLBDirectory = "/data/vascrem/testoutput/TestFlowThroughNonSymetricCollapse/"+i+"/HemeLBFluid/"
    #     mHemeLBPath ="/home/vascrem/hemelb-dev/"
    #     print "GmyUnstructuredGridReader for ", i
    #     GmyUnstructuredGridReader = "python " +mHemeLBPath+ "Tools/hemeTools/converters/GmyUnstructuredGridReader.py " + mHemeLBDirectory + "config.xml >nul"
    #     #subprocess.call(GmyUnstructuredGridReader, shell=True)
    #     print "GenerateFlowVtus for ", i
    #     GenerateFlowVtus = "python " +mHemeLBPath+ "Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory + "results/Extracted/surface-pressure.xtr " + mHemeLBDirectory + "results/Extracted/wholegeometry-velocity.xtr " + mHemeLBDirectory + "results/Extracted/surface-traction.xtr >nul"
    #     subprocess.call(GenerateFlowVtus, shell=True)


    # print "Done"
