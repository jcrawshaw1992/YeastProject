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



working_directory = '/Users/jcrawshaw/docker-polnet-master/NewMesh' 

def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")



def vtu2stl(iter):

    print "  Convert vtu to stl    "
    # Read the VTU file from disk
    vtu_reader = vtk.vtkXMLUnstructuredGridReader()
    vtu_reader.SetFileName(working_directory + '/config.vtk')

    extract_surface_filter = vtk.vtkDataSetSurfaceFilter()
    extract_surface_filter.AddInputConnection(vtu_reader.GetOutputPort())

    # Write out the data in unstructured grid format
    stl_writer = vtk.vtkSTLWriter()
    stl_writer.SetInputConnection(extract_surface_filter.GetOutputPort())
    stl_writer.SetFileName(working_directory + '/config.stl')
    stl_writer.Write()


if __name__=="__main__":

    vtu2stl(0)

