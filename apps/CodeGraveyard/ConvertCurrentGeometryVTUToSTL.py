import subprocess
import vtk
import shutil
import os
from xml.etree import ElementTree
import glob
import numpy as np
import time
import pdb
import string
import math
import sys
sys.path.insert(1, '/Applications/ParaView-5.4.1-822-g597adef982.app/Contents/Python/')
from paraview.simple import *



if __name__=="__main__":

     # Define arguments to be parsed
    parser = ArgumentParser(description='Run a vascular remodelling simulation')
    parser.add_argument('--num_iterations', default=Iterations, type=int, help='Number of Hemelb/Chaste iterations to be run (optional, default is 5).')
    parser.add_argument('--compute_radii', dest='compute_radii', action='store_true', help='Use VMTK to compute axis radii.')
    parser.add_argument('--output_postfix', dest='output_postfix', default='', help='This string will be added to ChasteWorkingDirectory to get the output folder.')
    parser.add_argument('--div_threshold', dest='div_threshold', default=1e10, help='This specifies the length that edges will divide at. (Defaults to no division, but 6e-4 is good for mm meshes')
    parser.add_argument('--mesh_scale', dest='mesh_scale', default=1e-3, help='This specifies what to scale the mesh by so that all distances are in meters (defaults to mm).')
    args = parser.parse_args()
   

   

def vtu2stl(file,WorkingDirectory , OutPutDirectory, iter):

    # print " Convert vtu to stl "
    # Read the VTU file from disk
    vtu_reader = vtk.vtkXMLUnstructuredGridReader()
    vtu_reader.SetFileName(WorkingDirectory +file)

    extract_surface_filter = vtk.vtkDataSetSurfaceFilter()
    extract_surface_filter.AddInputConnection(vtu_reader.GetOutputPort())
    # Write out the data in unstructured grid format
    stl_writer = vtk.vtkSTLWriter()
    stl_writer.SetInputConnection(extract_surface_filter.GetOutputPort())
    stl_writer.SetFileName(OutPutDirectory + 'mesh_' +str(iter)+'.stl')
    stl_writer.Write()

