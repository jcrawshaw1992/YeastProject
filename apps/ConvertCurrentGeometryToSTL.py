import vtk
import os
import math
from argparse import ArgumentParser

if __name__=="__main__":

     # Define arguments to be parsed
    parser = ArgumentParser(description='Converting the lastest chaste vtu into an stl that can be read into the remeshing tool')
    parser.add_argument('-ChasteOutput', dest='vtuDirectory',type=str, help=' Need the outputdirectory from which to sourse the latest .vtu to be defined')
    parser.add_argument('-stlOutput',dest='stlDirectory' , type=str,  default='/testouput/', help=' Where to save the stl, good to have it in an easy to reach location')
    args = parser.parse_args()
    # print args.vtuDirectory
    # LargestFile =0
    # for file in os.listdir(args.vtuDirectory):
    #         if file.startswith("mesh") and file.endswith(".vtu"):
    #             Index = file.find('_',1,len(file)) +1
    #             FileNumber = int(file[Index:-4])
    #             if (FileNumber> LargestFile):
    #                 LargestFile = FileNumber
    # # print LargestFile

    # vtuFile = args.vtuDirectory + 'mesh_'+str(LargestFile)+'.vtu'
    vtuFile = args.vtuDirectory + 'config.vtu'
    print vtuFile
    print " Convert vtu to stl "
    # Read the VTU file from disk
    vtu_reader = vtk.vtkXMLUnstructuredGridReader()
    vtu_reader.SetFileName(vtuFile)
    extract_surface_filter = vtk.vtkDataSetSurfaceFilter()
    extract_surface_filter.AddInputConnection(vtu_reader.GetOutputPort())
    # Write out the data in unstructured grid format
    stl_writer = vtk.vtkSTLWriter()
    stl_writer.SetInputConnection(extract_surface_filter.GetOutputPort())
    stl_writer.SetFileName(args.stlDirectory)
    stl_writer.Write()