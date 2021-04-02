import vtk
from argparse import ArgumentParser

if __name__=="__main__":

    parser = ArgumentParser(description='Get the radius file ')
    
    parser.add_argument('-Directory', type=str, help='Need to supply a input folder')
    args = parser.parse_args()    
    vtuFile  = args.Directory + "config.vtu"
    stlFile  = args.Directory + "config2.stl"
    # _______________________
    #  Convert vtu to stl
    #  _______________________
    print "  Convert vtu to stl    "
    # Read the VTU file from disk
    vtu_reader = vtk.vtkXMLUnstructuredGridReader()
    vtu_reader.SetFileName(vtuFile)

    extract_surface_filter = vtk.vtkDataSetSurfaceFilter()
    extract_surface_filter.AddInputConnection(vtu_reader.GetOutputPort())

    # Write out the data in unstructured grid format
    stl_writer = vtk.vtkSTLWriter()
    stl_writer.SetInputConnection(extract_surface_filter.GetOutputPort())
    stl_writer.SetFileName(stlFile)
    stl_writer.Write()
    print "  Doneies    "



 