#!/usr/bin/env python
# from https://gist.githubusercontent.com/thewtex/8263132/raw/079eca1603ed750d1707516ef51eed33e4934b44/vtk-unstructuredgrid-to-stl.py
"""Convert PolyData in .vtk files to STL files."""

import sys
import vtk

def convert(InputFile, OutputFile): 
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(InputFile)

    surface_filter = vtk.vtkDataSetSurfaceFilter()
    surface_filter.SetInputConnection(reader.GetOutputPort())

    triangle_filter = vtk.vtkTriangleFilter()
    triangle_filter.SetInputConnection(surface_filter.GetOutputPort())

    writer = vtk.vtkSTLWriter()
    writer.SetFileName(OutputFile)
    writer.SetInputConnection(triangle_filter.GetOutputPort())
    writer.SetFileTypeToBinary()
    writer.Write()

if __name__ == '__main__':
    main()