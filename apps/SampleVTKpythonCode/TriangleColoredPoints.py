import vtk

def get_program_parameters():
    import argparse
    description = 'Generate a triangle with colored points and write it to a .vtp file.'
    epilogue = '''
   '''
    parser = argparse.ArgumentParser(description=description, epilog=epilogue)
    parser.add_argument('filename', help='A required vtp filename.', nargs='?',
                        const='TestTriangleColoredPoints.vtp',
                        type=str, default='TestTriangleColoredPoints.vtp')
    args = parser.parse_args()
    return args.filename


def main():
    colors = vtk.vtkNamedColors()

    filename = 'TrialCenterlines.vtp'#get_program_parameters()

    # setup points and vertices
    Points = vtk.vtkPoints()
    Vertices = vtk.vtkCellArray()

    id = Points.InsertNextPoint(1.0, 0.0, 0.0)
    Vertices.InsertNextCell(1)
    Vertices.InsertCellPoint(id)
    id = Points.InsertNextPoint(0.0, 0.0, 0.0)
    Vertices.InsertNextCell(1)
    Vertices.InsertCellPoint(id)
    id = Points.InsertNextPoint(0.0, 1.0, 0.0)
    Vertices.InsertNextCell(1)
    Vertices.InsertCellPoint(id)

    # # setup colors
    # Colors = vtk.vtkUnsignedCharArray()
    # Colors.SetNumberOfComponents(3)
    # Colors.SetName("Colors")
    # Colors.InsertNextTuple3(*colors.GetColor3ub('Red'))
    # Colors.InsertNextTuple3(*colors.GetColor3ub('LimeGreen'))
    # Colors.InsertNextTuple3(*colors.GetColor3ub('Blue'))

     # setup colors
    Radii = vtk.vtkDoubleArray()
    Radii.SetNumberOfComponents(1)
    Radii.SetName("Radius")
    Radii.InsertNextTuple1(1)
    Radii.InsertNextTuple1(2)
    Radii.InsertNextTuple1(3)


    # ssert point_data.GetArrayName(0) == 'Radius', "VTP file doesn't contain array Radius"
    #     radii = point_data.GetArray(0)
    #     radii2 = point_data.GetArray(1)
        

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(Points)
    polydata.SetVerts(Vertices)
    # polydata.GetPointData().SetScalars(Colors)
    polydata.GetPointData().Setdouble(Radii)
    polydata.Modified()

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(polydata)
    writer.Write()


if __name__ == '__main__':
    main()