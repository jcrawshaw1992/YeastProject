import os
import vtk
import argparse
import pdb
import numpy
import meshio

def convertFile(vtkFile, outdir):
    print os.path.isfile(vtkFile)
    if os.path.isfile(vtkFile):
        basename = os.path.basename(vtkFile)
        print("Copying file:", basename)
        basename = os.path.splitext(basename)[0]

        outfile = os.path.join(outdir, basename+".vtu")

        reader = vtk.vtkPolyDataReader() 
   
        reader.SetFileName(outfile)
        reader.Update()
        output = reader.GetOutput()

        Nodes = []
        ElementList = []
        # Add the first Node so this isnt an emtpy list
        Nodes.append(output.GetCell(0).GetPoints().GetPoint(0))
        
        # loop over all of the elements in the vtk polydata file 
        # For each element, the polydata has recorded the location of each node,
        # but not the node index. As such here we loop over the elements 
        # and select out the node, which are put in the Node list. 
        # As I loop over the elements I will also record which node indices are in each element
        for i in range(output.GetNumberOfCells()):
            Element = []
            pts = output.GetCell(i).GetPoints()  
            np_pts = numpy.array([pts.GetPoint(i) for i in range(pts.GetNumberOfPoints())]) 
            
            # Saving each of the nodes in this element
            for j in [0,1,2]:
                AddNodeToList = 1 
                # want to add the Node to the Nodes list, But need to check if it is in the list yet
                # Do this by looping over list and seeing if it is in there
                for k in Nodes: 
                    if(k == pts.GetPoint(j)):
                        AddNodeToList =0
                if (AddNodeToList ==1):
                    # print 'Adding Node'
                        # Filling in the array storing Node locations
                    Nodes.append(pts.GetPoint(j))
                NodeIndex = Nodes.index(pts.GetPoint(j))
                Element.append(NodeIndex)
                # This works to fill in Nodes
            ElementList.append(Element)
        pdb.set_trace()

        # Save the nodes as points and elements for the meshio writer
        points = numpy.array(Nodes)
        elements = {
        "triangle": numpy.array(ElementList
        )
        }    

        meshio.write_points_cells(
        "NewPlexus2.vtu",
        points,
        elements,
        # Optionally provide extra data on points, cells, etc.
        # point_data=point_data,
        # cell_data=cell_data,
        # field_data=field_data
        )
 


        print 'Finished'

def convertFiles(indir, outdir):
    files = os.listdir(indir)

    files = [ os.path.join(indir,f) for f in files if f.endswith('.vtk') ]

    ret = 0
    print("In:", indir)
    print("Out:", outdir)
    for f in files:
        ret += convertFile(f, outdir)
    print("Successfully converted %d out of %d files." % (ret, len(files)))

def run():
    inf = "/Users/jcrawshaw/docker-polnet-master/ScalledMesh/"
    # out = "/Users/jcrawshaw/docker-polnet-master/NewMesh/"
    # convertFiles(args.indir, args.outdir)
    print "B"
    convertFiles(inf,inf)

def GenerateAVTUTriangle():

    points = numpy.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [2.0, 0.0, 0.0]
        ])
    elements = {
    "triangle": numpy.array([
        [0, 1, 2], [1, 2, 3]
        ])
    }      
    meshio.write_points_cells(
        "foo.vtu",
        points,
        elements,
        # Optionally provide extra data on points, cells, etc.
        # point_data=point_data,
        # cell_data=cell_data,
        # field_data=field_data
        )

    pdb.set_trace()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="VTK to STL converter")
    parser.add_argument('-ifile', help="Path to input directory.")
    args = parser.parse_args()    

    # ' ------- Setting up args ------- '

    vtkFile  = "/Users/jcrawshaw/docker-polnet-master/ScalledMesh/PlexusMesh.vtk"# inargs.ifile #+ 'results_from_time_0/'
    OutPutDirectory = "/Users/jcrawshaw/docker-polnet-master/ScalledMesh/"
    print "A"
    convertFile(vtkFile, OutPutDirectory)


    # GenerateAVTUTriangle()
    #https://stackoverflow.com/questions/45717843/writing-vtk-unstrucutredgrid-to-file-in-python
    # run
    


    print 'Finito'
 


















    # Grave Yard

    # # We'll create the building blocks of polydata including data attributes.
    # cube    = vtk.vtkPolyData()
    # points  = vtk.vtkPoints()
    # polys   = vtk.vtkCellArray()
    # scalars = vtk.vtkFloatArray()
    

    # # Load the point, cell, and data attributes.
    # for i in range(8):
    #     points.InsertPoint(i, x[i])
    # for i in range(6):
    #     polys.InsertNextCell( mkVtkIdList(pts[i]) )
    # for i in range(8):
    #     scalars.InsertTuple1(i,i)
    # pdb.set_trace()

    # writer = vtk.vtkUnstructuredGridWriter()
    # unstructuredGrid.SetCells(points, polys)
    # # writer.SetInputConnection(cube)
    # writer.SetFileName("/Users/jcrawshaw/docker-polnet-master/NewMesh/Cube.vtu")

    # # We now assign the pieces to the vtkPolyData.
    # cube.SetPoints(points)
    # del points
    # cube.SetPolys(polys)
    # del polys
    # cube.GetPointData().SetScalars(scalars)
    # del scalars

    # # Now we'll look at it.
    # cubeMapper = vtk.vtkPolyDataMapper()
    # if vtk.VTK_MAJOR_VERSION <= 5:
    #     cubeMapper.SetInput(cube)
    # else:
    #     cubeMapper.SetInputData(cube)
    # cubeMapper.SetScalarRange(0,7)
    # cubeActor = vtk.vtkActor()
    # cubeActor.SetMapper(cubeMapper)

    # # The usual rendering stuff.
    # camera = vtk.vtkCamera()
    # camera.SetPosition(1,1,1)
    # camera.SetFocalPoint(0,0,0)

    # renderer = vtk.vtkRenderer()
    # renWin   = vtk.vtkRenderWindow()
    # renWin.AddRenderer(renderer)

    # iren = vtk.vtkRenderWindowInteractor()
    # iren.SetRenderWindow(renWin)

    # renderer.AddActor(cubeActor)
    # renderer.SetActiveCamera(camera)
    # renderer.ResetCamera()
    # renderer.SetBackground(1,1,1)

    # renWin.SetSize(300,300)

    # # interact with data
    # renWin.Render()
    # iren.Start()

    # # Clean up
    # del cube
    # del cubeMapper
    # del cubeActor
    # del camera
    # del renderer
    # del renWin
    # del iren