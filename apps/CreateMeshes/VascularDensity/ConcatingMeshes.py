import math
import stl
from stl import mesh
import numpy


#
import subprocess
import vtk
import os
import pdb
import string
import math as m

import array as A
# from vtk import vtkXMLUnstructuredGridWriter as Writer, VTK_VERSION
# from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
import numpy as np
import meshio



import math

from os import path



# docker pull pymesh/pymesh:py3.7-slim
# sudo docker cp /Volumes/Hardrive/Projects/FSI/NetworkDensity/Collectedstls/combined.stl  pymesh/pymesh:/~



# find the max dimensions, so we can know the bounding box, getting the height,
# width, length (because these are the step size)...
def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")




def Remove(Folder):
    rm = 'rm ' +Folder
    subprocess.call(rm, shell=True)

def Remesh(Directory,Mesh, NewMesh, server):

    PreRemesh  = Directory+'CombindMesh.off'
    Remeshed =  Directory+'CombindMesh_Remeshed.off'
    Remeshedstl =  Directory+NewMesh
    print NewMesh
    # command = 'vmtksurfaceremeshing -ifile '+Directory+Mesh +' -iterations 5 -edgelength 0.01 -elementsizemode "edgelength" -ofile ' + Remeshed
    # subprocess.call(command, shell=True)

    command = 'meshio-convert  ' + Directory+Mesh + ' ' +PreRemesh
    subprocess.call(command, shell=True)
    if path.exists(PreRemesh)==0:
        pause()
    server = 0

    # if server ==1:
    #     CGALRemeshingCommand = "(cd  /home/vascrem/CGAL-5.0.2/Polygon_mesh_processing/examples/Polygon_mesh_processing/;./isotropic_remeshing_ForChaste -input " + PreRemesh  + " -output " + Remeshed + " -target_edge_length 0.015 -iterations 5 > null)"
    #     # SystemOutput = std::system(CGALRemeshingCommand.c_str()); // system only takes char *
    #     subprocess.call(CGALRemeshingCommand, shell=True)
    # else:
    #     CGALRemeshingCommand = "(cd  ~/Documents/CGAL-5.0.2/examples/Polygon_mesh_processing/;./isotropic_remeshing_ForChaste -input " + PreRemesh  + " -output " + Remeshed + " -target_edge_length 0.015 -iterations 5 > null)"
    #     subprocess.call(CGALRemeshingCommand, shell=True)
    # if path.exists(Remeshed)==0:

    command = 'vmtksurfaceremeshing -ifile '+Directory+Mesh +' -iterations 5 -edgelength 0.01 -elementsizemode "edgelength" -ofile ' + Remeshedstl
    subprocess.call(command, shell=True)

    # command = 'meshio-convert  ' + Remeshed + ' ' +Remeshedstl
    # subprocess.call(command, shell=True)
    # Remove(Remeshed)
    # Remove(PreRemesh)


def ScaleAndClippMeshes(Directory,List , CreateUpperRegion, CreateMiddelRegion, CreateBottomlRegion, SingleHoneycomb,TopHoneycomb, BottomHoneycomb):
    from vtk.util import numpy_support
    import ConvertVTKtoSTL
    import vtk
    import clip
    import vmtk
    from vmtk import pypes
    if CreateUpperRegion ==1:
        for i in List:
            # ----  Clip Edges -------------# 
            # print "Jess is good"
            vtkMesh = Directory +'MeshClipped'+str(i)+'.vtk'
            Clip1 = Directory +'NewClippedMesh_'+str(i)+'.vtk'
    
            clip.clip_surface_with_plane(vtkMesh,(0,-0.7,0), (0,1,0), Clip1)

            Outputstl = Directory+"Mesh_PreScale"+str(i)+".stl"
            ConvertVTKtoSTL.convert(Clip1, Outputstl)

            Outputvtu = Directory+"Mesh_PreScale"+str(i)+".vtu"
            command = 'meshio-convert  ' + Outputstl + ' ' +Outputvtu
            subprocess.call(command, shell=True)

            ScaledMesh = Directory+"Mesh_ScaledAndClipped"+str(i)+".vtu"
            Scale = 'vmtkmeshscaling -ifile '+ Outputvtu + ' -scale 0.20  --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ ScaledMesh
            subprocess.call(Scale, shell=True)

            stlMesh = Directory +'UpperSection'+str(i)+'.stl '
            command = 'meshio-convert  ' + ScaledMesh + ' ' +stlMesh
            subprocess.call(command, shell=True)

            Remove(Clip1)
            Remove(Outputvtu)
            Remove(ScaledMesh)
            Remove(Outputstl)

    if CreateMiddelRegion ==1:
        # ----  Clip Edges -------------# 
        # print "Jess is good"
        vtkMesh = Directory +'MeshClipped1.0.vtk'
        Clipped1 = Directory +'CenteralRegion.vtk'

        clip.clip_surface_with_plane(vtkMesh,(0,0.7,0), (0,-1,0), Clipped1)
        clip.clip_surface_with_plane(Clipped1,(0,-0.7,0), (0,1,0), Clipped1)

        Outputstl = Directory+"Mesh_PreScale.stl"
        ConvertVTKtoSTL.convert(Clipped1, Outputstl)

        Outputvtu = Directory+"Mesh_PreScale.vtu"
        command = 'meshio-convert  ' + Outputstl + ' ' +Outputvtu
        subprocess.call(command, shell=True)

        ScaledMesh = Directory+"Mesh_ScaledAndClipped.vtu"
        Scale = 'vmtkmeshscaling -ifile '+ Outputvtu + ' -scale 0.20  --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ ScaledMesh
        subprocess.call(Scale, shell=True)

        stlMesh = Directory +'MiddelSection.stl '
        command = 'meshio-convert  ' + ScaledMesh + ' ' +stlMesh
        subprocess.call(command, shell=True)

        Remove(Clipped1)
        Remove(Outputvtu)
        Remove(ScaledMesh)
        Remove(Outputstl)


    if CreateBottomlRegion ==1:
        # ----  Clip Edges -------------# 
        # print "Jess is good"
        vtkMesh = Directory +'MeshClipped1.0.vtk'
        Clipped1 = Directory +'CenteralRegion.vtk'

        clip.clip_surface_with_plane(vtkMesh,(0,-0.7,0), (0,-1,0), Clipped1)

        Outputstl = Directory+"Mesh_PreScale.stl"
        ConvertVTKtoSTL.convert(Clipped1, Outputstl)

        Outputvtu = Directory+"Mesh_PreScale.vtu"
        command = 'meshio-convert  ' + Outputstl + ' ' +Outputvtu
        subprocess.call(command, shell=True)

        ScaledMesh = Directory+"Mesh_ScaledAndClipped.vtu"
        Scale = 'vmtkmeshscaling -ifile '+ Outputvtu + ' -scale 0.20  --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ ScaledMesh
        subprocess.call(Scale, shell=True)

        stlMesh = Directory +'BottomSection.stl '
        command = 'meshio-convert  ' + ScaledMesh + ' ' +stlMesh
        subprocess.call(command, shell=True)

        Remove(Clipped1)
        Remove(Outputvtu)
        Remove(ScaledMesh)
        Remove(Outputstl)

    
    if SingleHoneycomb ==1:
        # ----  Clip Edges -------------# 
        # print "Jess is good"
        vtkMesh = Directory +'MeshClipped1.0.vtk'
        Clipped1 = Directory +'CenteralRegion.vtk'

        clip.clip_surface_with_plane(vtkMesh,(0,0.7,0), (0,-1,0), Clipped1)
        clip.clip_surface_with_plane(Clipped1,(0,-0.7,0), (0,1,0), Clipped1)
        clip.clip_surface_with_plane(Clipped1,(3.900000050663948,0,0), (-1,0,0), Clipped1)

        Outputstl = Directory+"Mesh_PreScale.stl"
        ConvertVTKtoSTL.convert(Clipped1, Outputstl)

        Outputvtu = Directory+"Mesh_PreScale.vtu"
        command = 'meshio-convert  ' + Outputstl + ' ' +Outputvtu
        subprocess.call(command, shell=True)

        ScaledMesh = Directory+"Mesh_ScaledAndClipped.vtu"
        Scale = 'vmtkmeshscaling -ifile '+ Outputvtu + ' -scale 0.20  --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ ScaledMesh
        subprocess.call(Scale, shell=True)

        stlMesh = Directory +'Honeycomb.stl '
        command = 'meshio-convert  ' + ScaledMesh + ' ' +stlMesh
        subprocess.call(command, shell=True)

    if TopHoneycomb ==1:
        # ----  Clip Edges -------------# 
        # print "Jess is good"
        vtkMesh = Directory +'MeshClipped1.0.vtk'
        Clipped1 = Directory +'CenteralRegion.vtk'

        clip.clip_surface_with_plane(vtkMesh,(0,0.7,0), (0,1,0), Clipped1)
        clip.clip_surface_with_plane(Clipped1,(3.900000050663948,0,0), (-1,0,0), Clipped1)

        Outputstl = Directory+"Mesh_PreScale.stl"
        ConvertVTKtoSTL.convert(Clipped1, Outputstl)

        Outputvtu = Directory+"Mesh_PreScale.vtu"
        command = 'meshio-convert  ' + Outputstl + ' ' +Outputvtu
        subprocess.call(command, shell=True)

        ScaledMesh = Directory+"Mesh_ScaledAndClipped.vtu"
        Scale = 'vmtkmeshscaling -ifile '+ Outputvtu + ' -scale 0.20  --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ ScaledMesh
        subprocess.call(Scale, shell=True)

        stlMesh = Directory +'TopHoneycomb.stl '
        command = 'meshio-convert  ' + ScaledMesh + ' ' +stlMesh
        subprocess.call(command, shell=True)

        Remove(Clipped1)
        Remove(Outputvtu)
        Remove(ScaledMesh)
        Remove(Outputstl)

    if BottomHoneycomb ==1:
        # ----  Clip Edges -------------# 
        # print "Jess is good"
        vtkMesh = Directory +'MeshClipped1.0.vtk'
        Clipped1 = Directory +'CenteralRegion.vtk'

        clip.clip_surface_with_plane(vtkMesh,(0,-0.7,0), (0,-1,0), Clipped1)
        clip.clip_surface_with_plane(Clipped1,(3.900000050663948,0,0), (-1,0,0), Clipped1)

        Outputstl = Directory+"Mesh_PreScale.stl"
        ConvertVTKtoSTL.convert(Clipped1, Outputstl)

        Outputvtu = Directory+"Mesh_PreScale.vtu"
        command = 'meshio-convert  ' + Outputstl + ' ' +Outputvtu
        subprocess.call(command, shell=True)

        ScaledMesh = Directory+"Mesh_ScaledAndClipped.vtu"
        Scale = 'vmtkmeshscaling -ifile '+ Outputvtu + ' -scale 0.20  --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ ScaledMesh
        subprocess.call(Scale, shell=True)

        stlMesh = Directory +'BottomHoneycomb.stl '
        command = 'meshio-convert  ' + ScaledMesh + ' ' +stlMesh
        subprocess.call(command, shell=True)

        Remove(Clipped1)
        Remove(Outputvtu)
        Remove(ScaledMesh)
        Remove(Outputstl)


def find_mins_maxs(obj):
    minx = obj.x.min()
    maxx = obj.x.max()
    miny = obj.y.min()
    maxy = obj.y.max()
    minz = obj.z.min()
    maxz = obj.z.max()
    return minx, maxx, miny, maxy, minz, maxz


def translate(_solid, step, padding, multiplier, axis):
    if 'x' == axis:
        items = 0, 3, 6
    elif 'y' == axis:
        items = 1, 4, 7
    elif 'z' == axis:
        items = 2, 5, 8
    else:
        raise RuntimeError('Unknown axis %r, expected x, y or z' % axis)

    # _solid.points.shape == [:, ((x, y, z), (x, y, z), (x, y, z))]
    _solid.points[:, items] += (step * multiplier) + (padding * multiplier)


def copy_obj(obj, dims, num_rows, num_cols, num_layers):
    w, l, h = dims
    copies = []
    for layer in range(num_layers):
        for row in range(num_rows):
            for col in range(num_cols):
                # skip the position where original being copied is
                if row == 0 and col == 0 and layer == 0:
                    continue
                _copy = mesh.Mesh(obj.data.copy())
                # translate(_copy, l, 0., row, 'y')
                # pad the space between objects by 10% of the dimension being
                # translated
                if col != 0:
                    translate(_copy, w, 0, col, 'x')
                if row != 0:
                    translate(_copy, l, 0., row, 'y')
                if layer != 0:
                    translate(_copy, h, 0, layer, 'z')
                copies.append(_copy)
    return copies


def CreateBottom():
    # Using an existing stl file:
    Stl1 ='/Volumes/Hardrive/Projects/FSI/NetworkDensity/Collectedstls/mesh_10.stl'
    Stl2 ='/Volumes/Hardrive/Projects/FSI/NetworkDensity/Collectedstls/mesh_10.stl'
    main_body = mesh.Mesh.from_file(Stl1)

    # rotate along Y
    # main_body.rotate([0.0, 0.5, 0.0], math.radians(90))

    minx, maxx, miny, maxy, minz, maxz = find_mins_maxs(main_body)
    w1 = maxx - minx
    l1 = maxy - miny
    h1 = maxz - minz
    # copies = copy_obj(main_body, (w1, l1, h1), 1, 1, 1) # THis will give how many copies of the mesh you want in each direction 

    # I wanted to add another related STL to the final STL
    twist_lock = mesh.Mesh.from_file(Stl2)
    BottomSection = mesh.Mesh.from_file(Stl2)
    minx, maxx, miny, maxy, minz, maxz = find_mins_maxs(twist_lock)
    w2 = maxx - minx
    l2 = maxy - miny
    h2 = maxz - minz

    # translate(_solid, step, padding, multiplier, axis):
    translate(twist_lock, w1-0.04, 0, 1, 'x')
    translate(BottomSection, -1*l2+0.08 , 0, 1, 'y')
    combined = mesh.Mesh(numpy.concatenate([main_body.data, twist_lock.data]))
                                                                   
    combined.save('/Volumes/Hardrive/Projects/FSI/NetworkDensity/Collectedstls/Bottom.stl', mode=stl.Mode.ASCII)  # save as ASCII



def ConcateMesh(Directory,NewDirectory, collapse, NumberOfExtraLevels, NumberOfCollums):

   
    Stl1 = Directory+'UpperSection'+str(collapse)+'.stl'
    Stl2 = Directory+'MiddelSection.stl'
    Stl3 = Directory+'BottomSection.stl'
    Stl4 = Directory+'UpperSection1.0.stl'
    
    main_body = mesh.Mesh.from_file(Stl1)
    MiddelSection = mesh.Mesh.from_file(Stl2)
    BottomSection = mesh.Mesh.from_file(Stl3)
    UpperSection = mesh.Mesh.from_file(Stl4)
    
    minx, maxx, miny, maxy, minz, maxz = find_mins_maxs(MiddelSection)
    w2 = maxx - minx
    l2 = maxy - miny
    h2 = maxz - minz

    # # translate
    translate(BottomSection, -NumberOfExtraLevels*l2 , 0, 1, 'y')
    translate(BottomSection, -w2+0.047 , 0, 1, 'x')
  
  
    # if NumberOfCollums>1:
    UpperCopies = copy_obj(UpperSection, (w2-0.047, -1*l2-1, h2), 1, NumberOfCollums, 1)
    BottomCopies = copy_obj(BottomSection, (w2-0.047, 0,0), 1, NumberOfCollums+1, 1)
    MiddelCopies = copy_obj(MiddelSection, (w2-0.047, -1*l2, h2), NumberOfExtraLevels+1, NumberOfCollums, 1)

    combined = mesh.Mesh(numpy.concatenate([main_body.data]+
                                [copy.data for copy in BottomCopies]+
                                [copy.data for copy in UpperCopies]+
                                [copy.data for copy in MiddelCopies]))
    combined.save(NewDirectory+ 'combined.stl', mode=stl.Mode.ASCII)  # save as ASCII
    # print 

def DensifyMesh(Directory, NewDirectory,collapse, NumberOfExtraLevels, NumberOfCollums):
    Stl1 = Directory+'UpperSection'+str(collapse)+'.stl'
    
    Stl2 = Directory+'MiddelSection.stl'
    Stl3 = Directory+'BottomSection.stl'
    Stl4 = Directory+'UpperSection1.0.stl'
    Stl4 = Directory+'TopHoneycomb.stl'
    Stl3 = Directory+'BottomHoneycomb.stl'
    Stl5 = Directory+'Honeycomb.stl'
    # Stl4 = Directory+'BottomHoneycomb.stl'
    

    main_body = mesh.Mesh.from_file(Stl1)
    MiddelSection = mesh.Mesh.from_file(Stl2)
    BottomSection = mesh.Mesh.from_file(Stl3)
    UpperSection = mesh.Mesh.from_file(Stl4)
    Honeycomb = mesh.Mesh.from_file(Stl5)

    minx, maxx, miny, maxy, minz, maxz = find_mins_maxs(MiddelSection)
    w2 = maxx - minx
    l2 = maxy - miny
    h2 = maxz - minz

    MiddelCopies = copy_obj(MiddelSection, (w2-0.047, -1*l2, h2), NumberOfExtraLevels+1, 1, 1)
    translate(BottomSection, -(NumberOfExtraLevels)*l2 , 0, 1, 'y')
    translate(BottomSection, (-w2+0.047)/2 , 0, 1, 'x')
    translate(UpperSection, (w2-0.047)/2, 0, 1, 'x')


    # if NumberOfCollums>1:
    UpperCopies = copy_obj(UpperSection, ( (w2-0.047)/2 , -1*l2-1, h2), 1, NumberOfCollums+1, 1)
    BottomCopies = copy_obj(BottomSection, ((w2-0.047)/2, 0,0), 1, NumberOfCollums+3, 1)


    # # translate
    minx, maxx, miny, maxy, minz, maxz = find_mins_maxs(Honeycomb)
    w2 = maxx - minx
    l2 = maxy - miny
    h2 = maxz - minz
    translate(Honeycomb, 2*w2-(0.047), 0, 1, 'x')
    HoneyCombCopies = copy_obj(Honeycomb, (w2-(0.047)/2, -1*l2,0), NumberOfExtraLevels+1, NumberOfCollums, 1)

    combined = mesh.Mesh(numpy.concatenate([main_body.data,Honeycomb.data]+
                                [copy.data for copy in BottomCopies]+
                                [copy.data for copy in UpperCopies]+
                                [copy.data for copy in HoneyCombCopies]+
                                [copy.data for copy in MiddelCopies]))
    combined.save(NewDirectory+ 'NewMeshcombined.stl', mode=stl.Mode.ASCII)  # save as ASCII


if __name__=="__main__":
    # Using an existing stl file:
    Remesh('/data/vascrem/MeshCollection/IdealisedNetwork/','ScalledMesh_NewNew.stl','RemeshedMatlabMesh.stl', 1)

    # def Remesh(Directory,Mesh, NewMesh, server):

    
    # First create a bunch of mesh bits to stick together
    # Directory = '/Volumes/Hardrive/Projects/FSI/NetworkDensity/Meshes/'
    # LevelDirectory = '/Volumes/Hardrive/Projects/FSI/NetworkDensity/Meshes/ConcatedMeshes/'

    # Directory = '/data/vascrem/MeshCollection/IdealisedNetwork/VascularDensity/'
    # if path.isdir('/data/vascrem/MeshCollection/IdealisedNetwork/VascularDensity/IncreasingMeshDensity/')==0:
    #     os.mkdir('/data/vascrem/MeshCollection/IdealisedNetwork/VascularDensity/IncreasingMeshDensity/')
    # List = [ '0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0']
    
    # # ScaleAndClippMeshes(Directory,List , 0, 0, 0, 0, 1, 1)
    
    # server=1
    # collapse = 0.4
    # NumberOfExtraLevels =2
    # NumberOfCollums =2
    
    # AdditionalLevels = [1,2,3,4,5]
    # List = [ '0.3']
    # for j in AdditionalLevels:
    #     for i in List:
    #         collapse = i
    #         LevelDirectory = Directory+'IncreasingMeshDensity/'+str(j)+'/'
    #         if path.isdir(LevelDirectory)==0:
    #             os.mkdir(LevelDirectory)
    #         DensifyMesh(Directory, LevelDirectory,collapse, j, j)
    #         NewMesh = 'mesh_'+str(collapse)+'.stl'
    #         Remesh(LevelDirectory,'combined.stl', NewMesh, server)





    # List = [ '0.3']
    # AdditionalLevels = [1,2,3,4,5]
    # for j in AdditionalLevels:
    #     for i in List:
    #         collapse = i
    #         LevelDirectory = Directory+'Levels_'+str(j+3)+'/'
    #         if path.isdir(LevelDirectory)==0:
    #             os.mkdir(LevelDirectory)
    #         NumberOfExtraLevels = j
    #         # ConcateMesh(Directory, LevelDirectory,collapse, NumberOfExtraLevels, 1)
    #         NewMesh = 'mesh_'+str(collapse)+'.stl'
    #         Remesh(LevelDirectory,'combined.stl', NewMesh, server)


