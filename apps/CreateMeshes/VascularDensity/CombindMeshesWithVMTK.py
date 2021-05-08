import math
import stl
from stl import mesh
import numpy

import subprocess
import vtk
import os
import pdb
import string
import math as m
import vmtk
from vmtk import pypes
import array as A
# from vtk import vtkXMLUnstructuredGridWriter as Writer, VTK_VERSION
# from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
import numpy as np
import meshio
from vtk.util import numpy_support
import vtk
import clip
import math



# find the max dimensions, so we can know the bounding box, getting the height,
# width, length (because these are the step size)...
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
    twist_lock2 = mesh.Mesh.from_file(Stl2)
    minx, maxx, miny, maxy, minz, maxz = find_mins_maxs(twist_lock)
    w2 = maxx - minx
    l2 = maxy - miny
    h2 = maxz - minz

    # translate(_solid, step, padding, multiplier, axis):
    translate(twist_lock, w1-0.04, 0, 1, 'x')
    translate(twist_lock2, -1*l2+0.08 , 0, 1, 'y')
    combined = mesh.Mesh(numpy.concatenate([main_body.data, twist_lock.data]))
                                                                   
    combined.save('/Volumes/Hardrive/Projects/FSI/NetworkDensity/Collectedstls/Bottom.stl', mode=stl.Mode.ASCII)  # save as ASCII


if __name__=="__main__":
    # Using an existing stl file:


    # New to get two meshes to start with and translate one  -- think they will need to be vtks

    # Mesh1 = '/data/vascrem/MeshCollection/IdealisedNetwork/VascularDensity/mesh_0.vtk'
    # Mesh2 =  '/data/vascrem/MeshCollection/IdealisedNetwork/VascularDensity/mesh_10.vtk'
    # TranslatedMesh =  '/data/vascrem/MeshCollection/IdealisedNetwork/VascularDensity/Trans.stl'
    
    # MergedVTK = '/data/vascrem/MeshCollection/IdealisedNetwork/VascularDensity/Merged.vtp'
    # Stl1  = '/data/vascrem/MeshCollection/IdealisedNetwork/VascularDensity/mesh_0.stl'
    # Stl2 =  '/data/vascrem/MeshCollection/IdealisedNetwork/VascularDensity/mesh_10.stl'

    # Mesh1vtu = '/data/vascrem/MeshCollection/IdealisedNetwork/VascularDensity/mesh_0.vtu'
    # Mesh2vtu =  '/data/vascrem/MeshCollection/IdealisedNetwork/VascularDensity/mesh_10.vtu'


    # main_body = mesh.Mesh.from_file(Stl1)

    # # rotate along Y
    # # main_body.rotate([0.0, 0.5, 0.0], math.radians(90))

    # minx, maxx, miny, maxy, minz, maxz = find_mins_maxs(main_body)
    # w1 = maxx - minx
    # l1 = maxy - miny
    # h1 = maxz - minz
    # # copies = copy_obj(main_body, (w1, l1, h1), 1, 1, 1) # THis will give how many copies of the mesh you want in each direction 

    # # I wanted to add another related STL to the final STL
    # twist_lock = mesh.Mesh.from_file(Stl2)
    # twist_lock2 = mesh.Mesh.from_file(Stl2)
    # # twist_lock3 = mesh.Mesh.from_file(Stl3)
    # # twist_lock4 = mesh.Mesh.from_file(Stl3)
    # minx, maxx, miny, maxy, minz, maxz = find_mins_maxs(twist_lock)
    # w2 = maxx - minx
    # l2 = maxy - miny
    # h2 = maxz - minz

    # # # translate(_solid, step, padding, multiplier, axis):
    # translate(twist_lock, w1-0.04, 0, 1, 'x')
    # twist_lock.save(TranslatedMesh, mode=stl.Mode.ASCII)  # save as ASCII


    # TranslatedMesh =  '/data/vascrem/MeshCollection/IdealisedNetwork/VascularDensity/Trans.stl'
    # Translatedvtp = '/data/vascrem/MeshCollection/IdealisedNetwork/VascularDensity/Trans.vtp'

    # # MeshioConvert  = 'meshio-convert '+ TranslatedMesh + ' '+ Translatedvtp
    # subprocess.call(MeshioConvert, shell=True)

    # Mesh1vtp = '/data/vascrem/MeshCollection/IdealisedNetwork/VascularDensity/mesh_0.vtp'
    # MeshioConvert  = 'meshio-convert '+ Mesh1 + ' '+ Mesh1vtp
    # subprocess.call(MeshioConvert, shell=True)

    # vmtkmeshmerge = 'vmtkmeshmerge -mesh1 '+ Mesh1vtp + ' -mesh2 '+ Translatedvtp +' -ofile ' + MergedVTK
    # subprocess.call(vmtkmeshmerge, shell=True)
    TranslatedVTK = '/data/vascrem/MeshCollection/IdealisedNetwork/VascularDensity/mesh_10.stl'
    Combined =  '/data/vascrem/MeshCollection/IdealisedNetwork/VascularDensity/combined.stl'
    CombinedClean=  '/data/vascrem/MeshCollection/IdealisedNetwork/VascularDensity/AAcombinedClean.stl'
    boolian ='vmtksurfacebooleanoperation -operation "difference" -ifile '+  Combined +' -ofile '+ CombinedClean +' -i2file ' + TranslatedVTK
    subprocess.call(boolian, shell=True)


    # Transform = 'vmtksurfacereferencesystemtransform -ifile '+Mesh2+ ' -referenceorigin  1 2 3 -referencesystems  '+Mesh2+ ' -ofile '+ TranslatedMesh
    # subprocess.call(Transform, shell=True)
        #   ScaledMesh = AngleDirectory+"Mesh_Scaled"+str(i)+".vtu"
        # # # ---- Interpolate the points in the centerlines file, this will reduce the refinment needed in the centerline modeller -------------# 
        # Scale = 'vmtkmeshscaling -ifile '+ Outputvtu + ' -scale 0.20  --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ ScaledMesh
        # subprocess.call(Scale, shell=True)

        # ScaledMeshstl = AngleDirectory+"ScaledMesh."+str(i)+".stl"
        # convert = 'meshio-convert '+ ScaledMesh +'  '+ ScaledMeshstl
        # subprocess.call(convert, shell=True)
        # print "Done one"
  



    # translate(twist_lock2, -1*l2+0.08 , 0, 1, 'y')
    # translate(twist_lock3, -1*l2+0.08 , 0, 1, 'y')
    # translate(twist_lock4, -2*l2+2*0.08 , 0, 1, 'y')
    # copies2 = copy_obj(twist_lock, (w2, -1*l2+0.08, h2), 2, 1, 1)
    # combined = mesh.Mesh(numpy.concatenate([main_body.data, twist_lock.data,twist_lock3.data]))
    
                                        
                             

    # combined.save('/Volumes/Hardrive/Projects/FSI/NetworkDensity/Collectedstls/combined.stl', mode=stl.Mode.BINARY)  # save as ASCII
