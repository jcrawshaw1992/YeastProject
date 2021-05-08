#!/usr/bin/env python
#
import os, sys
import clip

import subprocess
import vtk
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
import math
import ConvertVTKtoSTL

# from vtk.numpy_interface import dataset_adapter as dsa



# ---------------------------------
  
def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")
    


if __name__=="__main__":
    
    Density = ['1','2','3','4','5']
    List = [ '0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0']
    for j in Density:
        Directory = '/data/vascrem/MeshCollection/IdealisedNetwork/VascularDensity/IncreasingMeshDensity/'+str(j)+'/'
        for i in List:
            STLFile= Directory+'mesh_'+str(i)+'.stl'
            ScalledSTLFile= Directory+'ScalledMesh_'+str(i)+'.stl'
            VTUFile= Directory+'mesh.vtu'
            command = 'meshio-convert  ' + STLFile + ' ' +VTUFile
            # subprocess.call(command, shell=True)
            ScalledVTU = Directory+'ScalledMesh_'+str(i)+'.vtu'
            # Scalling section 
            # if j == '1':
            #     Scale = 'vmtkmeshscaling -ifile '+ VTUFile + ' -scale 0.77 --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ ScalledVTU
            #     subprocess.call(Scale, shell=True)
            #     # 1.465526907489212
            # if j =='2':
            #     Scale = 'vmtkmeshscaling -ifile '+ VTUFile + ' -scale 0.625 --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ ScalledVTU
            #     subprocess.call(Scale, shell=True)
            #     # 1.4494820904073549
            # if j =='3':
            #     Scale = 'vmtkmeshscaling -ifile '+ VTUFile + ' -scale 0.52 --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ ScalledVTU
            #     subprocess.call(Scale, shell=True)
            #     # 1.4378934675422443
            # if j =='4':
            #     Scale = 'vmtkmeshscaling -ifile '+ VTUFile + ' -scale 0.45 --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ ScalledVTU
            #     subprocess.call(Scale, shell=True)
            #     # 1.4285466576576455
            # if j =='5':
            #     Scale = 'vmtkmeshscaling -ifile '+ VTUFile + ' -scale 0.397 --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ ScalledVTU
            #     subprocess.call(Scale, shell=True)
            #     # 1.4285466576576455


            # Translation section 

    List = [ '0','1','2','3','4','5','6','7','8','9','1']
    for j in Density:
        Directory = '/data/vascrem/MeshCollection/IdealisedNetwork/VascularDensity/IncreasingMeshDensity/'+str(j)+'/'
        for i in List:
            
            VTKFIle = Directory+'ScalledMesh_'+i+'.vtk'
            Clipped_Mesh =  Directory+'ClippedMesh_'+i+'.vtk'
            # if j == '1':
            #     clip.clip_surface_with_plane(VTKFIle,(1.465526907489212,0,0), (-1,0,0), Clipped_Mesh)
            # if j =='2':
            #     clip.clip_surface_with_plane(VTKFIle,(1.4494820904073549,0,0), (-1,0,0), Clipped_Mesh)
            #     # 1.4494820904073549
            # if j =='3':
            #     clip.clip_surface_with_plane(VTKFIle,(1.4378934675422443,0,0), (-1,0,0), Clipped_Mesh)
            #     # 1.4378934675422443
            # if j =='4':
            #     clip.clip_surface_with_plane(VTKFIle,(1.4285466576576455,0,0), (-1,0,0), Clipped_Mesh)
            #     # 1.4285466576576455
            # if j =='5':
            #     clip.clip_surface_with_plane(VTKFIle,(1.4285466576576455,0,0), (-1,0,0), Clipped_Mesh)
            #     # 1.4285466576576455

            Meshstl =  Directory+'mesh.'+i+'.stl'
            ConvertVTKtoSTL.convert(Clipped_Mesh, Meshstl)



            # command = 'meshio-convert  ' + ScalledVTU + ' ' +ScalledSTLFile
            # subprocess.call(command, shell=True)

            #### 4 
            # Scale = 'vmtkmeshscaling -ifile '+ VTKFile + ' -scale 0.45 --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ VTKFile
            # subprocess.call(Scale, shell=True)

            # # #### 3
            # Scale = 'vmtkmeshscaling -ifile '+ VTKFile + ' -scale 0.52 --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ VTKFile
            # subprocess.call(Scale, shell=True)

            # #### 2
            # Scale = 'vmtkmeshscaling -ifile '+ VTKFile + ' -scale 0.625 --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ VTKFile
            # subprocess.call(Scale, shell=True)

            #### 1
            # 


            #### 5
            # Scale = 'vmtkmeshscaling -ifile '+ VTKFile + ' -scale 0.397 --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ VTKFile
            # subprocess.call(Scale, shell=True)




    # vtkMesh = Directory +'Mesh.vtk'
    # Clipped = Directory +'NewClippedMesh_.vtk'

    # command = 'meshio-convert  ' + STLFile + ' ' +vtkMesh
    # subprocess.call(command, shell=True)

    # command = 'meshio-info '  +vtkMesh
    # subprocess.call(command, shell=True)

    # clip.clip_surface_with_plane(vtkMesh,(1.48,0,0), (-1,0,0), Clipped)


    # 1.48
    
    # Clipped_Mesh =Directory+'Clippedcombined.vtk'
    # ----  Clip Edges -------------# 
    print "Jess is good"
    # print Clipped_Mesh
    # clip.clip_surface_with_plane(VTKFile,(1.1452452423654245,0,0), (1,0,0), Clipped_Mesh)
    # clip.clip_surface_with_plane(Clipped_Mesh,(1.8001460600894323,0,0), (-1,0,0), Clipped_Mesh)
    # clip.clip_surface_with_plane(Clipped_Mesh,(1.7985872167796197,-0.7,-0.06), (0,1,0), Clipped_Mesh)

    # print 'about to convert'
    # # ----  Convert files vtk to stl :)   -------------#    

    # Outputstl = Directory+'Clipped.stl'
    # ConvertVTKtoSTL.convert(Clipped_Mesh, Outputstl)


    # Outputvtu = AngleDirectory+"Mesh_PreScale"+str(i)+".vtu"
    # command = 'meshio-convert  ' + Outputstl + ' ' +Outputvtu
    # subprocess.call(command, shell=True)

    # ScaledMesh = AngleDirectory+"Mesh_Scaled"+str(i)+".vtu"
    # # # ---- Interpolate the points in the centerlines file, this will reduce the refinment needed in the centerline modeller -------------# 
        # Scale = 'vmtkmeshscaling -ifile '+ Outputvtu + ' -scale 0.20  --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ ScaledMesh
        # subprocess.call(Scale, shell=True)

        # ScaledMeshstl = AngleDirectory+"ScaledMesh."+str(i)+".stl"
        # convert = 'meshio-convert '+ ScaledMesh +'  '+ ScaledMeshstl
        # subprocess.call(convert, shell=True)
        # print "Done one"
    


        # print 'Scale'
        # # ----  Scale files  -------------# 
        # Scale = 'vmtkmeshscaling -ifile '+ Outputvtu + ' -scale 0.20  --pipe vmtkmeshwriter -entityidsarray CellEntityIds -ofile '+ ScaledMesh
        # subprocess.call(Scale, shell=True)

        # # ----  Convert files  -------------# 
        # convert = 'meshio-convert '+ ScaledMesh +'  '+ ScaledMeshstl
        # subprocess.call(convert, shell=True)


# print " Finished CLipping  "