
#include "Exception.hpp"
#include <boost/filesystem.hpp>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning
#include <vtkXMLUnstructuredGridWriter.h>
  #include <vtkTriangle.h>
  #include <vtkPoints.h>
  #include <vtkTetra.h>
  #include <vtkGeometryFilter.h>
  #include <vtkTriangleFilter.h>
  #include <vtkSTLWriter.h>
  #include <vtkCleanPolyData.h>
  #include <vtkIdList.h>
  #include <vtkCellArray.h>
  #include <vtkDataSetTriangleFilter.h>
 #include <vtkExtractUnstructuredGrid.h>
 #include <vtkDataSetSurfaceFilter.h>
 //#include <dolfin.h>
 #include "MultiFormatMeshWriter.hpp"
 
 template<unsigned DIM>
 MultiFormatMeshWriter<DIM>::MultiFormatMeshWriter()
     : mpVtkMesh(),
       mpMesh(),
       mFilepath(),
          mOutputFormat(MeshFormat::VTU)
    {
    
    }
    
    template<unsigned DIM>
    boost::shared_ptr<MultiFormatMeshWriter<DIM> > MultiFormatMeshWriter<DIM>::Create()
    {
        MAKE_PTR(MultiFormatMeshWriter<DIM>, pSelf);
        return pSelf;
  }
  
  template<unsigned DIM>
  MultiFormatMeshWriter<DIM>::~MultiFormatMeshWriter()
  {
  
  }
  
  template<unsigned DIM>
  void MultiFormatMeshWriter<DIM>::SetMesh(vtkSmartPointer<vtkUnstructuredGrid> pMesh)
    {
        mpVtkMesh = pMesh;
    }
    
    template<unsigned DIM>
    void MultiFormatMeshWriter<DIM>::SetMesh(boost::shared_ptr<DiscreteContinuumMesh<DIM> > pMesh)
    {
        mpMesh = pMesh;
    }
    
  template<unsigned DIM>
  void MultiFormatMeshWriter<DIM>::SetFilename(const std::string& filename)
  {
      mFilepath = filename;
  }
  
  template<unsigned DIM>
  void MultiFormatMeshWriter<DIM>::SetOutputFormat(MeshFormat::Value outputFormat)
  {
      mOutputFormat = outputFormat;
 }
 
 template<unsigned DIM>
 void MultiFormatMeshWriter<DIM>::Write()
 {
     if(mFilepath == "")
     {
         EXCEPTION("Output file not specified for mesh writer");
     }
 
     if(mOutputFormat == MeshFormat::VTU or mOutputFormat == MeshFormat::STL)
     {
         // If there is a DiscreteContinuum mesh convert it to vtk format first
         if(mpMesh)
         {
             vtkSmartPointer<vtkPoints> p_vtk_points = vtkSmartPointer<vtkPoints>::New();
             std::vector<c_vector<double, DIM> > node_locations = mpMesh->GetNodeLocations();
             p_vtk_points->SetNumberOfPoints(node_locations.size());
             for(unsigned idx=0; idx<node_locations.size(); idx++)
             {
                 if(DIM==3)
                 {
                     p_vtk_points->InsertPoint(idx, node_locations[idx][0], node_locations[idx][1], node_locations[idx][2]);
                 }
                 else
                 {
                     p_vtk_points->InsertPoint(idx, node_locations[idx][0], node_locations[idx][1], 0.0);
                 }
             }
             mpVtkMesh->SetPoints(p_vtk_points);
 
             // Add vtk tets or triangles
             std::vector<std::vector<unsigned> > element_connectivity =  mpMesh->GetConnectivity();
             unsigned num_elements = element_connectivity.size();
             mpVtkMesh->Allocate(num_elements, num_elements);
 
             for(unsigned idx=0; idx<num_elements; idx++)
             {
                 if(DIM==3)
                 {
                     vtkSmartPointer<vtkTetra> p_vtk_element = vtkSmartPointer<vtkTetra>::New();
                     unsigned num_nodes = element_connectivity[idx].size();
                     for(unsigned jdx=0; jdx<num_nodes; jdx++)
                     {
                         p_vtk_element->GetPointIds()->SetId(jdx, element_connectivity[idx][jdx]);
                     }
                     mpVtkMesh->InsertNextCell(p_vtk_element->GetCellType(), p_vtk_element->GetPointIds());
                 }
                 else
                 {
                     vtkSmartPointer<vtkTriangle> p_vtk_element = vtkSmartPointer<vtkTriangle>::New();
                     unsigned num_nodes = element_connectivity[idx].size();
                     for(unsigned jdx=0; jdx<num_nodes; jdx++)
                     {
                         p_vtk_element->GetPointIds()->SetId(jdx, element_connectivity[idx][jdx]);
                     }
                     mpVtkMesh->InsertNextCell(p_vtk_element->GetCellType(), p_vtk_element->GetPointIds());
                 }
             }
         }
 
         if(!mpVtkMesh)
         {
             EXCEPTION("No mesh has been set, cannot do write.");
         }
 
         if(mOutputFormat == MeshFormat::VTU)
         {
             vtkSmartPointer<vtkXMLUnstructuredGridWriter> p_writer1 = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
             p_writer1->SetFileName((mFilepath + ".vtu").c_str());
             #if VTK_MAJOR_VERSION <= 5
                 p_writer1->SetInput(mpVtkMesh);
             #else
                 p_writer1->SetInputData(mpVtkMesh);
             #endif
             p_writer1->Write();
         }
         else
         {
             vtkSmartPointer<vtkGeometryFilter> p_geom_filter = vtkSmartPointer<vtkGeometryFilter>::New();
             #if VTK_MAJOR_VERSION <= 5
                 p_geom_filter->SetInput(mpVtkMesh);
             #else
                 p_geom_filter->SetInputData(mpVtkMesh);
             #endif
 
             p_geom_filter->Update();
 
             vtkSmartPointer<vtkTriangleFilter> p_tri_filter = vtkSmartPointer<vtkTriangleFilter>::New();
             p_tri_filter->SetInputConnection(p_geom_filter->GetOutputPort());
 
             vtkSmartPointer<vtkCleanPolyData> p_clean_filter = vtkSmartPointer<vtkCleanPolyData>::New();
             p_clean_filter->SetInputConnection(p_tri_filter->GetOutputPort());
             p_clean_filter->Update();
 
             vtkSmartPointer<vtkSTLWriter> p_writer1 = vtkSmartPointer<vtkSTLWriter>::New();
             p_writer1->SetFileName((mFilepath + ".stl").c_str());
 
             #if VTK_MAJOR_VERSION <= 5
                 p_writer1->SetInput(p_clean_filter->GetOutput());
             #else
                 p_writer1->SetInputData(p_clean_filter->GetOutput());
             #endif
             p_writer1->SetFileTypeToASCII();
             p_writer1->Write();
         }
     }
 //    else if (mOutputFormat == MeshFormat::DOLFIN)
 //    {
 //        // Use DiscreteContinuum mesh directly
 //        if(mpMesh)
 //        {
 //            dolfin::MeshEditor editor;
 //            dolfin::Mesh dolfin_mesh;
 //            editor.open(dolfin_mesh, DIM, DIM);
 //
 //            std::vector<std::vector<double> > node_locations = mpMesh->GetNodeLocations();
 //            std::vector<std::vector<unsigned> > element_connectivity =  mpMesh->GetConnectivity();
 //            editor.init_vertices(node_locations.size());
 //            editor.init_cells(element_connectivity.size());
 //
 //            for(unsigned idx=0; idx<node_locations.size(); idx++)
 //            {
 //                if(DIM==2)
 //                {
 //                    editor.add_vertex(idx, node_locations[idx][0],
 //                                      node_locations[idx][1]);
 //                }
 //                else
 //                {
 //                    editor.add_vertex(idx, node_locations[idx][0],
 //                                      node_locations[idx][1],
 //                                      node_locations[idx][2]);
 //                }
 //
 //            }
 //
 //            for(unsigned idx=0; idx<element_connectivity.size(); idx++)
 //            {
 //                if(DIM==2)
 //                {
 //                    editor.add_cell(idx, element_connectivity[idx][0],
 //                                    element_connectivity[idx][1],
 //                                    element_connectivity[idx][2]);
 //                }
 //                else
 //                {
 //                    editor.add_cell(idx, element_connectivity[idx][0],
 //                                    element_connectivity[idx][1],
 //                                    element_connectivity[idx][2],
 //                                    element_connectivity[idx][3]);
 //                }
 //            }
 //
 //            editor.close();
 //            dolfin_mesh.init();
 //
 //            // Write the mesh to file
 //            dolfin::File mesh_file (mFilepath + ".xml");
 //            mesh_file << (dolfin_mesh);
 //        }
 //
  //    }
  }
  
  // Explicit instantiation
  template class MultiFormatMeshWriter<2>;
  template class MultiFormatMeshWriter<3>;