#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

#include <cxxtest/TestSuite.h>
#include <cstdio>
#include <ctime>
#include <cmath>


#include "VtkMeshWriter.hpp"
#include "Debug.hpp"
#include "CommandLineArguments.hpp"
#include "PetscSetupAndFinalize.hpp"


#include "GenericMeshReader.hpp"
#include "VtkMeshReader.hpp"

#include <vtkXMLPolyDataReader.h>
#include <vtkPLYWriter.h>


#include <vtkAppendFilter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersion.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <vtkPolyDataReader.h>


#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
// #include <Python.h>
#include <iostream>
#include <array>

#include "AbstractCellBasedTestSuite.hpp"


class TestRunFlowInPipe : public AbstractCellBasedTestSuite
{


public:

	void TestMeshReaderAndWriter() throw (Exception)
	{


        std::string inputFileName = "mesh/test/data/CrapMesh.vtk";
        std::string outputFileName = "mesh/test/data/GeneatedVTU.vtu";

;
        // vtkSmartPointer<vtkPolyDataReader> Reader = vtkSmartPointer<vtkPolyDataReader>::New();
vtkPolyDataReader *vpdr = vtkPolyDataReader::New();

        // vtkUnstructuredGrid* reader = vtkSmartPointer<vtkPolyDataReader>::New();
        
        // vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
        // Reader->SetFileName(inputFileName.c_str());
        // Reader->Update();
        // // output = reader.GetOutput()

        // vtkSmartPointer<vtkPLYWriter> writer = vtkSmartPointer<vtkPLYWriter>::New();
        // writer->SetFileName(outputFileName.c_str());
        // writer->SetInputConnection(reader->GetOutputPort());
        // writer->Update();





//         TRACE("A")
//         vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
//         TRACE("B")
//         sphereSource->SetCenter(0.0, 0.0, 0.0);
//         sphereSource->SetRadius(5.0);
//         sphereSource->Update();

//         TRACE("C")
//         // Combine the two data sets
//         vtkSmartPointer<vtkAppendFilter> appendFilter =
//         vtkSmartPointer<vtkAppendFilter>::New();
//         TRACE("D")
//         #if VTK_MAJOR_VERSION <= 5
//         TRACE("ETRACE("A")")
//         appendFilter->AddInput(sphereSource->GetOutput());
//         #else
//         TRACE("E")
//         appendFilter->AddInputData(sphereSource->GetOutput());
//         #endif
//         appendFilter->Update();
//             TRACE("G")
//         vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
//     vtkSmartPointer<vtkUnstructuredGrid>::New();
//         unstructuredGrid->ShallowCopy(appendFilter->GetOutput());
// TRACE("H")
//         // Write the unstructured grid
//         vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
//         vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
//         TRACE("I")
//         writer->SetFileName("UnstructuredGrid.vtu");
//         #if VTK_MAJOR_VERSION <= 5
//         writer->SetInput(unstructuredGrid);
//         #else
//         writer->SetInputData(unstructuredGrid);
//         #endif
//         writer->Write();



        //  std::shared_ptr<AbstractMeshReader<2,3> > p_mesh_reader = GenericMeshReader<2,3>("mesh/test/data/CrapMesh.vtk");
            //We have to make a mesh so that we can get the node connectivity list back
            // DistributedTetrahedralMesh<2,3> mesh;
        //     mesh.ConstructFromMeshReader(*p_mesh_reader);

        // VtkMeshReader<2,3> mesh_reader("mesh/test/data/CrapMesh.vtk");
        // VtkMeshReader<2,2> vtk_reader("mesh/test/data/CrapMesh.vtk");

        // TetrahedralMesh<2,3> mesh;
        // mesh.ConstructFromMeshReader(mesh_reader);
        // VtkMeshWriter<2,3> writer("TestVtkMeshWriter", "CrapMesh", false);


	



		// VtkMeshWriter<2,3> mesh_writer(output_dir, "config", false);
     	// MutableMesh<2,3>* p_mesh = &(dynamic_cast<MeshBasedCellPopulation<2,3>*>(&(p_simulator->rGetCellPopulation()))->rGetMesh());
		// mesh_writer.WriteFilesUsingMesh(*p_mesh2);

	}

};

#endif /*TESTRELAXATION_HPP_*/


