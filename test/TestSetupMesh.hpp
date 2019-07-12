#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

//  You might like to fix up the boundaries here

#include <cxxtest/TestSuite.h>
#include <cstdio>
#include <ctime>
#include <cmath>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "HoneycombMeshGenerator.hpp"

#include "CylindricalHoneycombMeshGenerator.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "AppliedForceModifier.hpp"
#include "SpringLengthModifier.hpp"

#include "CellIdWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"

#include "CellsGenerator.hpp"
// #include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "AppliedForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "SmartPointers.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VtkMeshWriter.hpp"
#include "CommandLineArguments.hpp"
#include "MembraneStiffnessForce.hpp"
#include "PetscSetupAndFinalize.hpp"


#include "HemodynamicTraction.hpp"


#include "XmlTools.hpp"
using namespace xsd::cxx::tree;

class TestSetupFlowInPipe : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:

    void TestSetupRunAndSavePipe() throw (Exception)
    {


// I want a regular mesh 

        int N_D = 20;
        int N_Z = 1.5 * N_D;
        double mesh_scale = 1e-3;
        Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z,1.5e-3, 12e-3);
        MutableMesh<2,3>* p_mesh = generator.GetMesh();
        p_mesh->Translate(-6.0e-3*unit_vector<double>(3,2));
        // std::string output_dir = "/";
         std::string output_dir = "CylindricalMesh/";
        // TRACE(output_dir);
        VtkMeshWriter<2,3> mesh_writer(output_dir, "config", false);
        p_mesh->Scale(1.0/mesh_scale,1.0/mesh_scale,1.0/mesh_scale); // so distances are back in original scale
        mesh_writer.WriteFilesUsingMesh(*p_mesh);
    
    }

       void offTestSetupRunAndSavePipe15() throw (Exception)
    {


// I want a regular mesh 

        int N_D = 50;
        int N_Z = 1.5 * N_D;
        double mesh_scale = 1e-3;
        Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, 0.15e-3, 1e-3);
        MutableMesh<2,3>* p_mesh = generator.GetMesh();
        // p_mesh->Translate(-6.0e-3*unit_vector<double>(3,2));
        std::string output_dir = "LargeCylinder/";
        // std::string output_dir = "CylindricalMesh/";
        // TRACE(output_dir);
        VtkMeshWriter<2,3> mesh_writer(output_dir, "config", false);
        p_mesh->Scale(1.0/mesh_scale,1.0/mesh_scale,1.0/mesh_scale); // so distances are back in original scale
        mesh_writer.WriteFilesUsingMesh(*p_mesh);
    
    }
     void offTestSetupRunAndSavePipe2() throw (Exception)
    {


// I want a regular mesh 

        int N_D = 50;
        int N_Z = 1.5 * N_D;
        double mesh_scale = 1e-3;
        Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, 0.2e-3, 1e-3);
        MutableMesh<2,3>* p_mesh = generator.GetMesh();
        // p_mesh->Translate(-6.0e-3*unit_vector<double>(3,2));
        std::string output_dir = "Growing_2/";
        // std::string output_dir = "CylindricalMesh/";
        // TRACE(output_dir);
        VtkMeshWriter<2,3> mesh_writer(output_dir, "config", false);
        p_mesh->Scale(1.0/mesh_scale,1.0/mesh_scale,1.0/mesh_scale); // so distances are back in original scale
        mesh_writer.WriteFilesUsingMesh(*p_mesh);
    
    }

   


};

#endif /*TESTRELAXATION_HPP_*/

