
#ifndef TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_CORRECTEDDRAG_HPP_
#define TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_CORRECTEDDRAG_HPP_

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <cstdio>
#include <ctime>
#include "AbstractCellBasedTestSuite.hpp"
#include "CellsGenerator.hpp"
#include "Debug.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"

#include "OffLatticeSimulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"

#include "MutableMesh.hpp"
#include "HistoryDepMutableMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "RemeshingModifier.hpp"
#include "VtkMeshReader.hpp"



class SurfaceAreaForceOnCylinder : public AbstractCellBasedTestSuite
{

public:

    void TestAreaFroceDragCorrectedEqui() throw(Exception)
    {
        clock_t t = clock();

        std::string working_directory = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/ShrunkPlexusWithLongInlets/" ;
		
        std::string mesh_file =  "/Users/jcrawshaw/Documents/testoutput/TestingRemeshingModifier/OriginalMesh.vtu";//working_directory +"SetUpData/TestOriginal_WithStretchedInlets.vtu" ;//"SetUpData/config.vtu";
		

        // Honeycomb3DCylinderMeshGenerator generator(20, 20, 1e-3,10e-3);
        // MutableMesh<2, 3>* p_mesh = generator.GetMesh();
        // HistoryDepMutableMesh<2, 3>* p_mesh2 = new HistoryDepMutableMesh<2, 3>;
        // p_mesh2 =dynamic_cast<HistoryDepMutableMesh<2, 3>*>(p_mesh);
 
        std::string output_directory = "TestingRemeshingModifier";
        // PRINT_VARIABLE(p_mesh2->GetNumNodes())
        
        VtkMeshReader<2,3> mesh_reader(mesh_file);
        HistoryDepMutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        PRINT_VARIABLE(mesh.GetNumNodes())
        // mesh.DeleteMesh();
      
        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2,3> cell_population(mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        

        // Start the remeshing
        double startime = 0;
        cell_population.SetChasteOutputDirectory(output_directory, startime );
        // cell_population.RemeshGeometry();
        cell_population.ExecuteHistoryDependentRemeshing();
        TRACE("In main code")


    // VtkMeshReader<2,3> mesh_reader2("/Users/jcrawshaw/Documents/testoutput/TestingRemeshingModifier/RemeshedGeometry.vtu");
    // HistoryDepMutableMesh<2,3> New_mesh2;
    // New_mesh2.ConstructFromMeshReader(mesh_reader2);
    // 
    // mesh.DeleteMesh();
    // mesh.AssignNewMesh(&New_mesh2);
    // New_mesh2.~HistoryDepMutableMesh();

    
//         // // Set up cell-based simulation
            OffLatticeSimulation<2, 3> simulator(cell_population);
            //         TRACE("Sent into off lattice sim")
            simulator.SetOutputDirectory(output_directory + "/WithremeshedMesh");


            simulator.SetEndTime(0.1); //(M_TIME_FOR_SIMULATION);
            simulator.SetDt(0.01);
            simulator.SetSamplingTimestepMultiple(1);
            simulator.SetUpdateCellPopulationRule(false); // No remeshing.
            TRACE("About to solve ")
            simulator.Solve();
            TRACE("Solved")

//         //       
//         // boost::shared_ptr<RemeshingModifier<2, 3> > p_Mesh_modifier(new RemeshingModifier<2, 3>());
//         // p_Mesh_modifier->SaveInitalConditions(cell_population);
//         // simulator.AddSimulationModifier(p_Mesh_modifier);

        


//     //     // p_Mesh_modifier->SetChasteOutputDirectory(output_directory, startime );


//     //     // // p_Mesh_modifier->RemeshGeometry();
//     //     // // p_Mesh_modifier->MappingAdaptedMeshToInitalGeometry(cell_population, *p_mesh);
//     //     // // std::map<unsigned, c_vector<double, 3> > RemeshedInitalConditions = p_Mesh_modifier->GetInitalConditions();

//     //     // // Need to think about what I do now with the cell populations 

    //     // // To reset before looping: this is usually done by the SetUp and TearDown methods
    //     // SimulationTime::Instance()->Destroy();
    //     // SimulationTime::Instance()->SetStartTime(0.0);
    //     t = clock() - t;
    //     std::cout<<"simulation time: "<<t/CLOCKS_PER_SEC<<" seconds"<<std::endl;

    }
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/
