
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

#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"

#include "HistoryDepMutableMesh.hpp"
#include "MutableMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "RemeshingModifier.hpp"
#include "VtkMeshReader.hpp"

#include "OutwardsPressure.hpp"

class SurfaceAreaForceOnCylinder : public AbstractCellBasedTestSuite
{

public:
    void TestAreaFroceDragCorrectedEqui() throw(Exception)
    {
        clock_t t = clock();

        std::string working_directory = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/ShrunkPlexusWithLongInlets/";
        std::string mesh_file = "/Users/jcrawshaw/Documents/testoutput/TestingRemeshingModifier/OriginalMesh.vtu"; //working_directory +"SetUpData/TestOriginal_WithStretchedInlets.vtu" ;//"SetUpData/config.vtu";


        std::string output_directory = "TestingRemeshingModifier";

        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        HistoryDepMutableMesh<2, 3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        
        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);

        // Start the remeshing
        double startime = 0;
        cell_population.SetChasteOutputDirectory(output_directory, startime);
        cell_population.ExecuteHistoryDependentRemeshing();
        TRACE("In main code")


        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory + "/WithremeshedMesh");
        simulator.SetEndTime(0.1);
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetUpdateCellPopulationRule(false);

        /*
        -----------------------------
        Constant Pressure ballance 
        ----------------------------
        */
        double P_blood = 0.2133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.1466542; // Pa == 1.1000e-05 mmHg

        double TransmuralPressure = P_blood - P_tissue;

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure(TransmuralPressure);
        simulator.AddForce(p_ForceOut);

        TRACE("About to solve ")
        simulator.Solve();
        TRACE("Solved")

        t = clock() - t;
        std::cout<<"simulation time: "<<t/CLOCKS_PER_SEC<<" seconds"<<std::endl;
    }
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/



        // VtkMeshReader<2,3> mesh_reader2("/Users/jcrawshaw/Documents/testoutput/TestingRemeshingModifier/RemeshedGeometry.vtu");
        // HistoryDepMutableMesh<2,3> New_mesh2;
        // New_mesh2.ConstructFromMeshReader(mesh_reader2);
        //
        // mesh.DeleteMesh();
        // mesh.AssignNewMesh(&New_mesh2);
        // New_mesh2.~HistoryDepMutableMesh();


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