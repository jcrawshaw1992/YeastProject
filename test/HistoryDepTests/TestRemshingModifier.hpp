
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
#include "MembraneBendingForce.hpp"

#include "FixedRegionBoundaryCondition.hpp"

#include "RemeshingTriggerModifier.hpp"

class SurfaceAreaForceOnCylinder : public AbstractCellBasedTestSuite
{

public:
void TestRemeshingOnV() throw(Exception)
    {
        clock_t t = clock();
        std::string mesh_file = "/Users/jcrawshaw/Documents/Projects/MeshMatlab/Meshes/V2D.vtu"; 
        std::string output_directory = "RemeshingComparison/V2D_SharpCorners";
        
        // Read in plexus
        VtkMeshReader<2,3> mesh_reader(mesh_file);
        HistoryDepMutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
 
        // // Create cells
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
        // cell_population.SetInitialAnlgesAcrossMembrane();
         cell_population.SaveInitalConditions();
        cell_population.SaveInitalConditionsAtTime0();
        cell_population.MarkBoundaryNodes();

        vmtksurfaceremeshing -ifile /Users/jcrawshaw/Documents/testoutput/RemeshingComparison/V2D_SharpCorners/CurrentMesh.stl -iterations  10   -area 0.0005 -maxarea 0.004 -ofile /Users/jcrawshaw/Documents/testoutput/RemeshingComparison/V2D_SharpCorners/Remeshed.stl
        cell_population.SetRelativePath(output_directory, startime);
        cell_population.SetTargetRemeshingElementArea(0.00005);
        cell_population.SetPrintRemeshedIC(1);
        cell_population.SetTargetRemeshingIterations(10);
        
        cell_population.ExecuteHistoryDependentRemeshing();

        OffLatticeSimulation<2, 3> simulator(cell_population);
        // simulator.SetOutputDirectory(output_directory);

        // simulator.SetEndTime(5);
        // simulator.SetDt(0.01);
        // simulator.SetSamplingTimestepMultiple(1);
        // simulator.SetUpdateCellPopulationRule(false);

        // std::vector<c_vector<double,3> > boundary_plane_points;
        // std::vector<c_vector<double,3> > boundary_plane_normals;

 
      
        // boost::shared_ptr<BoundariesModifier<2, 3> > p_Boundary_modifier(new BoundariesModifier<2, 3>());
        // p_Boundary_modifier->CreateBoundaryNodes(cell_population,boundary_plane_normals, boundary_plane_points);
        // p_Boundary_modifier->SetupSolve(cell_population,output_directory );
        // std::map<unsigned, c_vector<unsigned, 2> > NearestNodesMap = p_Boundary_modifier->GetNeighbouringNodesMap();


        // /*
        // -----------------------------
        // Compressive tissue pressure
        // ----------------------------
        // */        
        // double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        // double P_tissue = 0.001466542; // Pa == 1.1000e-05 mmHg

        // boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        // p_ForceOut->SetPressure((P_blood-P_tissue));
        // p_ForceOut->SetNearestNeighboursMap(NearestNodesMap);
        // simulator.AddForce(p_ForceOut);

        // /*
        // -----------------------------
        // Bending Force
        // ----------------------------
        // */
        // double membrane_constant = 1e-16;
        // boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        // p_membrane_force->SetMembraneStiffness(membrane_constant, 10, 10);       
        // simulator.AddForce(p_membrane_force);

        // /*       
        // -----------------------------
        // Boundaries 
        // ----------------------------
        // */

        // //Create a plane boundary to represent the inlet and pass them to the simulation
        // c_vector<long double, 3> Boundary1 = Create_c_vector(0, 0, 0.01);
        // c_vector<long double, 3> Normal1 = Create_c_vector(0, 0, -1);

        // c_vector<long double, 3> Boundary2 = Create_c_vector(0, 0, 0);
        // c_vector<long double, 3> Normal2 = Create_c_vector(0, 0, 1);
        
        // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary1, Normal1, 1));
        // simulator.AddCellPopulationBoundaryCondition(p_condition_1);

        // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary2, Normal2, 1));
        // simulator.AddCellPopulationBoundaryCondition(p_condition_2);

        // boost::shared_ptr<RemeshingTriggerModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerModifier<2, 3>());
        // p_Mesh_modifier->SetRemeshingInterval(30);
        // simulator.AddSimulationModifier(p_Mesh_modifier);


        // TRACE("About to solve ")
        // simulator.Solve();
        // TRACE("Solved")

        // SimulationTime::Instance()->Destroy();
        // SimulationTime::Instance()->SetStartTime(0.0);
        // t = clock() - t;
        // std::cout << "simulation time: " << t / CLOCKS_PER_SEC << " seconds" << std::endl;
    }

    void OffTestAreaFroceDragCorrectedEqui() throw(Exception)
    {
        clock_t t = clock();

        std::string working_directory = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/ShrunkPlexusWithLongInlets/";
        std::string mesh_file = "/Users/jcrawshaw/Documents/testoutput/TestingRemeshingModifier/OriginalMesh.vtu"; 
        
        std::string output_directory = "TestingChasteRemeshing/second";

        Honeycomb3DCylinderMeshGenerator generator(40, 40,10, 10);
        MutableMesh<2, 3>* p_mesh = generator.GetMesh();
        HistoryDepMutableMesh<2, 3>*  mesh= static_cast<HistoryDepMutableMesh<2, 3>*>(p_mesh);


        // // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh->GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(*mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);

        // Start the remeshing
        double startime = 0;
        cell_population.SetChasteOutputDirectory(output_directory, startime);
        
        cell_population.ExecuteHistoryDependentRemeshing();

        OffLatticeSimulation<2, 3> simulator(cell_population);
        // simulator.SetOutputDirectory(output_directory);

        // simulator.SetEndTime(5);
        // simulator.SetDt(0.01);
        // simulator.SetSamplingTimestepMultiple(1);
        // simulator.SetUpdateCellPopulationRule(false);

        // std::vector<c_vector<double,3> > boundary_plane_points;
        // std::vector<c_vector<double,3> > boundary_plane_normals;

 
      
        // boost::shared_ptr<BoundariesModifier<2, 3> > p_Boundary_modifier(new BoundariesModifier<2, 3>());
        // p_Boundary_modifier->CreateBoundaryNodes(cell_population,boundary_plane_normals, boundary_plane_points);
        // p_Boundary_modifier->SetupSolve(cell_population,output_directory );
        // std::map<unsigned, c_vector<unsigned, 2> > NearestNodesMap = p_Boundary_modifier->GetNeighbouringNodesMap();


        // /*
        // -----------------------------
        // Compressive tissue pressure
        // ----------------------------
        // */        
        // double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        // double P_tissue = 0.001466542; // Pa == 1.1000e-05 mmHg

        // boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        // p_ForceOut->SetPressure((P_blood-P_tissue));
        // p_ForceOut->SetNearestNeighboursMap(NearestNodesMap);
        // simulator.AddForce(p_ForceOut);

        // /*
        // -----------------------------
        // Bending Force
        // ----------------------------
        // */
        // double membrane_constant = 1e-16;
        // boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        // p_membrane_force->SetMembraneStiffness(membrane_constant, 10, 10);       
        // simulator.AddForce(p_membrane_force);

        // /*       
        // -----------------------------
        // Boundaries 
        // ----------------------------
        // */

        // //Create a plane boundary to represent the inlet and pass them to the simulation
        // c_vector<long double, 3> Boundary1 = Create_c_vector(0, 0, 0.01);
        // c_vector<long double, 3> Normal1 = Create_c_vector(0, 0, -1);

        // c_vector<long double, 3> Boundary2 = Create_c_vector(0, 0, 0);
        // c_vector<long double, 3> Normal2 = Create_c_vector(0, 0, 1);
        
        // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary1, Normal1, 1));
        // simulator.AddCellPopulationBoundaryCondition(p_condition_1);

        // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary2, Normal2, 1));
        // simulator.AddCellPopulationBoundaryCondition(p_condition_2);

        // boost::shared_ptr<RemeshingTriggerModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerModifier<2, 3>());
        // p_Mesh_modifier->SetRemeshingInterval(30);
        // simulator.AddSimulationModifier(p_Mesh_modifier);


        // TRACE("About to solve ")
        // simulator.Solve();
        // TRACE("Solved")

        // SimulationTime::Instance()->Destroy();
        // SimulationTime::Instance()->SetStartTime(0.0);
        // t = clock() - t;
        // std::cout << "simulation time: " << t / CLOCKS_PER_SEC << " seconds" << std::endl;
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
