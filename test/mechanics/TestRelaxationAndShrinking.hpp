#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

#include <cxxtest/TestSuite.h>

#include <cstdio>
#include <ctime>
#include <cmath>

#include "AppliedForceModifier.hpp"
#include "SpringLengthModifier.hpp"
#include "OffLatticeSimulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
// #include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "LinearSpringWithRestLengthDependentSpringConstantsForce.hpp"
#include "LinearSpringWithVariableSpringConstantsForce.hpp"
#include "MembraneStiffnessForce.hpp"
#include "AppliedForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "MeshBasedCellPopulation.hpp"
// #include "FixedRegionBoundaryCondition.hpp"
#include "SmartPointers.hpp"
#include "CellLabel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VtkMeshWriter.hpp"
#include "InitialConditionsGenerator.hpp"
#include "Debug.hpp"

#include "PetscSetupAndFinalize.hpp"



class TestRelaxation : public AbstractCellBasedTestSuite
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

//     void  noTestCylinderSteadyStateWorking() throw (Exception)
//     {
//         double applied_pressure = 5.0;
//         double spring_constant = 15.0;

//         VtkMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/cylinder_validation/GenerateCylinder/cyl_20_26.vtu"); // Also 40_52
//         //TrianglesMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/open_ended_box");

//         MutableMesh<2,3> mesh;
//         mesh.ConstructFromMeshReader(mesh_reader);

//         // Create cells
//         MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
//         std::vector<CellPtr> cells;
//         CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
//         cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

//         // Create a cell population
//         MeshBasedCellPopulation<2,3> cell_population(mesh, cells);
//         cell_population.SetWriteVtkAsPoints(true);
//         cell_population.SetOutputMeshInVtk(true);
//         cell_population.CalculateRestLengths();

//         // Stop the simulation loading the forces from file so providing manually
//         for (AbstractCellPopulation<2,3>::Iterator cell_iter = cell_population.Begin();
//                  cell_iter != cell_population.End();
//                  ++cell_iter)
//         {
//             double voronoi_cell_area = cell_population.GetVolumeOfCell(*cell_iter);

//             c_vector<double, 3> location = cell_population.GetLocationOfCellCentre(*cell_iter);
//             location[2] =0.0;
//             location /= norm_2(location);

//             cell_iter->GetCellData()->SetItem("applied_force_x", applied_pressure * voronoi_cell_area * location[0]);
//             cell_iter->GetCellData()->SetItem("applied_force_y", applied_pressure * voronoi_cell_area * location[1]);
//             cell_iter->GetCellData()->SetItem("applied_force_z", applied_pressure * voronoi_cell_area * location[2]);
//         }

//         // Calculate rest lengths so that in equilibrium with applied force
//         InitialConditionsGenerator generator;
//         generator.SetInitialRestLengths(cell_population, spring_constant);

//         // Set up cell-based simulation
//         OffLatticeSimulation<2,3> simulator(cell_population);
//         simulator.SetOutputDirectory("TestCylinderSteadyState");
//         simulator.SetEndTime(2.0); //20.0
//         simulator.SetSamplingTimestepMultiple(12);
//         simulator.SetUpdateCellPopulationRule(false);

//         // Create a force law and pass it to the simulation
//         boost::shared_ptr<GeneralisedLinearSpringForce<2,3> > p_force(new GeneralisedLinearSpringForce<2,3>());
//         p_force->SetMeinekeSpringStiffness(spring_constant);
//         simulator.AddForce(p_force);

//         // Create a force law to apply the saved applied forces
//         boost::shared_ptr<AppliedForce<2,3> > p_applied_force(new AppliedForce<2,3>());
//         simulator.AddForce(p_applied_force);

//         // Fix the end points
//         // Create a plane boundary and pass them to the simulation points less than x=-4.9 wont move
//         boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_boundary_condition_1(new FixedRegionBoundaryCondition<2,3>(&cell_population,-4.9*unit_vector<double>(3,2), -1.0*unit_vector<double>(3,2), 10.0));
//         simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_1);

//         simulator.Solve();

//         //Check nothings moved
//         for (AbstractCellPopulation<2,3>::Iterator cell_iter = cell_population.Begin();
//             cell_iter != cell_population.End();
//             ++cell_iter)
//         {
//             c_vector<double,3> cell_location = cell_population.GetLocationOfCellCentre(*cell_iter);
//             cell_location(2)=0.0;
//             TS_ASSERT_DELTA(norm_2(cell_location), 1.5, 1e-5);

//         }
//     }

//     void  noTestCylinderSteadyStateNonUniformMesh() throw (Exception)
//     {
//         double applied_pressure = 5.0;
//         double spring_constant = 15.0;

//         VtkMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/cylinder_validation/GenerateCylinder/cyl_581_nodes.vtu");
//         MutableMesh<2,3> mesh;
//         mesh.ConstructFromMeshReader(mesh_reader);

//         // Create cells
//         MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
//         std::vector<CellPtr> cells;
//         CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
//         cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

//         // Create a cell population
//         MeshBasedCellPopulation<2,3> cell_population(mesh, cells);
//         cell_population.SetWriteVtkAsPoints(true);
//         cell_population.SetOutputMeshInVtk(true);
//         cell_population.CalculateRestLengths();

//         // Stop the simulation loading the forces from file so providing manually
//         for (AbstractCellPopulation<2,3>::Iterator cell_iter = cell_population.Begin();
//                  cell_iter != cell_population.End();
//                  ++cell_iter)
//         {
//             double voronoi_cell_area = cell_population.GetVolumeOfCell(*cell_iter);

//             c_vector<double, 3> location = cell_population.GetLocationOfCellCentre(*cell_iter);
//             location[2] =0.0;
//             location /= norm_2(location);

//             cell_iter->GetCellData()->SetItem("applied_force_x", applied_pressure * voronoi_cell_area * location[0]);
//             cell_iter->GetCellData()->SetItem("applied_force_y", applied_pressure * voronoi_cell_area * location[1]);
//             cell_iter->GetCellData()->SetItem("applied_force_z", applied_pressure * voronoi_cell_area * location[2]);
//         }

//         // Calculate rest lengths so that in equilibrium with applied force
//         InitialConditionsGenerator generator;
//         generator.SetInitialRestLengths(cell_population, spring_constant);

//         // Set up cell-based simulation
//         OffLatticeSimulation<2,3> simulator(cell_population);
//         simulator.SetOutputDirectory("TestCylinderSteadyState");
//         simulator.SetEndTime(20.0); //20.0
//         simulator.SetSamplingTimestepMultiple(12);
//         simulator.SetUpdateCellPopulationRule(false);

//         // Create a force law and pass it to the simulation
//         boost::shared_ptr<GeneralisedLinearSpringForce<2,3> > p_force(new GeneralisedLinearSpringForce<2,3>());
//         p_force->SetMeinekeSpringStiffness(spring_constant);
//         simulator.AddForce(p_force);

//         // Create a force law to apply the saved applied forces
//         boost::shared_ptr<AppliedForce<2,3> > p_applied_force(new AppliedForce<2,3>());
//         simulator.AddForce(p_applied_force);

//         // Fix the end points
//         // Create a plane boundary and pass them to the simulation points less than x=-4.9 wont move
//         boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_boundary_condition_1(new FixedRegionBoundaryCondition<2,3>(&cell_population,-4.9*unit_vector<double>(3,2), -1.0*unit_vector<double>(3,2), 10.0));
//         simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_1);

//         simulator.Solve();

//         //Check nothings moved
//         for (AbstractCellPopulation<2,3>::Iterator cell_iter = cell_population.Begin();
//             cell_iter != cell_population.End();
//             ++cell_iter)
//         {
//             c_vector<double,3> cell_location = cell_population.GetLocationOfCellCentre(*cell_iter);
//             cell_location(2)=0.0;
//             TS_ASSERT_DELTA(norm_2(cell_location), 1.5, 0.2); //Low tolerance as mesh not stationairy however if you dont recaluclate rest lengths this fails.

//         }
//     }


//     void  noTestCylinderApproximateSteadyState() throw (Exception)
//     {
//         double applied_pressure = 5.0;
//         double spring_constant = 15.0;

//         VtkMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/cylinder_validation/GenerateCylinder/cyl_581_nodes.vtu"); // Also 40_52

//         //VtkMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/cylinder_validation/GenerateCylinder/cyl_20_26.vtu"); // Also 40_52
//         //TrianglesMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/open_ended_box");

//         MutableMesh<2,3> mesh;
//         mesh.ConstructFromMeshReader(mesh_reader);

//         // Create cells
//         MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
//         std::vector<CellPtr> cells;
//                 CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;    // Jess updated this line because this version of chaste doesnt seem to have the old cell generator, might, i didnt really check ----CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
//         cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

//         // Create a cell population
//         MeshBasedCellPopulation<2,3> cell_population(mesh, cells);
//         cell_population.SetWriteVtkAsPoints(true);
//         cell_population.SetOutputMeshInVtk(true);
//         cell_population.CalculateRestLengths();

//         // Stop the simulation loading the forces from file so providing manually
//         for (AbstractCellPopulation<2,3>::Iterator cell_iter = cell_population.Begin();
//                  cell_iter != cell_population.End();
//                  ++cell_iter)
//         {
//             double voronoi_cell_area = cell_population.GetVolumeOfCell(*cell_iter);

//             c_vector<double, 3> location = cell_population.GetLocationOfCellCentre(*cell_iter);
//             location[2] =0.0;
//             location /= norm_2(location);

//             cell_iter->GetCellData()->SetItem("applied_force_x", applied_pressure * voronoi_cell_area * location[0]);
//             cell_iter->GetCellData()->SetItem("applied_force_y", applied_pressure * voronoi_cell_area * location[1]);
//             cell_iter->GetCellData()->SetItem("applied_force_z", applied_pressure * voronoi_cell_area * location[2]);
//         }

//         // Calculate rest lengths so that in equilibrium with applied force
//         InitialConditionsGenerator generator;
//         generator.SetSimpleInitialRestLengths(cell_population, spring_constant);

//         // Set up cell-based simulation
//         OffLatticeSimulation<2,3> simulator(cell_population);
//         simulator.SetOutputDirectory("TestCylinderApproximateSteadyState");
//         simulator.SetEndTime(2.0); //20.0
//         simulator.SetSamplingTimestepMultiple(12);
//         simulator.SetUpdateCellPopulationRule(false);

//         // Create a force law and pass it to the simulation
//         boost::shared_ptr<GeneralisedLinearSpringForce<2,3> > p_force(new GeneralisedLinearSpringForce<2,3>());
//         p_force->SetMeinekeSpringStiffness(spring_constant);
//         simulator.AddForce(p_force);

//         // Create a force law to apply the saved applied forces
//         boost::shared_ptr<AppliedForce<2,3> > p_applied_force(new AppliedForce<2,3>());
//         simulator.AddForce(p_applied_force);

//         // Fix the end points
//         // Create a plane boundary and pass them to the simulation points less than x=-4.9 wont move
//         boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_boundary_condition_1(new FixedRegionBoundaryCondition<2,3>(&cell_population,-4.9*unit_vector<double>(3,2), -1.0*unit_vector<double>(3,2), 10.0));
//         simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_1);

//         // Fix the end points
//         // Create a plane boundary and pass them to the simulation points greater than x=4.9 wont move
//         boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_boundary_condition_2(new FixedRegionBoundaryCondition<2,3>(&cell_population,4.9*unit_vector<double>(3,2), 1.0*unit_vector<double>(3,2), 10.0));
//         simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_2);

//         simulator.Solve();

//         //Check nothings moved
//         for (AbstractCellPopulation<2,3>::Iterator cell_iter = cell_population.Begin();
//             cell_iter != cell_population.End();
//             ++cell_iter)
//         {
//             c_vector<double,3> cell_location = cell_population.GetLocationOfCellCentre(*cell_iter);
//             cell_location(2)=0.0;
//             TS_ASSERT_DELTA(norm_2(cell_location), 1.5, 0.06);

//         }
//     }

//     /*
//      * Cells with z < 0 have non zero wall shear stress
//      * z>0 has ~zero wall shear stress
//      */
//     void noTestCylinderShrinking() throw (Exception)
//     {
//         double applied_pressure = 1.0;
//         double spring_constant = 15.0;
//         double low_wss = 0.01;
//         double high_wss = 1.0;
//         double wss_threshold = 0.1;

//         VtkMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/cylinder_validation/GenerateCylinder/cyl_581_nodes.vtu"); // 581
//         MutableMesh<2,3> mesh;
//         mesh.ConstructFromMeshReader(mesh_reader);

//         // Create cells
//         MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
//         std::vector<CellPtr> cells;
//         CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
//         cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);



//         // Create a cell population
//         MeshBasedCellPopulation<2,3> cell_population(mesh, cells);
//         cell_population.SetWriteVtkAsPoints(true);
//         cell_population.SetOutputMeshInVtk(true);
//         cell_population.CalculateRestLengths();

//         // Stop the simulation loading the forces from file so providing manually
//         for (AbstractCellPopulation<2,3>::Iterator cell_iter = cell_population.Begin();
//                  cell_iter != cell_population.End();
//                  ++cell_iter)
//         {
//             double voronoi_cell_area = cell_population.GetVolumeOfCell(*cell_iter);

//             c_vector<double, 3> location = cell_population.GetLocationOfCellCentre(*cell_iter);
//             location[2] =0.0;
//             location /= norm_2(location);

//             cell_iter->GetCellData()->SetItem("applied_force_x", applied_pressure * voronoi_cell_area * location[0]);
//             cell_iter->GetCellData()->SetItem("applied_force_y", applied_pressure * voronoi_cell_area * location[1]);
//             cell_iter->GetCellData()->SetItem("applied_force_z", applied_pressure * voronoi_cell_area * location[2]);


//             cell_iter->GetCellData()->SetItem("applied_shear_stress_x", 0.0);
//             cell_iter->GetCellData()->SetItem("applied_shear_stress_y", 0.0);
//             cell_iter->GetCellData()->SetItem("applied_shear_stress_z", high_wss);
//             cell_iter->GetCellData()->SetItem("applied_shear_stress_mag", high_wss);

//             if(cell_population.GetLocationOfCellCentre(*cell_iter)[2]>0)
//             {
//                 cell_iter->GetCellData()->SetItem("applied_shear_stress_z", low_wss);
//                 cell_iter->GetCellData()->SetItem("applied_shear_stress_mag", low_wss);
//             }
//         }

//         // Calculate rest lengths so that in equilibrium with applied force
//         InitialConditionsGenerator generator;
//         generator.SetApproximateInitialRestLengths(cell_population, spring_constant);

//         // Set up cell-based simulation
//         OffLatticeSimulation<2,3> simulator(cell_population);
//         simulator.SetOutputDirectory("TestCylinderRelaxationWithShrinking");
//         simulator.SetEndTime(10.0); //20.0
//         simulator.SetSamplingTimestepMultiple(12);
//         simulator.SetUpdateCellPopulationRule(false);

//         boost::shared_ptr<SpringLengthModifier<2,3> > p_length_modifier(new SpringLengthModifier<2,3>());
//         p_length_modifier->SetShearThreshold(wss_threshold);
//         p_length_modifier->SetReductionFactor(0.001);
//         p_length_modifier->SetMaxDivisions(500);
//         simulator.AddSimulationModifier(p_length_modifier);

//         // Create a force law and pass it to the simulation
//         boost::shared_ptr<GeneralisedLinearSpringForce<2,3> > p_force(new GeneralisedLinearSpringForce<2,3>());
//         p_force->SetMeinekeSpringStiffness(spring_constant);
//         simulator.AddForce(p_force);

//         // Create a force law to apply the saved applied forces
//         boost::shared_ptr<AppliedForce<2,3> > p_applied_force(new AppliedForce<2,3>());
//         simulator.AddForce(p_applied_force);

//         // Fix the end points
//         // Create a plane boundary and pass them to the simulation points less than x=-4.9 wont move
// //        boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_boundary_condition_1(new FixedRegionBoundaryCondition<2,3>(&cell_population,-4.9*unit_vector<double>(3,2), -1.0*unit_vector<double>(3,2), 10.0));
// //        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_1);

//         // Fix the end points
//         // Create a plane boundary and pass them to the simulation points greater than x=4.9 wont move
// //        boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_boundary_condition_2(new FixedRegionBoundaryCondition<2,3>(&cell_population,4.9*unit_vector<double>(3,2), 1.0*unit_vector<double>(3,2), 10.0));
// //        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_2);

//         simulator.Solve();

//         //Check nothings moved
//         for (AbstractCellPopulation<2,3>::Iterator cell_iter = cell_population.Begin();
//             cell_iter != cell_population.End();
//             ++cell_iter)
//         {
//             c_vector<double,3> cell_location = cell_population.GetLocationOfCellCentre(*cell_iter);
//             cell_location(2)=0.0;
//             //TS_ASSERT_DELTA(norm_2(cell_location), 1.5, 0.06);

//         }
//     }

    /*
     * Loads the network mesh and pre calculated tractions
     */
    void TestExampleNetwork() throw (Exception)
    {
       double spring_constant = 1.0; //0.2
       double membrane_constant = 1e-12;
       double drag_coeficient = 1.0; //0.1
       double wss_threshold = 0.025;

       VtkMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/network_example/0min-mask-no-frame-iolet-extensions_corrected_tubed_openends_remeshed.vtu");
       MutableMesh<2,3> mesh;
       mesh.ConstructFromMeshReader(mesh_reader);
       mesh.Scale(1e-6, 1e-6, 1e-6);

       // Create cells
       MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
       std::vector<CellPtr> cells;
       CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
       cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

       // Create a cell population
       MeshBasedCellPopulation<2,3> cell_population(mesh, cells);
       cell_population.SetWriteVtkAsPoints(true);
       cell_population.SetOutputMeshInVtk(true);
       cell_population.CalculateRestLengths();
       cell_population.SetDampingConstantNormal(drag_coeficient);

       // Use an Applied Force modifier to read in the Flow details
       boost::shared_ptr<AppliedForceModifier<2,3> > p_force_modifier(new AppliedForceModifier<2,3>());
       p_force_modifier->SetResetTractionsOnCells(true,"projects/VascularRemodelling/test/data/network_example/surface-tractions.xtr");
       TRACE("Before Load");
       p_force_modifier->LoadTractionFromFile();
       TRACE("Before Cell Data");
       p_force_modifier->UpdateCellData(cell_population);
       TRACE("After Cell Data");

       // Set up cell-based simulation
       OffLatticeSimulation<2,3> simulator(cell_population);
       simulator.SetOutputDirectory("TestNetworkRelaxationWithShrinking");
       simulator.SetEndTime(10.0); //20.0
       simulator.SetSamplingTimestepMultiple(12);
       simulator.SetUpdateCellPopulationRule(false);

       boost::shared_ptr<SpringLengthModifier<2,3> > p_length_modifier(new SpringLengthModifier<2,3>());
       p_length_modifier->SetShearThreshold(wss_threshold);
       p_length_modifier->SetReductionFactor(0.01);
       p_length_modifier->SetMaxDivisions(90);
       simulator.AddSimulationModifier(p_length_modifier);

       // Create a force law and pass it to the simulation
       boost::shared_ptr<GeneralisedLinearSpringForce<2,3> > p_force(new GeneralisedLinearSpringForce<2,3>());
       p_force->SetMeinekeSpringStiffness(spring_constant);
       simulator.AddForce(p_force);

       boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
       p_membrane_force->SetupInitialMembrane(mesh);
       p_membrane_force->SetMembraneStiffness(membrane_constant);
       simulator.AddForce(p_membrane_force);

       simulator.Solve();
    }


//    void noTestSquare() throw (Exception)
//  {
//      // This data file is in mm
//      TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/square_in_3d");
//      MutableMesh<2,3> mesh;
//      mesh.ConstructFromMeshReader(mesh_reader);
//
//      // Create cells
//      MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
//      std::vector<CellPtr> cells;
//      CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
//      cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);
//
//      // Stop the simulation loading the forces from file so providing manually
//      for (unsigned i=0; i<cells.size(); i++)
//      {
//          cells[i]->GetCellData()->SetItem("applied_shear_stress_x", 0.0);
//          cells[i]->GetCellData()->SetItem("applied_shear_stress_y", 0.0);
//          cells[i]->GetCellData()->SetItem("applied_shear_stress_z", 0.5);
//          cells[i]->GetCellData()->SetItem("applied_shear_stress_mag", 0.5);
//      }
//
//
//      // Create a cell population
//      MeshBasedCellPopulation<2,3> cell_population(mesh, cells);
//      cell_population.SetWriteVtkAsPoints(true);
//      cell_population.SetOutputMeshInVtk(true);
//      cell_population.CalculateRestLengths();
//
//      // Set up cell-based simulation
//      AppliedForceOffLatticeSimulation<2,3> simulator(cell_population,"", 10.0,1.0);
//      simulator.SetOutputDirectory("TestSquareRelaxation");
//      simulator.SetEndTime(20.0); //20.0
//      simulator.SetSamplingTimestepMultiple(12);
//      simulator.SetUpdateCellPopulationRule(false);
//
//      // Stop the simulation loading the forces from file
//      simulator.SetResetTractionsOnCells(false,"");
//
//
//
//      // Create a force law and pass it to the simulation
//      boost::shared_ptr<LinearSpringWithRestLengthDependentSpringConstantsForce<2,3> > p_force(new LinearSpringWithRestLengthDependentSpringConstantsForce<2,3>());
//      p_force->SetCutOffLength(1.5);
//      simulator.AddForce(p_force);
//
//      // Create a force law to load in the applied forces
//      //boost::shared_ptr<AppliedForce<2,3> > p_applied_force(new AppliedForce<2,3>());
//      //simulator.AddForce(p_applied_force);
//
////        // Create two plane boundaries and pass them to the simulation
////        c_vector<double,3> point;
////        c_vector<double,3> normal;
////
////        point(0) = 0.001;
////        point(1) = 0.001;
////        point(2) = 0.0;
////
////        normal(0) = -1.0;
////        normal(1) = 0.0;
////        normal(2) = 0.0;
////
////        boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_fixed_condition_1(new FixedRegionBoundaryCondition<2,3>(&cell_population,point,normal));
////        simulator.AddCellPopulationBoundaryCondition(p_fixed_condition_1);
//
////        normal(0) = 0.0;
////        normal(1) = -1.0;
////
////        boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_fixed_condition_2(new FixedRegionBoundaryCondition<2,3>(&cell_population,point,normal));
////        simulator.AddCellPopulationBoundaryCondition(p_fixed_condition_2);
//
//
//      simulator.Solve();
//  }

//    void noTestPipe() throw (Exception)
//  {
//      // This data file is in mm
//      VtkMeshReader<2,3> mesh_reader("projects/ozzy/test/data/straight_vessel.vtu");
//      MutableMesh<2,3> mesh;
//      mesh.ConstructFromMeshReader(mesh_reader);
//      mesh.Scale(1e-3,1e-3,1e-3); // so distances are in m
//
//      // Create cells
//      MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
//      std::vector<CellPtr> cells;
//      CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
//      cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);
//
//      boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
//      //Select some cells to label
//      for (unsigned i = 0; i<mesh.GetNumNodes(); i++)
//      {
//          c_vector<double,3> position = mesh.GetNode(i)->rGetLocation();
//          c_vector<double,3> point;
//          point(0)=0.5*(-0.0288726- 0.0296406502875)+0.001;
//          point(1)=0.5*(-0.260683 - 0.256549081113)+0.0001;
//          point(2)=0.5*(-0.313772 - 0.314901550959);
//
//          double radius = 0.0012;
//
//          if (norm_2(point-position)<radius)
//          {
//              cells[i]->AddCellProperty(p_label);
//              //PRINT_VARIABLE(i);
//          }
//      }
//
//      // Create a cell population
//      MeshBasedCellPopulation<2,3> cell_population(mesh, cells);
//      cell_population.SetWriteVtkAsPoints(true);
//      cell_population.SetOutputMeshInVtk(true);
//      cell_population.SetOutputCellProliferativeTypes(true);
//      cell_population.SetOutputCellMutationStates(true);
//
//      cell_population.CalculateRestLengths();
//      cell_population.SetDampingConstantNormal(1e5); // to rescale parameters for length in m
//      cell_population.SetDampingConstantMutant(1e5); // to rescale parameters for length in m
//
//      // Set up cell-based simulation
//      AppliedForceOffLatticeSimulation<2,3> simulator(cell_population);
//      simulator.SetOutputDirectory("TestPipeRelaxation");
//      simulator.SetEndTime(1.0); //20.0
//      //simulator.SetDt(0.0001); //20.0
//      simulator.SetSamplingTimestepMultiple(12);
//      simulator.SetUpdateCellPopulationRule(false);
//
//
//      // Create a force law and pass it to the simulation
//      boost::shared_ptr<LinearSpringWithRestLengthDependentSpringConstantsForce<2,3> > p_force(new LinearSpringWithRestLengthDependentSpringConstantsForce<2,3>());
//      //p_force->SetCutOffLength(1.5);
//      //p_force->SetMeinekeSpringStiffness(15e5*1e-5); // to rescale parameters for length in m and to balance with applied traction
//      p_force->SetMeinekeSpringStiffness(80e5*1e-5); // to rescale parameters for length in m and to balance with applied traction
//      simulator.AddForce(p_force);
//
//        // Create a force law to load in the applied forces
//        boost::shared_ptr<AppliedForce<2,3> > p_applied_force(new AppliedForce<2,3>());
//        simulator.AddForce(p_applied_force);
//
////        // Create a plane boundary and pass them to the simulation
//        c_vector<double,3> point;
//        c_vector<double,3> normal;
//
////        // This inlet position is taken from the Hemelb input XML file
////        point(0) = 0.0226;
////        point(1) = 0.0582;
////        point(2) = 0.0624;
////
////        normal(0) = 0.3465;
////        normal(1) = -0.4535;
////        normal(2) = 0.8211;
//
//        point(0)=-0.0288726;
//        point(1)=-0.260683;
//        point(2)=-0.313772;
//
//        normal(0)=-0.232271;
//        normal(1)= 0.95232;
//        normal(2)=-0.197833;
//
////
//      boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_inlet_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population,point,normal));
//      simulator.AddCellPopulationBoundaryCondition(p_inlet_condition);
//
//
////        point(0) = 0.0182;
////        point(1) = 0.0598;
////        point(2) = 0.0552;
////
////        normal(0) = -0.7705;
////        normal(1) = -0.2612;
////        normal(2) = -0.5813;
//
//      point(0) = -0.0296406502875;
//      point(1) = -0.256549081113;
//      point(2) = -0.314901550959;
//
//      normal(0) = 0.250133733542;
//      normal(1) = -0.949421892861;
//      normal(2) = 0.189818820722;
//
//      boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_outlet_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population,point,normal));
//      simulator.AddCellPopulationBoundaryCondition(p_outlet_condition);
//
//      TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(),1148u);
//      TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumNodes(), 1148u);
//      simulator.Solve();
//      TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumNodes(),2053u);
//      TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(),2053u);
//
//      VtkMeshWriter<2,3> mesh_writer("TestPipeRelaxation", "end_mesh", false);
//      mesh_writer.WriteFilesUsingMesh(dynamic_cast<MeshBasedCellPopulation<2,3>*>(&(simulator.rGetCellPopulation()))->rGetMesh());
//
//      TRACE("Before Save");
//      CellBasedSimulationArchiver<2,AppliedForceOffLatticeSimulation<2,3>, 3>::Save(&simulator);
//      TRACE("After Save");
//
//  }

//  void noTestPipeSaveAndLoad() throw (Exception)
//  {
//      // 0-1 hours from above simulation
//      {
//          // This data file is in mm
//          VtkMeshReader<2,3> mesh_reader("projects/ozzy/test/data/pipe.vtu");
//          MutableMesh<2,3> mesh;
//          mesh.ConstructFromMeshReader(mesh_reader);
//          mesh.Scale(1e-3,1e-3,1e-3); // so distances are in m
//
//          // Create cells
//          MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
//          std::vector<CellPtr> cells;
//          CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
//          cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);
//
//          // Create a cell population
//          MeshBasedCellPopulation<2,3> cell_population(mesh, cells);
//          cell_population.SetWriteVtkAsPoints(true);
//          cell_population.SetOutputMeshInVtk(true);
//          cell_population.CalculateRestLengths();
//          cell_population.SetDampingConstantNormal(1e5); // to rescale parameters for length in m
//
//          // Set up cell-based simulation
//          AppliedForceOffLatticeSimulation<2,3> simulator(cell_population,"");
//          simulator.SetOutputDirectory("TestPipeRelaxationSaveAndLoad");
//          simulator.SetEndTime(1.0); //20.0
//          simulator.SetSamplingTimestepMultiple(12);
//          simulator.SetUpdateCellPopulationRule(false);
//
//
//          // Create a force law and pass it to the simulation
//          boost::shared_ptr<LinearSpringWithRestLengthDependentSpringConstantsForce<2,3> > p_force(new LinearSpringWithRestLengthDependentSpringConstantsForce<2,3>());
//          //p_force->SetCutOffLength(1.5);
//          p_force->SetMeinekeSpringStiffness(15e5*1e-5); // to rescale parameters for length in m and to balance with applied traction
//          simulator.AddForce(p_force);
//
//          // Create a force law to load in the applied forces
//          boost::shared_ptr<AppliedForce<2,3> > p_applied_force(new AppliedForce<2,3>());
//          simulator.AddForce(p_applied_force);
//
//          // Create a plane boundary and pass them to the simulation
//          c_vector<double,3> point;
//          c_vector<double,3> normal;
//
//          // This inlet position is taken from the Hemelb input XML file
//          point(0) = 0.0226;
//          point(1) = 0.0582;
//          point(2) = 0.0624;
//
//          normal(0) = 0.3465;
//          normal(1) = -0.4535;
//          normal(2) = 0.8211;
//
//          boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_inlet_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population,point,normal));
//          simulator.AddCellPopulationBoundaryCondition(p_inlet_condition);
//
//          // This outlet position is taken from the Hemelb input XML file
//      //      point(0) = 0.0182;
//      //      point(1) = 0.0598;
//      //      point(2) = 0.0552;
//      //
//      //      normal(0) = -0.7705;
//      //      normal(1) = -0.2612;
//      //      normal(2) = -0.5813;
//      //
//      //      boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_outlet_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population,point,normal));
//      //      simulator.AddCellPopulationBoundaryCondition(p_outlet_condition);
//
//          simulator.Solve();
//
//          //Save simulation
//          TRACE("Before Save");
//          CellBasedSimulationArchiver<2,AppliedForceOffLatticeSimulation<2,3>, 3>::Save(&simulator);
//          TRACE("After Save");
//
//      }
//      // 0-1 hours from above simulation
//      {
//          TRACE("Before Load");
//          AppliedForceOffLatticeSimulation<2,3>* p_simulator = CellBasedSimulationArchiver<2, AppliedForceOffLatticeSimulation<2,3>, 3 >::Load("TestPipeRelaxationSaveAndLoad", 1.0);
//          TRACE("After Load");
//          p_simulator->SetEndTime(2.0);
//          p_simulator->SetResetTractionsOnCells(false);
//          p_simulator->Solve();
//
//          delete p_simulator;
//      }
//
//  }

};

#endif /*TESTRELAXATION_HPP_*/

