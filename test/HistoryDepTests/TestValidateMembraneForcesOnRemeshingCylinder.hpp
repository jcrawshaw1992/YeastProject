#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>
#include "RadialForce.hpp"
// #include <cstdio>
// #include <ctime>
// #include <cmath>
// #include <vector> 

#include "Debug.hpp"
#include "VtkMeshReader.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "SmartPointers.hpp"
 
#include "AbstractCellBasedTestSuite.hpp"
#include "OffLatticeSimulation.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "HistoryDepMutableMesh.hpp"

// #include "AppliedForceModifier.hpp"
#include "RemeshingTriggerOnHeteroMeshModifier.hpp"
#include "FixedRegionBoundaryConditionWithRemeshing.hpp"
#include "MembraneDeformationForce.hpp"
// #include "OutwardsPressureWithBreaks.hpp"
#include "OutwardsPressure.hpp"
#include "MembraneStiffnessForce.hpp"
#include "HemeLBForce.hpp"


class TestRemeshing  : public AbstractCellBasedTestSuite
{
public:

 void TestMembraneForcesOnCylinderWithRemeshing() throw (Exception)
    {
        double startime =0;
        double EndTime = 120;    

        std::string output_dir = "TestRemeshing2";
       
        double scale = 1e3;
        double Length = 10e-6 * scale/2;
        double Radius = 1e-6 * scale;

        Honeycomb3DCylinderMeshGenerator generator(40, 30, Radius, Length);
        MutableMesh<2, 3>* p_mesh = generator.GetMesh();
        HistoryDepMutableMesh<2, 3>* mesh= static_cast<HistoryDepMutableMesh<2, 3>*>(p_mesh);
    
       // Create the cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh->GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(*mesh, cells);
        cell_population.SetChasteOutputDirectory(output_dir, startime);
        // cell_population.SetInitialAnlgesAcrossMembrane(); // Dont worry about this for now, I think there is something moff
        cell_population.SetRelativePath(output_dir, startime);
        cell_population.SetTargetRemeshingEdgeLength(0.2e-6 * scale); 
        // cell_population.EdgeLengthVariable(1.2); 
        cell_population.SetPrintRemeshedIC(0);
        cell_population.SetTargetRemeshingIterations(3);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        cell_population.SetRemeshingSoftwear("VMTK");
        cell_population.SetOperatingSystem("mac");
        // Set population to output all data to results files
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2,3> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetDt(0.005);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(EndTime);

        /*
        -----------------------------
        RemeshingTriggerOnHeteroMeshModifier
        ----------------------------
        */  
        boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnHeteroMeshModifier<2, 3>());
        // p_Mesh_modifier->SetMembraneStrength(3);
        std::map<double, c_vector<long double, 4> > GrowthMaps;
        GrowthMaps[1] = Create_c_vector(pow(10, -7), pow(10, -7), pow(10, -10), 0);
        //Strength , hetro, stepsize, setupsolve
        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 0, 100, 1);

        p_Mesh_modifier->SetRemeshingInterval(100);// I have turned this off because I need to know what will happen without remeshing, and then with remeshing
        simulator.AddSimulationModifier(p_Mesh_modifier);

        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */
       
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        boost::shared_ptr<RadialForce> p_ForceOut(new RadialForce());
        p_ForceOut->SetPressure((P_blood - P_tissue));// needs to be negative for server ?? 
        simulator.AddForce(p_ForceOut);

        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        simulator.AddForce(p_shear_force);




        /*
        -----------------------------
        Boundary conditions
        ----------------------------
        */


        c_vector<double, 3> PlaneNormal1 = Create_c_vector(0, 0, 1);
        c_vector<double, 3> Point1 = Create_c_vector(0, 0, 0.001);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(0, 0, -1);

        c_vector<double, 3> Point2 = Create_c_vector(0, 0, 0.004);

        std::vector<c_vector<double, 3> > boundary_plane_points;
        std::vector<c_vector<double, 3> > boundary_plane_normals;

        boundary_plane_points.push_back(Point1);
        boundary_plane_normals.push_back(PlaneNormal1);

        boundary_plane_points.push_back(Point2);
        boundary_plane_normals.push_back(PlaneNormal2);

        for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
            boost::shared_ptr<FixedRegionBoundaryConditionWithRemeshing<2, 3> > p_condition(new FixedRegionBoundaryConditionWithRemeshing<2, 3>(&cell_population, boundary_plane_points[boundary_id], -boundary_plane_normals[boundary_id], 0.005));
            simulator.AddCellPopulationBoundaryCondition(p_condition);
        }

        simulator.Solve();
        CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(&simulator);
}


//   void offTestMembraneForcesOnCylinderWithoutRemeshing() throw (Exception)
//     {
//         double startime =0;
//         double EndTime = 20;    
//         std::string output_dir = "TestMembraneForcesWithRemeshing/Cylinder/WithoutRemeshing/";
       
//         double scale = 1e3;
//         double Length = 50e-6 * scale;
//         double Radius = 1e-6 * scale;

//         Honeycomb3DCylinderMeshGenerator generator(40, 100, Radius, Length);
//         MutableMesh<2, 3>* p_mesh = generator.GetMesh();
//         HistoryDepMutableMesh<2, 3>* mesh= static_cast<HistoryDepMutableMesh<2, 3>*>(p_mesh);
    
//         // Create the cells 
//         MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
//         std::vector<CellPtr> cells;
//         CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
//         cells_generator.GenerateBasicRandom(cells, mesh->GetNumNodes(), p_differentiated_type);

//         // Create a cell population
//         HistoryDepMeshBasedCellPopulation<2, 3> cell_population(*mesh, cells);
//         cell_population.SetChasteOutputDirectory(output_dir, startime);
//         // cell_population.SetInitialAnlgesAcrossMembrane(); // Dont worry about this for now, I think there is something moff
//         cell_population.SetRelativePath(output_dir, startime);
//         cell_population.SetTargetRemeshingEdgeLength(0.65* scale); 
//         cell_population.EdgeLengthVariable(1.2); 
//         cell_population.SetPrintRemeshedIC(0);
//         cell_population.SetTargetRemeshingIterations(2);
//         cell_population.SetWriteVtkAsPoints(true);
//         cell_population.SetOutputMeshInVtk(true);
//         cell_population.SetRemeshingSoftwear("CGAL");
//         // Set population to output all data to results files
//         cell_population.AddCellWriter<CellProliferativeTypesWriter>();

//         // Set up cell-based simulation
//         OffLatticeSimulation<2,3> simulator(cell_population);
//         simulator.SetOutputDirectory(output_dir);
//         simulator.SetSamplingTimestepMultiple(400);
//         simulator.SetDt(0.005);
//         simulator.SetUpdateCellPopulationRule(false);
//         simulator.SetEndTime(EndTime);

//         /*
//         -----------------------------
//         RemeshingTriggerOnHeteroMeshModifier
//         ----------------------------
//         */  
//         boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnHeteroMeshModifier<2, 3>());
//         // p_Mesh_modifier->SetMembraneStrength(3);
//         std::map<double, c_vector<long double, 4> > GrowthMaps;
//         GrowthMaps[1] = Create_c_vector(pow(10, -7), pow(10, -7), pow(10, -10), 0);
//         //Strength , hetro, stepsize, setupsolve
//         p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 0, 100, 1);

//         // p_Mesh_modifier->SetRemeshingInterval(900);// I have turned this off because I need to know what will happen without remeshing, and then with remeshing
//         simulator.AddSimulationModifier(p_Mesh_modifier);

//         /*
//         -----------------------------
//         Constant Compressive tissue pressure
//         ----------------------------
//         */
       
//         double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
//         double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

//         boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
//         p_ForceOut->SetPressure((P_blood - P_tissue));// needs to be negative for server ?? 
//         simulator.AddForce(p_ForceOut);

//         /*
//         -----------------------------
//         Membrane forces
//         ----------------------------
//         */
//         boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
//         simulator.AddForce(p_shear_force);




//         /*
//         -----------------------------
//         Boundary conditions
//         ----------------------------
//         */


//         c_vector<double, 3> PlaneNormal1 = Create_c_vector(0, 0, 1);
//         c_vector<double, 3> Point1 = Create_c_vector(0, 0, 1e-6 * scale);

//         c_vector<double, 3> PlaneNormal2 = Create_c_vector(0, 0, -1);

//         c_vector<double, 3> Point2 = Create_c_vector(0, 0, 49e-6 * scale);

//         std::vector<c_vector<double, 3> > boundary_plane_points;
//         std::vector<c_vector<double, 3> > boundary_plane_normals;

//         boundary_plane_points.push_back(Point1);
//         boundary_plane_normals.push_back(PlaneNormal1);

//         boundary_plane_points.push_back(Point2);
//         boundary_plane_normals.push_back(PlaneNormal2);

//         for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
//         {
//             boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(&cell_population, boundary_plane_points[boundary_id], -boundary_plane_normals[boundary_id], 0.5));
//             simulator.AddCellPopulationBoundaryCondition(p_condition);
//         }

//         simulator.Solve();
//         CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(&simulator);
// }




};




#endif /*TESTRELAXATION_HPP_*/



// %%%% Cylinder
//   void TestMembraneForcesOnCylinderWithoutRemeshing() throw (Exception)
//     {
//         double startime =0;
//         double EndTime = 20;    
//         std::string output_dir = "TestRemeshing/";
       
//         double scale = 1e3;
//         double Length = 50e-6 * scale;
//         double Radius = 1e-6 * scale;

//         Honeycomb3DCylinderMeshGenerator generator(40, 100, Radius, Length);
//         MutableMesh<2, 3>* p_mesh = generator.GetMesh();
//         HistoryDepMutableMesh<2, 3>* mesh= static_cast<HistoryDepMutableMesh<2, 3>*>(p_mesh);
    
//         // Create the cells 
//         MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
//         std::vector<CellPtr> cells;
//         CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
//         cells_generator.GenerateBasicRandom(cells, mesh->GetNumNodes(), p_differentiated_type);

//         // Create a cell population
//         HistoryDepMeshBasedCellPopulation<2, 3> cell_population(*mesh, cells);
//         cell_population.SetChasteOutputDirectory(output_dir, startime);
//         // cell_population.SetInitialAnlgesAcrossMembrane(); // Dont worry about this for now, I think there is something moff
//         cell_population.SetRelativePath(output_dir, startime);
//         cell_population.SetTargetRemeshingEdgeLength(0.65* scale); 
//         cell_population.EdgeLengthVariable(1.2); 
//         cell_population.SetPrintRemeshedIC(0);
//         cell_population.SetTargetRemeshingIterations(2);
//         cell_population.SetWriteVtkAsPoints(true);
//         cell_population.SetOutputMeshInVtk(true);
//         cell_population.SetRemeshingSoftwear("CGAL");
//         // Set population to output all data to results files
//         cell_population.AddCellWriter<CellProliferativeTypesWriter>();

//         // Set up cell-based simulation
//         OffLatticeSimulation<2,3> simulator(cell_population);
//         simulator.SetOutputDirectory(output_dir);
//         simulator.SetSamplingTimestepMultiple(400);
//         simulator.SetDt(0.005);
//         simulator.SetUpdateCellPopulationRule(false);
//         simulator.SetEndTime(EndTime);

//         /*
//         -----------------------------
//         RemeshingTriggerOnHeteroMeshModifier
//         ----------------------------
//         */  
//         boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnHeteroMeshModifier<2, 3>());
//         // p_Mesh_modifier->SetMembraneStrength(3);
//         std::map<double, c_vector<long double, 4> > GrowthMaps;
//         GrowthMaps[1] = Create_c_vector(pow(10, -7), pow(10, -7), pow(10, -10), 0);
//         //Strength , hetro, stepsize, setupsolve
//         p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 0, 100, 1);

//         // p_Mesh_modifier->SetRemeshingInterval(900);// I have turned this off because I need to know what will happen without remeshing, and then with remeshing
//         simulator.AddSimulationModifier(p_Mesh_modifier);

//         /*
//         -----------------------------
//         Constant Compressive tissue pressure
//         ----------------------------
//         */
       
//         double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
//         double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

//         boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
//         p_ForceOut->SetPressure((P_blood - P_tissue));// needs to be negative for server ?? 
//         simulator.AddForce(p_ForceOut);

//         /*
//         -----------------------------
//         Membrane forces
//         ----------------------------
//         */
//         boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
//         simulator.AddForce(p_shear_force);




//         /*
//         -----------------------------
//         Boundary conditions
//         ----------------------------
//         */


//         c_vector<double, 3> PlaneNormal1 = Create_c_vector(0, 0, 1);
//         c_vector<double, 3> Point1 = Create_c_vector(0, 0, 1e-6 * scale);

//         c_vector<double, 3> PlaneNormal2 = Create_c_vector(0, 0, -1);

//         c_vector<double, 3> Point2 = Create_c_vector(0, 0, 49e-6 * scale);

//         std::vector<c_vector<double, 3> > boundary_plane_points;
//         std::vector<c_vector<double, 3> > boundary_plane_normals;

//         boundary_plane_points.push_back(Point1);
//         boundary_plane_normals.push_back(PlaneNormal1);

//         boundary_plane_points.push_back(Point2);
//         boundary_plane_normals.push_back(PlaneNormal2);

//         for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
//         {
//             boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(&cell_population, boundary_plane_points[boundary_id], -boundary_plane_normals[boundary_id], 0.5));
//             simulator.AddCellPopulationBoundaryCondition(p_condition);
//         }

//         simulator.Solve();
//         CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(&simulator);
// }


