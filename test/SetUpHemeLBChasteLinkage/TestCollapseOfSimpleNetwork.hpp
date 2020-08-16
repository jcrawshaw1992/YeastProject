#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

//  You might like to fix up the boundaries here

// 1FriendshopBetrayalIs
// jhj.maclean@gmail.com

#include <cxxtest/TestSuite.h>
#include <cstdio>
#include <ctime>
#include <cmath>

#include "Debug.hpp"

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "OffLatticeSimulation.hpp"

#include "CellIdWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"

#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"

#include "AbstractCellBasedTestSuite.hpp"
// #include "MeshBasedCellPopulation.hpp"

#include "FixedRegionBoundaryCondition.hpp"
#include "SmartPointers.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VtkMeshWriter.hpp"
// #include "CommandLineArguments.hpp"
// #include "MembraneStiffnessForce.hpp"
#include "PetscSetupAndFinalize.hpp"

// #include "PressureForce.hpp"
#include "XmlTools.hpp"

#include "UblasCustomFunctions.hpp"

#include "VtkMeshReader.hpp"

// #include "CellMutationStatesWriter.hpp"
#include "OutwardsPressure.hpp"

#include "MembraneDeformationForce.hpp"

#include "RemeshingTriggerModifier.hpp"
#include "RemeshingTriggerOnHeteroMeshModifier.hpp"
#include "HistoryDepMutableMesh.hpp"
#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "ElementQualityOutputModifier.hpp"
#include "HemeLBForce.hpp"




using namespace xsd::cxx::tree;

class TestRemeshing  : public AbstractCellBasedTestSuite
{
public:



void TestDeformationGrowthTOEqui() throw (Exception)
   {        
        // // /Users/jcrawshaw/docker-polnet-master/IdealiseNetworks/3_by_4/
        // std::string mesh_file =  "/Users/jcrawshaw/docker-polnet-master/IdealiseNetworks/3_by_4/IdealNetwork.vtu";//"projects/VascularRemodelling/test/data/bifurcation_cut/config.vtu";//
        // std::string output_dir = "TestFlowThroughSimpleNetwork/GrowthToEquilibrium/";
        
        // double mesh_scale = 1e-1;  //1e-3; 
        // double startime = 0;
       
        // // Read in the plexus 
        // VtkMeshReader<2,3> mesh_reader(mesh_file);
        // HistoryDepMutableMesh<2,3> mesh;
        // mesh.ConstructFromMeshReader(mesh_reader);
        // mesh.Scale(mesh_scale,mesh_scale,mesh_scale); // so distances are in mm
      
        // MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        // std::vector<CellPtr> cells;
        // CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        // cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // // Create a cell population
        // HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
        // cell_population.SetChasteOutputDirectory(output_dir, startime);
        // cell_population.SetInitialAnlgesAcrossMembrane();
        // cell_population.SetRelativePath(output_dir, startime);
        // // cell_population.SetTargetRemeshingEdgeLength(1e-3);
        // cell_population.SetPrintRemeshedIC(1);
        // // cell_population.SetTargetRemeshingIterations(5);
        // cell_population.SetWriteVtkAsPoints(true);
        // cell_population.SetOutputMeshInVtk(true);
        // // Set population to output all data to results files
        // cell_population.AddCellWriter<CellProliferativeTypesWriter>();

        // // Set up cell-based simulation
        // OffLatticeSimulation<2,3> simulator(cell_population);
        // simulator.SetOutputDirectory(output_dir);
        // simulator.SetSamplingTimestepMultiple(100);
        // simulator.SetDt(0.001);
        // simulator.SetUpdateCellPopulationRule(false);
        // simulator.SetEndTime(100);

        // /*
        //   -----------------------------
        //   RemeshingTriggerOnHeteroMeshModifier
        //   ----------------------------
        // */  

        // boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnHeteroMeshModifier<2, 3>());
        // simulator.AddSimulationModifier(p_Mesh_modifier);

        // /*
        // -----------------------------
        // Constant Force
        // ----------------------------
        // */        

        // // double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        // // double P_tissue = 0.001466542; // Pa == 1.1000e-05 mmHg

        // boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        // p_ForceOut->SetPressure(0.00035);
        // simulator.AddForce(p_ForceOut);

        // /*
        // -----------------------------
        // SMembrane forces
        // ----------------------------
        // */

        // boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        // simulator.AddForce(p_shear_force);

        // // Create a plane boundary to represent the inlets/outlets and pass them to the simulation
        // std::vector<c_vector<double,3> > boundary_plane_points;
        // std::vector<c_vector<double,3> > boundary_plane_normals;

        // // TRACE("These BCs need fixing, they no longer cover the plexus -- but this is not the current problem ")
        // boundary_plane_points.push_back(Create_c_vector(0.5,0,0));
        // boundary_plane_normals.push_back(Create_c_vector (1,0,0));

        // boundary_plane_points.push_back(Create_c_vector(0.5,1.4,0));
        // boundary_plane_normals.push_back(Create_c_vector (1,0,0));

        // boundary_plane_points.push_back(Create_c_vector(0.5,-1.4,0));
        // boundary_plane_normals.push_back(Create_c_vector (1,0,0));


        // boundary_plane_points.push_back(Create_c_vector(7.4,0,0));
        // boundary_plane_normals.push_back(Create_c_vector (-1,0,0));

        // boundary_plane_points.push_back(Create_c_vector(7.4,1.4,0));
        // boundary_plane_normals.push_back(Create_c_vector (-1,0,0));

        // boundary_plane_points.push_back(Create_c_vector(7.4,-1.4,0));
        // boundary_plane_normals.push_back(Create_c_vector (-1,0,0));

        // for(unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        // {
        //     boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population, boundary_plane_points[boundary_id],-boundary_plane_normals[boundary_id],0.05));
        //     simulator.AddCellPopulationBoundaryCondition(p_condition);
        // }
     	// simulator.Solve();
        // CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);

    }


// /Users/jcrawshaw/petsc-3.7.7/include/petscsys.h:152:6: error: "PETSc was


// void OffTestDeformationOfHoneyCombNetworkHetroAfterEqui() throw (Exception)
//    {        
//         std::string mesh_file = "/Users/jcrawshaw/docker-polnet-master/IdealiseNetworks/Mesh.vtu";
//         std::string output_dir = "FlowThroughCollapsingNetwork/";
        
//         double mesh_scale = 1e-1;  //1e-3; 
//         double startime = 0;
       
//         // Read in the plexus 
//         VtkMeshReader<2,3> mesh_reader(mesh_file);
//         HistoryDepMutableMesh<2,3> mesh;
//         mesh.ConstructFromMeshReader(mesh_reader);
//         mesh.Scale(mesh_scale,mesh_scale,mesh_scale); // so distances are in mm
//       // Create the cells 
      
//         MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
//         std::vector<CellPtr> cells;
//         CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
//         cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

//         // Create a cell population
//         HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
//         cell_population.SetChasteOutputDirectory(output_dir, startime);
//         cell_population.SetInitialAnlgesAcrossMembrane();
//         cell_population.SetRelativePath(output_dir, startime);
//         cell_population.SetTargetRemeshingEdgeLength(1e-2);
//         cell_population.SetPrintRemeshedIC(1);
//         cell_population.SetTargetRemeshingIterations(2);
//         cell_population.SetWriteVtkAsPoints(true);
//         cell_population.SetOutputMeshInVtk(true);
//         // Set population to output all data to results files
//         cell_population.AddCellWriter<CellProliferativeTypesWriter>();

//         // Set up cell-based simulation
//         OffLatticeSimulation<2,3> simulator(cell_population);
//         simulator.SetOutputDirectory(output_dir);
//         simulator.SetSamplingTimestepMultiple(100);
//         simulator.SetDt(0.01);
//         simulator.SetUpdateCellPopulationRule(false);
//         simulator.SetEndTime(60);

//         /*
//           -----------------------------
//           RemeshingTriggerOnHeteroMeshModifier
//           ----------------------------
//         */  
//         boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnHeteroMeshModifier<2, 3>());
//         // p_Mesh_modifier->SetRemeshingInterval(100); //(1000);
// // 
//           std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
//           //         KA,     Kalpha           Ks                                 Kb
//           GrowthMaps[2] = Create_c_vector(0, pow(10, -6.8124), pow(10, -7), 0);
//           GrowthMaps[1] = Create_c_vector(0, pow(10, -4), pow(10, -5), 0);

          
//           c_vector<double, 3> UpperPlanePoint = Create_c_vector(0.32,0,0);
//           c_vector<double, 3> LowerPlanePoint = Create_c_vector(0.36,0,0);

//           c_vector<double, 3> UpperPlaneNormal = Create_c_vector(1,0,0);
//           c_vector<double, 3> LowerPlaneNormal = Create_c_vector(-1,0,0);

//           p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);
//           p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 2, 1,100, 1); 
//           simulator.AddSimulationModifier(p_Mesh_modifier);

  
//         // boost::shared_ptr<ElementQualityOutputModifier<2, 3> > p_ElementWriter_modifier(new ElementQualityOutputModifier<2, 3>());
//         // // p_ElementWriter_modifier->SetElementMetricsFileName("/Users/jcrawshaw/Documents/testoutput/"+output_dir+"/ElementMetricsFile.txt");
//         // p_ElementWriter_modifier->SetElementMetricsFileName("/Users/jcrawshaw/Documents/testoutput/"+output_dir);
//         // p_ElementWriter_modifier->SetWriteoutInterval(200);
//         // simulator.AddSimulationModifier(p_ElementWriter_modifier);

//         /*
//         -----------------------------
//         HemeLB Force
//         ----------------------------
//         */        
//         c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0);
//         c_vector<double, 3> Point1 = Create_c_vector(0.1*mesh_scale,0,0);

//         c_vector<double, 3> PlaneNormal2 = Create_c_vector(-1,0,0);
//         c_vector<double, 3> Point2 = Create_c_vector(6.7*mesh_scale,0,0);

//         boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
//         p_ForceOut->Inlets(PlaneNormal1, Point1, 0.00025, "Inlet");
//         p_ForceOut->Inlets(PlaneNormal2, Point2, 0.00005, "Outlet");
//         p_ForceOut->SetUpHemeLBConfiguration(output_dir, mesh, cell_population);
//         simulator.AddForce(p_ForceOut);



//         // double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
//         // double P_tissue = 0.001466542; // Pa == 1.1000e-05 mmHg

//         // boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
//         // p_ForceOut->SetPressure((P_blood-P_tissue));
//         // simulator.AddForce(p_ForceOut);
//         /*
//         -----------------------------
//         SMembrane forces
//         ----------------------------
//          */

//         boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
//         simulator.AddForce(p_shear_force);

//         // Create a plane boundary to represent the inlets/outlets and pass them to the simulation

//         std::vector<c_vector<double,3> > boundary_plane_points;
//         std::vector<c_vector<double,3> > boundary_plane_normals;

//         // TRACE("These BCs need fixing, they no longer cover the plexus -- but this is not the current problem ")
//         boundary_plane_points.push_back(Create_c_vector(0,0,0));
//         boundary_plane_normals.push_back(Create_c_vector(1,0,0));

//         boundary_plane_points.push_back(Create_c_vector(7,0,0));
//         boundary_plane_normals.push_back(Create_c_vector(-1,0,0));

//         for(unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
//         {
//             boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population, boundary_plane_points[boundary_id],-boundary_plane_normals[boundary_id],0.05));
//             simulator.AddCellPopulationBoundaryCondition(p_condition);
//         }


//      	simulator.Solve();
//     }





 };



#endif /*TESTRELAXATION_HPP_*/