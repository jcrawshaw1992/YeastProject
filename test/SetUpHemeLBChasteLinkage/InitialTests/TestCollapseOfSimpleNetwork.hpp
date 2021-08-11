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
// #include "XmlTools.hpp"

#include "UblasCustomFunctions.hpp"

#include "VtkMeshReader.hpp"

// #include "CellMutationStatesWriter.hpp"
#include "OutwardsPressure.hpp"

#include "MembraneDeformationForce.hpp"

#include "RemeshingTriggerModifier.hpp"
#include "RemeshingTriggerOnHeteroMeshModifier.hpp"
#include "HistoryDepMutableMesh.hpp"
#include "HistoryDepMeshBasedCellPopulation.hpp"
// #include "ElementQualityOutputModifier.hpp" -- there are some issues here 
#include "HemeLBForce.hpp"

#include "Honeycomb3DCylinderMeshGenerator.hpp"




// using namespace xsd::cxx::tree;

class TestRemeshing  : public AbstractCellBasedTestSuite
{
public:



void offTestDeformationGrowthTOEqui() throw (Exception)
   {        
        // /Users/jcrawshaw/docker-polnet-master/IdealiseNetworks/3_by_4/
        std::string mesh_file =  "/Users/jcrawshaw/docker-polnet-master/IdealiseNetworks/3_by_4/IdealNetwork.vtu";//"projects/VascularRemodelling/test/data/bifurcation_cut/config.vtu";//
        std::string output_dir = "TestFlowThroughSimpleNetwork/GrowthToEquilibrium/";
        
        double mesh_scale = 1e-1;  //1e-3; 
        double startime = 0;
       
        // Read in the plexus 
        VtkMeshReader<2,3> mesh_reader(mesh_file);
        HistoryDepMutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Scale(mesh_scale,mesh_scale,mesh_scale); // so distances are in mm
      
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
        cell_population.SetChasteOutputDirectory(output_dir, startime);
        cell_population.SetInitialAnlgesAcrossMembrane();
        cell_population.SetRelativePath(output_dir, startime);
        // cell_population.SetTargetRemeshingEdgeLength(1e-3);
        cell_population.SetPrintRemeshedIC(1);
        // cell_population.SetTargetRemeshingIterations(5);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        // Set population to output all data to results files
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2,3> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetDt(0.01);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(1);

        /*
          -----------------------------
          RemeshingTriggerOnHeteroMeshModifier
          ----------------------------
        */  

        boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnHeteroMeshModifier<2, 3>());
        simulator.AddSimulationModifier(p_Mesh_modifier);

        // /*
        // -----------------------------
        // Constant Force
        // ----------------------------
        // */        

        // double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        // double P_tissue = 0.001466542; // Pa == 1.1000e-05 mmHg

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure(-0.00035);
        simulator.AddForce(p_ForceOut);

        // /*
        // -----------------------------
        // SMembrane forces
        // ----------------------------
        // */

        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        simulator.AddForce(p_shear_force);

        // Create a plane boundary to represent the inlets/outlets and pass them to the simulation
        std::vector<c_vector<double,3> > boundary_plane_points;
        std::vector<c_vector<double,3> > boundary_plane_normals;

        // TRACE("These BCs need fixing, they no longer cover the plexus -- but this is not the current problem ")
        boundary_plane_points.push_back(Create_c_vector(0.5,0,0));
        boundary_plane_normals.push_back(Create_c_vector (1,0,0));

        boundary_plane_points.push_back(Create_c_vector(0.5,1.4,0));
        boundary_plane_normals.push_back(Create_c_vector (1,0,0));

        boundary_plane_points.push_back(Create_c_vector(0.5,-1.4,0));
        boundary_plane_normals.push_back(Create_c_vector (1,0,0));


        boundary_plane_points.push_back(Create_c_vector(7.4,0,0));
        boundary_plane_normals.push_back(Create_c_vector (-1,0,0));

        boundary_plane_points.push_back(Create_c_vector(7.4,1.4,0));
        boundary_plane_normals.push_back(Create_c_vector (-1,0,0));

        boundary_plane_points.push_back(Create_c_vector(7.4,-1.4,0));
        boundary_plane_normals.push_back(Create_c_vector (-1,0,0));

        for(unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
            boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population, boundary_plane_points[boundary_id],-boundary_plane_normals[boundary_id],0.05));
            simulator.AddCellPopulationBoundaryCondition(p_condition);
        }
      	simulator.Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);

    }


// /Users/jcrawshaw/petsc-3.7.7/include/petscsys.h:152:6: error: "PETSc was


void TestDeformationOfHoneyCombNetworkHetroAfterEqui() throw (Exception)
   {  


        /*
          -----------------------------
          Testing the fluid structure interaction on a basic cylinder 

          If this code works everything is fine
          ----------------------------
        */  
 
       double startime = 0;
       
        double scale = 1e3;
        double Length = 50e-6 * scale;
        double Radius = 1e-6 * scale;

        double Inlet =  0;
        double Outlet = Length;

        
        std::string output_dir = "TestHemeLBForceOnCylinder/";
    
        Honeycomb3DCylinderMeshGenerator generator(40, 100, Radius, Length);
        MutableMesh<2, 3>* p_mesh = generator.GetMesh();
        HistoryDepMutableMesh<2, 3>* mesh= static_cast<HistoryDepMutableMesh<2, 3>*>(p_mesh);
      
      
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh->GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(*mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2,3> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetDt(0.01);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(60);

        /*
          -----------------------------
          RemeshingTriggerOnHeteroMeshModifier
          ----------------------------
        */  
        boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnHeteroMeshModifier<2, 3>());
        std::map<double, c_vector<long double, 4> > GrowthMaps; 
        GrowthMaps[2] = Create_c_vector(0, pow(10, -6.8124), pow(10, -7), 0);
        GrowthMaps[1] = Create_c_vector(0, pow(10, -4), pow(10, -5), 0);
        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 2, 1,100, 1); 
        simulator.AddSimulationModifier(p_Mesh_modifier);

        /*
        -----------------------------
        HemeLB Force
        ----------------------------
        */        
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point1 = Create_c_vector(Inlet+ 0.001,0,0);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point2 = Create_c_vector(Outlet-0.001,0,0);

        boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
        TRACE("Set inlets and outlets :) ")
        p_ForceOut->Inlets(PlaneNormal1, Point1, 0.00025, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, 0.00005, "Outlet");
        TRACE("GOING in for the SetUpHemeLBConfiguration")
        p_ForceOut->SetUpHemeLBConfiguration(output_dir,  cell_population);
        simulator.AddForce(p_ForceOut);


        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        simulator.AddForce(p_shear_force);

        // Create a plane boundary to represent the inlets/outlets and pass them to the simulation

        std::vector<c_vector<double,3> > boundary_plane_points;
        std::vector<c_vector<double,3> > boundary_plane_normals;

        // // TRACE("These BCs need fixing, they no longer cover the plexus -- but this is not the current problem ")
        // boundary_plane_points.push_back(Create_c_vector(0,0,0));
        // boundary_plane_normals.push_back(Create_c_vector(1,0,0));

        // boundary_plane_points.push_back(Create_c_vector(7,0,0));
        // boundary_plane_normals.push_back(Create_c_vector(-1,0,0));

      //   for(unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
      //   {
      //       boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population, boundary_plane_points[boundary_id],-boundary_plane_normals[boundary_id],0.05));
      //       simulator.AddCellPopulationBoundaryCondition(p_condition);
      //   }


     	// simulator.Solve();
    }

 };



#endif /*TESTRELAXATION_HPP_*/