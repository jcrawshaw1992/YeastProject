#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>
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
#include "FixedRegionBoundaryCondition.hpp"
#include "MembraneDeformationForce.hpp"
#include "OutwardsPressure.hpp"
#include "HemeLBForce.hpp"


class TestRemeshing  : public AbstractCellBasedTestSuite
{
public:
 void OffTestGrowToEquiCylinder() throw (Exception)
    {
  
        double EndTime = 50;
        double scale = 1e3;
        double Length = 50e-6 * scale;
        double Radius = 1e-6 * scale;

        unsigned N_D = 40;
        unsigned N_Z = 100;
        
        std::string output_dir = "TestHemeLBChaste/";

        // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
        // in um will be too large and break chaste without carefull playing with or a tiny time step

        Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
        MutableMesh<2, 3>* p_mesh = generator.GetMesh();
        HistoryDepMutableMesh<2, 3>* mesh= static_cast<HistoryDepMutableMesh<2, 3>*>(p_mesh);
   
        // Create the cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh->GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(*mesh, cells);
        // cell_population.SetChasteOutputDirectory(output_dir, 0);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        // Set population to output all data to results files
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
  

        // Set up cell-based simulation
        OffLatticeSimulation<2,3> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetDt(0.02);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(EndTime);
        
        /*
        -----------------------------
        RemeshingTriggerOnHeteroMeshModifier
        ----------------------------
        */  
        boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnHeteroMeshModifier<2, 3>());
        p_Mesh_modifier->SetMembraneStrength(5);
        simulator.AddSimulationModifier(p_Mesh_modifier);

        /*
        -----------------------------
        HemeLB Force
        ----------------------------
        */        
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(0,0,1);
        c_vector<double, 3> Point1 = Create_c_vector(0,0,1e-6 * scale);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(0,0,-1);
        c_vector<double, 3> Point2 = Create_c_vector(0,0,49e-6 * scale);

        double P_blood = 0.002133152*1.5; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542*1.5; // Pa == 1.5000e-05 mmHg

        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */        

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure((P_blood-P_tissue)/1);
        p_ForceOut->SetRadiusThreshold(5*Radius);
        simulator.AddForce(p_ForceOut);
 
        /*
        -----------------------------
        Boundary conditions
        ----------------------------
        */
        std::vector<c_vector<double,3> > boundary_plane_points;
        std::vector<c_vector<double,3> > boundary_plane_normals;

        boundary_plane_points.push_back(Point1);
        boundary_plane_normals.push_back(PlaneNormal1);

        boundary_plane_points.push_back(Point2);
        boundary_plane_normals.push_back(PlaneNormal2);

        for(unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
          boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population, boundary_plane_points[boundary_id],boundary_plane_normals[boundary_id],0.5));
          simulator.AddCellPopulationBoundaryCondition(p_condition);
        }
      
     	  simulator.Solve();
        CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(&simulator);
}

 void offTestCollapsingCylinder() throw (Exception)
    {        
        std::string output_dir = "TestHemeLBChaste/";
        double scale = 1e3; double EndTime = 50;

        // Load and fix any settings in the simulator 
        OffLatticeSimulation<2,3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2,3>, 3 >::Load(output_dir, EndTime);  

        /* Update the ouput directory for the population  */ 
        static_cast<HistoryDepMeshBasedCellPopulation<2,3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2,3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);
    
        /* Remove the constant pressure force   */ 
        // p_simulator->RemoveForce(0); // TRACE("RemoveForce will only work with the edit I made in OffLatticeSimulation.cpp line 69" )
        p_simulator->RemoveAllForces();
        p_simulator->SetEndTime(EndTime+800);
        p_simulator->SetSamplingTimestepMultiple(100);
        p_simulator->SetDt(0.002);
        p_simulator->SetOutputDirectory(output_dir+"AddingHemeLBForceToArchivedSimulation/");

        /*
        -----------------------------
        Add the HemeLB Force
        ----------------------------
        */        
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(0,0,1);
        c_vector<double, 3> Point1 = Create_c_vector(0,0,1e-6 * scale);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(0,0,-1);
        c_vector<double, 3> Point2 = Create_c_vector(0,0,49e-6 * scale);

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = (0.002133152 - 0.001466542)*1.02; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later 
        double OutletPressure = (0.002133152 - 0.001466542)*(0.98);

        boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(EndTime);
        p_ForceOut->SetFluidSolidIterations(10000);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir+"AddingHemeLBForceToArchivedSimulation/", p_simulator->rGetCellPopulation());
        p_simulator->AddForce(p_ForceOut);


        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_simulator->AddForce(p_shear_force);

        /* 
        -----------------------------
        Edit  RemeshingTriggerOnHeteroMeshModifier
        ----------------------------
         */  

        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2,3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnHeteroMeshModifier<2, 3> >(*iter);
      
        // Upstream 
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0,0,16e-6* scale);
        c_vector<double, 3> UpperPlaneNormal = Create_c_vector(0,0,1);
        // Down stream 
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0,0,34e-6 * scale);
        c_vector<double, 3> LowerPlaneNormal = Create_c_vector(0,0,-1);
        p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);

     	  p_simulator->Solve();
        CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(p_simulator);
  }


 void TestGrowToEquiIdealNetwork() throw (Exception)
    {
  
        double EndTime = 50;
        double scale = 1e-2;        
        std::string output_dir = "TestHemeLBOnNetwork/";
        std::string mesh_file = "/home/vascrem/MeshCollection/3_by_4/Mesh.vtu";
          
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Scale(scale,scale,scale);

       // Create the cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
  
        // Set up cell-based simulation
        OffLatticeSimulation<2,3> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetDt(0.02);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(EndTime);
        
        /*
        -----------------------------
        RemeshingTriggerOnHeteroMeshModifier
        ----------------------------
        */  
        boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnHeteroMeshModifier<2, 3>());
        p_Mesh_modifier->SetMembraneStrength(5);
        simulator.AddSimulationModifier(p_Mesh_modifier);

        /*
        -----------------------------
        HemeLB Force
        ----------------------------
        */        
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point1 = Create_c_vector(0.004,0,0);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point2 = Create_c_vector(0.004,-0.014,0);

        c_vector<double, 3> PlaneNormal3 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point3 = Create_c_vector(0.004,0.014,0);

        c_vector<double, 3> PlaneNormal4 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point4 = Create_c_vector(0.075,0.014,0);

        c_vector<double, 3> PlaneNormal5 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point5 = Create_c_vector(0.075,0,0);

        c_vector<double, 3> PlaneNormal6 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point6 = Create_c_vector(0.075,-0.014,0);

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */        

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure((P_blood-P_tissue)/1);
        // p_ForceOut->SetRadiusThreshold(0.002);
        p_ForceOut->SetInitialPosition(cell_population); 
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
        std::vector<c_vector<double,3> > boundary_plane_points;
        std::vector<c_vector<double,3> > boundary_plane_normals;

        boundary_plane_points.push_back(Point1);
        boundary_plane_normals.push_back(PlaneNormal1);

        boundary_plane_points.push_back(Point2);
        boundary_plane_normals.push_back(PlaneNormal2);
        
        Boundary_plane_points.push_back(Point3);
        boundary_plane_normals.push_back(PlaneNormal3);
        
        Boundary_plane_points.push_back(Point4);
        boundary_plane_normals.push_back(PlaneNormal4);
        
        Boundary_plane_points.push_back(Point5);
        boundary_plane_normals.push_back(PlaneNormal5);

        
        Boundary_plane_points.push_back(Point6);
        boundary_plane_normals.push_back(PlaneNormal6);

        for(unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
          boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population, boundary_plane_points[boundary_id],boundary_plane_normals[boundary_id],0.01));
          simulator.AddCellPopulationBoundaryCondition(p_condition);
        }
      
     	  simulator.Solve();
        CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(&simulator);
}

//  void offTestCollapsingIdeaNework() throw (Exception)
//     {        
//         std::string output_dir = "TestHemeLBChasteIdealNetwork/";
      
//         double scale = 1e3; double EndTime = 100;

//         // Load and fix any settings in the simulator 
//         OffLatticeSimulation<2,3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2,3>, 3 >::Load(output_dir, EndTime);  

//         /* Update the ouput directory for the population  */ 
//         static_cast<HistoryDepMeshBasedCellPopulation<2,3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
//         static_cast<HistoryDepMeshBasedCellPopulation<2,3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);
    
//         /* Remove the constant pressure force   */ 
//         // p_simulator->RemoveForce(0); // TRACE("RemoveForce will only work with the edit I made in OffLatticeSimulation.cpp line 69" )
//         p_simulator->SetEndTime(EndTime+800);
//         p_simulator->SetSamplingTimestepMultiple(200);
//         p_simulator->SetDt(0.002);
//         p_simulator->SetOutputDirectory(output_dir+"TestCollapseWithConstantForce/");

//         /* 
//         -----------------------------
//         Edit  RemeshingTriggerOnHeteroMeshModifier
//         ----------------------------
//          */  

//         std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
//         boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2,3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnHeteroMeshModifier<2, 3> >(*iter);
      

//         // Upstream 
//         c_vector<double, 3> UpperPlanePoint = Create_c_vector(1.8049900160683177,0.7,0)* scale;
//         c_vector<double, 3> UpperPlaneNormal = Create_c_vector(1,0,0);
//         // Down stream 
//         c_vector<double, 3> LowerPlanePoint = Create_c_vector(2.5456314376506,0.7,0)* scale;
//         c_vector<double, 3> LowerPlaneNormal = Create_c_vector(-1,0,0);
//         p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);

  /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        // boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        // simulator.AddForce(p_shear_force);

//      	  p_simulator->Solve();
//          CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(p_simulator);
//   }


};




#endif /*TESTRELAXATION_HPP_*/

