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
 
void TestGrowToEquiIdealNetwork() throw (Exception)
    {
        double EndTime = 50;
        double scale = 1e-2;        
        std::string output_dir =  "CollapsingNetwork_Test2/";
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
        simulator.SetDt(0.00001);
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
        Constant Compressive tissue pressure
        ----------------------------
        */    
       
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure((P_blood-P_tissue)*0.00001); //*0.000000001);
        p_ForceOut->SetInitialPosition(cell_population); 
        simulator.AddForce(p_ForceOut);



        /*
        -----------------------------
        HemeLB Force
        ----------------------------
        */        
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point1 = Create_c_vector(0.005,0,0);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point2 = Create_c_vector(0.005,-0.014,0);

        c_vector<double, 3> PlaneNormal3 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point3 = Create_c_vector(0.005,0.014,0);

        c_vector<double, 3> PlaneNormal4 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point4 = Create_c_vector(0.072,0.014,0);

        c_vector<double, 3> PlaneNormal5 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point5 = Create_c_vector(0.072,0,0);

        c_vector<double, 3> PlaneNormal6 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point6 = Create_c_vector(0.072,-0.014,0);

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
        
        boundary_plane_points.push_back(Point3);
        boundary_plane_normals.push_back(PlaneNormal3);
        
        boundary_plane_points.push_back(Point4);
        boundary_plane_normals.push_back(PlaneNormal4);
        
        boundary_plane_points.push_back(Point5);
        boundary_plane_normals.push_back(PlaneNormal5);

        boundary_plane_points.push_back(Point6);
        boundary_plane_normals.push_back(PlaneNormal6);

        for(unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
          boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population, boundary_plane_points[boundary_id],boundary_plane_normals[boundary_id],0.01));
          simulator.AddCellPopulationBoundaryCondition(p_condition);
        }
      
     	  simulator.Solve();
        CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(&simulator);
}

void offTestCollapsingIdeaNework() throw (Exception)
    {        
        std::string output_dir = "TestHemeLBOnNetwork/";
      
        double scale = 1e-2; double EndTime = 50;

        // Load and fix any settings in the simulator 
        OffLatticeSimulation<2,3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2,3>, 3 >::Load(output_dir, EndTime);  

        /* Update the ouput directory for the population  */ 
        static_cast<HistoryDepMeshBasedCellPopulation<2,3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2,3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);
    
        /* Remove the constant pressure force   */ 
        // p_simulator->RemoveForce(0); // TRACE("RemoveForce will only work with the edit I made in OffLatticeSimulation.cpp line 69" )
        p_simulator->RemoveAllForces();  
        p_simulator->SetEndTime(EndTime+800);
        p_simulator->SetSamplingTimestepMultiple(200);
        p_simulator->SetDt(0.002);
        p_simulator->SetOutputDirectory(output_dir+"AddingHemeLBForceToArchivedSimulation_2/");

            
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point1 = Create_c_vector(0.0043,0,0);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point2 = Create_c_vector(0.0043,-0.014,0);

        c_vector<double, 3> PlaneNormal3 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point3 = Create_c_vector(0.0043,0.014,0);

        c_vector<double, 3> PlaneNormal4 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point4 = Create_c_vector(0.0745,0.014,0);

        c_vector<double, 3> PlaneNormal5 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point5 = Create_c_vector(0.0745,0,0);

        c_vector<double, 3> PlaneNormal6 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point6 = Create_c_vector(0.0745,-0.014,0);


        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = 0.05*(P_blood - P_tissue)*1.02; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later 
        double OutletPressure = 0.05*(P_blood - P_tissue)*(0.98);

        boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal3, Point3, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal4, Point4, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal5, Point5, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal6, Point6, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(EndTime);
        p_ForceOut->SetFluidSolidIterations(5000);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir+"AddingHemeLBForceToArchivedSimulation_2/", p_simulator->rGetCellPopulation(),0);
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
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0.036,-0.014,0);
        c_vector<double, 3> UpperPlaneNormal = Create_c_vector(1,0,0);
        // Down stream 
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0.042,-0.014,0);
        c_vector<double, 3> LowerPlaneNormal = Create_c_vector(-1,0,0);
        p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);

     	  p_simulator->Solve();
         CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(p_simulator);
  }
void  offTestMiddelBranchCollapsing() throw (Exception)
    {        
        std::string output_dir = "TestHemeLBOnNetwork/";
      
        double scale = 1e-2; double EndTime = 50;

        // Load and fix any settings in the simulator 
        OffLatticeSimulation<2,3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2,3>, 3 >::Load(output_dir, EndTime);  

        /* Update the ouput directory for the population  */ 
        static_cast<HistoryDepMeshBasedCellPopulation<2,3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2,3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);
    
        /* Remove the constant pressure force   */ 
        // p_simulator->RemoveForce(0); // TRACE("RemoveForce will only work with the edit I made in OffLatticeSimulation.cpp line 69" )
        p_simulator->RemoveAllForces();  
        p_simulator->SetEndTime(EndTime+800);
        p_simulator->SetSamplingTimestepMultiple(200);
        p_simulator->SetDt(0.002);
        p_simulator->SetOutputDirectory(output_dir+"CollapsingMiddelBranch/");

            
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point1 = Create_c_vector(0.0043,0,0);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point2 = Create_c_vector(0.0043,-0.014,0);

        c_vector<double, 3> PlaneNormal3 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point3 = Create_c_vector(0.0043,0.014,0);

        c_vector<double, 3> PlaneNormal4 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point4 = Create_c_vector(0.0745,0.014,0);

        c_vector<double, 3> PlaneNormal5 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point5 = Create_c_vector(0.0745,0,0);

        c_vector<double, 3> PlaneNormal6 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point6 = Create_c_vector(0.0745,-0.014,0);


        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = 0.04*(P_blood - P_tissue)*1.02; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later 
        double OutletPressure = 0.04*(P_blood - P_tissue)*(0.98);

        boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal3, Point3, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal4, Point4, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal5, Point5, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal6, Point6, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(EndTime);
        p_ForceOut->SetFluidSolidIterations(5000);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir+"CollapsingMiddelBranch/", p_simulator->rGetCellPopulation());
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
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0.036,0.00,0);
        c_vector<double, 3> UpperPlaneNormal = Create_c_vector(1,0,0);
        // Down stream 
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0.042,0.00,0);
        c_vector<double, 3> LowerPlaneNormal = Create_c_vector(-1,0,0);
        p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);

     	  p_simulator->Solve();
         CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(p_simulator);
  }


 void offTestMiddelUpperCollapsing() throw (Exception)
    {        
        std::string output_dir = "TestHemeLBOnNetwork/";
      
        double scale = 1e-2; double EndTime = 50;

        // Load and fix any settings in the simulator 
        OffLatticeSimulation<2,3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2,3>, 3 >::Load(output_dir, EndTime);  

        /* Update the ouput directory for the population  */ 
        static_cast<HistoryDepMeshBasedCellPopulation<2,3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2,3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);
    
        /* Remove the constant pressure force   */ 
        // p_simulator->RemoveForce(0); // TRACE("RemoveForce will only work with the edit I made in OffLatticeSimulation.cpp line 69" )
        p_simulator->RemoveAllForces();  
        p_simulator->SetEndTime(EndTime+400);
        p_simulator->SetSamplingTimestepMultiple(1400);
        p_simulator->SetDt(0.002);
        p_simulator->SetOutputDirectory(output_dir+"CollapsingUpperBranch/");

            
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point1 = Create_c_vector(0.0043,0,0);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point2 = Create_c_vector(0.0043,-0.014,0);

        c_vector<double, 3> PlaneNormal3 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point3 = Create_c_vector(0.0043,0.014,0);

        c_vector<double, 3> PlaneNormal4 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point4 = Create_c_vector(0.0745,0.014,0);

        c_vector<double, 3> PlaneNormal5 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point5 = Create_c_vector(0.0745,0,0);

        c_vector<double, 3> PlaneNormal6 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point6 = Create_c_vector(0.0745,-0.014,0);


        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = 0.04*(P_blood - P_tissue)*1.02; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later 
        double OutletPressure = 0.04*(P_blood - P_tissue)*(0.98);

        boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal3, Point3, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal4, Point4, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal5, Point5, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal6, Point6, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(EndTime);
        p_ForceOut->SetFluidSolidIterations(5000);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir+"CollapsingUpperBranch/", p_simulator->rGetCellPopulation());
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
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0.036,0.014,0);
        c_vector<double, 3> UpperPlaneNormal = Create_c_vector(1,0,0);
        // Down stream 
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0.042,0.014,0);
        c_vector<double, 3> LowerPlaneNormal = Create_c_vector(-1,0,0);
        p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);

     	  p_simulator->Solve();
         CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(p_simulator);
  }


 void offTestMiddelLowerCollapsing() throw (Exception)
    {        
        std::string output_dir = "TestHemeLBOnNetwork/";
      
        double scale = 1e-2; double EndTime = 50;

        // Load and fix any settings in the simulator 
        OffLatticeSimulation<2,3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2,3>, 3 >::Load(output_dir, EndTime);  

        /* Update the ouput directory for the population  */ 
        static_cast<HistoryDepMeshBasedCellPopulation<2,3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2,3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);
    
        /* Remove the constant pressure force   */ 
        // p_simulator->RemoveForce(0); // TRACE("RemoveForce will only work with the edit I made in OffLatticeSimulation.cpp line 69" )
        p_simulator->RemoveAllForces();  
        p_simulator->SetEndTime(EndTime+200);
        p_simulator->SetSamplingTimestepMultiple(1400);
        p_simulator->SetDt(0.002);
        p_simulator->SetOutputDirectory(output_dir+"CollapsingLowerBranch/");

            
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point1 = Create_c_vector(0.0043,0,0);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point2 = Create_c_vector(0.0043,-0.014,0);

        c_vector<double, 3> PlaneNormal3 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point3 = Create_c_vector(0.0043,0.014,0);

        c_vector<double, 3> PlaneNormal4 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point4 = Create_c_vector(0.0745,0.014,0);

        c_vector<double, 3> PlaneNormal5 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point5 = Create_c_vector(0.0745,0,0);

        c_vector<double, 3> PlaneNormal6 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point6 = Create_c_vector(0.0745,-0.014,0);


        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = 0.04*(P_blood - P_tissue)*1.02; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later 
        double OutletPressure = 0.04*(P_blood - P_tissue)*(0.98);

        boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal3, Point3, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal4, Point4, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal5, Point5, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal6, Point6, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(EndTime);
        p_ForceOut->SetFluidSolidIterations(5000);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir+"CollapsingLowerBranch/", p_simulator->rGetCellPopulation());
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
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0.036,-0.014,0);
        c_vector<double, 3> UpperPlaneNormal = Create_c_vector(1,0,0);
        // Down stream 
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0.042,-0.014,0);
        c_vector<double, 3> LowerPlaneNormal = Create_c_vector(-1,0,0);
        p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);

     	  p_simulator->Solve();
         CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(p_simulator);
  }


void offTestCollapseWithoutHemeLB() throw (Exception)
    {        
        std::string output_dir = "TestHemeLBOnNetwork/";
      
        double scale = 1e-2; double EndTime = 50;

        // Load and fix any settings in the simulator 
        OffLatticeSimulation<2,3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2,3>, 3 >::Load(output_dir, EndTime);  

        /* Update the ouput directory for the population  */ 
        static_cast<HistoryDepMeshBasedCellPopulation<2,3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2,3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);
    
        /* Remove the constant pressure force   */ 
        // p_simulator->RemoveForce(0); // TRACE("RemoveForce will only work with the edit I made in OffLatticeSimulation.cpp line 69" )
        // p_simulator->RemoveAllForces();  
        p_simulator->SetEndTime(EndTime+800);
        p_simulator->SetSamplingTimestepMultiple(200);
        p_simulator->SetDt(0.1);
        p_simulator->SetOutputDirectory(output_dir+"MiddleCollapse_HemeLBFree/");

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
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0.036,0.00,0);
        c_vector<double, 3> UpperPlaneNormal = Create_c_vector(1,0,0);
        // Down stream 
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0.042,0.00,0);
        c_vector<double, 3> LowerPlaneNormal = Create_c_vector(-1,0,0);
        p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);

     	  p_simulator->Solve();
         CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(p_simulator);
  }

};




#endif /*TESTRELAXATION_HPP_*/

