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
// #include "OutwardsPressureWithBreaks.hpp"
#include "OutwardsPressure.hpp"
#include "MembraneStiffnessForce.hpp"
#include "HemeLBForce.hpp"


class TestRemeshing  : public AbstractCellBasedTestSuite
{
public:

void TestHemeLBOnNetwork9() throw (Exception)
    {         

        double EndTime = 10;
        double scale = 1e-2;        
        std::string output_dir = "TestFlowThroughNonSymetricCollapse/9/";
        std::string mesh_file =   "/home/vascrem/MeshCollection/IdealisedNetwork/IdealMeshWIthCentralCollapse/Nonsymmetric/Mesh9.vtu";
          
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
        simulator.SetDt(0.0002);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(1);
        
        /*
        -----------------------------
        HemeLB Force
        ----------------------------
        */      
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point1 = Create_c_vector(0.005,0,0);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point2 = Create_c_vector(0.005,-0.014,0);

        c_vector<double, 3> PlaneNormal3 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point3 = Create_c_vector(0.072,0,0);

        c_vector<double, 3> PlaneNormal4 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point4 = Create_c_vector(0.072,-0.014,0);

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = (P_blood - P_tissue)*1.02; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later 
        double OutletPressure = (P_blood - P_tissue)*(0.98);

        boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal3, Point3, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal4, Point4, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(0);
        p_ForceOut->SetFluidSolidIterations(5000);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir, cell_population,1);
  }





 void TestHemeLBOnNetwork3() throw (Exception)
    {         

        double EndTime = 10;
        double scale = 1e-2;        
        std::string output_dir = "TestFlowThroughNonSymetricCollapse/3/";
        std::string mesh_file =   "/home/vascrem/MeshCollection/IdealisedNetwork/IdealMeshWIthCentralCollapse/Nonsymmetric/Mesh3.vtu";
          
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
        simulator.SetDt(0.0002);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(1);
        
        /*
        -----------------------------
        HemeLB Force
        ----------------------------
        */      
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point1 = Create_c_vector(0.005,0,0);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point2 = Create_c_vector(0.005,-0.014,0);

        c_vector<double, 3> PlaneNormal3 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point3 = Create_c_vector(0.072,0,0);

        c_vector<double, 3> PlaneNormal4 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point4 = Create_c_vector(0.072,-0.014,0);

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = (P_blood - P_tissue)*1.02; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later 
        double OutletPressure = (P_blood - P_tissue)*(0.98);

        boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal3, Point3, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal4, Point4, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(0);
        p_ForceOut->SetFluidSolidIterations(5000);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir, cell_population,1);
  }


void TestHemeLBOnNetwork1() throw (Exception)
    {         

        double EndTime = 10;
        double scale = 1e-2;        
        std::string output_dir = "TestFlowThroughNonSymetricCollapse/1/";
        std::string mesh_file =   "/home/vascrem/MeshCollection/IdealisedNetwork/IdealMeshWIthCentralCollapse/Nonsymmetric/Mesh1.vtu";
          
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
        simulator.SetDt(0.0002);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(1);
        
        /*
        -----------------------------
        HemeLB Force
        ----------------------------
        */      
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point1 = Create_c_vector(0.005,0,0);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point2 = Create_c_vector(0.005,-0.014,0);

        c_vector<double, 3> PlaneNormal3 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point3 = Create_c_vector(0.072,0,0);

        c_vector<double, 3> PlaneNormal4 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point4 = Create_c_vector(0.072,-0.014,0);

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = (P_blood - P_tissue)*1.02; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later 
        double OutletPressure = (P_blood - P_tissue)*(0.98);

        boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal3, Point3, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal4, Point4, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(0);
        p_ForceOut->SetFluidSolidIterations(5000);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir, cell_population,1);
  }

void TestHemeLBOnNetwork2() throw (Exception)
    {         

        double EndTime = 10;
        double scale = 1e-2;        
        std::string output_dir = "TestFlowThroughNonSymetricCollapse/2/";
        std::string mesh_file =   "/home/vascrem/MeshCollection/IdealisedNetwork/IdealMeshWIthCentralCollapse/Nonsymmetric/Mesh2.vtu";
          
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
        simulator.SetDt(0.0002);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(1);
        
        /*
        -----------------------------
        HemeLB Force
        ----------------------------
        */      
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point1 = Create_c_vector(0.005,0,0);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point2 = Create_c_vector(0.005,-0.014,0);

        c_vector<double, 3> PlaneNormal3 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point3 = Create_c_vector(0.072,0,0);

        c_vector<double, 3> PlaneNormal4 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point4 = Create_c_vector(0.072,-0.014,0);

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = (P_blood - P_tissue)*1.02; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later 
        double OutletPressure = (P_blood - P_tissue)*(0.98);

        boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal3, Point3, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal4, Point4, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(0);
        p_ForceOut->SetFluidSolidIterations(5000);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir, cell_population,1);
  }



void TestHemeLBOnNetwork4() throw (Exception)
    {         

        double EndTime = 10;
        double scale = 1e-2;        
        std::string output_dir = "TestFlowThroughNonSymetricCollapse/4/";
        std::string mesh_file =   "/home/vascrem/MeshCollection/IdealisedNetwork/IdealMeshWIthCentralCollapse/Nonsymmetric/Mesh4.vtu";
          
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
        simulator.SetDt(0.0002);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(1);
        
        /*
        -----------------------------
        HemeLB Force
        ----------------------------
        */      
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point1 = Create_c_vector(0.005,0,0);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point2 = Create_c_vector(0.005,-0.014,0);

        c_vector<double, 3> PlaneNormal3 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point3 = Create_c_vector(0.072,0,0);

        c_vector<double, 3> PlaneNormal4 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point4 = Create_c_vector(0.072,-0.014,0);

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = (P_blood - P_tissue)*1.02; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later 
        double OutletPressure = (P_blood - P_tissue)*(0.98);

        boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal3, Point3, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal4, Point4, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(0);
        p_ForceOut->SetFluidSolidIterations(5000);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir, cell_population,1);
  }


void TestHemeLBOnNetwork5() throw (Exception)
    {         

        double EndTime = 10;
        double scale = 1e-2;        
        std::string output_dir = "TestFlowThroughNonSymetricCollapse/5/";
        std::string mesh_file =   "/home/vascrem/MeshCollection/IdealisedNetwork/IdealMeshWIthCentralCollapse/Nonsymmetric/Mesh5.vtu";
          
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
        simulator.SetDt(0.0002);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(1);
        
        /*
        -----------------------------
        HemeLB Force
        ----------------------------
        */      
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point1 = Create_c_vector(0.005,0,0);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point2 = Create_c_vector(0.005,-0.014,0);

        c_vector<double, 3> PlaneNormal3 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point3 = Create_c_vector(0.072,0,0);

        c_vector<double, 3> PlaneNormal4 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point4 = Create_c_vector(0.072,-0.014,0);

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = (P_blood - P_tissue)*1.02; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later 
        double OutletPressure = (P_blood - P_tissue)*(0.98);

        boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal3, Point3, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal4, Point4, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(0);
        p_ForceOut->SetFluidSolidIterations(5000);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir, cell_population,1);
  }

void TestHemeLBOnNetwork6() throw (Exception)
    {         

        double EndTime = 10;
        double scale = 1e-2;        
        std::string output_dir = "TestFlowThroughNonSymetricCollapse/6/";
        std::string mesh_file =   "/home/vascrem/MeshCollection/IdealisedNetwork/IdealMeshWIthCentralCollapse/Nonsymmetric/Mesh6.vtu";
          
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
        simulator.SetDt(0.0002);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(1);
        
        /*
        -----------------------------
        HemeLB Force
        ----------------------------
        */      
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point1 = Create_c_vector(0.005,0,0);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point2 = Create_c_vector(0.005,-0.014,0);

        c_vector<double, 3> PlaneNormal3 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point3 = Create_c_vector(0.072,0,0);

        c_vector<double, 3> PlaneNormal4 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point4 = Create_c_vector(0.072,-0.014,0);

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = (P_blood - P_tissue)*1.02; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later 
        double OutletPressure = (P_blood - P_tissue)*(0.98);

        boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal3, Point3, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal4, Point4, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(0);
        p_ForceOut->SetFluidSolidIterations(5000);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir, cell_population,1);
  }


void TestHemeLBOnNetwork7() throw (Exception)
    {         

        double EndTime = 10;
        double scale = 1e-2;        
        std::string output_dir = "TestFlowThroughNonSymetricCollapse/7/";
        std::string mesh_file =   "/home/vascrem/MeshCollection/IdealisedNetwork/IdealMeshWIthCentralCollapse/Nonsymmetric/Mesh7.vtu";
          
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
        simulator.SetDt(0.0002);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(1);
        
        /*
        -----------------------------
        HemeLB Force
        ----------------------------
        */      
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point1 = Create_c_vector(0.005,0,0);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point2 = Create_c_vector(0.005,-0.014,0);

        c_vector<double, 3> PlaneNormal3 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point3 = Create_c_vector(0.072,0,0);

        c_vector<double, 3> PlaneNormal4 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point4 = Create_c_vector(0.072,-0.014,0);

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = (P_blood - P_tissue)*1.02; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later 
        double OutletPressure = (P_blood - P_tissue)*(0.98);

        boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal3, Point3, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal4, Point4, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(0);
        p_ForceOut->SetFluidSolidIterations(5000);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir, cell_population,1);
  }


void TestHemeLBOnNetwork8() throw (Exception)
    {         

        double EndTime = 10;
        double scale = 1e-2;        
        std::string output_dir = "TestFlowThroughNonSymetricCollapse/8/";
        std::string mesh_file =   "/home/vascrem/MeshCollection/IdealisedNetwork/IdealMeshWIthCentralCollapse/Nonsymmetric/Mesh8.vtu";
          
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
        simulator.SetDt(0.0002);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(1);
        
        /*
        -----------------------------
        HemeLB Force
        ----------------------------
        */      
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point1 = Create_c_vector(0.005,0,0);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point2 = Create_c_vector(0.005,-0.014,0);

        c_vector<double, 3> PlaneNormal3 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point3 = Create_c_vector(0.072,0,0);

        c_vector<double, 3> PlaneNormal4 = Create_c_vector(-1,0,0);
        c_vector<double, 3> Point4 = Create_c_vector(0.072,-0.014,0);

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = (P_blood - P_tissue)*1.02; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later 
        double OutletPressure = (P_blood - P_tissue)*(0.98);

        boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal3, Point3, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal4, Point4, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(0);
        p_ForceOut->SetFluidSolidIterations(5000);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir, cell_population,1);
  }



};




#endif /*TESTRELAXATION_HPP_*/
