#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>

#include "Debug.hpp"
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
#include <iostream>
#include <sstream>
#include <iomanip>

// #include "VtkMeshReader.hpp"

// #include "AppliedForceModifier.hpp"
#include "RemeshingTriggerOnHeteroMeshModifier.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "MembraneDeformationForce.hpp"
// #include "OutwardsPressure.hpp"

#include "HemeLBForce.hpp"

class TestRemeshing  : public AbstractCellBasedTestSuite
{
public:

    void TestHemeLBInBasicCylinder() throw (Exception)
    {
        double scale = 1e3;
        double Length = 50e-6 * scale;
        double Radius = 1e-6 * scale;

        unsigned N_D = 40;
        unsigned N_Z = 100;
        
        std::string output_dir = "TestRunning/";
    
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
        simulator.SetEndTime(20);

      
        /*
        -----------------------------
        Add the HemeLB Force
        ----------------------------
        */        
        TRACE("Here to add HemeLB forces")
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(0,0,1);
        c_vector<double, 3> Point1 = Create_c_vector(0,0,1e-6 * scale);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(0,0,-1);
        c_vector<double, 3> Point2 = Create_c_vector(0,0,49e-6 * scale);

        // double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        // double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = (0.002133152 - 0.001466542)*1.2; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later 
        double OutletPressure = (0.002133152 - 0.001466542)*(0.98);

        boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(0);
        p_ForceOut->SetFluidSolidIterations(2000);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir, simulator.rGetCellPopulation());
        simulator.AddForce(p_ForceOut);

        /*
        -----------------------------
        RemeshingTriggerOnHeteroMeshModifier
        ----------------------------
        */  
        boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnHeteroMeshModifier<2, 3>());
        p_Mesh_modifier->SetMembraneStrength(10);
        simulator.AddSimulationModifier(p_Mesh_modifier);

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
        //Create a plane boundary to represent the inlet and pass them to the simulation
        c_vector<long double, 3> Boundary1 = Create_c_vector(0, 0, 1e-6 * scale);
        c_vector<long double, 3> Normal1 = Create_c_vector(0, 0, 1);

        c_vector<long double, 3> Boundary2 = Create_c_vector(0,0,49e-6 * scale);
        c_vector<long double, 3> Normal2 = Create_c_vector(0, 0, -1);
        
        boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary1, Normal1, 1));
        simulator.AddCellPopulationBoundaryCondition(p_condition_1);

        boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary2, Normal2, 1));
        simulator.AddCellPopulationBoundaryCondition(p_condition_2);


        simulator.Solve();
        
    
  }


//  void TestHemeLBReader() throw (Exception)
//     {
//       // boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
//       // p_ForceOut->LoadTractionFromFile();

//       double var = 1.0000000001;
//       std::stringstream ss;
//       ss << setprecision(16) << var;
//       std::string str;
//       ss >> str;
//       std::cout << str << std::endl;
  // std::string update_xml_file = "python projects/VascularRemodelling/apps/update_xml_file.py -period 100 -directory /data/vascrem/testoutput/TestHemeLBForce/HemeLBFluid -InitalConditions 0.0003 -ConvergenceTermination true -AveragePressure 0.004"; 
  //     int SystemOutput = std::system(update_xml_file.c_str());

  // scp /Users/jcrawshaw/Chaste/projects/VascularRemodelling/src/MembraneForces/HemeLBForce.cpp vascrem@josborne.science.unimelb.edu.au:/home/vascrem/Chaste/projects/VascularRemodelling/src/MembraneForces/HemeLBForce.cpp

    

  //  scp vascrem@josborne.science.unimelb.edu.au:surface-pressure_526.vtu /Users/jcrawshaw/


//     }




};




#endif /*TESTRELAXATION_HPP_*/
