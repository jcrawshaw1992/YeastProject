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

#include "HemeLBForce.hpp"


class TestRemeshing  : public AbstractCellBasedTestSuite
{
public:
 void TestHemeLBReader() throw (Exception)
    {
      boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
      p_ForceOut->LoadTractionFromFile();
    }
  //   void OffTestHemeLBInBasicCylinder() throw (Exception)
  //   {
  //       double scale = 1e3;
  //       double Length = 50e-6 * scale;
  //       double Radius = 1e-6 * scale;

  //       unsigned N_D = 40;
  //       unsigned N_Z = 100;
        
  //       std::string output_dir = "TestHemeLBForce/";
    
  //       Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
  //       MutableMesh<2, 3>* p_mesh = generator.GetMesh();
  //       HistoryDepMutableMesh<2, 3>* mesh= static_cast<HistoryDepMutableMesh<2, 3>*>(p_mesh);
   
  //       // Create the cells 
  //       MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
  //       std::vector<CellPtr> cells;
  //       CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
  //       cells_generator.GenerateBasicRandom(cells, mesh->GetNumNodes(), p_differentiated_type);

  //       // Create a cell population
  //       HistoryDepMeshBasedCellPopulation<2, 3> cell_population(*mesh, cells);
  //       // cell_population.SetChasteOutputDirectory(output_dir, 0);
  //       cell_population.SetWriteVtkAsPoints(true);
  //       cell_population.SetOutputMeshInVtk(true);
  //       // Set population to output all data to results files
  //       cell_population.AddCellWriter<CellProliferativeTypesWriter>();

  //       // Set up cell-based simulation
  //       OffLatticeSimulation<2,3> simulator(cell_population);
  //       simulator.SetOutputDirectory(output_dir); 

  //       /*
  //       -----------------------------
  //       Add the HemeLB Force
  //       ----------------------------
  //       */        
  //       TRACE("Here to add HemeLB forces")
  //       c_vector<double, 3> PlaneNormal1 = Create_c_vector(0,0,1);
  //       c_vector<double, 3> Point1 = Create_c_vector(0,0,1e-6 * scale);

  //       c_vector<double, 3> PlaneNormal2 = Create_c_vector(0,0,-1);
  //       c_vector<double, 3> Point2 = Create_c_vector(0,0,49e-6 * scale);

  //       double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
  //       double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

  //       double InletPressure = 1.5*(0.002133152 - 0.001466542)*1.2; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later 
  //       double OutletPressure = 1.5*(0.002133152 - 0.001466542)*(0.98);

  //       boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
  //       // p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
  //       // p_ForceOut->Inlets(PlaneNormal2, Point2, OutletPressure, "Outlet");
  //       p_ForceOut->SetStartTime(0);
  //       // p_ForceOut->SetFluidSolidIterations(100);
  //       p_ForceOut->SetUpHemeLBConfiguration(output_dir, simulator.rGetCellPopulation());
  //       // simulator.AddForce(p_ForceOut);
    
  // }





};




#endif /*TESTRELAXATION_HPP_*/

