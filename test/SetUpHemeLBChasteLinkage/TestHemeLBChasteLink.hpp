#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

//  You might like to fix up the boundaries here

#include <cxxtest/TestSuite.h>
#include <cstdio>
#include <ctime>
#include <cmath>

#include "Debug.hpp"

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "XmlTools.hpp"
#include "UblasCustomFunctions.hpp"
#include "VtkMeshReader.hpp"
#include "SmartPointers.hpp"
 #include "AbstractCellBasedTestSuite.hpp"

#include "OffLatticeSimulation.hpp"
#include "CellIdWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VtkMeshWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "HistoryDepMutableMesh.hpp"
#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "HemeLBForce.hpp"


#include <boost/process.hpp> 
#include <string> 
#include <vector> 


#include <cstdio>
#include <iostream>
#include <memory>
#include <string>
#include <array>
#include <utility>
#include <fstream>

#include <iomanip>

#include <sstream> 

#include "vtkXMLUnstructuredGridReader.h"
#include <vtkCellCenters.h>

#include <vtkPolyData.h>

#include <vtkCellLocator.h>

#include "AppliedForceModifier.hpp"

using namespace xsd::cxx::tree;

class TestRemeshing  : public AbstractCellBasedTestSuite
{
public:





void offTestRUnningtings() throw (Exception)
    {
    
        std::string file = "/Users/jcrawshaw/Documents/testoutput/TestHemeLBChasteLinkage/HemeLBFluid/results/Extracted/surface-pressure_3348.vtu";
    
        
      vtkSmartPointer<vtkXMLUnstructuredGridReader> Reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
      Reader->SetFileName(file.c_str());
      Reader->Update();


        vtkSmartPointer<vtkCellLocator> cellLocator = 
        vtkSmartPointer<vtkCellLocator>::New();
        cellLocator->SetDataSet(Reader->GetOutput());
        cellLocator->BuildLocator();



      vtkUnstructuredGrid* Output = Reader->GetOutput();
      vtkPointData* point_data = Reader->GetOutput()->GetPointData();
      double NumberOfDataPoints = Reader->GetOutput()->GetNumberOfPoints();
      PRINT_VARIABLE(NumberOfDataPoints)


      // Write all of the coordinates of the points in the vtkPolyData to the console.
      // for(int i = 0; i < Output->GetNumberOfPoints(); i++)
      // {
      //   double p[3];
      //   Output->GetPoint(i,p);
      //   // // This is identical to:
      //   // // polydata->GetPoints()->GetPoint(i,p);
      //   std::cout << "Point " << i << " : (" << p[0] << " " << p[1] << " " << p[2] << ")" << std::endl;
      //   c_vector<double,3> Point = p;
      //   PRINT_VECTOR(Point)
      // }

      // // 
      std::cout << *point_data << std::endl;
      // std::cout << *Output << std::endl;

      vtkCellData *cellData = Output->GetCellData();


  vtkIdTypeArray* CellLocations = Output->GetCellLocationsArray();

  vtkIdType *locs = Output->GetCellLocationsArray()->GetPointer(2);

  vtkIdType maxid = Output->GetCellLocationsArray()->GetMaxId();
  PRINT_2_VARIABLES(maxid, locs)
  // std::cout <<Output->GetCellLocationsArray() << std::endl;
  //  vtkPointData* point_data2 = CellLocations->GetPointData();


      // std::cout << *CellLocations<< std::endl;
      //  std::cout << CellLocations->GetDataType()<< std::endl;
      // std::cout << CellLocations.vtkDataArray() << std::endl;

      // std::cout << *cellData << std::endl;
      // std::cout << *Output<< std::endl;

        vtkDataArray* data = cellData->GetArray(0);
        // std::cout << data->GetAttribute() << std::endl;
        // cout << "name " << data->GetName() << endl;
        //  cout << "Data" << *data << endl;
        PRINT_VARIABLE(NumberOfDataPoints)
        for (int j = 0; j < data->GetNumberOfTuples(); j++)
        {
          c_vector<double,3> Location = Create_c_vector(0,0,0);
          double value = data->GetTuple1(j);
          for (int i = CellLocations->GetTuple1(j); i<CellLocations->GetTuple1(j )+9; i++)
          {
             double p[3];
             Output->GetPoint(i,p);
             for (int k=0; k<3; k++)
             {Location[k]+=p[k]; }
            

          }
          Location/=3;
          // PRINT_VECTOR(Location)
            
            int comp;
          //  std::cout << CellLocations->GetTuple1(j )<< std::endl;
        }
  





    }



 void TestHemeLBChasteLinkage() throw (Exception)
    {

        double mesh_scale = 1e-3; 
        double scale = 1e3;
        double Length = 20e-6 * scale; //12e-3; //12e-3
        double Radius = 5e-6 * scale;

  
        unsigned N_D = 40;
        unsigned N_Z = 20;//N_D * 3;
        
        std::string output_dir = "TestHemeLBChasteLinkage/";
    
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
        simulator.SetSamplingTimestepMultiple(3);
        simulator.SetDt(0.02);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(4);

        /*
        -----------------------------
        HemeLB Force
        ----------------------------
        */        
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(0,0,1);
        c_vector<double, 3> Point1 = Create_c_vector(0,0,1e-6 * scale);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(0,0,-1);
        c_vector<double, 3> Point2 = Create_c_vector(0,0,19e-6 * scale);

        boost::shared_ptr<HemeLBForce<2,3>> p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, 0.0225, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, 0.02, "Outlet");
        p_ForceOut->SetUpHemeLBConfiguration(output_dir, *mesh, cell_population);
        simulator.AddForce(p_ForceOut);



     	  simulator.Solve();
}



};




#endif /*TESTRELAXATION_HPP_*/

