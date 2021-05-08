#ifndef TestBendingForceParameterSweep_HPP_
#define TestBendingForceParameterSweep_HPP_


#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

// #include <boost/filesystem.hpp>// Needed here to avoid serialization errors (on Boost<1.37)
// #include <ctime>
// #include <boost/algorithm/string/predicate.hpp>
// #include <cstdlib>
// // #include <direct.h>
// #include <boost/filesystem.hpp>
// #include <mpi.h>
// #include <unistd.h>  
// #include <fstream>

// #include <iostream>

// #include <cstdio>
// #include <memory>
// #include <string>
// #include <array>
// #include <utility>

#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h> // Need to read the centerlines file

#include <algorithm>
#include <vector>
#include <cmath>
#include "VtkMeshWriter.hpp" // Need a mesh writer -- writing out the current mesh as a vtu so HemeLB can read
#include <sstream>
#include <iomanip>

// #include <map>
#include "vtkXMLUnstructuredGridReader.h"
#include <vtkPolyData.h>
// #include "MeshBasedCellPopulation.hpp"
#include "Debug.hpp"
// #include <math.h>

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "PetscTools.hpp"
#include "XdrFileReader.hpp"

#include "UblasCustomFunctions.hpp"
#include "AbstractCellBasedTestSuite.hpp"
// #include "HemeLBForce.hpp"



class GetShearStress : public AbstractCellBasedTestSuite
{
    public:


    void TestGettingShearStress() throw(Exception)
    {
        std::vector < c_vector<double,3> > mAppliedTractions;
        std::vector < c_vector<double,3> > mAppliedTangentTractions;
        std::vector < c_vector<double,3> > mAppliedPosition;
        std::vector <double > mAppliedPressure;
    

    // LoadTractionFromVTKFile()
      std::string file = "/Volumes/Hardrive/Projects/FSI/UpperBranchCollapse/CollectedResults/surface-traction_5.5.vtu";

      vtkSmartPointer<vtkXMLUnstructuredGridReader> Reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
      Reader->SetFileName(file.c_str());
      Reader->Update();

      std::cout<< Reader->GetOutput()->GetPointData()->GetArray("traction")<<std::endl;
      vtkUnstructuredGrid* Output = Reader->GetOutput();
      vtkCellData *cellData = Output->GetCellData();
      vtkDataArray* TractionData = cellData->GetArray(0);
      std::cout<<  TractionData <<std::endl;
      

// vtkXMLUnstructuredGridReader

    //   vtkUnstructuredGrid* Output = Reader->GetOutput();
    //   vtkPointData* point_data = Reader->GetOutput()->GetPointData();

    //   double NumberOfDataPoints = Reader->GetOutput()->GetNumberOfPoints();

//       // Write all of the coordinates of the points in the vtkPolyData to the console.
//       for(int i = 0; i < Output->GetNumberOfPoints(); i++)
//       {
//         double p[3];
//         Output->GetPoint(i,p);
//         c_vector<double,3> Point;
//         Point[0]= p[0]; Point[1]= p[1]; Point[2]= p[2];
//         mAppliedPosition.push_back(Point*1e3);
//       }



//       PRINT_VARIABLE(NumberOfDataPoints)
//     //   std::cout << *point_data << std::endl;
//       std::cout << *Output << std::endl;
//     //   vtkPointData* point_data = Reader->GetOutput()->GetPointData();


//     // //   This will get the fluid property at each point -- still need the corrds 
//       vtkCellData *cellData = Output->GetCellData();

//       vtkDataSet *DataSet = Reader->GetOutputAsDataSet();

//       std::cout << *cellData << std::endl;
//       std::cout << *DataSet << std::endl;

      
//     //   for (int i = 0; i < cellData->GetNumberOfArrays(); i++)
//     //   {
//           PRINT_VARIABLE(cellData->GetNumberOfArrays());
//           vtkDataArray* TractionData = cellData->GetArray(0);
//           cout << "name " << TractionData->GetName() << endl;
//           PRINT_VARIABLE(TractionData->GetNumberOfTuples());

//           for (int j = 0; j < TractionData->GetNumberOfTuples(); j++)
//           {
//             //   double Tractions = TractionData->GetTuple1(j);
//             //  PRINT_VARIABLE(Tractions)
//             //   cout << "  value " << j << "th is " << value << endl;
//             //   mAppliedPressure.push_back(value);

        //   }
    // //   }
    // PRINT_2_VARIABLES(data->GetNumberOfTuples(), Output->GetNumberOfPoints())
    // assert(mAppliedPosition.size() == Output->GetNumberOfPoints());
    // assert(mAppliedTractions.size() == number_fluid_sites);
    // assert(mAppliedTangentTractions.size() == number_fluid_sites);


 }


 
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/