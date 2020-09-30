#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

// Must be included before other cell_based headers
// #include "PetscSetupAndFinalize.hpp"
#include <boost/filesystem.hpp>
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "Debug.hpp"
#include "SmartPointers.hpp"
#include "UblasCustomFunctions.hpp"
// #include <mpi.h>
#include <array>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <unistd.h>
#include <utility>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h> // Need to read the centerlines file

//-----

#include <boost/algorithm/string/predicate.hpp>
#include <boost/serialization/base_object.hpp>
#include <cstdlib>
#include <ctime>
#include "ChasteSerialization.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>
#include <vtkPolyData.h>
#include "VtkMeshWriter.hpp" // Need a mesh writer -- writing out the current mesh as a vtu so HemeLB can read

#include <boost/serialization/base_object.hpp>

#include <map>
#include "UblasCustomFunctions.hpp"
#include "XdrFileReader.hpp"

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:
   void Test0PressureHemeLBI2nlet() throw(Exception)
    {
        std::string CenterlinesFile = "/data/vascrem/testoutput/TestHemeLBOnNetwork/CollapsingMiddelBranch/HemeLB_results_from_time_0/Centerlines_20.vtp";
        /* Read the centerlines file and get the max radius */
        vtkSmartPointer<vtkXMLPolyDataReader> Reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        Reader->SetFileName(CenterlinesFile.c_str());
        Reader->Update();
        // std::cout << *Reader << std::endl;
        vtkPointData* point_data = Reader->GetOutput()->GetPointData();
        vtkPoints* Points = Reader->GetOutput()->GetPoints();
        double NumberOfDataPoints = Reader->GetOutput()->GetNumberOfPoints();
        // assert(point_data->HasArray(mRadiusdataName.c_str())); // No Radi point data
        // std::cout<<*point_data<<std::endl;
        //  std::cout<< *(Reader->GetOutput()->GetPoints()->GetPoint(0))<<std::endl;
        // std::cout<< *(Reader->GetOutput()->GetPoints())<<std::endl;
        vtkDataArray* p_scalars = point_data->GetArray("Radius");

        assert(p_scalars->GetNumberOfComponents() == 1); // Radi are scalars, so should only have one component for each data point, otherwise there is a problem

        std::vector<double> XValue;
        double MinRadius = 1000000;
        double Radius_0 = 0.0020;
        for (unsigned i = 0; i < NumberOfDataPoints; i++)
        {
            double* data = p_scalars->GetTuple(i); //RadiVector.push_back(*data);
            double Radius_t = *data;
            double RadiusChange = (Radius_t - Radius_0) / Radius_0 * 100; //-> Lets go for a 10% decrease -- change this descrese percent later and need to have some understanding of what the last radius was

            if (abs(RadiusChange) > 30 & * data> 0)
            {
              PRINT_2_VARIABLES(Radius_0,Radius_t)
              PRINT_2_VARIABLES(RadiusChange, *(Reader->GetOutput()->GetPoints()->GetPoint(i)))
                // Want to add a 0 pressure boundary here
                XValue.push_back(*(Reader->GetOutput()->GetPoints()->GetPoint(i)));
            }
            if (*data<MinRadius& * data> 0)
            {
                MinRadius = *data;
            }
        }
        double mRadius = MinRadius * 0.001;

        PRINT_VARIABLE(XValue.size())
        // determine if I need aditional 0 inlets ect
        if (XValue.size() > 0)
        {
            PRINT_VECTOR(XValue)
            // Additional inlets to be added
            if (XValue.size() == 1)
            {
                // Just one inlet
            }
            if (XValue.size() > 1)
            {
                //Going to have a start and end to the inlets -- if there are more than one region needing cutting off, im not sure what to do, but will deal with later
            }
        }
    }


};

#endif /*TESTRELAXATION_HPP_*/
