#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

#include <math.h>
// #include <vtkCellData.h>
#include <vtkCleanPolyData.h>
#include <vtkDataObject.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>  // I think one of these requires mex.h, which I dont know where it is, and i shouldnt need it  
#include <vtkSmartPointer.h>
#include <vtkTubeFilter.h>
#include <vtkVersion.h>
#include <vtkXMLPolyDataWriter.h>

#include "AbstractCellBasedTestSuite.hpp"

#include <fstream>

#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <stdlib.h>     /* atof */
#include <iostream> 
#include <string> 

#include "Debug.hpp"

// #define acces_2d(matrix, rowIndex, columnIndex, numRows) (matrix[(columnIndex) * (numRows) + (rowIndex)])

// using namespace xsd::cxx::tree;
#define convert(T, V) static constexpr T const T##_##V(V)
class TestRemeshing  : public AbstractCellBasedTestSuite
{
    #define convert(T, V) static constexpr T const T##_##V(V)
public:

     void WritePolyDataToDisk(const vtkSmartPointer<vtkPolyData> polyData, const std::string& fileName)
    {
        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(polyData);

        // Optional - set the mode. The default is binary.
        //writer->SetDataModeToBinary();
        //writer->SetDataModeToAscii();

        writer->Write();
    }

    void TestGenerateIdealMesh() throw (Exception)
    {
         std::string CenterlineEdgeFile; std::string CenterlinePointsFile; std::string output_file;double RadiValue;double Collapse;
         double LeftBound;  double RightBound;
        if (CommandLineArguments::Instance()->OptionExists("-CenterlineEdges"))
        {
            CenterlineEdgeFile = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-CenterlineEdges"); 
            CenterlinePointsFile = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-CenterlinePoints");  
            output_file = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-ofile");
            RadiValue = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-Radius");
            Collapse = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-Collapse");
            LeftBound = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-LeftBound");
            RightBound = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-RightBound");
                                                                    
        }else
        {
            CenterlineEdgeFile = "/Users/jcrawshaw/docker-polnet-master/IdealiseNetworks/3_by_4/CenterlineEdges.txt";
            CenterlinePointsFile = "/Users/jcrawshaw/docker-polnet-master/IdealiseNetworks/3_by_4/CenterlinePoints.txt";
            output_file = "/Users/jcrawshaw/docker-polnet-master/IdealiseNetworks/3_by_4/TestCenterlines.vtp";
            RadiValue =0.2;
            Collapse = 0.5;
            LeftBound =3.5; 
            RightBound = 4.5;
            
        }  
        int SPACE_DIM =2;
        
        // std::vector< std::vector<double>> CenterlinesPoints;
        // Get the points from the points.txt and create the vtk points
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        std::string line;
        ifstream PointsFile (CenterlinePointsFile);
        int isblank;
        std::vector<double> CollapsedRegion;
        double counter = 0;
        if (PointsFile.is_open())
        {
            while ( getline (PointsFile,line) )
            {
                if (SPACE_DIM ==3)
                {
                    // std::cout<<line<< std::endl;
                    points->InsertNextPoint(atof(&line[1]), atof(&line[4]), atof(&line[7]));
                }else if (SPACE_DIM ==2)
                {             
                    // std::cout<<line<< std::endl;
                    int i =1;
                    if(line[1]==' ') 
                    {
                        i =2;
                    }
                    double X = atof(&line[i]);
                    // Need to find the next corrdinate point on this line
                    while (line[i]!=' ')
                    {
                        // Skip 
                        i +=1;
                    }
                    i +=1;
                    double Y = atof(&line[i]);
                    points->InsertNextPoint(X,Y ,0);
                    if (Y ==0  &&  X >LeftBound && X<RightBound)
                    {
                        CollapsedRegion.push_back(counter);
                    }
                    counter +=1;
                }
            }
            PointsFile.close();
        }else 
        {
            cout << "Unable to open points file"; 
        }
        PointsFile.close();

        // Create a polydata object and add the points to it.
        vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
        polydata->SetPoints(points);

        // Create the array containing edge lengths
        vtkSmartPointer<vtkDoubleArray> lengths_array = vtkSmartPointer<vtkDoubleArray>::New();
        lengths_array->SetNumberOfComponents(1);
        lengths_array->SetName("EdgeLengths");

        // Allocate memory and add the definition of edges to the polydata object
        polydata->Allocate();

        vtkIdType edge_ends[2];    
        ifstream EdgesFile (CenterlineEdgeFile);
        if (EdgesFile.is_open())
        {
            while ( getline (EdgesFile,line) )
            { 
                int i =1;
                if(line[1]==' ') 
                {
                    i =2;
                }
                edge_ends[0] = atof(&line[i]); 
                while (line[i]!=' ')
                {
                    // Skip 
                    i +=1;
                }
                i +=1;
                edge_ends[1] = atof(&line[i]); 
 
                polydata->InsertNextCell(VTK_LINE, 2, edge_ends);
                double edge_length = 1;//DistanceTwoVertices(vertices, num_vertices, edge_ends) * pixel_to_um;
                lengths_array->InsertNextValue(edge_length);
            }
            EdgesFile.close();
        }else 
        {
            cout << "Unable to open edges file"; 
        }

          // Create and fill the array containing the plexus radii at each vertex.
        vtkSmartPointer<vtkDoubleArray> radii_array = vtkSmartPointer<vtkDoubleArray>::New();
        radii_array->SetNumberOfComponents(1);
        radii_array->SetName("Radius");
        
        for (double vertex_id = 0; vertex_id < points->GetNumberOfPoints(); vertex_id++)
        {
            if (std::find(CollapsedRegion.begin(), CollapsedRegion.end(), vertex_id) != CollapsedRegion.end()) /* v does contain x */
            {
              radii_array->InsertNextValue(RadiValue*Collapse);
            }else
            {
                radii_array->InsertNextValue(RadiValue);
            }
        }
        polydata->GetPointData()->AddArray(radii_array);
        polydata->GetPointData()->SetActiveScalars("Radius");

        std::string output_filename(output_file);
        WritePolyDataToDisk(polydata, output_filename);
    }
};
#endif /*TESTRELAXATION_HPP_*/

