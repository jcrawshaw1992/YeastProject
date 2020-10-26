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
         std::string CenterlineEdgeFile; std::string CenterlinePointsFile; std::string output_file;double RadiValue;
        if (CommandLineArguments::Instance()->OptionExists("-CenterlineEdges"))
        {
            CenterlineEdgeFile = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-CenterlineEdges"); 
            CenterlinePointsFile = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-CenterlinePoints");  
            output_file = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-ofile");
            RadiValue = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-Radius");
        }else
        {
            CenterlineEdgeFile = "/Users/jcrawshaw/docker-polnet-master/IdealiseNetworks/CenterlineEdges.txt";
            CenterlinePointsFile = "/Users/jcrawshaw/docker-polnet-master/IdealiseNetworks/CenterlinePoints.txt";
            output_file = "/Users/jcrawshaw/docker-polnet-master/IdealiseNetworks/TestCenterlines.vtp";
            RadiValue =0.2;
        }
        int SPACE_DIM =2;
        
        // std::vector< std::vector<double>> CenterlinesPoints;
        // Get the points from the points.txt and create the vtk points
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        std::string line;
        ifstream PointsFile (CenterlinePointsFile);
        int isblank;
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
                    // PRINT_2_VARIABLES(X,Y)
                    // std::cout<<" 0: ("<<line[0]<<") 1: ("<<line[1] <<") 2: ("<<line[2]<< ") 3: ("<<line[3]<<") 4: ("<<line[4]<<") 5: ("<<line[5]<<") 6: ("<<line[6]<<") 7: ("<<line[7]<<") 8: ("<<line[8]<< ") 9: ("<<line[9]<< ") 10: ("<<line[10]<< ") 11: ("<<line[11]<< ") 12: ("<<line[12]<< ") 13: ("<<line[13]<< ") 14: ("<<line[14]<< ") 15: ("<<line[15]<< ") 16: ("<<line[16]<< std::endl;
                    // std::cout<<"\n\n"<< std::endl;
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
       
                // PRINT_2_VARIABLES( edge_ends[0] ,  edge_ends[1] )
                //  std::cout<<" 0: ("<<line[0]<<") 1: ("<<line[1] <<") 2: ("<<line[2]<< ") 3: ("<<line[3]<<") 4: ("<<line[4]<<") 5: ("<<line[5]<<") 6: ("<<line[6]<<") 7: ("<<line[7]<<") 8: ("<<line[8]<< std::endl;
                //  std::cout<<line<< std::endl;
                // std::cout<<"\n\n"<< std::endl;
 
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
        
        for (unsigned vertex_id = 0; vertex_id < points->GetNumberOfPoints(); vertex_id++)
        {
            

            // if (vertex_id == points->GetNumberOfPoints()-2  || vertex_id == points->GetNumberOfPoints()-1)
            // {
            //  radii_array->InsertNextValue(RadiValue); 
            // }else if (vertex_id == points->GetNumberOfPoints()-3 || vertex_id == points->GetNumberOfPoints()-4)
            // {
            //  radii_array->InsertNextValue(RadiValue); 
            // }
            // else{
                radii_array->InsertNextValue(RadiValue);
            // }
        }
        polydata->GetPointData()->AddArray(radii_array);
        polydata->GetPointData()->SetActiveScalars("Radius");

        std::string output_filename(output_file);
        WritePolyDataToDisk(polydata, output_filename);
    }
};
#endif /*TESTRELAXATION_HPP_*/


 // Create and fill a point container
        // vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        // for (unsigned vertex = 0; vertex < NumberOfPoints; vertex++)
        // {
        //     points->InsertNextPoint(CenterlinesPoints[vertex][0],
        //                             CenterlinesPoints[vertex][1],
        //                             CenterlinesPoints[vertex][2]);
        //                             PRINT_3_VARIABLES(CenterlinesPoints[vertex][0],
        //                             CenterlinesPoints[vertex][1],
        //                             CenterlinesPoints[vertex][2])
        // }

         // vtkIdType edge_ends[2];
        // for (unsigned edge_id = 0; edge_id < num_edges; edge_id++)
        // {
        //     // The definition of edges comes from MATLAB and and the vertices are one-based indexed. Offset them.
        //     edge_ends[0] = 0;//acces_2d(edges, edge_id, 0, num_edges) - 1;
        //     edge_ends[1] = 1;//acces_2d(edges, edge_id, 1, num_edges) - 1;

        //     polydata->InsertNextCell(VTK_LINE, 2, edge_ends);
        //     double edge_length = 1;//DistanceTwoVertices(vertices, NumberOfPoints, edge_ends) * pixel_to_um;
        //     lengths_array->InsertNextValue(edge_length);

        //     // acces_2d(edge_length_radii, edge_id, 0, num_edges) = edge_length;
        //     // acces_2d(edge_length_radii, edge_id, 1, num_edges) = radii_fudge_factor * pixel_to_um * (radii[edge_ends[0]] + radii[edge_ends[1]]) / 2;
        // }
 
     
