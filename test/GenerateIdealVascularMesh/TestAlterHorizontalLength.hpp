#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

// Must be included before other cell_based headers
#include "PetscSetupAndFinalize.hpp"

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"

#include "Debug.hpp"
#include "SmartPointers.hpp"
#include "VtkMeshReader.hpp"
#include "VtkMeshWriter.hpp"

#include "AbstractCellBasedTestSuite.hpp"

#include "MeshBasedCellPopulation.hpp"

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:
    void TestParametersOnBifucation() throw(Exception)
    {

        int Collapsing[11] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        for (int j = 0; j < 11; ++j)
        {
            std::stringstream out;
            out << Collapsing[j];
            std::string MeshNumber = out.str();
        
            std::string mesh_file = "/Volumes/Hardrive/Projects/FSI/VaryingEdgeLength/CollectedVtus/mesh_" + MeshNumber + ".vtu";

            double Scalling[1] = {1};//0.2, 0.4, 0.6, 0.8, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3 };
            for (int i = 0; i < 1; ++i)
            {
                
                VtkMeshReader<2, 3> mesh_reader(mesh_file);
                MutableMesh<2, 3> mesh;
                mesh.ConstructFromMeshReader(mesh_reader);

                double Scale = Scalling[i];

                double A = 0.9325189791142847;
                double B = 1.8109551126806118;
                double C = 2.6272663534520184;
                double D = 3.475951713466752;
                double E = 4.334920411905742;
                double F = 5.190067179472811; //5.2055084892613595;
                double G = 6.011974880051098; //6.005186130930892;
                double H = 6.866473248571025; //6.8893658760407375;

                // First round
                // for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                //      node_iter != mesh.GetNodeIteratorEnd();
                //      ++node_iter)
                // {
                //     c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                //     c_vector<double, 3> NewLocation = InitalLocation;

                //     if (InitalLocation[0] < A) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                //     {

                //         NewLocation[0] = NewLocation[0] * Scale;

                //         (node_iter)->rGetModifiableLocation() = NewLocation;
                //     }
                //     else
                //     {

                //         NewLocation[0] += A * (Scale - 1);
                //         (node_iter)->rGetModifiableLocation() = NewLocation;
                //     }
                // }

                // B += A * (Scale - 1);
                // C += A * (Scale - 1);
                // D += A * (Scale - 1);
                // E += A * (Scale - 1);
                // F += A * (Scale - 1);
                // G += A * (Scale - 1);
                // H += A * (Scale - 1);

                // // Second
                // for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                //      node_iter != mesh.GetNodeIteratorEnd();
                //      ++node_iter)
                // {
                //     c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                //     c_vector<double, 3> NewLocation  = InitalLocation;

                //     double UpperBound = B;
                //     double LowerBound = C;
                //     if (InitalLocation[0] >= UpperBound && InitalLocation[0] <= LowerBound) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                //     {
                //         NewLocation[0] = (NewLocation[0]) * Scale - (UpperBound) * (Scale - 1);
                //         (node_iter)->rGetModifiableLocation() = NewLocation;
                //     }
                //     else if (InitalLocation[0] > LowerBound)
                //     {

                //         NewLocation[0] += (LowerBound - UpperBound) * (Scale - 1);
                //         (node_iter)->rGetModifiableLocation() = NewLocation;
                //     }
                // }

                // D += (C - B) * (Scale - 1);
                // E += (C - B) * (Scale - 1);
                // F += (C - B) * (Scale - 1);
                // G += (C - B) * (Scale - 1);
                // H += (C - B) * (Scale - 1);
                // // third
                // for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                //      node_iter != mesh.GetNodeIteratorEnd();
                //      ++node_iter)
                // {
                //     c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                //     c_vector<double, 3> NewLocation = InitalLocation;

                //     double UpperBound = D;
                //     double LowerBound = E;
                //     if (InitalLocation[0] >= UpperBound && InitalLocation[0] <= LowerBound) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                //     {

                //         NewLocation[0] = (NewLocation[0]) * Scale - (UpperBound) * (Scale - 1);
                //         (node_iter)->rGetModifiableLocation() = NewLocation;
                //     }
                //     else if (InitalLocation[0] > LowerBound)
                //     {

                //         NewLocation[0] += (LowerBound - UpperBound) * (Scale - 1);
                //         (node_iter)->rGetModifiableLocation() = NewLocation;
                //     }
                // }

                // F += (E - D) * (Scale - 1);
                // G += (E - D) * (Scale - 1);
                // H += (E - D) * (Scale - 1);

                // // Fourth
                // for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                //      node_iter != mesh.GetNodeIteratorEnd();
                //      ++node_iter)
                // {
                //     c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                //     c_vector<double, 3> NewLocation = InitalLocation;

                //     double UpperBound = F;
                //     double LowerBound = G;
                //     if (InitalLocation[0] >= UpperBound && InitalLocation[0] <= LowerBound) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                //     {
                //         NewLocation[0] = (NewLocation[0]) * Scale - (UpperBound) * (Scale - 1);
                //         (node_iter)->rGetModifiableLocation() = NewLocation;
                //     }
                //     else if (InitalLocation[0] > LowerBound)
                //     {

                //         NewLocation[0] += (LowerBound - UpperBound) * (Scale - 1);
                //         (node_iter)->rGetModifiableLocation() = NewLocation;
                //     }
                // }
                // H += (G - F) * (Scale - 1);

                // Fitfh
                // for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                //      node_iter != mesh.GetNodeIteratorEnd();
                //      ++node_iter)
                // {
                //     c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                //     c_vector<double, 3> NewLocation = InitalLocation;

                //     double UpperBound = H;
                //     if (InitalLocation[0] > UpperBound) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                //     {
                //         NewLocation[0] = NewLocation[0] * Scale - UpperBound * (Scale - 1);
                //         (node_iter)->rGetModifiableLocation() = NewLocation;
                //     }
                // }

                //Scale back to 0.2
                for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                    node_iter != mesh.GetNodeIteratorEnd(); ++node_iter)
                {
                    c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                    c_vector<double, 3> NewLocation = InitalLocation *0.2;

                    (node_iter)->rGetModifiableLocation() = NewLocation;
              
                }

                std::string OutputMeshFile = "mesh_" + MeshNumber; 
                    
                std::stringstream out;
                out <<  Scale;
                std::string Scal = out.str();

                std::string output_dir = "VariableEdgeLength/Scalling"+Scal+"/";
                VtkMeshWriter<2, 3> mesh_writer1(output_dir, OutputMeshFile, false);
                mesh_writer1.WriteFilesUsingMesh(mesh);


                // Convert the .vtu into an .stl 
                std::string vtuFile = "/Users/jcrawshaw/Documents/testoutput/VariableEdgeLength/Scalling"+Scal+"/mesh_"+MeshNumber+".vtu";
                std::string stlFile = "/Users/jcrawshaw/Documents/testoutput/VariableEdgeLength/Scalling"+Scal+"/mesh_"+MeshNumber+".stl";
                std::string vtu2stlCommand = "meshio-convert " + vtuFile + " " + stlFile + " > null";
                int SystemOutput = std::system(vtu2stlCommand.c_str()); // system only takes char *.

            }
        }
    }
};

#endif /*TESTRELAXATION_HPP_*/
