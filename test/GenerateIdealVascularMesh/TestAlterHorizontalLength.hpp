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
    void offTestScaleRightAngle() throw(Exception)
    {

        // int Collapsing[11] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        // std::string Collapsing[11] = { "0", "0.1227", "0.2248", "0.3178", "0.417", "0.5124", "0.6119", "0.708", "0.8059", "0.9032", "1.0" };
        
        std::string Collapsing[11] = { "0.3178", "0.417", "0.5124", "0.6119", "0.708", "0.8059"};
        
        // std::string Collapsing[1] = {"0.1227"};
        for (int j = 0; j < 11; ++j)
        {

            double Angle = 2.2;
            // std::stringstream out;
            // out << Collapsing[j];
            // std::string MeshNumber = out.str();
            std::string MeshNumber = Collapsing[j];

            // Convert the .stl into a .vtn
            // std::string stlFile = "/Volumes/Hardrive/VASCREM_server/AngleVariation_3X3Network/AngleCollected/PI_6/mesh_" + MeshNumber + ".stl";
            std::string stlFile = "/Volumes/Hardrive/VASCREM_server/MeshCollection/IdealisedNetwork/AngleVariation_3X3NetworkExtended/PI_2.2/ScaledMesh." + MeshNumber + ".stl";

            std::string mesh_file = "/Volumes/Hardrive/VASCREM_server/MeshCollection/IdealisedNetwork/AngleVariation_3X3NetworkExtended/PI_2.2/ScaledMesh." + MeshNumber + ".vtu";
            std::string MeshioCommand = "meshio-convert " + stlFile + " " + mesh_file + " > null";
            int SystemOutput = std::system(MeshioCommand.c_str()); // system only takes char *.

            double Scalling[15] = { 1, 0.2, 0.4, 0.6, 0.8, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3 };
            for (int i = 0; i < 15; ++i)
            {

                VtkMeshReader<2, 3> mesh_reader(mesh_file);
                MutableMesh<2, 3> mesh;
                mesh.ConstructFromMeshReader(mesh_reader);

                double Scale = Scalling[i];

                double A;
                double B;
                double C;
                double D;
                double E;
                double F;
                double G;
                double H;

                if (Angle == 2.2)
                {

                    A = 0.25855930139587047;
                    B = 0.3544873077828776;
                    C = 0.5749512125739644;
                    D = 0.67;
                    E = 0.8917585890837827;
                    F = 0.9878080133140323;
                    G = 1.211095011657897;
                    H = 1.3065911826148642;
                }
                else
                {

                    A = 0.9325189791142847;
                    B = 1.8109551126806118;
                    C = 2.6272663534520184;
                    D = 3.475951713466752;
                    E = 4.334920411905742;
                    F = 5.190067179472811;
                    G = 6.011974880051098;
                    H = 6.866473248571025;
                }

                // First round
                for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                     node_iter != mesh.GetNodeIteratorEnd();
                     ++node_iter)
                {
                    c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                    c_vector<double, 3> NewLocation = InitalLocation;

                    if (InitalLocation[0] < A) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                    {

                        NewLocation[0] = NewLocation[0] * Scale;

                        (node_iter)->rGetModifiableLocation() = NewLocation;
                    }
                    else
                    {

                        NewLocation[0] += A * (Scale - 1);
                        (node_iter)->rGetModifiableLocation() = NewLocation;
                    }
                }

                B += A * (Scale - 1);
                C += A * (Scale - 1);
                D += A * (Scale - 1);
                E += A * (Scale - 1);
                F += A * (Scale - 1);
                G += A * (Scale - 1);
                H += A * (Scale - 1);

                // Second
                for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                     node_iter != mesh.GetNodeIteratorEnd();
                     ++node_iter)
                {
                    c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                    c_vector<double, 3> NewLocation = InitalLocation;

                    double UpperBound = B;
                    double LowerBound = C;
                    if (InitalLocation[0] >= UpperBound && InitalLocation[0] <= LowerBound) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                    {
                        NewLocation[0] = (NewLocation[0]) * Scale - (UpperBound) * (Scale - 1);
                        (node_iter)->rGetModifiableLocation() = NewLocation;
                    }
                    else if (InitalLocation[0] > LowerBound)
                    {

                        NewLocation[0] += (LowerBound - UpperBound) * (Scale - 1);
                        (node_iter)->rGetModifiableLocation() = NewLocation;
                    }
                }

                D += (C - B) * (Scale - 1);
                E += (C - B) * (Scale - 1);
                F += (C - B) * (Scale - 1);
                G += (C - B) * (Scale - 1);
                H += (C - B) * (Scale - 1);
                // third
                for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                     node_iter != mesh.GetNodeIteratorEnd();
                     ++node_iter)
                {
                    c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                    c_vector<double, 3> NewLocation = InitalLocation;

                    double UpperBound = D;
                    double LowerBound = E;
                    if (InitalLocation[0] >= UpperBound && InitalLocation[0] <= LowerBound) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                    {

                        NewLocation[0] = (NewLocation[0]) * Scale - (UpperBound) * (Scale - 1);
                        (node_iter)->rGetModifiableLocation() = NewLocation;
                    }
                    else if (InitalLocation[0] > LowerBound)
                    {

                        NewLocation[0] += (LowerBound - UpperBound) * (Scale - 1);
                        (node_iter)->rGetModifiableLocation() = NewLocation;
                    }
                }

                F += (E - D) * (Scale - 1);
                G += (E - D) * (Scale - 1);
                H += (E - D) * (Scale - 1);

                // Fourth
                for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                     node_iter != mesh.GetNodeIteratorEnd();
                     ++node_iter)
                {
                    c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                    c_vector<double, 3> NewLocation = InitalLocation;

                    double UpperBound = F;
                    double LowerBound = G;
                    if (InitalLocation[0] >= UpperBound && InitalLocation[0] <= LowerBound) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                    {
                        NewLocation[0] = (NewLocation[0]) * Scale - (UpperBound) * (Scale - 1);
                        (node_iter)->rGetModifiableLocation() = NewLocation;
                    }
                    else if (InitalLocation[0] > LowerBound)
                    {

                        NewLocation[0] += (LowerBound - UpperBound) * (Scale - 1);
                        (node_iter)->rGetModifiableLocation() = NewLocation;
                    }
                }
                H += (G - F) * (Scale - 1);

                // Fitfh
                for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                     node_iter != mesh.GetNodeIteratorEnd();
                     ++node_iter)
                {
                    c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                    c_vector<double, 3> NewLocation = InitalLocation;

                    double UpperBound = H;
                    if (InitalLocation[0] > UpperBound) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                    {
                        NewLocation[0] = NewLocation[0] * Scale - UpperBound * (Scale - 1);
                        (node_iter)->rGetModifiableLocation() = NewLocation;
                    }
                }

                // //Scale back to 0.2
                // for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                //     node_iter != mesh.GetNodeIteratorEnd(); ++node_iter)
                // {
                //     c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                //     c_vector<double, 3> NewLocation = InitalLocation *0.2;

                //     (node_iter)->rGetModifiableLocation() = NewLocation;

                // }

                std::string OutputMeshFile = "mesh_" + MeshNumber;

                std::stringstream out;
                out << Scale;
                std::string Scal = out.str();

                std::string output_dir = "VariableEdgeAndAngle/Scalling" + Scal + "/";
                VtkMeshWriter<2, 3> mesh_writer1(output_dir, OutputMeshFile, false);
                mesh_writer1.WriteFilesUsingMesh(mesh);

                // Convert the .vtu into an .stl
                std::string vtuFile = "/Users/jcrawshaw/Documents/testoutput/VariableEdgeAndAngle/Scalling" + Scal + "/mesh_" + MeshNumber + ".vtu";
                std::string stlFile = "/Users/jcrawshaw/Documents/testoutput/VariableEdgeAndAngle/Scalling" + Scal + "/mesh_" + MeshNumber + ".stl";
                std::string vtu2stlCommand = "meshio-convert " + vtuFile + " " + stlFile + " > null";
                int SystemOutput = std::system(vtu2stlCommand.c_str()); // system only takes char *.
            }
        }
    }

    void TestScaleRightOver3() throw(Exception)
    {

        // int Collapsing[11] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        // std::string Collapsing[11] = { "0", "0.1227", "0.2248", "0.3178", "0.417", "0.5124", "0.6119", "0.708", "0.8059", "0.9032", "1.0" };
        std::string Collapsing[6] = { "0.3178", "0.417", "0.5124", "0.6119", "0.708", "0.8059"};
       
        std::string Angle[4] = { "5", "3", "2.2", "6"};

        for (int j = 0; j < 6; ++j)
        {
            // std::stringstream out;
            // out << Collapsing[j];
            // std::string MeshNumber = out.str();
            std::string MeshNumber = Collapsing[j];

            for (int k = 0; k < 4; ++k)
            {
                std::string SelectedAngle = Angle[k];

                // Convert the .stl into a .vtn
                // std::string stlFile = "/Volumes/Hardrive/VASCREM_server/AngleVariation_3X3Network/AngleCollected/PI_6/mesh_" + MeshNumber + ".stl";
                std::string stlFile = "/Volumes/Hardrive/VASCREM_server/MeshCollection/IdealisedNetwork/AngleVariation_3X3NetworkExtended/PI_" + SelectedAngle + "/ScaledMesh." + MeshNumber + ".stl";

                std::string mesh_file = "/Volumes/Hardrive/VASCREM_server/MeshCollection/IdealisedNetwork/AngleVariation_3X3NetworkExtended/PI_" + SelectedAngle + "/ScaledMesh." + MeshNumber + ".vtu";
                std::string MeshioCommand = "meshio-convert " + stlFile + " " + mesh_file + " > null";
                // int SystemOutput = std::system(MeshioCommand.c_str()); // system only takes char *.

                double Scalling[10] = {0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2 };
                for (int i = 0; i < 10; ++i)
                {

                    VtkMeshReader<2, 3> mesh_reader(mesh_file);
                    MutableMesh<2, 3> mesh;
                    mesh.ConstructFromMeshReader(mesh_reader);

                    double Scale = Scalling[i];

                    double A;
                    double B;
                    double C;
                    double D;
                    double E;
                    double F;
                    double G;
                    double H;

                    if (SelectedAngle == "3")
                    {
                        A = 0.2237704726996297;
                        B = 0.34860911720103555;
                        C = 0.555683235673123;
                        D = 0.6794005694932105;
                        E =  0.8805250776721328;
                        F = 1.006400658645822;
                        G = 1.2150057186948247;
                        H = 1.3340147413729737;
                    }
                    else if (SelectedAngle == "5")
                    {
                        A = 0.1473897786924556;
                        B = 0.3634977388552934;
                        C = 0.4918463320958366;
                        D = 0.7171314297331874;//0.7202346867380669;
                        E = 0.8429869429456978;//0.849051504162333;
                        F = 1.0670553131580023;
                        G = 1.198363036164769;
                        H = 1.4174248685185353;
                    }
                    else if (SelectedAngle == "2.2")
                    {

                        A = 0.25855930139587047;
                        B = 0.3544873077828776;
                        C = 0.5749512125739644;
                        D = 0.67;
                        E = 0.8917585890837827;
                        F = 0.9878080133140323;
                        G = 1.211095011657897;
                        H = 1.3065911826148642;
                    }
                    else if (SelectedAngle == "6")
                    {

                        A = 0.11227749655811817; // 0.11639389678044389;
                        B = 0.3701005853263882;
                        C = 0.4633189069381215;
                        D = 0.7343337694628858;

                        E = 0.8242707003858315;

                        F = 1.096257629452491;

                        G = 1.1864628656760259;
                        H = 1.4560118307576777;
                    }
                    else
                    {
                        A = 0.9325189791142847;
                        B = 1.8109551126806118;
                        C = 2.6272663534520184;
                        D = 3.475951713466752;
                        E = 4.334920411905742;
                        F = 5.190067179472811;
                        G = 6.011974880051098;
                        H = 6.866473248571025;
                    }
                    
                    
                    
                    double TrueScalling = (E-D); // Treat (E-D) as the 1 that will know be scalled. 

                    double UnitLength =1;

                    double ToScaleToUnitLength = 1/TrueScalling;

                    double scaling = Scale *ToScaleToUnitLength;
                    // First round
                    for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                         node_iter != mesh.GetNodeIteratorEnd();
                         ++node_iter)
                    {
                        c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                        c_vector<double, 3> NewLocation = InitalLocation;

                        if (InitalLocation[0] < A) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                        {

                            NewLocation[0] = NewLocation[0] * Scale;

                            (node_iter)->rGetModifiableLocation() = NewLocation;
                        }
                        else
                        {

                            NewLocation[0] += A * (Scale - 1);
                            (node_iter)->rGetModifiableLocation() = NewLocation;
                        }
                    }

                    B += A * (Scale - 1);
                    C += A * (Scale - 1);
                    D += A * (Scale - 1);
                    E += A * (Scale - 1);
                    F += A * (Scale - 1);
                    G += A * (Scale - 1);
                    H += A * (Scale - 1);

                    // Second
                    for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                         node_iter != mesh.GetNodeIteratorEnd();
                         ++node_iter)
                    {
                        c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                        c_vector<double, 3> NewLocation = InitalLocation;

                        double UpperBound = B;
                        double LowerBound = C;
                        if (InitalLocation[0] >= UpperBound && InitalLocation[0] <= LowerBound) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                        {
                            NewLocation[0] = (NewLocation[0]) * Scale - (UpperBound) * (Scale - 1);
                            (node_iter)->rGetModifiableLocation() = NewLocation;
                        }
                        else if (InitalLocation[0] > LowerBound)
                        {

                            NewLocation[0] += (LowerBound - UpperBound) * (Scale - 1);
                            (node_iter)->rGetModifiableLocation() = NewLocation;
                        }
                    }

                    D += (C - B) * (Scale - 1);
                    E += (C - B) * (Scale - 1);
                    F += (C - B) * (Scale - 1);
                    G += (C - B) * (Scale - 1);
                    H += (C - B) * (Scale - 1);
                    // third
                    for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                         node_iter != mesh.GetNodeIteratorEnd();
                         ++node_iter)
                    {
                        c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                        c_vector<double, 3> NewLocation = InitalLocation;

                        double UpperBound = D;
                        double LowerBound = E;
                        if (InitalLocation[0] >= UpperBound && InitalLocation[0] <= LowerBound) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                        {

                            NewLocation[0] = (NewLocation[0]) * Scale - (UpperBound) * (Scale - 1);
                            (node_iter)->rGetModifiableLocation() = NewLocation;
                        }
                        else if (InitalLocation[0] > LowerBound)
                        {

                            NewLocation[0] += (LowerBound - UpperBound) * (Scale - 1);
                            (node_iter)->rGetModifiableLocation() = NewLocation;
                        }
                    }

                    F += (E - D) * (Scale - 1);
                    G += (E - D) * (Scale - 1);
                    H += (E - D) * (Scale - 1);

                    // Fourth
                    for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                         node_iter != mesh.GetNodeIteratorEnd();
                         ++node_iter)
                    {
                        c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                        c_vector<double, 3> NewLocation = InitalLocation;

                        double UpperBound = F;
                        double LowerBound = G;
                        if (InitalLocation[0] >= UpperBound && InitalLocation[0] <= LowerBound) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                        {
                            NewLocation[0] = (NewLocation[0]) * Scale - (UpperBound) * (Scale - 1);
                            (node_iter)->rGetModifiableLocation() = NewLocation;
                        }
                        else if (InitalLocation[0] > LowerBound)
                        {

                            NewLocation[0] += (LowerBound - UpperBound) * (Scale - 1);
                            (node_iter)->rGetModifiableLocation() = NewLocation;
                        }
                    }
                    H += (G - F) * (Scale - 1);

                    // Fitfh
                    for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                         node_iter != mesh.GetNodeIteratorEnd();
                         ++node_iter)
                    {
                        c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                        c_vector<double, 3> NewLocation = InitalLocation;

                        double UpperBound = H;
                        if (InitalLocation[0] > UpperBound) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                        {
                            NewLocation[0] = NewLocation[0] * Scale - UpperBound * (Scale - 1);
                            (node_iter)->rGetModifiableLocation() = NewLocation;
                        }
                    }

                    // //Scale back to 0.2
                    // for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                    //     node_iter != mesh.GetNodeIteratorEnd(); ++node_iter)
                    // {
                    //     c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                    //     c_vector<double, 3> NewLocation = InitalLocation *0.2;

                    //     (node_iter)->rGetModifiableLocation() = NewLocation;

                    // }

                    std::string OutputMeshFile = "mesh_" + MeshNumber;

                    std::stringstream out;
                    out << Scale;
                    std::string Scal = out.str();

                    std::string output_dir = "VariableEdgeAndAngle2nd/PI_" + SelectedAngle + "/Scalling" + Scal + "/";
                    VtkMeshWriter<2, 3> mesh_writer1(output_dir, OutputMeshFile, false);
                    mesh_writer1.WriteFilesUsingMesh(mesh);

                    // Convert the .vtu into an .stl
                    std::string vtuFile = "/Users/jcrawshaw/Documents/testoutput/VariableEdgeAndAngle2nd/PI_" + SelectedAngle + "/Scalling" + Scal + "/mesh_" + MeshNumber + ".vtu";
                    std::string stlFile = "/Users/jcrawshaw/Documents/testoutput/VariableEdgeAndAngle2nd/PI_" + SelectedAngle + "/Scalling" + Scal + "/mesh_" + MeshNumber + ".stl";
                    std::string vtu2stlCommand = "meshio-convert " + vtuFile + " " + stlFile + " > null";
                    int SystemOutput = std::system(vtu2stlCommand.c_str()); // system only takes char *.
                }
            }
        }
    }

    void offTestScaleRightOver4() throw(Exception)
    {

        // int Collapsing[11] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        // std::string Collapsing[11] = {"0","0.1227","0.2248", "0.3178", "0.417", "0.5124", "0.6119", "0.708", "0.8059", "0.9032", "1.0"};
        // std::string Collapsing[11] = {"0","0.1227","0.2248", "0.3178", "0.417", "0.5124", "0.6119", "0.708", "0.8059", "0.9032", "1.0"};
        std::string Collapsing[11] = {  "4", "5", "6", "7", "8" };

        std::string Angle[1] = { "4" };
        for (int k = 0; k < 1; ++k)
        {
            std::string SelectedAngle = Angle[k];

            for (int j = 0; j < 11; ++j)
            {
                // std::stringstream out;
                // out << Collapsing[j];
                // std::string MeshNumber = out.str();
                std::string MeshNumber = Collapsing[j];

                // Convert the .stl into a .vtn
                // std::string stlFile = "/Volumes/Hardrive/VASCREM_server/AngleVariation_3X3Network/AngleCollected/PI_6/mesh_" + MeshNumber + ".stl";
                // std::string stlFile = "/Volumes/Hardrive/VASCREM_server/MeshCollection/IdealisedNetwork/AngleVariation_3X3NetworkExtended/PI_"+SelectedAngle+"/ScaledMesh." + MeshNumber + ".stl";

                std::string mesh_file = "/Volumes/Hardrive/VASCREM_server/MeshCollection/IdealisedNetwork/AngleVariation_3X3NetworkExtended/PI_" + SelectedAngle + "/ScaledMesh." + MeshNumber + ".vtu";
                // std::string MeshioCommand = "meshio-convert " + stlFile + " " + mesh_file + " > null";
                // int SystemOutput = std::system(MeshioCommand.c_str()); // system only takes char *.

                double Scalling[15] = { 1, 0.2, 0.4, 0.6, 0.8, 1.2, 1.4, 1.6, 1.8, 2};
                for (int i = 0; i < 15; ++i)
                {

                    VtkMeshReader<2, 3> mesh_reader(mesh_file);
                    MutableMesh<2, 3> mesh;
                    mesh.ConstructFromMeshReader(mesh_reader);

                    double Scale = Scalling[i];

                    double A;
                    double B;
                    double C;
                    double D;
                    double E;
                    double F;
                    double G;
                    double H;

                    if (SelectedAngle == "3")
                    {
                        A = 0.2237704726996297;
                        B = 0.34860911720103555;
                        C = 0.555683235673123;
                        D = 0.6901191872103832;
                        E = 0.8769862837326823;

                        F = 1.006400658645822;
                        G = 1.2150057186948247;
                        H = 1.3340147413729737;
                    }
                    else if (SelectedAngle == "5")
                    {
                        A = 0.22845383082074686;
                        B = 0.34784642224668405;
                        C = 0.5629143113002606;
                        D = 0.6904865871118324;
                        E = 0.8761501353326354;
                        F = 1.0048333682604662;
                        G = 1.2094830291058754;
                        H = 1.3335130430796431;
                    }
                    else if (SelectedAngle == "2.2")
                    {

                        A = 0.25855930139587047;
                        B = 0.3544873077828776;
                        C = 0.5749512125739644;
                        D = 0.67;
                        E = 0.8917585890837827;
                        F = 0.9878080133140323;
                        G = 1.211095011657897;
                        H = 1.3065911826148642;
                    }
                    else if (SelectedAngle == "6")
                    {

                        A = 0.11639389678044389;
                        B = 0.36848021505958845;
                        C = 0.4577377613362114;
                        D = 0.7617086720001398;
                        E = 0.8038645877330738;
                        F = 1.0962260802115176;
                        G = 1.198273046717881;
                        H = 1.448964047735432;
                    }
                    else if (SelectedAngle == "4")
                    {
                        A = 0.9325189791142847/0.2;
                        B = 1.8109551126806118/0.2;
                        C = 2.6272663534520184/0.2;
                        D = 3.475951713466752/0.2;
                        E = 4.334920411905742/0.2;
                        F = 5.190067179472811/0.2;
                        G = 6.011974880051098/0.2;
                        H = 6.866473248571025/0.2;
                    }

                    // First round
                    for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                         node_iter != mesh.GetNodeIteratorEnd();
                         ++node_iter)
                    {
                        c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                        c_vector<double, 3> NewLocation = InitalLocation;

                        if (InitalLocation[0] < A) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                        {

                            NewLocation[0] = NewLocation[0] * Scale;

                            (node_iter)->rGetModifiableLocation() = NewLocation;
                        }
                        else
                        {

                            NewLocation[0] += A * (Scale - 1);
                            (node_iter)->rGetModifiableLocation() = NewLocation;
                        }
                    }

                    B += A * (Scale - 1);
                    C += A * (Scale - 1);
                    D += A * (Scale - 1);
                    E += A * (Scale - 1);
                    F += A * (Scale - 1);
                    G += A * (Scale - 1);
                    H += A * (Scale - 1);

                    // Second
                    for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                         node_iter != mesh.GetNodeIteratorEnd();
                         ++node_iter)
                    {
                        c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                        c_vector<double, 3> NewLocation = InitalLocation;

                        double UpperBound = B;
                        double LowerBound = C;
                        if (InitalLocation[0] >= UpperBound && InitalLocation[0] <= LowerBound) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                        {
                            NewLocation[0] = (NewLocation[0]) * Scale - (UpperBound) * (Scale - 1);
                            (node_iter)->rGetModifiableLocation() = NewLocation;
                        }
                        else if (InitalLocation[0] > LowerBound)
                        {

                            NewLocation[0] += (LowerBound - UpperBound) * (Scale - 1);
                            (node_iter)->rGetModifiableLocation() = NewLocation;
                        }
                    }

                    D += (C - B) * (Scale - 1);
                    E += (C - B) * (Scale - 1);
                    F += (C - B) * (Scale - 1);
                    G += (C - B) * (Scale - 1);
                    H += (C - B) * (Scale - 1);
                    // third
                    for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                         node_iter != mesh.GetNodeIteratorEnd();
                         ++node_iter)
                    {
                        c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                        c_vector<double, 3> NewLocation = InitalLocation;

                        double UpperBound = D;
                        double LowerBound = E;
                        if (InitalLocation[0] >= UpperBound && InitalLocation[0] <= LowerBound) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                        {

                            NewLocation[0] = (NewLocation[0]) * Scale - (UpperBound) * (Scale - 1);
                            (node_iter)->rGetModifiableLocation() = NewLocation;
                        }
                        else if (InitalLocation[0] > LowerBound)
                        {

                            NewLocation[0] += (LowerBound - UpperBound) * (Scale - 1);
                            (node_iter)->rGetModifiableLocation() = NewLocation;
                        }
                    }

                    F += (E - D) * (Scale - 1);
                    G += (E - D) * (Scale - 1);
                    H += (E - D) * (Scale - 1);

                    // Fourth
                    for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                         node_iter != mesh.GetNodeIteratorEnd();
                         ++node_iter)
                    {
                        c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                        c_vector<double, 3> NewLocation = InitalLocation;

                        double UpperBound = F;
                        double LowerBound = G;
                        if (InitalLocation[0] >= UpperBound && InitalLocation[0] <= LowerBound) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                        {
                            NewLocation[0] = (NewLocation[0]) * Scale - (UpperBound) * (Scale - 1);
                            (node_iter)->rGetModifiableLocation() = NewLocation;
                        }
                        else if (InitalLocation[0] > LowerBound)
                        {

                            NewLocation[0] += (LowerBound - UpperBound) * (Scale - 1);
                            (node_iter)->rGetModifiableLocation() = NewLocation;
                        }
                    }
                    H += (G - F) * (Scale - 1);

                    // Fitfh
                    for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                         node_iter != mesh.GetNodeIteratorEnd();
                         ++node_iter)
                    {
                        c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                        c_vector<double, 3> NewLocation = InitalLocation;

                        double UpperBound = H;
                        if (InitalLocation[0] > UpperBound) // &&  ((InitalLocation[1]<0.05&& InitalLocation[1]>-0.05) || (InitalLocation[1]<-0.24&& InitalLocation[1]>-0.32))  )
                        {
                            NewLocation[0] = NewLocation[0] * Scale - UpperBound * (Scale - 1);
                            (node_iter)->rGetModifiableLocation() = NewLocation;
                        }
                    }

                    // //Scale back to 0.2
                    // for (typename AbstractMesh<2, 3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                    //     node_iter != mesh.GetNodeIteratorEnd(); ++node_iter)
                    // {
                    //     c_vector<double, 3> InitalLocation = node_iter->rGetLocation();
                    //     c_vector<double, 3> NewLocation = InitalLocation *0.2;

                    //     (node_iter)->rGetModifiableLocation() = NewLocation;

                    // }

                    std::string OutputMeshFile = "mesh_" + MeshNumber;

                    std::stringstream out;
                    out << Scale;
                    std::string Scal = out.str();

                    std::string output_dir = "VariableEdgeAndAngle/PI_" + SelectedAngle + "/Scalling" + Scal + "/";
                    VtkMeshWriter<2, 3> mesh_writer1(output_dir, OutputMeshFile, false);
                    mesh_writer1.WriteFilesUsingMesh(mesh);

                    // Convert the .vtu into an .stl
                    std::string vtuFile = "/Users/jcrawshaw/Documents/testoutput/VariableEdgeAndAngle/PI_" + SelectedAngle + "/Scalling" + Scal + "/mesh_" + MeshNumber + ".vtu";
                    std::string stlFile = "/Users/jcrawshaw/Documents/testoutput/VariableEdgeAndAngle/PI_" + SelectedAngle + "/Scalling" + Scal + "/mesh_" + MeshNumber + ".stl";
                    std::string vtu2stlCommand = "meshio-convert " + vtuFile + " " + stlFile + " > null";
                    int SystemOutput = std::system(vtu2stlCommand.c_str()); // system only takes char *.
                }
            }
        }
    }
};

#endif /*TESTRELAXATION_HPP_*/
