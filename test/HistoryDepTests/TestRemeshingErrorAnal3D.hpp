#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

//  You might like to fix up the boundaries here

#include <cmath>
#include <cstdio>
#include <ctime>
#include <cxxtest/TestSuite.h>
#include "Debug.hpp"
// Must be included before other cell_based headers
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "CellIdWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellsGenerator.hpp"
#include "CommandLineArguments.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "VtkMeshWriter.hpp"
// #include "MembraneStiffnessForce.hpp"
#include "PetscSetupAndFinalize.hpp"
// #include "PressureForce.hpp"
#include "BoundariesModifier.hpp"
#include "CellMutationStatesWriter.hpp"
#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "HistoryDepMutableMesh.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "MembraneDeformationForce.hpp"
#include "MembraneForcesBasic.hpp"
#include "MembranePropertiesSecModifier.hpp"
#include "OutwardsPressure.hpp"
#include "RemeshingTriggerModifier.hpp"
#include "UblasCustomFunctions.hpp"
#include "VtkMeshReader.hpp"
// #include "XmlTools.hpp"

// using namespace xsd::cxx::tree;

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:


 void TestCreateCylinder() throw(Exception)
    {
        unsigned N_D = 80;
        unsigned N_Z = 100;

        double Length = 10;
        double Radius = 3;

        Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
        MutableMesh<2, 3>* p_mesh = generator.GetMesh();

        std::string output_directory = "RemeshingErrorAnal3D";
        VtkMeshWriter<2, 3> mesh_writer0(output_directory, "Cylinder", false);
        mesh_writer0.WriteFilesUsingMesh(*p_mesh);
    }



    void OffTestRemeshingCylinder_CourseInitial_VMTK() throw(Exception)
    {

        std::vector<double> RootMeanSquared;
        std::vector<double> RootSTDSquared;
        std::vector<double> EdgeL;
        std::vector<double> NumbOfNodes;

        double EdgeLength = 1;
        double Iter = 10;
        std::string output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder/VMTK/";
        double startime = 0;
        // I want to control the edge length, So i am just going to write a new cylinder based of the 'perfect' cylinder
        // I have in the based file here

        // Convert the .vtu into an .stl
        std::string stlfile = output_dir + "Cylinder.stl";

        // std::string vtu2offCommand = "meshio-convert " + output_dir + "config.vtu " + stlfile;
        // std::system(vtu2offCommand.c_str()); // system only takes char *.

        std::string Remeshedstl = output_dir + "RemeshedCylinder.stl";
        // std::string RemeshCommand = "vmtksurfaceremeshing -ifile " + stlfile + " -iterations " + std::to_string(20) + " -edgelength "+ std::to_string(1) +" -ofile " + Remeshedstl;
        // std::system(RemeshCommand.c_str());

        // std::string Remeshedvtu = output_dir + "CourseInital/RemeshedCylinder.vtu";
        // std::string stl2vtuCommand = " meshio-convert " + Remeshedstl + " " + Remeshedvtu;
        // std::system(stl2vtuCommand.c_str());

        while (EdgeLength > 0.01)
        {
            // Read in the new cylinder that is generated to my desire
            std::string mesh_file = "/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder/VMTK/CourseInital/RemeshedCylinder.vtu";
            VtkMeshReader<2, 3> mesh_reader(mesh_file);
            HistoryDepMutableMesh<2, 3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            std::stringstream out;
            out << EdgeLength;
            std::string mesh_size = out.str();
            output_dir = "RemeshingErrorAnal3D/Cylinder/VMTK/CourseInital/EdgeLength/" + mesh_size;

            // Create the cells
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

            // Create a cell population
            HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
            cell_population.SetChasteOutputDirectory(output_dir, startime);
            cell_population.SetInitialAnlgesAcrossMembrane();
            cell_population.SetTargetRemeshingEdgeLength(EdgeLength);
            cell_population.SetTargetRemeshingIterations(Iter);

            cell_population.SetRemeshingSoftwear("VMTK");
            cell_population.SetRelativePath(output_dir, startime);
            cell_population.SetPrintRemeshedIC(1);
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.SetOutputMeshInVtk(true);

            VtkMeshWriter<2, 3> mesh_writer0(output_dir, "InitalOriginalMesh", false);
            mesh_writer0.WriteFilesUsingMesh(mesh);

            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {

                c_vector<double, 3> InitalLocation = cell_population.GetNode(i)->rGetLocation();

                double R = sqrt(InitalLocation[0] * InitalLocation[0] + InitalLocation[1] * InitalLocation[1]);
                double Angle;
                // double Scalled_R;
                if (InitalLocation[0] >= 0)
                {
                    Angle = atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] <= 0)
                {
                    Angle = M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] >= 0)
                {
                    Angle = -M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }

                double X = R * cos(Angle) + 2 * sin(InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 2 * sin(Angle);

                c_vector<double, 3> DeformedLocation = Create_c_vector(X, Y, InitalLocation[2]);

                cell_population.GetNode(i)->rGetModifiableLocation() = DeformedLocation;
            }

            VtkMeshWriter<2, 3> mesh_writer1(output_dir, "DeformedOriginalMesh", false);
            mesh_writer1.WriteFilesUsingMesh(mesh);

            cell_population.ExecuteHistoryDependentRemeshing();

            VtkMeshWriter<2, 3> mesh_writer2(output_dir, "DeformedRemeshedMesh", false);
            mesh_writer2.WriteFilesUsingMesh(mesh);

            std::map<unsigned, c_vector<double, 3> > InitalNodePos = cell_population.GetInitalNodePositions();
            std::map<unsigned, c_vector<double, 3> > DeformedIC;
            std::vector<double> Error;

            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {
                c_vector<double, 3> InitalLocation = InitalNodePos[cell_population.GetNode(i)->GetIndex()];

                double R = sqrt(InitalLocation[0] * InitalLocation[0] + InitalLocation[1] * InitalLocation[1]);
                double Angle;
                // double Scalled_R;
                if (InitalLocation[0] >= 0)
                {
                    Angle = atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] <= 0)
                {
                    Angle = M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] >= 0)
                {
                    Angle = -M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }

                double X = R * cos(Angle) + 2 * sin(InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 2 * sin(Angle);

                c_vector<double, 3> DeformedLocation_IC = Create_c_vector(X, Y, InitalLocation[2]);

                DeformedIC[cell_population.GetNode(i)->GetIndex()] = DeformedLocation_IC;
                c_vector<double, 3> DeformedNode = cell_population.GetNode(i)->rGetLocation();

                Error.push_back(norm_2(DeformedNode - DeformedLocation_IC));

                // Now have a deformed place to put this and the error :)
            }
            double MeanError = 0;
            for (std::vector<double>::iterator iter = Error.begin(); iter != Error.end(); ++iter)
            {
                MeanError += (*iter) * (*iter);
            }
            MeanError = sqrt(MeanError);

            double rmsE = 0;
            for (std::vector<double>::iterator iter = Error.begin(); iter != Error.end(); ++iter)
            {
                rmsE += ((*iter) - MeanError) * ((*iter) - MeanError);
            }
            rmsE = sqrt(rmsE / Error.size());
            // Loop over nodes and deform the relevant inital poition and then test error -- write something u pin maths functions for the relevant erros
            VtkMeshWriter<2, 3> mesh_writer(output_dir, "DeformedOriginalMesh", false);
            mesh_writer.WriteFilesUsingMesh(mesh);

            RootMeanSquared.push_back(MeanError);
            RootSTDSquared.push_back(rmsE);
            NumbOfNodes.push_back(mesh.GetNumNodes());

            EdgeL.push_back(EdgeLength);
            EdgeLength = EdgeLength - 0.1;

            output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder/VMTK/CourseInital/EdgeLength/" + mesh_size;
            // Convert the 4 meshes I need off to stl
            std::string ConvertCommand = "meshio-convert " + output_dir + "/InitalOriginalMesh.vtu " + output_dir + "/InitalOriginalMesh.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
            ConvertCommand = "meshio-convert " + output_dir + "/DeformedRemeshedMesh.vtu " + output_dir + "/DeformedRemeshedMesh.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
            ConvertCommand = "meshio-convert " + output_dir + "/DeformedOriginalMesh.vtu " + output_dir + "/DeformedOriginalMesh.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
            ConvertCommand = "meshio-convert " + output_dir + "/results_from_time_0/NewInitalConfiguration1.vtu " + output_dir + "/NewInitalConfiguration.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
        }

        std::vector<double>::iterator j = RootSTDSquared.begin();
        std::vector<double>::iterator k = EdgeL.begin();
        std::vector<double>::iterator l = NumbOfNodes.begin();

        std::ofstream f("/Users/jcrawshaw/Documents/Projects/MeshMatlab/ErrorAna_3D_WavyDef_Cylinder_CourseInitial_VMTK.txt");
        for (std::vector<double>::iterator i = RootMeanSquared.begin(); i != RootMeanSquared.end(); ++j, ++k, ++l, ++i)
        {
            f << *k << " " << *i << " " << *j << " " << *l << '\n';
        }
    }

    void OffTestRemeshingCylinder_FineInitial_VMTK() throw(Exception)
    {

        std::vector<double> RootMeanSquared;
        std::vector<double> RootSTDSquared;
        std::vector<double> EdgeL;
        std::vector<double> NumbOfNodes;

        double EdgeLength = 1;
        double Iter = 10;
        std::string output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder/VMTK/";
        double startime = 0;
        // I want to control the edge length, So i am just going to write a new cylinder based of the 'perfect' cylinder
        // I have in the based file here

        // Convert the .vtu into an .stl
        std::string stlfile = output_dir + "Cylinder.stl";

        // std::string vtu2offCommand = "meshio-convert " + output_dir + "config.vtu " + stlfile;
        // std::system(vtu2offCommand.c_str()); // system only takes char *.

        std::string Remeshedstl = output_dir + "RemeshedCylinder.stl";
        // std::string RemeshCommand = "vmtksurfaceremeshing -ifile " + stlfile + " -iterations " + std::to_string(20) + " -edgelength "+ std::to_string(0.05) +" -ofile " + Remeshedstl;
        // std::system(RemeshCommand.c_str());

        std::string Remeshedvtu = output_dir + "FineInital/RemeshedCylinder.vtu";
        std::string stl2vtuCommand = " meshio-convert " + Remeshedstl + " " + Remeshedvtu;
        std::system(stl2vtuCommand.c_str());

        while (EdgeLength > 0.01)
        {
            // Read in the new cylinder that is generated to my desire
            std::string mesh_file = "/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder/VMTK/FineInital/RemeshedCylinder.vtu";
            VtkMeshReader<2, 3> mesh_reader(mesh_file);
            HistoryDepMutableMesh<2, 3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            std::stringstream out;
            out << EdgeLength;
            std::string mesh_size = out.str();
            output_dir = "RemeshingErrorAnal3D/Cylinder/VMTK/FineInital/EdgeLength/" + mesh_size;

            // Create the cells
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

            // Create a cell population
            HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
            cell_population.SetChasteOutputDirectory(output_dir, startime);
            cell_population.SetInitialAnlgesAcrossMembrane();
            cell_population.SetTargetRemeshingEdgeLength(EdgeLength);
            cell_population.SetTargetRemeshingIterations(Iter);
            cell_population.SetRemeshingSoftwear("VMTK");
            cell_population.SetRelativePath(output_dir, startime);
            cell_population.SetPrintRemeshedIC(1);
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.SetOutputMeshInVtk(true);

            VtkMeshWriter<2, 3> mesh_writer0(output_dir, "InitalOriginalMesh", false);
            mesh_writer0.WriteFilesUsingMesh(mesh);

            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {

                c_vector<double, 3> InitalLocation = cell_population.GetNode(i)->rGetLocation();

                double R = sqrt(InitalLocation[0] * InitalLocation[0] + InitalLocation[1] * InitalLocation[1]);
                double Angle;
                // double Scalled_R;
                if (InitalLocation[0] >= 0)
                {
                    Angle = atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] <= 0)
                {
                    Angle = M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] >= 0)
                {
                    Angle = -M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }

                double X = R * cos(Angle) + 2 * sin(InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 2 * sin(Angle);

                c_vector<double, 3> DeformedLocation = Create_c_vector(X, Y, InitalLocation[2]);

                cell_population.GetNode(i)->rGetModifiableLocation() = DeformedLocation;
            }

            VtkMeshWriter<2, 3> mesh_writer1(output_dir, "DeformedOriginalMesh", false);
            mesh_writer1.WriteFilesUsingMesh(mesh);

            cell_population.ExecuteHistoryDependentRemeshing();

            VtkMeshWriter<2, 3> mesh_writer2(output_dir, "DeformedRemeshedMesh", false);
            mesh_writer2.WriteFilesUsingMesh(mesh);

            std::map<unsigned, c_vector<double, 3> > InitalNodePos = cell_population.GetInitalNodePositions();
            std::map<unsigned, c_vector<double, 3> > DeformedIC;
            std::vector<double> Error;

            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {
                c_vector<double, 3> InitalLocation = InitalNodePos[cell_population.GetNode(i)->GetIndex()];

                double R = sqrt(InitalLocation[0] * InitalLocation[0] + InitalLocation[1] * InitalLocation[1]);
                double Angle;
                // double Scalled_R;
                if (InitalLocation[0] >= 0)
                {
                    Angle = atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] <= 0)
                {
                    Angle = M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] >= 0)
                {
                    Angle = -M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }

                double X = R * cos(Angle) + 2 * sin(InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 2 * sin(Angle);

                c_vector<double, 3> DeformedLocation_IC = Create_c_vector(X, Y, InitalLocation[2]);

                DeformedIC[cell_population.GetNode(i)->GetIndex()] = DeformedLocation_IC;
                c_vector<double, 3> DeformedNode = cell_population.GetNode(i)->rGetLocation();

                Error.push_back(norm_2(DeformedNode - DeformedLocation_IC));

                // Now have a deformed place to put this and the error :)
            }
            double MeanError = 0;
            for (std::vector<double>::iterator iter = Error.begin(); iter != Error.end(); ++iter)
            {
                MeanError += (*iter) * (*iter);
            }
            MeanError = sqrt(MeanError);

            double rmsE = 0;
            for (std::vector<double>::iterator iter = Error.begin(); iter != Error.end(); ++iter)
            {
                rmsE += ((*iter) - MeanError) * ((*iter) - MeanError);
            }
            rmsE = sqrt(rmsE / Error.size());
            // Loop over nodes and deform the relevant inital poition and then test error -- write something u pin maths functions for the relevant erros
            VtkMeshWriter<2, 3> mesh_writer(output_dir, "DeformedOriginalMesh", false);
            mesh_writer.WriteFilesUsingMesh(mesh);

            RootMeanSquared.push_back(MeanError);
            RootSTDSquared.push_back(rmsE);
            NumbOfNodes.push_back(mesh.GetNumNodes());

            EdgeL.push_back(EdgeLength);

            output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder/VMTK/FineInital/EdgeLength/" + mesh_size;
            // Convert the 4 meshes I need off to stl
            std::string ConvertCommand = "meshio-convert " + output_dir + "/InitalOriginalMesh.vtu " + output_dir + "/InitalOriginalMesh.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
            ConvertCommand = "meshio-convert " + output_dir + "/DeformedRemeshedMesh.vtu " + output_dir + "/DeformedRemeshedMesh.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
            ConvertCommand = "meshio-convert " + output_dir + "/DeformedOriginalMesh.vtu " + output_dir + "/DeformedOriginalMesh.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
            ConvertCommand = "meshio-convert " + output_dir + "/results_from_time_0/NewInitalConfiguration1.vtu " + output_dir + "/NewInitalConfiguration.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *

            EdgeLength = EdgeLength - 0.1;
        }

        std::vector<double>::iterator j = RootSTDSquared.begin();
        std::vector<double>::iterator k = EdgeL.begin();
        std::vector<double>::iterator l = NumbOfNodes.begin();

        std::ofstream f("/Users/jcrawshaw/Documents/Projects/MeshMatlab/ErrorAna_3D_WavyDef_Cylinder_FineInitial_VMTK.txt");
        for (std::vector<double>::iterator i = RootMeanSquared.begin(); i != RootMeanSquared.end(); ++j, ++k, ++l, ++i)
        {
            f << *k << " " << *i << " " << *j << " " << *l << '\n';
        }
    }

    void OffTestRemeshingCylinder_CourseRemeshedMesh_VMTK() throw(Exception)
    {
        std::vector<double> RootMeanSquared;
        std::vector<double> RootSTDSquared;
        std::vector<double> EdgeL;
        std::vector<double> NumbOfNodes;
        double EdgeLength = 1;
        double Iter = 10;
        double startime = 0;
        while (EdgeLength > 0.01)
        {
            std::stringstream out;
            out << EdgeLength;
            std::string mesh_size = out.str();
            std::string output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder/VMTK/CourseRemeshed/";

            // I want to control the edge length, So i am just going to write a new cylinder based of the 'perfect' cylinder
            // I have in the based file here

            // Convert the .vtu into an .stl
            std::string stlfile = output_dir + "Cylinder.stl";

            std::string vtu2offCommand = "meshio-convert ~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder/VMTK/Cylinder.vtu " + stlfile;
            std::system(vtu2offCommand.c_str()); // system only takes char *.

            std::string Remeshedstl = output_dir + "RemeshedCylinder.stl";
            std::string RemeshCommand = "vmtksurfaceremeshing -ifile " + stlfile + " -iterations " + std::to_string(20) + " -edgelength " + std::to_string(EdgeLength) + " -ofile " + Remeshedstl + " -elementsizemode edgelength";
            std::system(RemeshCommand.c_str());

            std::string Remeshedvtu = output_dir + "RemeshedCylinder.vtu";
            std::string stl2vtuCommand = " meshio-convert " + Remeshedstl + " " + Remeshedvtu;
            std::system(stl2vtuCommand.c_str());

            output_dir = "RemeshingErrorAnal3D/Cylinder/VMTK/CourseRemeshed/EdgeLength/" + mesh_size;

            // Read in the new cylinder that is generated to my desire
            std::string mesh_file = "/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder/VMTK/CourseRemeshed/RemeshedCylinder.vtu";
            VtkMeshReader<2, 3> mesh_reader(mesh_file);
            HistoryDepMutableMesh<2, 3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            // Create the cells
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

            // Create a cell population
            HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
            cell_population.SetChasteOutputDirectory(output_dir, startime);
            cell_population.SetInitialAnlgesAcrossMembrane();
            cell_population.SetTargetRemeshingEdgeLength(1);
            cell_population.SetTargetRemeshingIterations(Iter);

            cell_population.SetRemeshingSoftwear("VMTK");
            cell_population.SetRelativePath(output_dir, startime);
            cell_population.SetPrintRemeshedIC(1);
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.SetOutputMeshInVtk(true);

            VtkMeshWriter<2, 3> mesh_writer0(output_dir, "InitalOriginalMesh", false);
            mesh_writer0.WriteFilesUsingMesh(mesh);

            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {

                c_vector<double, 3> InitalLocation = cell_population.GetNode(i)->rGetLocation();

                double R = sqrt(InitalLocation[0] * InitalLocation[0] + InitalLocation[1] * InitalLocation[1]);
                double Angle;
                // double Scalled_R;
                if (InitalLocation[0] >= 0)
                {
                    Angle = atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] <= 0)
                {
                    Angle = M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] >= 0)
                {
                    Angle = -M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }

                double X = R * cos(Angle) + 2 * sin(InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 2 * sin(Angle);

                c_vector<double, 3> DeformedLocation = Create_c_vector(X, Y, InitalLocation[2]);

                cell_population.GetNode(i)->rGetModifiableLocation() = DeformedLocation;
            }

            VtkMeshWriter<2, 3> mesh_writer1(output_dir, "DeformedOriginalMesh", false);
            mesh_writer1.WriteFilesUsingMesh(mesh);

            cell_population.ExecuteHistoryDependentRemeshing();

            VtkMeshWriter<2, 3> mesh_writer2(output_dir, "DeformedRemeshedMesh", false);
            mesh_writer2.WriteFilesUsingMesh(mesh);

            std::map<unsigned, c_vector<double, 3> > InitalNodePos = cell_population.GetInitalNodePositions();
            std::map<unsigned, c_vector<double, 3> > DeformedIC;
            std::vector<double> Error;

            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {
                c_vector<double, 3> InitalLocation = InitalNodePos[cell_population.GetNode(i)->GetIndex()];

                double R = sqrt(InitalLocation[0] * InitalLocation[0] + InitalLocation[1] * InitalLocation[1]);
                double Angle;
                // double Scalled_R;
                if (InitalLocation[0] >= 0)
                {
                    Angle = atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] <= 0)
                {
                    Angle = M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] >= 0)
                {
                    Angle = -M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }

                double X = R * cos(Angle) + 2 * sin(InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 2 * sin(Angle);

                c_vector<double, 3> DeformedLocation_IC = Create_c_vector(X, Y, InitalLocation[2]);

                DeformedIC[cell_population.GetNode(i)->GetIndex()] = DeformedLocation_IC;
                c_vector<double, 3> DeformedNode = cell_population.GetNode(i)->rGetLocation();

                Error.push_back(norm_2(DeformedNode - DeformedLocation_IC));

                // Now have a deformed place to put this and the error :)
            }
            double MeanError = 0;
            for (std::vector<double>::iterator iter = Error.begin(); iter != Error.end(); ++iter)
            {
                MeanError += (*iter) * (*iter);
            }
            MeanError = sqrt(MeanError);

            double rmsE = 0;
            for (std::vector<double>::iterator iter = Error.begin(); iter != Error.end(); ++iter)
            {
                rmsE += ((*iter) - MeanError) * ((*iter) - MeanError);
            }
            rmsE = sqrt(rmsE / Error.size());
            // Loop over nodes and deform the relevant inital poition and then test error -- write something u pin maths functions for the relevant erros
            VtkMeshWriter<2, 3> mesh_writer(output_dir, "DeformedOriginalMesh", false);
            mesh_writer.WriteFilesUsingMesh(mesh);

            RootMeanSquared.push_back(MeanError);
            RootSTDSquared.push_back(rmsE);
            NumbOfNodes.push_back(mesh.GetNumNodes());

            EdgeL.push_back(EdgeLength);
            EdgeLength = EdgeLength - 0.1;

            output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder/VMTK/CourseRemeshedMesh/EdgeLength/" + mesh_size;
            // Convert the 4 meshes I need off to stl
            std::string ConvertCommand = "meshio-convert " + output_dir + "/InitalOriginalMesh.vtu " + output_dir + "/InitalOriginalMesh.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
            ConvertCommand = "meshio-convert " + output_dir + "/DeformedRemeshedMesh.vtu " + output_dir + "/DeformedRemeshedMesh.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
            ConvertCommand = "meshio-convert " + output_dir + "/DeformedOriginalMesh.vtu " + output_dir + "/DeformedOriginalMesh.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
            ConvertCommand = "meshio-convert " + output_dir + "/results_from_time_0/NewInitalConfiguration1.vtu " + output_dir + "/NewInitalConfiguration.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
        }

        std::vector<double>::iterator j = RootSTDSquared.begin();
        std::vector<double>::iterator k = EdgeL.begin();
        std::vector<double>::iterator l = NumbOfNodes.begin();

        std::ofstream f("/Users/jcrawshaw/Documents/Projects/MeshMatlab/ErrorAna_3D_WavyDef_Cylinder_CourseRemeshsed_VMTK.txt");
        for (std::vector<double>::iterator i = RootMeanSquared.begin(); i != RootMeanSquared.end(); ++j, ++k, ++l, ++i)
        {
            f << *k << " " << *i << " " << *j << " " << *l << '\n';
        }
    }

    void OffTestRemeshingCylinder_FineRemeshsed_VMTK() throw(Exception)
    {

        std::vector<double> RootMeanSquared;
        std::vector<double> RootSTDSquared;
        std::vector<double> EdgeL;
        std::vector<double> NumbOfNodes;
        double Iter = 10;
        double EdgeLength = 1;
        double startime = 0;
        while (EdgeLength > 0.01)
        {

            // I want to control the edge length, So i am just going to write a new cylinder based of the 'perfect' cylinder
            // I have in the based file here

            std::stringstream out;
            out << EdgeLength;
            std::string mesh_size = out.str();

            // Read in the new cylinder that is generated to my desire
            std::string mesh_file = "/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder/VMTK/CourseRemeshed/EdgeLength/" + mesh_size + "/InitalOriginalMesh.vtu";

            VtkMeshReader<2, 3> mesh_reader(mesh_file);
            HistoryDepMutableMesh<2, 3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            std::string output_dir = "RemeshingErrorAnal3D/Cylinder/VMTK/FineRemeshed/EdgeLength/" + mesh_size;

            // Create the cells
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

            // Create a cell population
            HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
            cell_population.SetChasteOutputDirectory(output_dir, startime);
            cell_population.SetInitialAnlgesAcrossMembrane();
            cell_population.SetTargetRemeshingEdgeLength(0.3);
            cell_population.SetTargetRemeshingIterations(Iter);
            cell_population.SetRemeshingSoftwear("VMTK");
            cell_population.SetRelativePath(output_dir, startime);
            cell_population.SetPrintRemeshedIC(1);
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.SetOutputMeshInVtk(true);

            VtkMeshWriter<2, 3> mesh_writer0(output_dir, "InitalOriginalMesh", false);
            mesh_writer0.WriteFilesUsingMesh(mesh);

            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {

                c_vector<double, 3> InitalLocation = cell_population.GetNode(i)->rGetLocation();

                double R = sqrt(InitalLocation[0] * InitalLocation[0] + InitalLocation[1] * InitalLocation[1]);
                double Angle;
                // double Scalled_R;
                if (InitalLocation[0] >= 0)
                {
                    Angle = atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] <= 0)
                {
                    Angle = M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] >= 0)
                {
                    Angle = -M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }

                double X = R * cos(Angle) + 2 * sin(InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 2 * sin(Angle);

                c_vector<double, 3> DeformedLocation = Create_c_vector(X, Y, InitalLocation[2]);

                cell_population.GetNode(i)->rGetModifiableLocation() = DeformedLocation;
            }

            VtkMeshWriter<2, 3> mesh_writer1(output_dir, "DeformedOriginalMesh", false);
            mesh_writer1.WriteFilesUsingMesh(mesh);

            cell_population.ExecuteHistoryDependentRemeshing();

            VtkMeshWriter<2, 3> mesh_writer2(output_dir, "DeformedRemeshedMesh", false);
            mesh_writer2.WriteFilesUsingMesh(mesh);

            std::map<unsigned, c_vector<double, 3> > InitalNodePos = cell_population.GetInitalNodePositions();
            std::map<unsigned, c_vector<double, 3> > DeformedIC;
            std::vector<double> Error;

            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {
                c_vector<double, 3> InitalLocation = InitalNodePos[cell_population.GetNode(i)->GetIndex()];

                double R = sqrt(InitalLocation[0] * InitalLocation[0] + InitalLocation[1] * InitalLocation[1]);
                double Angle;
                // double Scalled_R;
                if (InitalLocation[0] >= 0)
                {
                    Angle = atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] <= 0)
                {
                    Angle = M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] >= 0)
                {
                    Angle = -M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }

                double X = R * cos(Angle) + 2 * sin(InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 2 * sin(Angle);

                c_vector<double, 3> DeformedLocation_IC = Create_c_vector(X, Y, InitalLocation[2]);

                DeformedIC[cell_population.GetNode(i)->GetIndex()] = DeformedLocation_IC;
                c_vector<double, 3> DeformedNode = cell_population.GetNode(i)->rGetLocation();

                Error.push_back(norm_2(DeformedNode - DeformedLocation_IC));

                // Now have a deformed place to put this and the error :)
            }
            double MeanError = 0;
            for (std::vector<double>::iterator iter = Error.begin(); iter != Error.end(); ++iter)
            {
                MeanError += (*iter) * (*iter);
            }
            MeanError = sqrt(MeanError);

            double rmsE = 0;
            for (std::vector<double>::iterator iter = Error.begin(); iter != Error.end(); ++iter)
            {
                rmsE += ((*iter) - MeanError) * ((*iter) - MeanError);
            }
            rmsE = sqrt(rmsE / Error.size());
            // Loop over nodes and deform the relevant inital poition and then test error -- write something u pin maths functions for the relevant erros
            VtkMeshWriter<2, 3> mesh_writer(output_dir, "DeformedOriginalMesh", false);
            mesh_writer.WriteFilesUsingMesh(mesh);

            RootMeanSquared.push_back(MeanError);
            RootSTDSquared.push_back(rmsE);
            NumbOfNodes.push_back(mesh.GetNumNodes());

            EdgeL.push_back(EdgeLength);

            output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder/VMTK/FineInital/EdgeLength/" + mesh_size;
            // Convert the 4 meshes I need off to stl
            std::string ConvertCommand = "meshio-convert " + output_dir + "/InitalOriginalMesh.vtu " + output_dir + "/InitalOriginalMesh.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
            ConvertCommand = "meshio-convert " + output_dir + "/DeformedRemeshedMesh.vtu " + output_dir + "/DeformedRemeshedMesh.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
            ConvertCommand = "meshio-convert " + output_dir + "/DeformedOriginalMesh.vtu " + output_dir + "/DeformedOriginalMesh.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
            ConvertCommand = "meshio-convert " + output_dir + "/results_from_time_0/NewInitalConfiguration1.vtu " + output_dir + "/NewInitalConfiguration.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *

            EdgeLength = EdgeLength - 0.1;
        }

        std::vector<double>::iterator j = RootSTDSquared.begin();
        std::vector<double>::iterator k = EdgeL.begin();
        std::vector<double>::iterator l = NumbOfNodes.begin();

        std::ofstream f("/Users/jcrawshaw/Documents/Projects/MeshMatlab/ErrorAna_3D_WavyDef_Cylinder_FineInitial_VMTK.txt");
        for (std::vector<double>::iterator i = RootMeanSquared.begin(); i != RootMeanSquared.end(); ++j, ++k, ++l, ++i)
        {
            f << *k << " " << *i << " " << *j << " " << *l << '\n';
        }
    }

    void TestRemeshingCylinder_CourseInitial_CGAL() throw(Exception)
    {

        std::vector<double> RootMeanSquared;
        std::vector<double> RootSTDSquared;
        std::vector<double> EdgeL;
        std::vector<double> NumbOfNodes;

        double EdgeLength = 1;
        double Iter = 2;
        std::string output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder/";
        double startime = 0;
        // I want to control the edge length, So i am just going to write a new cylinder based of the 'perfect' cylinder
        // I have in the based file here

        // // Convert the .vtu into an .off
        std::string offfile = output_dir + "Cylinder.off";
        std::string vtu2offCommand = "meshio-convert " + output_dir + "config.vtu " + offfile;
        std::system(vtu2offCommand.c_str()); // system only takes char *.

        // Now excute the CGAL command to remesh the current geometry - not the input and output within this file have to be pre-set. I will explore if I can make this more neat later, should care.... dont care
        std::string CGALRemeshingCommand = "(cd  ~/Documents/CGAL-5.0.2/examples/Polygon_mesh_processing/;./isotropic_remeshing_ForChaste -input " + offfile + " -output " + output_dir + "RemeshedCylinder.off -target_edge_length " + std::to_string(0.5) + " -iterations " + std::to_string(Iter) + " )";
        std::system(CGALRemeshingCommand.c_str()); // system only takes char *

        // Now ned to convert from .off back to a .vtu
        std::string Remeshedvtu = output_dir + "CourseInital/RemeshedCylinder.vtu";
        std::string off2vtuCommand = "meshio-convert " + output_dir + "RemeshedCylinder.off " + Remeshedvtu;
        std::system(off2vtuCommand.c_str()); // system only takes char *

        while (EdgeLength > 0.05)
        {
            // Read in the new cylinder that is generated to my desire
            std::string mesh_file = "/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder/CourseInital/RemeshedCylinder.vtu";
            VtkMeshReader<2, 3> mesh_reader(mesh_file);
            HistoryDepMutableMesh<2, 3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            std::stringstream out;
            out << EdgeLength;
            std::string mesh_size = out.str();
            output_dir = "RemeshingErrorAnal3D/Cylinder/CourseInital/EdgeLength/" + mesh_size;

            // Create the cells
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

            // Create a cell population
            HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
            cell_population.SetChasteOutputDirectory(output_dir, startime);
            cell_population.SetInitialAnlgesAcrossMembrane();
            cell_population.SetTargetRemeshingEdgeLength(EdgeLength);
            cell_population.SetTargetRemeshingIterations(4);

            cell_population.SetRemeshingSoftwear("CGAL");
            cell_population.SetRelativePath(output_dir, startime);
            cell_population.SetPrintRemeshedIC(1);
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.SetOutputMeshInVtk(true);

            VtkMeshWriter<2, 3> mesh_writer0(output_dir, "InitalOriginalMesh", false);
            mesh_writer0.WriteFilesUsingMesh(mesh);

            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {

                c_vector<double, 3> InitalLocation = cell_population.GetNode(i)->rGetLocation();

                double R = sqrt(InitalLocation[0] * InitalLocation[0] + InitalLocation[1] * InitalLocation[1]);
                double Angle;
                // double Scalled_R;
                if (InitalLocation[0] >= 0)
                {
                    Angle = atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] <= 0)
                {
                    Angle = M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] >= 0)
                {
                    Angle = -M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }

                double X = R * cos(Angle) + 2 * sin(InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 2 * sin(Angle);

                c_vector<double, 3> DeformedLocation = Create_c_vector(X, Y, InitalLocation[2]);

                cell_population.GetNode(i)->rGetModifiableLocation() = DeformedLocation;
            }

            VtkMeshWriter<2, 3> mesh_writer1(output_dir, "DeformedOriginalMesh", false);
            mesh_writer1.WriteFilesUsingMesh(mesh);

            cell_population.ExecuteHistoryDependentRemeshing();

            VtkMeshWriter<2, 3> mesh_writer2(output_dir, "DeformedRemeshedMesh", false);
            mesh_writer2.WriteFilesUsingMesh(mesh);

            std::map<unsigned, c_vector<double, 3> > InitalNodePos = cell_population.GetInitalNodePositions();
            std::map<unsigned, c_vector<double, 3> > DeformedIC;
            std::vector<double> Error;

            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {
                c_vector<double, 3> InitalLocation = InitalNodePos[cell_population.GetNode(i)->GetIndex()];

                double R = sqrt(InitalLocation[0] * InitalLocation[0] + InitalLocation[1] * InitalLocation[1]);
                double Angle;
                // double Scalled_R;
                if (InitalLocation[0] >= 0)
                {
                    Angle = atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] <= 0)
                {
                    Angle = M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] >= 0)
                {
                    Angle = -M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }

                double X = R * cos(Angle) + 2 * sin(InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 2 * sin(Angle);

                c_vector<double, 3> DeformedLocation_IC = Create_c_vector(X, Y, InitalLocation[2]);

                DeformedIC[cell_population.GetNode(i)->GetIndex()] = DeformedLocation_IC;
                c_vector<double, 3> DeformedNode = cell_population.GetNode(i)->rGetLocation();

                Error.push_back(norm_2(DeformedNode - DeformedLocation_IC));

                // Now have a deformed place to put this and the error :)
            }
            double MeanError = 0;
            for (std::vector<double>::iterator iter = Error.begin(); iter != Error.end(); ++iter)
            {
                MeanError += (*iter) * (*iter);
            }
            MeanError = sqrt(MeanError);

            double rmsE = 0;
            for (std::vector<double>::iterator iter = Error.begin(); iter != Error.end(); ++iter)
            {
                rmsE += ((*iter) - MeanError) * ((*iter) - MeanError);
            }
            rmsE = sqrt(rmsE / Error.size());
            // Loop over nodes and deform the relevant inital poition and then test error -- write something u pin maths functions for the relevant erros
            VtkMeshWriter<2, 3> mesh_writer(output_dir, "DeformedOriginalMesh", false);
            mesh_writer.WriteFilesUsingMesh(mesh);

            RootMeanSquared.push_back(MeanError);
            RootSTDSquared.push_back(rmsE);
            NumbOfNodes.push_back(mesh.GetNumNodes());

            EdgeL.push_back(EdgeLength);
            EdgeLength = 0.9 * EdgeLength;

            output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder/CourseInital/EdgeLength/" + mesh_size;
            // Convert the 4 meshes I need off to stl
            std::string ConvertCommand = "meshio-convert " + output_dir + "/InitalOriginalMesh.vtu " + output_dir + "/InitalOriginalMesh.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
            ConvertCommand = "meshio-convert " + output_dir + "/DeformedRemeshedMesh.vtu " + output_dir + "/DeformedRemeshedMesh.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
            ConvertCommand = "meshio-convert " + output_dir + "/DeformedOriginalMesh.vtu " + output_dir + "/DeformedOriginalMesh.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
            ConvertCommand = "meshio-convert " + output_dir + "/results_from_time_0/NewInitalConfiguration1.vtu " + output_dir + "/NewInitalConfiguration.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *

            EdgeLength = 0.9 * EdgeLength;
        }

        std::vector<double>::iterator j = RootSTDSquared.begin();
        std::vector<double>::iterator k = EdgeL.begin();
        std::vector<double>::iterator l = NumbOfNodes.begin();

        std::ofstream f("/Users/jcrawshaw/Documents/Projects/MeshMatlab/ErrorAna_3D_WavyDef_Cylinder_CourseInitial_CGAL.txt");
        for (std::vector<double>::iterator i = RootMeanSquared.begin(); i != RootMeanSquared.end(); ++j, ++k, ++l, ++i)
        {
            f << *k << " " << *i << " " << *j << " " << *l << '\n';
        }
    }

    void TestRemeshingCylinder_FineInitial_CGAL() throw(Exception)
    {

        std::vector<double> RootMeanSquared;
        std::vector<double> RootSTDSquared;
        std::vector<double> EdgeL;
        std::vector<double> NumbOfNodes;

        double EdgeLength = 0.05;
        std::string output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder/";
        double startime = 0;
        // I want to control the edge length, So i am just going to write a new cylinder based of the 'perfect' cylinder
        // I have in the based file here

        // // Convert the .vtu into an .off
        std::string offfile = output_dir + "Cylinder.off";
        std::string vtu2offCommand = "meshio-convert " + output_dir + "config.vtu " + offfile;
        std::system(vtu2offCommand.c_str()); // system only takes char *.

        double Iter = 2;
        // Now excute the CGAL command to remesh the current geometry - not the input and output within this file have to be pre-set. I will explore if I can make this more neat later, should care.... dont care
        std::string CGALRemeshingCommand = "(cd  ~/Documents/CGAL-5.0.2/examples/Polygon_mesh_processing/;./isotropic_remeshing_ForChaste -input " + offfile + " -output " + output_dir + "RemeshedCylinder.off -target_edge_length " + std::to_string(0.01) + " -iterations " + std::to_string(Iter) + " )";
        std::system(CGALRemeshingCommand.c_str()); // system only takes char *

        // Now ned to convert from .off back to a .vtu
        std::string Remeshedvtu = output_dir + "FineInital/RemeshedCylinder.vtu";
        std::string off2vtuCommand = "meshio-convert " + output_dir + "RemeshedCylinder.off " + Remeshedvtu;
        std::system(off2vtuCommand.c_str()); // system only takes char *

        while (EdgeLength > 0.001)
        {
            // Read in the new cylinder that is generated to my desire
            std::string mesh_file = "/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder/FineInital/RemeshedCylinder.vtu";
            VtkMeshReader<2, 3> mesh_reader(mesh_file);
            HistoryDepMutableMesh<2, 3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            std::stringstream out;
            out << EdgeLength;
            std::string mesh_size = out.str();
            output_dir = "RemeshingErrorAnal3D/Cylinder/CourseInital/EdgeLength/" + mesh_size;

            // Create the cells
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

            // Create a cell population
            HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
            cell_population.SetChasteOutputDirectory(output_dir, startime);
            cell_population.SetInitialAnlgesAcrossMembrane();
            cell_population.SetRemeshingSoftwear("CGAL");
            cell_population.SetTargetRemeshingEdgeLength(EdgeLength);
            cell_population.SetRelativePath(output_dir, startime);
            cell_population.SetPrintRemeshedIC(1);
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.SetOutputMeshInVtk(true);

            VtkMeshWriter<2, 3> mesh_writer0(output_dir, "InitalOriginalMesh", false);
            mesh_writer0.WriteFilesUsingMesh(mesh);

            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {

                c_vector<double, 3> InitalLocation = cell_population.GetNode(i)->rGetLocation();

                double R = sqrt(InitalLocation[0] * InitalLocation[0] + InitalLocation[1] * InitalLocation[1]);
                double Angle;
                // double Scalled_R;
                if (InitalLocation[0] >= 0)
                {
                    Angle = atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] <= 0)
                {
                    Angle = M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] >= 0)
                {
                    Angle = -M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }

                double X = R * cos(Angle) + 2 * sin(InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 2 * sin(Angle);

                c_vector<double, 3> DeformedLocation = Create_c_vector(X, Y, InitalLocation[2]);

                cell_population.GetNode(i)->rGetModifiableLocation() = DeformedLocation;
            }

            VtkMeshWriter<2, 3> mesh_writer1(output_dir, "DeformedOriginalMesh", false);
            mesh_writer1.WriteFilesUsingMesh(mesh);

            cell_population.ExecuteHistoryDependentRemeshing();

            VtkMeshWriter<2, 3> mesh_writer2(output_dir, "DeformedRemeshedMesh", false);
            mesh_writer2.WriteFilesUsingMesh(mesh);

            std::map<unsigned, c_vector<double, 3> > InitalNodePos = cell_population.GetInitalNodePositions();
            std::map<unsigned, c_vector<double, 3> > DeformedIC;
            std::vector<double> Error;

            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {
                c_vector<double, 3> InitalLocation = InitalNodePos[cell_population.GetNode(i)->GetIndex()];

                double R = sqrt(InitalLocation[0] * InitalLocation[0] + InitalLocation[1] * InitalLocation[1]);
                double Angle;
                // double Scalled_R;
                if (InitalLocation[0] >= 0)
                {
                    Angle = atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] <= 0)
                {
                    Angle = M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] >= 0)
                {
                    Angle = -M_PI + atan(InitalLocation[1] / InitalLocation[0]);
                }

                double X = R * cos(Angle) + 2 * sin(InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 2 * sin(Angle);

                c_vector<double, 3> DeformedLocation_IC = Create_c_vector(X, Y, InitalLocation[2]);

                DeformedIC[cell_population.GetNode(i)->GetIndex()] = DeformedLocation_IC;
                c_vector<double, 3> DeformedNode = cell_population.GetNode(i)->rGetLocation();

                Error.push_back(norm_2(DeformedNode - DeformedLocation_IC));

                // Now have a deformed place to put this and the error :)
            }
            double MeanError = 0;
            for (std::vector<double>::iterator iter = Error.begin(); iter != Error.end(); ++iter)
            {
                MeanError += (*iter) * (*iter);
            }
            MeanError = sqrt(MeanError);

            double rmsE = 0;
            for (std::vector<double>::iterator iter = Error.begin(); iter != Error.end(); ++iter)
            {
                rmsE += ((*iter) - MeanError) * ((*iter) - MeanError);
            }
            rmsE = sqrt(rmsE / Error.size());
            // Loop over nodes and deform the relevant inital poition and then test error -- write something u pin maths functions for the relevant erros
            VtkMeshWriter<2, 3> mesh_writer(output_dir, "DeformedOriginalMesh", false);
            mesh_writer.WriteFilesUsingMesh(mesh);

            RootMeanSquared.push_back(MeanError);
            RootSTDSquared.push_back(rmsE);
            NumbOfNodes.push_back(mesh.GetNumNodes());

            EdgeL.push_back(EdgeLength);
            EdgeLength = 0.9 * EdgeLength;

            output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder/FineInital/EdgeLength/" + mesh_size;
            // Convert the 4 meshes I need off to stl
            std::string ConvertCommand = "meshio-convert " + output_dir + "/InitalOriginalMesh.vtu " + output_dir + "/InitalOriginalMesh.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
            ConvertCommand = "meshio-convert " + output_dir + "/DeformedRemeshedMesh.vtu " + output_dir + "/DeformedRemeshedMesh.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
            ConvertCommand = "meshio-convert " + output_dir + "/DeformedOriginalMesh.vtu " + output_dir + "/DeformedOriginalMesh.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *
            ConvertCommand = "meshio-convert " + output_dir + "/results_from_time_0/NewInitalConfiguration1.vtu " + output_dir + "/NewInitalConfiguration.stl";
            std::system(ConvertCommand.c_str()); // system only takes char *

            EdgeLength = 0.9 * EdgeLength;
        }

        std::vector<double>::iterator j = RootSTDSquared.begin();
        std::vector<double>::iterator k = EdgeL.begin();
        std::vector<double>::iterator l = NumbOfNodes.begin();

        std::ofstream f("/Users/jcrawshaw/Documents/Projects/MeshMatlab/ErrorAna_3D_WavyDef_Cylinder_FineInitial_CGAL.txt");
        for (std::vector<double>::iterator i = RootMeanSquared.begin(); i != RootMeanSquared.end(); ++j, ++k, ++l, ++i)
        {
            f << *k << " " << *i << " " << *j << " " << *l << '\n';
        }
    }
};

//        boundary_plane_points.push_back(Create_c_vector(0.523,0.734,-2.25));
//         boundary_plane_normals.push_back(Create_c_vector(0.999,-0.03557,0.0024));

//         boundary_plane_points.push_back(Create_c_vector(0.5971,0.45457,0.00166));
//         boundary_plane_normals.push_back(Create_c_vector(0.75169,0.659,-0.0038));

//         boundary_plane_points.push_back(Create_c_vector(0.78357,0.348,-0.00654));
//         boundary_plane_normals.push_back(Create_c_vector(0.1423,0.9849,-0.098));

//         boundary_plane_points.push_back(Create_c_vector(1.0688,0.40478,0.0029 ));
//         boundary_plane_normals.push_back(Create_c_vector( -0.72739,0.6858,-0.0218 ));

//         boundary_plane_points.push_back(Create_c_vector(1.08,0.732,0.004458));
//         boundary_plane_normals.push_back(Create_c_vector(-0.8,-0.59,0.0647));

//         boundary_plane_points.push_back(Create_c_vector(0.90644,1.00186,0.0029));
//         boundary_plane_normals.push_back(Create_c_vector(-0.96445,-0.5025,0.013317));
// // ---
//         boundary_plane_points.push_back(Create_c_vector(0.729,1.016095,0.00905));
//         boundary_plane_normals.push_back(Create_c_vector(0.6199170377,-0.78,0.04815));

#endif /*TESTRELAXATION_HPP_*/
