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

        double Length = 6;
        double Radius = 1;

        Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
        MutableMesh<2, 3>* p_mesh = generator.GetMesh();

        std::string output_directory = "RemeshingErrorAnal3D/Cylinder_WaveyDef/VMTK";
        VtkMeshWriter<2, 3> mesh_writer0(output_directory, "Cylinder", false);
        mesh_writer0.WriteFilesUsingMesh(*p_mesh);
    }

    void TestRemeshingCylinder_CourseInitial_VMTK() throw(Exception)
    {

        std::vector<double> RootMeanSquared;
        std::vector<double> RootSTDSquared;
        std::vector<double> EdgeL;
        std::vector<double> NumbOfNodes;

        double EdgeLength = 0.9;
        double Iter = 10;
        std::string output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/VMTK/";
        double startime = 0;
        // I want to control the edge length, So i am just going to write a new cylinder based of the 'perfect' cylinder
        // I have in the based file here

        // Convert the .vtu into an .stl
        std::string stlfile = output_dir + "Cylinder.stl";

        std::string vtu2offCommand = "meshio-convert " + output_dir + "Cylinder.vtu " + stlfile;
        std::system(vtu2offCommand.c_str()); // system only takes char *.

        std::string Remeshedstl = output_dir + "RemeshedCylinder.stl";
        std::string RemeshCommand = "vmtksurfaceremeshing -ifile " + stlfile + " -iterations " + std::to_string(20) + " -edgelength "+ std::to_string(0.9) +" -ofile " + Remeshedstl + " -elementsizemode edgelength";
        std::system(RemeshCommand.c_str());

        std::string Remeshedvtu = output_dir + "CourseInital/RemeshedCylinder.vtu";
        std::string stl2vtuCommand = " meshio-convert " + Remeshedstl + " " + Remeshedvtu;
        std::system(stl2vtuCommand.c_str());

        while (EdgeLength > 0.01)
        {
            // Read in the new cylinder that is generated to my desire
            std::string mesh_file = "/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/VMTK/CourseInital/RemeshedCylinder.vtu";
            VtkMeshReader<2, 3> mesh_reader(mesh_file);
            HistoryDepMutableMesh<2, 3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            std::stringstream out;
            out << EdgeLength;
            std::string mesh_size = out.str();
            output_dir = "RemeshingErrorAnal3D/Cylinder_WaveyDef/VMTK/CourseInital/EdgeLength/" + mesh_size;

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

                double X = R * cos(Angle) + sin(1.5*InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 1.5 * sin(Angle);

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

                double X = R * cos(Angle) + sin(1.5*InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 1.5 * sin(Angle);

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
            MeanError = sqrt(MeanError/ Error.size());

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

            output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/VMTK/CourseInital/EdgeLength/" + mesh_size;
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

        std::ofstream f("/Users/jcrawshaw/Documents/Projects/MeshMatlab/ErrorAnal_Cylinder_WaveyDef_LessCourseInitial.txt");
        for (std::vector<double>::iterator i = RootMeanSquared.begin(); i != RootMeanSquared.end(); ++j, ++k, ++l, ++i)
        {
            f << *k << " " << *i << " " << *j << " " << *l << '\n';
        }
    }

    void TestRemeshingCylinder_FineInitial_VMTK() throw(Exception)
    {
        std::vector<double> RootMeanSquared;
        std::vector<double> RootSTDSquared;
        std::vector<double> EdgeL;
        std::vector<double> NumbOfNodes;

        double EdgeLength = 2;
        double Iter = 10;
        std::string output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/VMTK/";
        double startime = 0;
        // I want to control the edge length, So i am just going to write a new cylinder based of the 'perfect' cylinder
        // I have in the based file here

        // Convert the .vtu into an .stl
        std::string stlfile = output_dir + "Cylinder.stl";

        std::string vtu2offCommand = "meshio-convert " + output_dir + "Cylinder.vtu " + stlfile;
        std::system(vtu2offCommand.c_str()); // system only takes char *.

        std::string Remeshedstl = output_dir + "RemeshedCylinder.stl";
        std::string RemeshCommand = "vmtksurfaceremeshing -ifile " + stlfile + " -iterations " + std::to_string(20) + " -edgelength "+ std::to_string(0.1) +" -ofile " + Remeshedstl + " -elementsizemode edgelength";
        std::system(RemeshCommand.c_str());

        std::string Remeshedvtu = output_dir + "FineInital/RemeshedCylinder.vtu";
        std::string stl2vtuCommand = " meshio-convert " + Remeshedstl + " " + Remeshedvtu;
        std::system(stl2vtuCommand.c_str());

        while (EdgeLength > 0.01)
        {
            // Read in the new cylinder that is generated to my desire
            std::string mesh_file = "/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/VMTK/FineInital/RemeshedCylinder.vtu";
            VtkMeshReader<2, 3> mesh_reader(mesh_file);
            HistoryDepMutableMesh<2, 3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            std::stringstream out;
            out << EdgeLength;
            std::string mesh_size = out.str();
            output_dir = "RemeshingErrorAnal3D/Cylinder_WaveyDef/VMTK/FineInital/EdgeLength/" + mesh_size;

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

                double X = R * cos(Angle) + sin(1.5*InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 1.5 * sin(Angle);

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

                double X = R * cos(Angle) + sin(1.5*InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 1.5 * sin(Angle);

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
            MeanError = sqrt(MeanError/ Error.size());

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

            output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/VMTK/FineInital/EdgeLength/" + mesh_size;
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




        std::ofstream f("/Users/jcrawshaw/Documents/Projects/MeshMatlab/ErrorAnal_Cylinder_WaveyDef_FineInitial2.txt");
        for (std::vector<double>::iterator i = RootMeanSquared.begin(); i != RootMeanSquared.end(); ++j, ++k, ++l, ++i)
        {
            f << *k << " " << *i << " " << *j << " " << *l << '\n';
        }
    }

    void TestRemeshingCylinder_CourseRemeshedMesh_VMTK() throw(Exception)
    {
        std::vector<double> RootMeanSquared;
        std::vector<double> RootSTDSquared;
        std::vector<double> EdgeL;
        std::vector<double> NumbOfNodes;
        double EdgeLength = 2;
        double Iter = 10;
        double startime = 0;
        while (EdgeLength > 0.01)
        {
            std::stringstream out;
            out << EdgeLength;
            std::string mesh_size = out.str();
            std::string output_dir;// = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/VMTK/CourseRemeshed/";

            // Read in the new cylinder that is generated to my desire
            std::string mesh_file = "/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/VMTK/CourseRemeshed/EdgeLength/" + mesh_size + "/InitalOriginalMesh.vtu";

            VtkMeshReader<2, 3> mesh_reader(mesh_file);
            HistoryDepMutableMesh<2, 3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            output_dir = "RemeshingErrorAnal3D/Cylinder_WaveyDef/VMTK/CourseRemeshed/EdgeLength/" + mesh_size;

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

                double X = R * cos(Angle) + sin(1.5*InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 1.5 * sin(Angle);

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

                double X = R * cos(Angle) + sin(1.5*InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 1.5 * sin(Angle);

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
            MeanError = sqrt(MeanError/ Error.size());

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

            output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/VMTK/CourseRemeshed/EdgeLength/" + mesh_size;
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

        std::ofstream f("/Users/jcrawshaw/Documents/Projects/MeshMatlab/ErrorAnal_Cylinder_WaveyDef_CourseRemeshsed2.txt");;
        for (std::vector<double>::iterator i = RootMeanSquared.begin(); i != RootMeanSquared.end(); ++j, ++k, ++l, ++i)
        {
            f << *k << " " << *i << " " << *j << " " << *l << '\n';
        }
    }

    void TestRemeshingCylinder_FineRemeshsed_VMTK() throw(Exception)
    {

        std::vector<double> RootMeanSquared;
        std::vector<double> RootSTDSquared;
        std::vector<double> EdgeL;
        std::vector<double> NumbOfNodes;
        double Iter = 10;
        double EdgeLength = 2;
        double startime = 0;
        while (EdgeLength > 0.01)
        {

            // I want to control the edge length, So i am just going to write a new cylinder based of the 'perfect' cylinder
            // I have in the based file here

            std::stringstream out;
            out << EdgeLength;
            std::string mesh_size = out.str();

            // Read in the new cylinder that is generated to my desire
            std::string mesh_file = "/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/VMTK/CourseRemeshed/EdgeLength/" + mesh_size + "/InitalOriginalMesh.vtu";

            VtkMeshReader<2, 3> mesh_reader(mesh_file);
            HistoryDepMutableMesh<2, 3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            std::string output_dir = "RemeshingErrorAnal3D/Cylinder_WaveyDef/VMTK/FineRemeshed/EdgeLength/" + mesh_size;

            // Create the cells
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

            // Create a cell population
            HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
            cell_population.SetChasteOutputDirectory(output_dir, startime);
            cell_population.SetInitialAnlgesAcrossMembrane();
            cell_population.SetTargetRemeshingEdgeLength(0.1);
            cell_population.SetTargetRemeshingIterations(Iter);
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

                double X = R * cos(Angle) + sin(1.5*InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 1.5 * sin(Angle);

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

                double X = R * cos(Angle) + sin(1.5*InitalLocation[2]); //+ InitalLocation[2]-1;
                double Y = R * 1.5 * sin(Angle);

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
            MeanError = sqrt(MeanError/ Error.size());

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

            output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/VMTK/FineInital/EdgeLength/" + mesh_size;
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

        std::ofstream f("/Users/jcrawshaw/Documents/Projects/MeshMatlab/ErrorAnal_Cylinder_WaveyDef_FineRemeshsed2.txt");
        for (std::vector<double>::iterator i = RootMeanSquared.begin(); i != RootMeanSquared.end(); ++j, ++k, ++l, ++i)
        {
            f << *k << " " << *i << " " << *j << " " << *l << '\n';
        }
    }
};

#endif /*TESTRELAXATION_HPP_*/
