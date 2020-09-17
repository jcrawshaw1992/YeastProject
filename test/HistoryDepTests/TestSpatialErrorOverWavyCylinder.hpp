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

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:

    void offTestCreateCylinder() throw(Exception)
    {
        unsigned N_D = 800;
        unsigned N_Z = 1000;

        double Length = 6;
        double Radius = 1;

        Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
        MutableMesh<2, 3>* p_mesh = generator.GetMesh();

        std::string output_directory = "RemeshingErrorAnal3D/Cylinder_WaveyDef/";
        VtkMeshWriter<2, 3> mesh_writer0(output_directory, "Cylinder", false);
        mesh_writer0.WriteFilesUsingMesh(*p_mesh);

        // Convert the .vtu
        std::string vtufile = "/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/Cylinder.vtu";
        std::string offfile = "/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/Cylinder.off";
        std::string stlfile = "/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/Cylinder.stl";

        std::string off2offCommand = "meshio-convert " + vtufile + offfile;
        std::system(off2offCommand.c_str()); // system only takes char *.

        std::string vtu2offCommand = "meshio-convert " + vtufile + stlfile;
        std::system(vtu2offCommand.c_str()); // system only takes char *.

        std::string MakeDirectory = "mkdir /Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/CoarseInitial";
        int SystemOutput = std::system(MakeDirectory.c_str()); // system only takes char *

        MakeDirectory = "mkdir /Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/FineInitial";
        SystemOutput = std::system(MakeDirectory.c_str()); // system only takes char *

        MakeDirectory = "mkdir /Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/CoarseRemeshed";
        SystemOutput = std::system(MakeDirectory.c_str()); // system only takes char *

        MakeDirectory = "mkdir /Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/FineRemeshed";
        SystemOutput = std::system(MakeDirectory.c_str()); // system only takes char *
    }

    void offTestRemeshingCylinder_CoarseInitial() throw(Exception)
    {

        std::vector<double> EdgeL;
        std::vector<double> Median;
        std::vector<double> Q25;
        std::vector<double> Q75;

        double EdgeLength = 0.4;
        double Iter = 10;
        std::string output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/CoarseInitial/";
        double startime = 0;
        // I want to control the edge length, So i am just going to write a new cylinder based of the 'perfect' cylinder
        // I have in the based file here
        std::string offfile = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/Cylinder.off";
        std::string CGALRemeshingCommand = "/Users/jcrawshaw/Documents/CGAL-5.0.2/examples/Polygon_mesh_processing/isotropic_remeshing_ForChaste -input " + offfile + " -output " + output_dir + "RemeshedCylinder.off -target_edge_length " + std::to_string(EdgeLength) + " -iterations " + std::to_string(10);
        int SystemOutput = std::system(CGALRemeshingCommand.c_str()); // system only takes char *

        std::string Remeshedvtu = output_dir + "RemeshedCylinder.vtu";
        std::string stl2vtuCommand = " meshio-convert " + output_dir + "RemeshedCylinder.off " + Remeshedvtu;
        std::system(stl2vtuCommand.c_str());
        EdgeLength = 0.5;
        while (EdgeLength > 0.03)
        {

            // Read in the new cylinder that is generated to my desire
            std::string mesh_file = "/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/CoarseInitial/RemeshedCylinder.vtu"; //Remeshedvtu;//"/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/CGAL/CourseInital/RemeshedCylinder.vtu";
            VtkMeshReader<2, 3> mesh_reader(mesh_file);
            HistoryDepMutableMesh<2, 3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            std::stringstream out;
            out << EdgeLength;
            std::string mesh_size = out.str();
            std::string output_directory = "RemeshingErrorAnal3D/Cylinder_WaveyDef/CoarseInitial/EdgeLength/" + mesh_size;

            // Create the cells
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

            // Create a cell population
            HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
            cell_population.SetChasteOutputDirectory(output_directory, startime);
            cell_population.SetInitialAnlgesAcrossMembrane();
            cell_population.SetTargetRemeshingEdgeLength(EdgeLength);
            cell_population.SetTargetRemeshingIterations(Iter);

            cell_population.SetRemeshingSoftwear("CGAL");
            cell_population.SetRelativePath(output_directory, startime);
            cell_population.SetPrintRemeshedIC(1);
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.SetOutputMeshInVtk(true);

            VtkMeshWriter<2, 3> mesh_writer0(output_directory, "InitalMesh", false);
            mesh_writer0.WriteFilesUsingMesh(mesh);
            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {
                c_vector<double, 3> InitalLocation = cell_population.GetNode(i)->rGetLocation();
                c_vector<double, 3> DeformedLocation = WavyDeformation(InitalLocation);
                cell_population.GetNode(i)->rGetModifiableLocation() = DeformedLocation;
            }

            VtkMeshWriter<2, 3> mesh_writer1(output_directory, "DeformedOriginalMesh", false);
            mesh_writer1.WriteFilesUsingMesh(mesh);

            // Remeshing should happen here --> This is where all the time is taken
            TRACE("Execute remeshing")
            cell_population.ExecuteHistoryDependentRemeshing();
            TRACE("Remeshing complete")
            VtkMeshWriter<2, 3> mesh_writer2(output_directory, "DeformedRemeshedMesh", false);
            mesh_writer2.WriteFilesUsingMesh(mesh);

            std::map<unsigned, c_vector<double, 3> > InitalNodePos = cell_population.GetInitalNodePositions();
            std::map<unsigned, c_vector<double, 3> > DeformedIC;
            std::vector<double> Error;

            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {
                c_vector<double, 3> InitalLocation = InitalNodePos[cell_population.GetNode(i)->GetIndex()];
                c_vector<double, 3> DeformedLocation_IC = WavyDeformation(InitalLocation);
                DeformedIC[cell_population.GetNode(i)->GetIndex()] = DeformedLocation_IC;
                c_vector<double, 3> DeformedNode = cell_population.GetNode(i)->rGetLocation();
                Error.push_back(norm_2(DeformedNode - DeformedLocation_IC));
            }

            VtkMeshWriter<2, 3> mesh_writer(output_directory, "DeformedOriginalMesh", false); // I think this is falsely labeled!
            mesh_writer.WriteFilesUsingMesh(mesh);

            std::vector<double> quartiles = Quantile<double>(Error, { 0.25, 0.5, 0.75 });
            typename std::vector<double>::iterator Iterator = quartiles.begin();
            double Quartile1 = *Iterator;
            std::advance(Iterator, 1);
            double Med = *Iterator;
            std::advance(Iterator, 1);
            double Quartile3 = *Iterator;

            Median.push_back(Med);
            Q25.push_back(Quartile1);
            Q75.push_back(Quartile3);

            EdgeL.push_back(EdgeLength);
            EdgeLength = EdgeLength - 0.05;
        }

        std::vector<double>::iterator a = Median.begin();
        std::vector<double>::iterator b = Q25.begin();
        std::vector<double>::iterator c = Q75.begin();
        std::ofstream f("/Users/jcrawshaw/Documents/Projects/MeshMatlab/Remeshing3D/RawSpatialErrorData/SpatialErrorAnal_3D_CoarseInitial.txt");
        for (std::vector<double>::iterator i = EdgeL.begin(); i != EdgeL.end(); ++a, ++b, ++c, ++i)
        {
            f << *i << " " << *b << " " << *a << " " << *c << '\n';
        }
    }

    void offTestRemeshingCylinder_FineInitial() throw(Exception)
    {

        std::vector<double> EdgeL;
        std::vector<double> Median;
        std::vector<double> Q25;
        std::vector<double> Q75;

        double EdgeLength = 0.05;
        double Iter = 10;
        std::string output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/FineInitial/";
        double startime = 0;
        // I want to control the edge length, So i am just going to write a new cylinder based of the 'perfect' cylinder
        // I have in the based file here
        std::string offfile = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/Cylinder.off";
        std::string CGALRemeshingCommand = "/Users/jcrawshaw/Documents/CGAL-5.0.2/examples/Polygon_mesh_processing/isotropic_remeshing_ForChaste -input " + offfile + " -output " + output_dir + "RemeshedCylinder.off -target_edge_length " + std::to_string(EdgeLength) + " -iterations " + std::to_string(10);
        int SystemOutput = std::system(CGALRemeshingCommand.c_str()); // system only takes char *

        std::string Remeshedvtu = output_dir + "RemeshedCylinder.vtu";
        std::string stl2vtuCommand = " meshio-convert " + output_dir + "RemeshedCylinder.off " + Remeshedvtu;
        std::system(stl2vtuCommand.c_str());
        EdgeLength = 0.5;

        while (EdgeLength > 0.03)
        {

            // Read in the new cylinder that is generated to my desire
            std::string mesh_file = "/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/FineInitial/RemeshedCylinder.vtu"; //Remeshedvtu;//"/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/CGAL/CourseInital/RemeshedCylinder.vtu";
            VtkMeshReader<2, 3> mesh_reader(mesh_file);
            HistoryDepMutableMesh<2, 3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            std::stringstream out;
            out << EdgeLength;
            std::string mesh_size = out.str();
            std::string output_directory = "RemeshingErrorAnal3D/Cylinder_WaveyDef/FineInitial/EdgeLength/" + mesh_size;

            // Create the cells
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

            // Create a cell population
            HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
            cell_population.SetChasteOutputDirectory(output_directory, startime);
            cell_population.SetInitialAnlgesAcrossMembrane();
            cell_population.SetTargetRemeshingEdgeLength(EdgeLength);
            cell_population.SetTargetRemeshingIterations(Iter);

            cell_population.SetRemeshingSoftwear("CGAL");
            cell_population.SetRelativePath(output_directory, startime);
            cell_population.SetPrintRemeshedIC(1);
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.SetOutputMeshInVtk(true);

            VtkMeshWriter<2, 3> mesh_writer0(output_directory, "InitalMesh", false);
            mesh_writer0.WriteFilesUsingMesh(mesh);
            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {
                c_vector<double, 3> InitalLocation = cell_population.GetNode(i)->rGetLocation();
                c_vector<double, 3> DeformedLocation = WavyDeformation(InitalLocation);
                cell_population.GetNode(i)->rGetModifiableLocation() = DeformedLocation;
            }

            VtkMeshWriter<2, 3> mesh_writer1(output_directory, "DeformedOriginalMesh", false);
            mesh_writer1.WriteFilesUsingMesh(mesh);

            // Remeshing should happen here --> This is where all the time is taken
            TRACE("Execute remeshing")
            cell_population.ExecuteHistoryDependentRemeshing();
            TRACE("Remeshing complete")
            VtkMeshWriter<2, 3> mesh_writer2(output_directory, "DeformedRemeshedMesh", false);
            mesh_writer2.WriteFilesUsingMesh(mesh);

            std::map<unsigned, c_vector<double, 3> > InitalNodePos = cell_population.GetInitalNodePositions();
            std::map<unsigned, c_vector<double, 3> > DeformedIC;
            std::vector<double> Error;

            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {
                c_vector<double, 3> InitalLocation = InitalNodePos[cell_population.GetNode(i)->GetIndex()];
                c_vector<double, 3> DeformedLocation_IC = WavyDeformation(InitalLocation);
                DeformedIC[cell_population.GetNode(i)->GetIndex()] = DeformedLocation_IC;
                c_vector<double, 3> DeformedNode = cell_population.GetNode(i)->rGetLocation();
                Error.push_back(norm_2(DeformedNode - DeformedLocation_IC));
            }

            VtkMeshWriter<2, 3> mesh_writer(output_directory, "DeformedOriginalMesh", false); // I think this is falsely labeled!
            mesh_writer.WriteFilesUsingMesh(mesh);

            std::vector<double> quartiles = Quantile<double>(Error, { 0.25, 0.5, 0.75 });
            typename std::vector<double>::iterator Iterator = quartiles.begin();
            double Quartile1 = *Iterator;
            std::advance(Iterator, 1);
            double Med = *Iterator;
            std::advance(Iterator, 1);
            double Quartile3 = *Iterator;

            Median.push_back(Med);
            Q25.push_back(Quartile1);
            Q75.push_back(Quartile3);

            EdgeL.push_back(EdgeLength);
            EdgeLength = EdgeLength - 0.05;
        }

        std::vector<double>::iterator a = Median.begin();
        std::vector<double>::iterator b = Q25.begin();
        std::vector<double>::iterator c = Q75.begin();
         
        std::ofstream f("/Users/jcrawshaw/Documents/Projects/MeshMatlab/Remeshing3D/RawSpatialErrorData/SpatialErrorAnal_3D_FineInitial.txt");
        for (std::vector<double>::iterator i = EdgeL.begin(); i != EdgeL.end(); ++a, ++b, ++c, ++i)
        {
            f << *i << " " << *b << " " << *a << " " << *c << '\n';
        }
    }

    void TestRemeshingCylinder_FineRemeshed() throw(Exception)
    {

        std::vector<double> EdgeL;
        std::vector<double> Median;
        std::vector<double> Q25;
        std::vector<double> Q75;

        double EdgeLength = 0.5;
        double Iter = 10;
        std::string output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/FineRemeshed/";
        double startime = 0;

        while (EdgeLength > 0.01)
        {
            std::stringstream out;
            out << EdgeLength;
            std::string mesh_size = out.str();

            std::string offfile = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/Cylinder.off";
            std::string CGALRemeshingCommand = "/Users/jcrawshaw/Documents/CGAL-5.0.2/examples/Polygon_mesh_processing/isotropic_remeshing_ForChaste -input " + offfile + " -output " + output_dir + "RemeshedCylinder.off -target_edge_length " + std::to_string(EdgeLength) + " -iterations " + std::to_string(10);
            int SystemOutput = std::system(CGALRemeshingCommand.c_str()); // system only takes char *

            std::string Remeshedvtu = output_dir + "RemeshedCylinder" + mesh_size + ".vtu";
            std::string stl2vtuCommand = " meshio-convert " + output_dir + "RemeshedCylinder.off " + Remeshedvtu;
            std::system(stl2vtuCommand.c_str());
            double RemeshedEdgeLength = 0.05;

            // Read in the new cylinder that is generated to my desire
            std::string mesh_file = "/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/FineRemeshed/RemeshedCylinder" + mesh_size + ".vtu"; //Remeshedvtu;//"/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/CGAL/CourseInital/RemeshedCylinder.vtu";
            VtkMeshReader<2, 3> mesh_reader(mesh_file);
            HistoryDepMutableMesh<2, 3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            std::string output_directory = "RemeshingErrorAnal3D/Cylinder_WaveyDef/FineRemeshed/EdgeLength/" + mesh_size + "/";

            // Create the cells
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

            // Create a cell population
            HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
            cell_population.SetChasteOutputDirectory(output_directory, startime);
            cell_population.SetInitialAnlgesAcrossMembrane();
            cell_population.SetTargetRemeshingEdgeLength(RemeshedEdgeLength);
            cell_population.SetTargetRemeshingIterations(Iter);

            cell_population.SetRemeshingSoftwear("CGAL");
            cell_population.SetRelativePath(output_directory, startime);
            cell_population.SetPrintRemeshedIC(1);
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.SetOutputMeshInVtk(true);

            VtkMeshWriter<2, 3> mesh_writer0(output_directory, "InitalMesh", false);
            mesh_writer0.WriteFilesUsingMesh(mesh);
            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {
                c_vector<double, 3> InitalLocation = cell_population.GetNode(i)->rGetLocation();
                c_vector<double, 3> DeformedLocation = WavyDeformation(InitalLocation);
                cell_population.GetNode(i)->rGetModifiableLocation() = DeformedLocation;
            }

            VtkMeshWriter<2, 3> mesh_writer1(output_directory, "DeformedOriginalMesh", false);
            mesh_writer1.WriteFilesUsingMesh(mesh);

            // Remeshing should happen here --> This is where all the time is taken
            TRACE("Execute remeshing")
            cell_population.ExecuteHistoryDependentRemeshing();
            TRACE("Remeshing complete")
            VtkMeshWriter<2, 3> mesh_writer2(output_directory, "DeformedRemeshedMesh", false);
            mesh_writer2.WriteFilesUsingMesh(mesh);

            std::map<unsigned, c_vector<double, 3> > InitalNodePos = cell_population.GetInitalNodePositions();
            std::map<unsigned, c_vector<double, 3> > DeformedIC;
            std::vector<double> Error;

            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {
                c_vector<double, 3> InitalLocation = InitalNodePos[cell_population.GetNode(i)->GetIndex()];
                c_vector<double, 3> DeformedLocation_IC = WavyDeformation(InitalLocation);
                DeformedIC[cell_population.GetNode(i)->GetIndex()] = DeformedLocation_IC;
                c_vector<double, 3> DeformedNode = cell_population.GetNode(i)->rGetLocation();
                Error.push_back(norm_2(DeformedNode - DeformedLocation_IC));
            }

            VtkMeshWriter<2, 3> mesh_writer(output_directory, "DeformedOriginalMesh", false); // I think this is falsely labeled!
            mesh_writer.WriteFilesUsingMesh(mesh);

            std::vector<double> quartiles = Quantile<double>(Error, { 0.25, 0.5, 0.75 });
            typename std::vector<double>::iterator Iterator = quartiles.begin();
            double Quartile1 = *Iterator;
            std::advance(Iterator, 1);
            double Med = *Iterator;
            std::advance(Iterator, 1);
            double Quartile3 = *Iterator;

            Median.push_back(Med);
            Q25.push_back(Quartile1);
            Q75.push_back(Quartile3);

            EdgeL.push_back(EdgeLength);
            EdgeLength = EdgeLength - 0.05;
        }

        std::vector<double>::iterator a = Median.begin();
        std::vector<double>::iterator b = Q25.begin();
        std::vector<double>::iterator c = Q75.begin();
                                                                                                        
        std::ofstream f("/Users/jcrawshaw/Documents/Projects/MeshMatlab/Remeshing3D/RawSpatialErrorData/SpatialErrorAnal_3D_FineRemeshed.txt");
        for (std::vector<double>::iterator i = EdgeL.begin(); i != EdgeL.end(); ++a, ++b, ++c, ++i)
        {
            f << *i << " " << *b << " " << *a << " " << *c << '\n';
        }
    }

    void TestRemeshingCylinder_CourseRemeshed() throw(Exception)
    {

        std::vector<double> EdgeL;
        std::vector<double> Median;
        std::vector<double> Q25;
        std::vector<double> Q75;

        double EdgeLength = 0.5; double RemeshedEdgeLength = 0.4;
        double Iter = 10;
        std::string output_dir = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/CoarseRemeshed/";
        double startime = 0;

        while (EdgeLength > 0.01)
        {
            std::stringstream out;
            out << EdgeLength;
            std::string mesh_size = out.str();

            std::string offfile = "~/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/Cylinder.off";
            std::string CGALRemeshingCommand = "/Users/jcrawshaw/Documents/CGAL-5.0.2/examples/Polygon_mesh_processing/isotropic_remeshing_ForChaste -input " + offfile + " -output " + output_dir + "RemeshedCylinder.off -target_edge_length " + std::to_string(EdgeLength) + " -iterations " + std::to_string(10);
            int SystemOutput = std::system(CGALRemeshingCommand.c_str()); // system only takes char *

            std::string Remeshedvtu = output_dir + "RemeshedCylinder" + mesh_size + ".vtu";
            std::string stl2vtuCommand = " meshio-convert " + output_dir + "RemeshedCylinder.off " + Remeshedvtu;
            std::system(stl2vtuCommand.c_str());

            // Read in the new cylinder that is generated to my desire
            std::string mesh_file = "/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/CoarseRemeshed/RemeshedCylinder" + mesh_size + ".vtu"; //Remeshedvtu;//"/Users/jcrawshaw/Documents/testoutput/RemeshingErrorAnal3D/Cylinder_WaveyDef/CGAL/CourseInital/RemeshedCylinder.vtu";
            VtkMeshReader<2, 3> mesh_reader(mesh_file);
            HistoryDepMutableMesh<2, 3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            std::string output_directory = "RemeshingErrorAnal3D/Cylinder_WaveyDef/CoarseRemeshed/EdgeLength/" + mesh_size + "/";

            // Create the cells
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

            // Create a cell population
            HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
            cell_population.SetChasteOutputDirectory(output_directory, startime);
            cell_population.SetInitialAnlgesAcrossMembrane();
            cell_population.SetTargetRemeshingEdgeLength(RemeshedEdgeLength);
            cell_population.SetTargetRemeshingIterations(Iter);

            cell_population.SetRemeshingSoftwear("CGAL");
            cell_population.SetRelativePath(output_directory, startime);
            cell_population.SetPrintRemeshedIC(1);
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.SetOutputMeshInVtk(true);

            VtkMeshWriter<2, 3> mesh_writer0(output_directory, "InitalMesh", false);
            mesh_writer0.WriteFilesUsingMesh(mesh);
            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {
                c_vector<double, 3> InitalLocation = cell_population.GetNode(i)->rGetLocation();
                c_vector<double, 3> DeformedLocation = WavyDeformation(InitalLocation);
                cell_population.GetNode(i)->rGetModifiableLocation() = DeformedLocation;
            }

            VtkMeshWriter<2, 3> mesh_writer1(output_directory, "DeformedOriginalMesh", false);
            mesh_writer1.WriteFilesUsingMesh(mesh);

            // Remeshing should happen here --> This is where all the time is taken
            TRACE("Execute remeshing")
            cell_population.ExecuteHistoryDependentRemeshing();
            TRACE("Remeshing complete")
            VtkMeshWriter<2, 3> mesh_writer2(output_directory, "DeformedRemeshedMesh", false);
            mesh_writer2.WriteFilesUsingMesh(mesh);

            std::map<unsigned, c_vector<double, 3> > InitalNodePos = cell_population.GetInitalNodePositions();
            std::map<unsigned, c_vector<double, 3> > DeformedIC;
            std::vector<double> Error;

            for (int i = 0; i < mesh.GetNumNodes(); i++)
            {
                c_vector<double, 3> InitalLocation = InitalNodePos[cell_population.GetNode(i)->GetIndex()];
                c_vector<double, 3> DeformedLocation_IC = WavyDeformation(InitalLocation);
                DeformedIC[cell_population.GetNode(i)->GetIndex()] = DeformedLocation_IC;
                c_vector<double, 3> DeformedNode = cell_population.GetNode(i)->rGetLocation();
                Error.push_back(norm_2(DeformedNode - DeformedLocation_IC));
            }

            VtkMeshWriter<2, 3> mesh_writer(output_directory, "DeformedOriginalMesh", false); // I think this is falsely labeled!
            mesh_writer.WriteFilesUsingMesh(mesh);

            std::vector<double> quartiles = Quantile<double>(Error, { 0.25, 0.5, 0.75 });
            typename std::vector<double>::iterator Iterator = quartiles.begin();
            double Quartile1 = *Iterator;
            std::advance(Iterator, 1);
            double Med = *Iterator;
            std::advance(Iterator, 1);
            double Quartile3 = *Iterator;

            Median.push_back(Med);
            Q25.push_back(Quartile1);
            Q75.push_back(Quartile3);

            EdgeL.push_back(EdgeLength);
            EdgeLength = EdgeLength - 0.05;
        }

        std::vector<double>::iterator a = Median.begin();
        std::vector<double>::iterator b = Q25.begin();
        std::vector<double>::iterator c = Q75.begin();
                                                                                                        
        std::ofstream f("/Users/jcrawshaw/Documents/Projects/MeshMatlab/Remeshing3D/RawSpatialErrorData/SpatialErrorAnal_3D_CoarseRemeshed.txt");
        for (std::vector<double>::iterator i = EdgeL.begin(); i != EdgeL.end(); ++a, ++b, ++c, ++i)
        {
            f << *i << " " << *b << " " << *a << " " << *c << '\n';
        }
    }


    template <typename T>
    static inline double Lerp(T v0, T v1, T t)
    {
        return (1 - t) * v0 + t * v1;
    }

    template <typename T>
    static inline std::vector<T> Quantile(const std::vector<T>& inData, const std::vector<T>& probs)
    {
        if (inData.empty())
        {
            return std::vector<T>();
        }

        if (1 == inData.size())
        {
            return std::vector<T>(1, inData[0]);
        }

        std::vector<T> data = inData;
        std::sort(data.begin(), data.end());
        std::vector<T> quantiles;

        for (size_t i = 0; i < probs.size(); ++i)
        {
            T poi = Lerp<T>(-0.5, data.size() - 0.5, probs[i]);

            size_t left = std::max(int64_t(std::floor(poi)), int64_t(0));
            size_t right = std::min(int64_t(std::ceil(poi)), int64_t(data.size() - 1));

            T datLeft = data.at(left);
            T datRight = data.at(right);

            T quantile = Lerp<T>(datLeft, datRight, poi - left);

            quantiles.push_back(quantile);
        }
        return quantiles;
    }

    c_vector<double, 3> WavyDeformation(c_vector<double, 3> InitalLocation)
    {
        double R = sqrt(InitalLocation[0] * InitalLocation[0] + InitalLocation[1] * InitalLocation[1]);
        double Angle; // double Scalled_R;
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

        // c_vector<double, 3> DeformedLocation = Create_c_vector(X, Y, sin(InitalLocation[2]);

        double X = R * cos(Angle) + 1 * sin(1.5 * InitalLocation[2]) - 0.5; //+ InitalLocation[2]-1;
        double Y = R * 1.5 * sin(Angle);

        c_vector<double, 3> DeformedLocation = Create_c_vector(X, Y, InitalLocation[2]);

        return DeformedLocation;
        // cell_population.GetNode(i)->rGetModifiableLocation() = DeformedLocation;
    }



};

#endif /*TESTRELAXATION_HPP_*/
