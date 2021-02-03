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

    void TestInstantDeformation() throw(Exception)
    {

        double RemeshedEdgeLength = 1e-1;
        double Iter = 2;
        std::string output_directory = "TestHemeLBOnNetwork/InstantGrowth";
        std::string mesh_file = "/Users/jcrawshaw/docker-polnet-master/IdealiseNetworks/SimpleNetwork.vtu";
        double startime = 0;

        // Read in the new cylinder that is generated to my desire
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
        cell_population.SetChasteOutputDirectory(output_directory, startime);
        // cell_population.SetInitialAnlgesAcrossMembrane();
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
                                               ///   pNode->
            c_vector<double, 3> InitalLocation = cell_population.GetNode(i)->rGetLocation();
        
            unsigned node_index = i;//rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<3>* p_node =  cell_population.GetNode(i);// rCellPopulation.GetNode(node_index);

            c_vector<double, 3> Normal = zero_vector<double>(3);

            std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();

            assert(containing_elements.size() > 0);
            for (std::set<unsigned>::iterator iter = containing_elements.begin();
                    iter != containing_elements.end();
                    ++iter)
            {
                Node<3>* pNode0 = mesh.GetNode(mesh.GetElement(*iter)->GetNodeGlobalIndex(0));
                Node<3>* pNode1 = mesh.GetNode(mesh.GetElement(*iter)->GetNodeGlobalIndex(1));
                Node<3>* pNode2 = mesh.GetNode(mesh.GetElement(*iter)->GetNodeGlobalIndex(2));

                c_vector<double, 3> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
                c_vector<double, 3> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3

                c_vector<double, 3> normalVector = VectorProduct(vector_12, vector_13);
                Normal += normalVector;
            }
            Normal /=norm_2(Normal);
        
            c_vector<double, 3> force = 10 *Normal;

            c_vector<double, 3> DeformedLocation = InitalLocation + Normal * 0.05;//RadialDeformation(InitalLocation);
            cell_population.GetNode(i)->rGetModifiableLocation() = DeformedLocation;
                
        }

        VtkMeshWriter<2, 3> mesh_writer1(output_directory, "DeformedOriginalMesh", false);
        mesh_writer1.WriteFilesUsingMesh(mesh);

        // Remeshing should happen here --> This is where all the time is taken
        TRACE("Execute remeshing")
        cell_population.ExecuteHistoryDependentRemeshing();
        TRACE("Remeshing complete")
        VtkMeshWriter<2, 3> mesh_writer2(output_directory, "DeformedAndRemeshedMesh", false);
        mesh_writer2.WriteFilesUsingMesh(mesh);

           
    }


    c_vector<double, 3> RadialDeformation(c_vector<double, 3> InitalLocation)
    {

        // Need the normal and the radial increase. which we shall say 2? 
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
