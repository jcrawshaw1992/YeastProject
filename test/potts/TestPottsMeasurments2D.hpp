
#ifndef TestPottsShearMigration2D_HPP_
#define TestPottsShearMigration2D_HPP_

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <cstdio>
#include <ctime>
// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "OnLatticeSimulation.hpp"

// Include the relevant cell writers
#include "CellAspectRatioWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellMajorAxisAngleWriter.hpp"
// #include "PottsBasedCellPopulation.hpp"
#include "ArbitraryWrappedAdhesionPottsUpdateRule.hpp"
#include "MeshBasedCellPopulation.hpp"
// New rule to make the cells wrap
#include "ArbitraryCurvatureConstraintPottsUpdateRule.hpp"

#include "PottsCellPropertiesModifier.hpp"

// #include "TractionDataLoader.hpp"
// #include "AdhesionPottsUpdateRule.hpp"
// #include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PottsMeshGenerator.hpp"
// #include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "LogFile.hpp"
#include "SmartPointers.hpp"
#include "Warnings.hpp"
// #include "MechanotaxisPottsUpdateRule.hpp"
#include "AspectRatioConstraintPottsUpdateRule.hpp"
#include "PottsMeshFromMutableMeshGenerator.hpp"

#include "VtkMeshReader.hpp"

#include "AbstractForce.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "RandomNumberGenerator.hpp"
// #include "ChemotaxisPottsUpdateRule.hpp"
#include "Debug.hpp"

#include "UblasCustomFunctions.hpp"

#include "PottsMeshFromMutableMeshGeneratorJess.hpp"

#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "MutableMesh.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"
// #include "ArbitrarySurfaceAreaConstraintPottsUpdateRule.hpp"

// #include "ArbitraryVolumeConstraintPottsUpdateRule.hpp"

#include "Honeycomb3DMeshGenerator.hpp"

// Needed for NodesOnlyMesh.
#include "PetscSetupAndFinalize.hpp"

#include <iostream>
#include <stdio.h>

// Adapted update rules

#include "ArbitraryPerimeterOnSurfacePottsUpdateRule.hpp"
#include "ArbitraryVolumeOnSurfacePottsUpdateRule.hpp"

#include "WrappedPottsBasedCellPopulation.hpp"

#include "PeriodicRectangleMeshGenerator.hpp"

class TestShearMigration : public AbstractCellBasedTestSuite
{
public:
    void TestPottsOnCylinderWithDeformation()
    {
        double N_D = 10;
        double N_Z = 20;
        double Width = 10;

        double Length = 20;

        PeriodicRectangleMeshGenerator generator(N_D, N_Z, Width, Length);
        MutableMesh<2, 3>* mutable_mesh = generator.GetMesh();
        std::vector<unsigned> BoundaryVector = generator.GetBoundaryVector();

        /*
		 * Setup Potts simulation
		 */

        PottsMeshFromMutableMeshGeneratorJess<3> potts_generator(*mutable_mesh);
        PottsArbitrarySurfaceIn3DMesh<3>* p_potts_mesh = potts_generator.GetMesh();
        p_potts_mesh->SetBoundaries(BoundaryVector);

        std::map<unsigned, unsigned> ElementPairing;

        //Select my nodes as I see fit

        // Cell 0   ==   * *
        unsigned i = 0;
        unsigned num_node = 21;
        std::vector<Node<3>*> element_nodes;
        std::vector<Node<3>*> element_nodes2;
        element_nodes.push_back(p_potts_mesh->GetNode(num_node));
        element_nodes2.push_back(p_potts_mesh->GetNode(num_node + 1));
        p_potts_mesh->AddElement(new PottsElement<3>(i, element_nodes));
        p_potts_mesh->AddElement(new PottsElement<3>(i + 1, element_nodes2));
        ElementPairing[i] = i + 1;
        ElementPairing[i + 1] = i;

        // //Cell 1   ==   * *
        //               *
        i = 2;
        num_node = 85;
        std::vector<Node<3>*> Aelement_nodes;
        std::vector<Node<3>*> Aelement_nodes2;
        Aelement_nodes.push_back(p_potts_mesh->GetNode(num_node));
        Aelement_nodes.push_back(p_potts_mesh->GetNode(75));
        Aelement_nodes2.push_back(p_potts_mesh->GetNode(num_node + 1));
        p_potts_mesh->AddElement(new PottsElement<3>(i, Aelement_nodes));
        p_potts_mesh->AddElement(new PottsElement<3>(i + 1, Aelement_nodes2));
        ElementPairing[i] = i + 1;
        ElementPairing[i + 1] = i;

        //Cell 2   ==   * *
        //             *
        i = 4;
        num_node = 53;
        std::vector<Node<3>*> Belement_nodes;
        std::vector<Node<3>*> Belement_nodes2;
        Belement_nodes.push_back(p_potts_mesh->GetNode(num_node));
        Belement_nodes.push_back(p_potts_mesh->GetNode(43));
        Belement_nodes2.push_back(p_potts_mesh->GetNode(num_node + 1));
        p_potts_mesh->AddElement(new PottsElement<3>(i, Belement_nodes));
        p_potts_mesh->AddElement(new PottsElement<3>(i + 1, Belement_nodes2));
        ElementPairing[i] = i + 1;
        ElementPairing[i + 1] = i;

        // Cell 3
        //
        i = 6;
        std::vector<Node<3>*> Celement_nodes;
        std::vector<Node<3>*> Celement_nodes2;
        Celement_nodes.push_back(p_potts_mesh->GetNode(3));
        Celement_nodes.push_back(p_potts_mesh->GetNode(4));
        Celement_nodes.push_back(p_potts_mesh->GetNode(5));
        Celement_nodes2.push_back(p_potts_mesh->GetNode(6));
        Celement_nodes.push_back(p_potts_mesh->GetNode(13));
        Celement_nodes.push_back(p_potts_mesh->GetNode(14));
        Celement_nodes2.push_back(p_potts_mesh->GetNode(15));
        Celement_nodes2.push_back(p_potts_mesh->GetNode(16));
        Celement_nodes2.push_back(p_potts_mesh->GetNode(17));

        Celement_nodes.push_back(p_potts_mesh->GetNode(23));
        Celement_nodes.push_back(p_potts_mesh->GetNode(24));
        Celement_nodes2.push_back(p_potts_mesh->GetNode(26));
        Celement_nodes2.push_back(p_potts_mesh->GetNode(27));
        Celement_nodes2.push_back(p_potts_mesh->GetNode(37));
        Celement_nodes2.push_back(p_potts_mesh->GetNode(47));

        p_potts_mesh->AddElement(new PottsElement<3>(i, Celement_nodes));
        p_potts_mesh->AddElement(new PottsElement<3>(i + 1, Celement_nodes2));
        ElementPairing[i] = i + 1;
        ElementPairing[i + 1] = i;

        // Cell 4
        //
        i = 8;
        std::vector<Node<3>*> Delement_nodes;
        std::vector<Node<3>*> Delement_nodes2;
        Delement_nodes.push_back(p_potts_mesh->GetNode(60));
        Delement_nodes.push_back(p_potts_mesh->GetNode(61));
        Delement_nodes.push_back(p_potts_mesh->GetNode(50));
        Delement_nodes.push_back(p_potts_mesh->GetNode(40));
        Delement_nodes.push_back(p_potts_mesh->GetNode(41));
        Delement_nodes.push_back(p_potts_mesh->GetNode(70));

        Delement_nodes2.push_back(p_potts_mesh->GetNode(79));
        Delement_nodes2.push_back(p_potts_mesh->GetNode(78));
        Delement_nodes2.push_back(p_potts_mesh->GetNode(69));
        Delement_nodes2.push_back(p_potts_mesh->GetNode(59));
        Delement_nodes2.push_back(p_potts_mesh->GetNode(89));
        Delement_nodes2.push_back(p_potts_mesh->GetNode(99));
        p_potts_mesh->AddElement(new PottsElement<3>(i, Delement_nodes));
        p_potts_mesh->AddElement(new PottsElement<3>(i + 1, Delement_nodes2));
        ElementPairing[i] = i + 1;
        ElementPairing[i + 1] = i;

        // Cell 5
        //
        i = 10; // | * *

        std::vector<Node<3>*> Eelement_nodes;
        std::vector<Node<3>*> Eelement_nodes2;
        Eelement_nodes.push_back(p_potts_mesh->GetNode(10));
        Eelement_nodes2.push_back(p_potts_mesh->GetNode(11));

        p_potts_mesh->AddElement(new PottsElement<3>(i, Eelement_nodes));
        p_potts_mesh->AddElement(new PottsElement<3>(i + 1, Eelement_nodes2));
        ElementPairing[i] = i + 1;
        ElementPairing[i + 1] = i;

        // Cell 6
        i = 12;
        std::vector<Node<3>*> Felement_nodes;
        std::vector<Node<3>*> Felement_nodes2;
        Felement_nodes.push_back(p_potts_mesh->GetNode(20));
        Felement_nodes2.push_back(p_potts_mesh->GetNode(39));

        p_potts_mesh->AddElement(new PottsElement<3>(i, Felement_nodes));
        p_potts_mesh->AddElement(new PottsElement<3>(i + 1, Felement_nodes2));
        ElementPairing[i] = i + 1;
        ElementPairing[i + 1] = i;

        // Cell 7
        i = 14;
        std::vector<Node<3>*> Gelement_nodes;
        std::vector<Node<3>*> Gelement_nodes2;
        Gelement_nodes.push_back(p_potts_mesh->GetNode(119));
        Gelement_nodes2.push_back(p_potts_mesh->GetNode(110));

        p_potts_mesh->AddElement(new PottsElement<3>(i, Gelement_nodes));
        p_potts_mesh->AddElement(new PottsElement<3>(i + 1, Gelement_nodes2));
        ElementPairing[i] = i + 1;
        ElementPairing[i + 1] = i;

        // Cell 8
        i = 16;
        std::vector<Node<3>*> Helement_nodes;
        std::vector<Node<3>*> Helement_nodes2;
        Helement_nodes.push_back(p_potts_mesh->GetNode(199));
        Helement_nodes2.push_back(p_potts_mesh->GetNode(198));

        p_potts_mesh->AddElement(new PottsElement<3>(i, Helement_nodes));
        p_potts_mesh->AddElement(new PottsElement<3>(i + 1, Helement_nodes2));
        ElementPairing[i] = i + 1;
        ElementPairing[i + 1] = i;

        // Cell 9
        i = 18;
        std::vector<Node<3>*> Jelement_nodes;
        std::vector<Node<3>*> Jelement_nodes2;
        Jelement_nodes.push_back(p_potts_mesh->GetNode(0));
        Jelement_nodes2.push_back(p_potts_mesh->GetNode(9));

        p_potts_mesh->AddElement(new PottsElement<3>(i, Jelement_nodes));
        p_potts_mesh->AddElement(new PottsElement<3>(i + 1, Jelement_nodes2));
        ElementPairing[i] = i + 1;
        ElementPairing[i + 1] = i;

        // Cell 10
        i = 20;
        std::vector<Node<3>*> Kelement_nodes;
        std::vector<Node<3>*> Kelement_nodes2;
        Kelement_nodes.push_back(p_potts_mesh->GetNode(180));
        Kelement_nodes2.push_back(p_potts_mesh->GetNode(190));

        p_potts_mesh->AddElement(new PottsElement<3>(i, Kelement_nodes));
        p_potts_mesh->AddElement(new PottsElement<3>(i + 1, Kelement_nodes2));
        ElementPairing[i] = i + 1;
        ElementPairing[i + 1] = i;

        // Cell 11
        i = 22;
        std::vector<Node<3>*> Lelement_nodes;
        std::vector<Node<3>*> Lelement_nodes2;
        Lelement_nodes.push_back(p_potts_mesh->GetNode(130));
        Lelement_nodes.push_back(p_potts_mesh->GetNode(131));
        Lelement_nodes.push_back(p_potts_mesh->GetNode(132));
        Lelement_nodes.push_back(p_potts_mesh->GetNode(133));
        Lelement_nodes.push_back(p_potts_mesh->GetNode(134));

        Lelement_nodes.push_back(p_potts_mesh->GetNode(140));
        Lelement_nodes.push_back(p_potts_mesh->GetNode(141));
        Lelement_nodes.push_back(p_potts_mesh->GetNode(142));
        Lelement_nodes.push_back(p_potts_mesh->GetNode(143));
        Lelement_nodes.push_back(p_potts_mesh->GetNode(144));
        Lelement_nodes.push_back(p_potts_mesh->GetNode(145));

        Lelement_nodes2.push_back(p_potts_mesh->GetNode(135));
        Lelement_nodes2.push_back(p_potts_mesh->GetNode(136));
        Lelement_nodes2.push_back(p_potts_mesh->GetNode(137));
        Lelement_nodes2.push_back(p_potts_mesh->GetNode(138));
        Lelement_nodes2.push_back(p_potts_mesh->GetNode(139));

        Lelement_nodes2.push_back(p_potts_mesh->GetNode(146));
        Lelement_nodes2.push_back(p_potts_mesh->GetNode(147));
        Lelement_nodes2.push_back(p_potts_mesh->GetNode(148));
        Lelement_nodes2.push_back(p_potts_mesh->GetNode(149));

        p_potts_mesh->AddElement(new PottsElement<3>(i, Lelement_nodes));
        p_potts_mesh->AddElement(new PottsElement<3>(i + 1, Lelement_nodes2));
        ElementPairing[i] = i + 1;
        ElementPairing[i + 1] = i;

        // // Randomly place the cell markers for seeding
        std::vector<CellPtr> potts_cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type); // Set the cell type for the seeds
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> potts_cells_generator; // Set the cell cylce for the seeds
        potts_cells_generator.GenerateBasicRandom(potts_cells, p_potts_mesh->GetNumElements(), p_diff_type); // Assigning the cell type and the cell cycle

        TRACE("Generate Potts cell population");
        // Create cell population linking potts mesh and cells
        WrappedPottsBasedCellPopulation<3> potts_population(*p_potts_mesh, potts_cells, ElementPairing);

        // // Add in all the cell writers
        potts_population.AddCellWriter<CellIdWriter>();
        potts_population.SetNumSweepsPerTimestep(1);

        // Set up Potts simulation
        OnLatticeSimulation<3> potts_simulator(potts_population);

        p_potts_mesh->UpdatePottsNodeLocationFromDelaunay();

        // I need to go through and fix up the center lines -- in this Test code I want the
        // center lines to be either at the center of a bunch of nodes .... or if there is only
        // two nodes in the cell pick one node to set at the center line and call it good -- this
        // will stop a cell that is wrapped about the bonadry having a center point that is ages away

        int counter = 0;
        for (std::map<unsigned, unsigned>::iterator it = ElementPairing.begin(); it != ElementPairing.end(); it++) // The way the map is constructed is such that it is doubled every second so i connected cell x cell y and i+1 connected cell y cell x
        {
            //  std::cout << it->first << '\t' << it->second << std::endl;
            if (counter % 2 == 0) // Only need to access one in every pair, from this we fix both
            {
                PottsElement<3>* p_element_1 = potts_population.GetElement(it->first);
                PottsElement<3>* p_element_2 = potts_population.GetElement(it->second);

                double CenterPoint = 0;

                // PRINT_VARIABLE(it->first)
                if (it->first == 8)
                {
                    CenterPoint = p_element_1->GetNodeLocation(0)[0];
                }
                else if (p_element_1->GetNumNodes() + p_element_2->GetNumNodes() != 2)
                {

                    // need to loop over and get the mid points
                    for (int i = 0; i < p_element_1->GetNumNodes(); ++i)
                    {
                        CenterPoint += p_element_1->GetNodeLocation(i)[0]; // ONly care about the x corrdinate
                    }
                    for (int i = 0; i < p_element_2->GetNumNodes(); ++i)
                    {
                        CenterPoint += p_element_2->GetNodeLocation(i)[0];
                    }

                    CenterPoint /= (p_element_1->GetNumNodes() + p_element_2->GetNumNodes());
                }
                else if (p_element_1->GetNumNodes() + p_element_2->GetNumNodes() == 2)
                {
                    CenterPoint = p_element_1->GetNodeLocation(0)[0];
                }
                CellPtr pCell1 = potts_population.GetCellUsingLocationIndex(it->first);
                CellPtr pCell2 = potts_population.GetCellUsingLocationIndex(it->second);

                pCell1->GetCellData()->SetItem("Center", CenterPoint);
                pCell2->GetCellData()->SetItem("Center", CenterPoint);
            }
            counter += 1;
        }

        /////////////////////////////

        CellPtr p_cell;
        CellPtr p_cell2;
        double Center;
        double Center2;

        // TS_ASSERT(1 + 1 > 1);

        // //Check cell perimeters
        double A = 0.596285;
        double B = 0.66666666666;
        double E = 0.5;

        // These ones all work
        //

        // Cell 0
        p_cell = potts_population.GetCellUsingLocationIndex(0);
        p_cell2 = potts_population.GetCellUsingLocationIndex(1);
        Center = p_cell->GetCellData()->GetItem("Center");
        Center2 = p_cell2->GetCellData()->GetItem("Center");
        TS_ASSERT(Center == Center2);

        double PerimeterOfCell = p_potts_mesh->GetPerimeterOfCoupledElements(0, 1, Center);
        double CalculatedPerim = 2 * B + 8 * A;
        TS_ASSERT_DELTA(PerimeterOfCell, CalculatedPerim, 0.1);

        // // // // Cell 1
        p_cell = potts_population.GetCellUsingLocationIndex(2);
        p_cell2 = potts_population.GetCellUsingLocationIndex(3);
        Center = p_cell->GetCellData()->GetItem("Center");
        Center2 = p_cell2->GetCellData()->GetItem("Center");
        TS_ASSERT(Center == Center2);

        double PerimeterOfCell2 = p_potts_mesh->GetPerimeterOfCoupledElements(2, 3, Center);
        double CalculatedPerim2 = 4 * B + 8 * A;
        TS_ASSERT_DELTA(PerimeterOfCell2, CalculatedPerim2, 0.1);

        // Cell 2
        p_cell = potts_population.GetCellUsingLocationIndex(4);
        p_cell2 = potts_population.GetCellUsingLocationIndex(5);
        Center = p_cell->GetCellData()->GetItem("Center");
        Center2 = p_cell2->GetCellData()->GetItem("Center");
        TS_ASSERT(Center == Center2);

        double PerimeterOfCell3 = p_potts_mesh->GetPerimeterOfCoupledElements(4, 5, Center);
        double CalculatedPerim3 = 4 * B + 10 * A;
        PRINT_2_VARIABLES(PerimeterOfCell3, CalculatedPerim3);
        TS_ASSERT_DELTA(PerimeterOfCell3, CalculatedPerim3, 0.1);

        // Cell 3

        p_cell = potts_population.GetCellUsingLocationIndex(6);
        p_cell2 = potts_population.GetCellUsingLocationIndex(7);
        Center = p_cell->GetCellData()->GetItem("Center");
        Center2 = p_cell2->GetCellData()->GetItem("Center");
        TS_ASSERT(Center == Center2);

        double PerimeterOfCell4 = p_potts_mesh->GetPerimeterOfCoupledElements(6, 7, Center);
        double CalculatedPerim4 = 11 * B + 20 * A + 8 * E;
        PRINT_2_VARIABLES(PerimeterOfCell4, CalculatedPerim4);
        TS_ASSERT_DELTA(PerimeterOfCell4, CalculatedPerim4, 0.1);

        // //   // Cell 5
        // //
        p_cell = potts_population.GetCellUsingLocationIndex(10);
        p_cell2 = potts_population.GetCellUsingLocationIndex(11);
        Center = p_cell->GetCellData()->GetItem("Center");
        Center2 = p_cell2->GetCellData()->GetItem("Center");
        TS_ASSERT(Center == Center2);
        double PerimeterOfCell5 = p_potts_mesh->GetPerimeterOfCoupledElements(10, 11, Center);
        double CalculatedPerim5 = 2 * B + 8 * A;
        PRINT_2_VARIABLES(PerimeterOfCell5, CalculatedPerim5);
        TS_ASSERT_DELTA(PerimeterOfCell5, CalculatedPerim5, 0.1);

        //  // Cell 6
        // //   //
        p_cell = potts_population.GetCellUsingLocationIndex(12);
        p_cell2 = potts_population.GetCellUsingLocationIndex(13);
        Center = p_cell->GetCellData()->GetItem("Center");
        Center2 = p_cell2->GetCellData()->GetItem("Center");
        TS_ASSERT(Center == Center2);
        double PerimeterOfCell6 = p_potts_mesh->GetPerimeterOfCoupledElements(12, 13, Center);
        double CalculatedPerim6 = 4 * B + 6 * A;
        PRINT_2_VARIABLES(PerimeterOfCell6, CalculatedPerim6);
        TS_ASSERT_DELTA(PerimeterOfCell6, CalculatedPerim6, 0.1);

        double EleIndex;

        //  // Cell 7
        //   //
        EleIndex = 7 * 2;
        p_cell = potts_population.GetCellUsingLocationIndex(EleIndex);
        p_cell2 = potts_population.GetCellUsingLocationIndex(EleIndex + 1);
        Center = p_cell->GetCellData()->GetItem("Center");
        Center2 = p_cell2->GetCellData()->GetItem("Center");
        TS_ASSERT(Center == Center2);
        double PerimeterOfCell7 = p_potts_mesh->GetPerimeterOfCoupledElements(EleIndex, EleIndex + 1, Center);
        double CalculatedPerim7 = 2 * B + 8 * A;
        PRINT_2_VARIABLES(PerimeterOfCell7, CalculatedPerim7);
        TS_ASSERT_DELTA(PerimeterOfCell7, CalculatedPerim7, 0.1);

        // Cell8
        //
        EleIndex = 8 * 2;
        TRACE("going to have to put some thought in to this one -- one is on the zip edge")
        p_cell = potts_population.GetCellUsingLocationIndex(EleIndex);
        p_cell2 = potts_population.GetCellUsingLocationIndex(EleIndex + 1);
        Center = p_cell->GetCellData()->GetItem("Center");
        Center2 = p_cell2->GetCellData()->GetItem("Center");
        TS_ASSERT(Center == Center2);
        double PerimeterOfCell8 = p_potts_mesh->GetPerimeterOfCoupledElements(EleIndex, EleIndex + 1, Center);
        double CalculatedPerim8 = B + 4 * A + 4 * E;
        PRINT_2_VARIABLES(PerimeterOfCell8, CalculatedPerim8);
        TS_ASSERT_DELTA(PerimeterOfCell8, CalculatedPerim8, 0.1);

        // //   // //  // Cell 4
        p_cell = potts_population.GetCellUsingLocationIndex(8);
        p_cell2 = potts_population.GetCellUsingLocationIndex(9);
        Center = p_cell->GetCellData()->GetItem("Center");
        Center2 = p_cell2->GetCellData()->GetItem("Center");
        TS_ASSERT(Center == Center2);

        double PerimeterOfCell5a = p_potts_mesh->GetPerimeterOfCoupledElements(8, 9, Center);
        double CalculatedPerim5a = 12 * B + 18 * A;
        PRINT_2_VARIABLES(PerimeterOfCell5a, CalculatedPerim5a);
        TS_ASSERT_DELTA(PerimeterOfCell5a, CalculatedPerim5a, 0.1);

        //Cell 9
        p_cell = potts_population.GetCellUsingLocationIndex(18);
        p_cell2 = potts_population.GetCellUsingLocationIndex(19);
        Center = p_cell->GetCellData()->GetItem("Center");
        Center2 = p_cell2->GetCellData()->GetItem("Center");
        TS_ASSERT(Center == Center2);

        double PerimeterOfCell9 = p_potts_mesh->GetPerimeterOfCoupledElements(18, 19, Center);
        double CalculatedPerim9 = B + 4 * A + 4 * E;
        PRINT_2_VARIABLES(PerimeterOfCell9, CalculatedPerim9);
        TS_ASSERT_DELTA(PerimeterOfCell9, CalculatedPerim9, 0.1);

        //Cell 10
        p_cell = potts_population.GetCellUsingLocationIndex(20);
        p_cell2 = potts_population.GetCellUsingLocationIndex(21);
        Center = p_cell->GetCellData()->GetItem("Center");
        Center2 = p_cell2->GetCellData()->GetItem("Center");
        TS_ASSERT(Center == Center2);

        double PerimeterOfCell10 = p_potts_mesh->GetPerimeterOfCoupledElements(20, 21, Center);
        double CalculatedPerim10 = 3 * B + 4 * A + 2 * E;
        PRINT_2_VARIABLES(PerimeterOfCell10, CalculatedPerim10);
        TS_ASSERT_DELTA(PerimeterOfCell10, CalculatedPerim10, 0.1);

        // // //Cell 11
        p_cell = potts_population.GetCellUsingLocationIndex(22);
        p_cell2 = potts_population.GetCellUsingLocationIndex(23);
        Center = p_cell->GetCellData()->GetItem("Center");
        Center2 = p_cell2->GetCellData()->GetItem("Center");
        TS_ASSERT(Center == Center2);

        double PerimeterOfCell11 = p_potts_mesh->GetPerimeterOfCoupledElements(22, 23, Center);
        double CalculatedPerim11 = 44 * A + 8 * B - 4 * B - 2 * A;
        PRINT_2_VARIABLES(PerimeterOfCell11, CalculatedPerim11);
        TS_ASSERT_DELTA(PerimeterOfCell11, CalculatedPerim11, 0.1);

        // potts_simulator.SetOutputDirectory("TestingPeriodicDomian/MeasurmentTests");
        // potts_simulator.SetEndTime(0.05);

        // potts_simulator.Solve();
        TRACE("Simulation Complete")
    }
};

#endif /*TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_*/
