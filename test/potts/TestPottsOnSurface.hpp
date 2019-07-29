
#ifndef TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_
#define TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellsGenerator.hpp"
#include "OnLatticeSimulation.hpp"
#include "OffLatticeSimulation.hpp"
// #include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "ArbitraryVolumeConstraintPottsUpdateRule.hpp"
#include "ArbitraryAdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PottsMeshGenerator.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "Warnings.hpp"
#include "LogFile.hpp"
#include "SmartPointers.hpp"
#include "ChemotaxisPottsUpdateRule.hpp"
#include "AspectRatioConstraintPottsUpdateRule.hpp"
#include "PottsMeshFromMutableMeshGenerator.hpp"
#include "VtkMeshReader.hpp"
#include "RandomNumberGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "AbstractForce.hpp"
#include "Debug.hpp"

#include "PottsArbitraryMeshFromMutableMeshGenerator.hpp"

#include "FixedG1GenerationalCellCycleModel.hpp"
// #include "PottsArbitrarySurfaceIn3DMesh.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"



// Needed for NodesOnlyMesh.
#include "PetscSetupAndFinalize.hpp"

class TestPottsOnSurface : public AbstractCellBasedTestSuite
{
public:

    void TestPottsOn2DIn3DDisk() throw(Exception)
    {
    	TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/disk_in_3d");
    	MutableMesh<2,3> original_delaunay_mesh;
    	original_delaunay_mesh.ConstructFromMeshReader(mesh_reader);

		PottsArbitraryMeshFromMutableMeshGenerator<3> potts_generator(original_delaunay_mesh);
		PottsArbitrarySurfaceIn3DMesh<3>* p_potts_mesh = potts_generator.GetMesh();

        // Make two triangular elements out of these nodes
        //What is this 
        std::vector<Node<3>*> element_nodes;
        element_nodes.push_back(p_potts_mesh->GetNode(100));

        p_potts_mesh->AddElement(new PottsElement<3>(0, element_nodes));


        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_potts_mesh->GetNumElements(), p_diff_type);

        // Create cell population
        PottsBasedCellPopulation<3> cell_population(*p_potts_mesh, cells);
        //cell_population.SetTemperature(1e-9);

        // Set up cell-based simulation
        OnLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("PottsOn3dDisc");
        simulator.SetDt(0.1);
        simulator.SetEndTime(100);


        // Create update rules and pass to the simulation
        MAKE_PTR(ArbitraryVolumeConstraintPottsUpdateRule<3>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(2.5);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(400);
        // simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);


        MAKE_PTR(ArbitraryAdhesionPottsUpdateRule<3>, p_adhesion_update_rule);
        p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02);
        p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16 * 10);
        // simulator.AddPottsUpdateRule(p_adhesion_update_rule);
        simulator.AddUpdateRule(p_adhesion_update_rule);

//        MAKE_PTR(AspectRatioConstraintPottsUpdateRule<3>, p_aspect_ratio_update_rule);
//        p_aspect_ratio_update_rule->SetTargetAspectRatio(1);
//        p_aspect_ratio_update_rule->SetAspectRatioEnergyParameter(1);
//        simulator.AddPottsUpdateRule(p_aspect_ratio_update_rule);

//        MAKE_PTR(ChemotaxisPottsUpdateRule<2>, p_chemotaxis_update_rule);
//        simulator.AddPottsUpdateRule(p_chemotaxis_update_rule);

        // Run simulation
        simulator.Solve();

        original_delaunay_mesh.Scale(1.5,3.0,1.0);
        p_potts_mesh->UpdatePottsNodeLocationFromDelaunay();

        simulator.SetEndTime(200);
        simulator.Solve();

    }
};

#endif /*TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_*/
