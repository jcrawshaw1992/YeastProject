
#ifndef TestPottsModelOn2DSquare_HPP_
#define TestPottsModelOn2DSquare_HPP_

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <cstdio>
#include <ctime>


// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellsGenerator.hpp"
#include "OnLatticeSimulation.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "CellIdWriter.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "ArbitraryVolumeConstraintPottsUpdateRule.hpp"
#include "ArbitraryAdhesionPottsUpdateRule.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PottsMeshGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "Warnings.hpp"
#include "LogFile.hpp"
#include "SmartPointers.hpp"
#include "AspectRatioConstraintPottsUpdateRule.hpp"
#include "PottsMeshFromMutableMeshGenerator.hpp"
#include "RandomNumberGenerator.hpp"
#include "Debug.hpp"
#include "UblasCustomFunctions.hpp"

#include "PottsMeshFromMutableMeshGeneratorJess.hpp"

#include "PottsArbitrarySurfaceIn3DMesh.hpp"


#include "ArbitrarySurfaceAreaConstraintPottsUpdateRule.hpp"
#include "Honeycomb3DMeshGenerator.hpp"

// Needed for NodesOnlyMesh.
#include "PetscSetupAndFinalize.hpp"

#include <iostream>
#include <stdio.h>

class TestPottsOnCylinder: public AbstractCellBasedTestSuite
{
public:

double ComputeMeshTotalArea(TetrahedralMesh<2,3>& mutable_mesh)
    {
    	double total_area = 0.;
    	for (TetrahedralMesh<2,3>::ElementIterator element_iter = mutable_mesh.GetElementIteratorBegin();
			 element_iter != mutable_mesh.GetElementIteratorEnd();
			 ++element_iter)
    	{
    		const c_vector<double, 3> AC = element_iter->GetNodeLocation(1) - element_iter->GetNodeLocation(0);
    		const c_vector<double, 3> AB = element_iter->GetNodeLocation(2) - element_iter->GetNodeLocation(0);
    		total_area += 0.5 * norm_2(VectorProduct(AC, AB));
    
    	}

    	PRINT_VARIABLE(total_area);
    	return total_area;
    }
	


    void TestPottsOnCylinderWithDeformation()
	{

            double N_D =10;
            double N_Z =10;
            Honeycomb3DMeshGenerator generator(N_D, N_Z, 1, 1);
            MutableMesh<2, 3>* mutable_mesh = generator.GetMesh();
            mutable_mesh->Translate( unit_vector<double>(3, 2));
        /*
		 * Setup Potts simulation
		 */
		
        PottsMeshFromMutableMeshGeneratorJess<3> potts_generator(* mutable_mesh);
		PottsArbitrarySurfaceIn3DMesh<3>* p_potts_mesh = potts_generator.GetMesh();
	
		

		unsigned max_num_cells = 10;  
		unsigned num_cells_seeded = 5;
		double target_volume = ComputeMeshTotalArea(* mutable_mesh) / max_num_cells;
        unsigned num_node;
		for (unsigned i=0; i<num_cells_seeded; i++) 
		{
			// Pick a node randomly avoiding boundary nodes
			do
			{
				num_node = RandomNumberGenerator::Instance()->randMod(p_potts_mesh->GetNumNodes());  // Selects a random node 

			}
			while(p_potts_mesh->GetNode(num_node)->IsBoundaryNode()); // XXX not sure what this does, Im not sure it does anything 
			PRINT_VARIABLE(i);

			std::vector<Node<3>*> element_nodes;
			element_nodes.push_back(p_potts_mesh->GetNode(num_node));
			p_potts_mesh->AddElement(new PottsElement<3>(i, element_nodes));
			std::cout << "Seeded at "<< num_node << std::endl;
		}

        // Randomly place the cell markers for seeding
		std::vector<CellPtr> potts_cells;
		MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type); // Set the cell type for the seeds 
		CellsGenerator<FixedG1GenerationalCellCycleModel, 3> potts_cells_generator; // Set the cell cylce for the seeds 
		potts_cells_generator.GenerateBasicRandom(potts_cells, p_potts_mesh->GetNumElements(), p_diff_type); // Assigning the cell type and the cell cycle 

        TRACE("Generate Potts cell population");
		// Create cell population linking potts mesh and cells
		PottsBasedCellPopulation<3> potts_population(*p_potts_mesh, potts_cells);
		potts_population.AddCellWriter<CellIdWriter>();
        potts_population.SetNumSweepsPerTimestep(1);
		// Set up Potts simulation
		OnLatticeSimulation<3> potts_simulator(potts_population);
		// potts_simulator.SetOutputDirectory("TestPottsOn2DSquare_Potts");
		// potts_simulator.SetDt(0.1);
		// potts_simulator.SetSamplingTimestepMultiple(1);


        // // Reduce temperature to avoid randomly removing cells before they grow (parameters will need proper fitting)
		// potts_population.SetTemperature(1e-9);

        /*
         * Simulate cell seeding from randomly placed cell markers
         */
		// MAKE_PTR(ArbitraryVolumeConstraintPottsUpdateRule<3>, p_volume_constraint_update_rule);
		// PRINT_VARIABLE(target_volume);
		// p_volume_constraint_update_rule->SetMatureCellTargetVolume(target_volume);
		// p_volume_constraint_update_rule->SetDeformationEnergyParameter(4000000);
		// potts_simulator.AddUpdateRule(p_volume_constraint_update_rule);

		// MAKE_PTR(ArbitraryAdhesionPottsUpdateRule<3>, p_adhesion_update_rule);
		// p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02 / 1000);
		// p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16 / 1000);
		// potts_simulator.AddUpdateRule(p_adhesion_update_rule);


        MAKE_PTR(ArbitrarySurfaceAreaConstraintPottsUpdateRule<3>, p_area_constraint_update_rule);
        p_area_constraint_update_rule->SetTargetSurfaceArea (150e-6); // (150e-6); // P (is 150e-12 in paper but value was found wrong and corrected for 150e-6)
        p_area_constraint_update_rule->SetSurfaceAreaEnergyParameter(0.016e18); //(0.016e18); // lambda_p
        

        // Check the RemoveInternalEdgesFunction is working, this will throw an error if not working 
        p_area_constraint_update_rule->TestRemoveInternalEdgesFunction();
        p_area_constraint_update_rule->FindElementNeighbours(* mutable_mesh);
		p_area_constraint_update_rule->CalculateLatticeSitePerimeters(* mutable_mesh);

        
        // Check the correct lattice perimeter is returned for the single lattice site Potts cell 
        unsigned currentNodeIndex = num_node;
        std::set<unsigned> containing_elements = potts_population.GetNode(currentNodeIndex)->rGetContainingElementIndices();
        unsigned current_element = (*containing_elements.begin());
        PottsElement<3>* pCurrentElement = potts_population.rGetMesh().GetElement(current_element);
        double LatticePerimeter = p_area_constraint_update_rule->GetSurfaceAreaOfElement(current_element, pCurrentElement);
        assert(LatticePerimeter - 0.34641 < 1e-5); // Check the perimiter is correct

        // potts_simulator.AddUpdateRule(p_area_constraint_update_rule);






// // Need to have a look at what this is doing, to asserts need to be DIM ==2 
//         MAKE_PTR(AspectRatioConstraintPottsUpdateRule<3>, p_aspect_ratio_update_rule);
//         p_aspect_ratio_update_rule->SetTargetAspectRatio(4);
//         p_aspect_ratio_update_rule->SetAspectRatioEnergyParameter(1e-7);
//         potts_simulator.AddUpdateRule(p_aspect_ratio_update_rule);


    






		// // Set up the inital seeding 
		// unsigned Confluence_time = 1;
        // potts_simulator.SetEndTime(Confluence_time);
        // potts_simulator.SetOutputDirectory("TestPottsOn2dSquare_Potts_Setup");
		// p_potts_mesh->UpdatePottsNodeLocationFromDelaunay();
		

		
		// // p_potts_mesh->UpdatePottsVolumeAndArea(* mutable_mesh);
		// TRACE("solve");
        // potts_simulator.Solve();




        





		// // Change the potts outdirectory now that the Potts population should have achieved confluence 
        // potts_simulator.SetOutputDirectory("TestPottsOnVesselGrowth_Potts");

        // Add mechanotaxis to Potts simulation for timesteps to follow
        // MAKE_PTR(MechanotaxisPottsUpdateRule<3>, p_mechanotaxis_update_rule);
        // p_mechanotaxis_update_rule->SetTractionCorrelationParameter(1e-8);
        // potts_simulator.AddUpdateRule(p_mechanotaxis_update_rule);

	}
};

#endif /*TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_*/




// Checked that i have a proper mesh 



        // std::vector<CellPtr> cells;
        // CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        // cells_generator.GenerateBasicRandom(cells, mutable_mesh->GetNumNodes(), p_diff_type);

        // // Create a cell population. Links the Mesh and the cells
        // MeshBasedCellPopulation<2,3> cell_population(* mutable_mesh, cells);

        // cell_population.SetWriteVtkAsPoints(true);
        // cell_population.SetOutputMeshInVtk(true);
        // // cell_population.SetOutputCellProliferativeTypes(true);
        // cell_population.CalculateRestLengths();
		// // cell_population.SetOutputCellIdData(true);
        // //cell_population.SetDampingConstantNormal(0.1);//e5);
        // //cell_population.SetDampingConstantMutant(0.1);//e5);


        // // Set up cell-based simulation
        // SimulationTime::Instance()->Destroy();
        // SimulationTime::Instance()->SetStartTime(0);
        // TRACE("NewSimulation");
        // OffLatticeSimulation<2,3> simulator(cell_population);
        // simulator.SetOutputDirectory("TestPottsOnVesselGrowth_Growth");
        // simulator.SetDt(0.002);
        // simulator.SetSamplingTimestepMultiple(500);
        // simulator.SetEndTime(13);
        // simulator.SetUpdateCellPopulationRule(false); // No remeshing.
        // simulator.Solve();
