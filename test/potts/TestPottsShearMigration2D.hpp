
#ifndef TestPottsShearMigration2D_HPP_
#define TestPottsShearMigration2D_HPP_

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <cstdio>
#include <ctime>
// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "CellsGenerator.hpp"
#include "OnLatticeSimulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"

// Include the relevant cell writers
#include "CellIdWriter.hpp"
#include "CellAspectRatioWriter.hpp"
#include "CellMajorAxisAngleWriter.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "ArbitraryAdhesionPottsUpdateRule.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
// New rule to make the cells wrap
#include "ArbitraryCurvatureConstraintPottsUpdateRule.hpp"
#include "AppliedForceModifier.hpp"

#include "PottsCellPropertiesModifier.hpp"

// #include "TractionDataLoader.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PottsMeshGenerator.hpp"
// #include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "Warnings.hpp"
#include "LogFile.hpp"
#include "SmartPointers.hpp"
#include "MechanotaxisPottsUpdateRule.hpp"
#include "AspectRatioConstraintPottsUpdateRule.hpp"
#include "PottsMeshFromMutableMeshGenerator.hpp"
// #include "VtkMeshReader.hpp"
#include "RandomNumberGenerator.hpp"
// #include "GeneralisedLinearSpringForce.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "AbstractForce.hpp"
// #include "ChemotaxisPottsUpdateRule.hpp"
#include "Debug.hpp"

#include "UblasCustomFunctions.hpp"

#include "PottsMeshFromMutableMeshGeneratorJess.hpp"

#include "PottsArbitrarySurfaceIn3DMesh.hpp"
#include "ArbitrarySurfaceAreaConstraintPottsUpdateRule.hpp"
#include "ArbitraryVolumeConstraintPottsUpdateRule.hpp"
#include "Honeycomb3DMeshGenerator.hpp"

// Needed for NodesOnlyMesh.
#include "PetscSetupAndFinalize.hpp"

#include <iostream>
#include <stdio.h>

class TestShearMigration: public AbstractCellBasedTestSuite
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

            // Regular square --hexagonal lattie
            double scale = 1e-3; 
             double N_D = 35;
             double N_Z = 50;
			 double Width = 35; //5e-3
			 double Length = 50; //30e-3

            Honeycomb3DMeshGenerator generator(N_D, N_Z,Width , Length);
            MutableMesh<2, 3>* mutable_mesh = generator.GetMesh();
            // mutable_mesh->Translate( unit_vector<double>(3, 2));
             TRACE("DINDINDINDIND");
        
        /*
		 * Setup Potts simulation
		 */
		
        PottsMeshFromMutableMeshGeneratorJess<3> potts_generator(* mutable_mesh);
        // PottsMeshFromMutableMeshGeneratorJess<2> potts_generator(* mutable_mesh);
		PottsArbitrarySurfaceIn3DMesh<3>* p_potts_mesh = potts_generator.GetMesh();
	
		unsigned max_num_cells =2;  
		unsigned num_cells_seeded = 3;
		double CellScale =1e6;//1e3;


		double LongAxis = 10/2 ;//* 1e-6 * CellScale ;//*scale ;
		double ShortAxis = 5.0/2.0 ;//*1e-6 * CellScale ;//*scale ;
		double CellArea = M_PI * LongAxis  * ShortAxis ;

		double AxisRatio = pow(LongAxis -ShortAxis,2)/ pow(LongAxis +ShortAxis,2);//  (LongAxis -ShortAxis)* (LongAxis -ShortAxis) /(LongAxis +ShortAxis)*(LongAxis +ShortAxis); 
		double CellPerimeter = M_PI * (LongAxis +ShortAxis)* (3 *AxisRatio * 1/(sqrt(-3*AxisRatio+4) +10) +1) ;   

		for (unsigned i=0; i<num_cells_seeded; i++) 
		{
			// Pick a node randomly avoiding boundary nodes
			unsigned num_node;
			do
			{
				double yPos =Length;
				while (yPos > Length/4)
				{
				  num_node = RandomNumberGenerator::Instance()->randMod(p_potts_mesh->GetNumNodes());  // Selects a random node 
				// Check nodes are on the left side 
				  yPos = mutable_mesh->GetNode(num_node)->rGetLocation()[1];
				}


			}			
			
			while(p_potts_mesh->GetNode(num_node)->IsBoundaryNode()); // XXX not sure what this does, Im not sure it does anything 
			PRINT_VARIABLE(i);

			std::vector<Node<3>*> element_nodes;
			element_nodes.push_back(p_potts_mesh->GetNode(num_node));
			p_potts_mesh->AddElement(new PottsElement<3>(i, element_nodes));
			std::cout << "Seeded at "<< num_node << std::endl;
		}

        // // Randomly place the cell markers for seeding
		std::vector<CellPtr> potts_cells;
		MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type); // Set the cell type for the seeds 
		CellsGenerator<FixedG1GenerationalCellCycleModel, 3> potts_cells_generator; // Set the cell cylce for the seeds 
		potts_cells_generator.GenerateBasicRandom(potts_cells, p_potts_mesh->GetNumElements(), p_diff_type); // Assigning the cell type and the cell cycle 

        TRACE("Generate Potts cell population");
		// Create cell population linking potts mesh and cells
		PottsBasedCellPopulation<3> potts_population(*p_potts_mesh, potts_cells);

		// // Add in all the cell writers 
		potts_population.AddCellWriter<CellIdWriter>();
        potts_population.SetNumSweepsPerTimestep(1);

		// Set up Potts simulation
		OnLatticeSimulation<3> potts_simulator(potts_population);
		potts_simulator.SetDt(0.1);
		potts_simulator.SetSamplingTimestepMultiple(2);
		
        // // // Reduce temperature to avoid randomly removing cells before they grow (parameters will need proper fitting)
		potts_population.SetTemperature(1e-9);

        // /*
        //  * Simulate cell seeding from randomly placed cell markers
        //  */

		MAKE_PTR(ArbitraryVolumeConstraintPottsUpdateRule<3>, p_volume_constraint_update_rule);
		p_volume_constraint_update_rule->SetMatureCellTargetVolume(CellArea);
		p_volume_constraint_update_rule->SetDeformationEnergyParameter(400);//(4000); //r(4000000000);
		potts_simulator.AddUpdateRule(p_volume_constraint_update_rule);

		MAKE_PTR(ArbitraryAdhesionPottsUpdateRule<3>, p_adhesion_update_rule);
		p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02/1000);
		p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16/1000);
		potts_simulator.AddUpdateRule(p_adhesion_update_rule);

		MAKE_PTR(ArbitrarySurfaceAreaConstraintPottsUpdateRule<3>, p_area_constraint_update_rule);
        p_area_constraint_update_rule->SetSurfaceAreaEnergyParameter(0.03e4);
		p_area_constraint_update_rule->SetTargetSurfaceArea(CellPerimeter); 
		potts_simulator.AddUpdateRule(p_area_constraint_update_rule);

		// std::string traction_file = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/FlowFolder/results/Extracted/surface-tractions.xtr";

		// Traction needs to be read in before the update, because the traction is reapplied at each lattive in the update step 
		// p_potts_mesh->TractionDataLoader(traction_file);
		p_potts_mesh->UpdatePottsNodeLocationFromDelaunay();
		double WallShearStress = 2e-3;//
		p_potts_mesh->SetConstantWallShearStress(WallShearStress);



		/*
        *  Add mechanotaxis to Potts simulation for timesteps to follow
        */

         MAKE_PTR(MechanotaxisPottsUpdateRule<3>, p_mechanotaxis_update_rule);
        // p_mechanotaxis_update_rule->SetTractionCorrelationParameter(30e2); //Parameter works for merks update rule 
		// p_mechanotaxis_update_rule->SetTractionCorrelationParameter(20e3);// Works for James update rule 
		 p_mechanotaxis_update_rule->SetTractionCorrelationParameter(10e5);// Works for ---  update rule 
        potts_simulator.AddUpdateRule(p_mechanotaxis_update_rule);




		/*
		* ADDITIONAL UPDATE RULES 
		*/

        MAKE_PTR(PottsCellPropertiesModifier<3>, p_modifier);
        potts_simulator.AddSimulationModifier(p_modifier);

        // // potts_simulator.SetStartTime(Confluence_time)
        potts_simulator.SetOutputDirectory("PottsPresentationMovies/ShearMigration");
        potts_simulator.SetEndTime(200);

        potts_simulator.Solve();




        // MAKE_PTR(ArbitraryCurvatureConstraintPottsUpdateRule<3>, p_Curvature_constraint_update_rule);
        // p_Curvature_constraint_update_rule->SetCurvatureEnergyParameter(0.016e18); //(0.016e18); // lambda_p
        // p_Curvature_constraint_update_rule->SetTargetCurvature(20993);// 7993) // p_area_constraint_update_rule->SetTargetSurfaceArea (150e-6); // (150e-6); // P (is 150e-12 in paper but value was found wrong and corrected for 150e-6)
        // potts_simulator.AddUpdateRule(p_Curvature_constraint_update_rule);

       



    	// // Need to have a look at what this is doing, to asserts need to be DIM ==2 
    	//         MAKE_PTR(AspectRatioConstraintPottsUpdateRule<3>, p_aspect_ratio_update_rule);
    	//         p_aspect_ratio_update_rule->SetTargetAspectRatio(4);
    	//         p_aspect_ratio_update_rule->SetAspectRatioEnergyParameter(1e-7);
    	//         potts_simulator.AddUpdateRule(p_aspect_ratio_update_rule);


        
	
	
	
	
		
		
		/*
			*  Add mechanotaxis to Potts simulation for timesteps to follow
			*/
			// MAKE_PTR(MechanotaxisPottsUpdateRule<3>, p_mechanotaxis_update_rule);
			// p_mechanotaxis_update_rule->SetTractionCorrelationParameter(0.01e20);//1e-8);
			// potts_simulator.AddUpdateRule(p_mechanotaxis_update_rule);

    


	
		// /*
        // -----------------------------
        // Inital Seeding
        // ----------------------------
        // */
		// double Confluence_time = 10;
        // potts_simulator.SetOutputDirectory("TestPottsOnSquare_Potts_Setup");
		// potts_simulator.Solve();
		// potts_simulator.SetOutputDirectory("TestPottsOnVesselGrowth_Potts");
        // potts_simulator.SetSamplingTimestepMultiple(1);
		// unsigned num_iterations = 60;




	}
};

#endif /*TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_*/


