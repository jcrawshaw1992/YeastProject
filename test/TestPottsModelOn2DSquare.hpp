
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
// #include "FixedRegionBoundaryCondition.hpp"
// #include "PlaneBasedCellKiller.hpp"
// #include "AbstractForce.hpp"
// #include "ChemotaxisPottsUpdateRule.hpp"
#include "Debug.hpp"

#include "UblasCustomFunctions.hpp"

#include "PottsMeshFromMutableMeshGeneratorJess.hpp"

#include "PottsArbitrarySurfaceIn3DMesh.hpp"

// #include "HoneycombMeshGenerator.hpp"


#include "Honeycomb3DCylinderMeshGenerator.hpp"

#include "ArbitrarySurfaceAreaConstraintPottsUpdateRule.hpp"
#include "ArbitraryVolumeConstraintPottsUpdateRule.hpp"


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

// Regualr square --hexagonal lattie
            // double N_D = 50;
            // double N_Z = 50;
            // Honeycomb3DMeshGenerator generator(N_D, N_Z, 1,1);
            // MutableMesh<2, 3>* mutable_mesh = generator.GetMesh();
            // mutable_mesh->Translate( unit_vector<double>(3, 2));



// Randon cylinder
		//   unsigned N = 581;
        //   std::stringstream out;
        //     out << N;
        //     std::string mesh_size = out.str();
        //     std::string mesh_file = "projects/EMBC2018/test/data/cyl_" + mesh_size + "_nodes.vtu";
        //     // std::string output_directory = "CylinderValidation/Random/" + mesh_size;

        //     // This data file is in mm
        //     VtkMeshReader<2,3> mesh_reader(mesh_file);
        //     MutableMesh<2,3> mutable_mesh;
        //     mutable_mesh.ConstructFromMeshReader(mesh_reader);
        //     double scaling = 1e-3;  // so distances are in m


// Regular cylinder 


// scale everything up 

double scale = 1e1;
PRINT_VARIABLE(scale);

			double NumberOfCells =4;
			unsigned N_D = 10;
			unsigned N_Z = 20;//N_D*1.5;
			double Radius = 5e-4 *scale;
			double Length = 30e-4 *scale;//2*96.6e-4 * scale; //12e-3;
			double trans = -Length/2;
			
			

			// MutableMesh<2, 3>* p_mesh = p_mesh_base;
			Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
			// Honeycomb3DCylinderMeshGenerator generator(N_D[N_D_index], N_Z, 1.5e-3, Length);
			MutableMesh<2, 3>* mutable_mesh = generator.GetMesh();
			mutable_mesh->Translate(trans * unit_vector<double>(3, 2));


        /*
		 * Setup Potts simulation
		 */
		
        // PottsMeshFromMutableMeshGeneratorJess<3> potts_generator(mutable_mesh);
        PottsMeshFromMutableMeshGeneratorJess<3> potts_generator(* mutable_mesh);
		PottsArbitrarySurfaceIn3DMesh<3>* p_potts_mesh = potts_generator.GetMesh();
	
		

		unsigned max_num_cells =5;  
		unsigned num_cells_seeded = 5;
		double target_volume = ComputeMeshTotalArea(* mutable_mesh) / max_num_cells;
		
		// Need to write something that stops a cell being doubled over on a node
		// double target_volume = ComputeMeshTotalArea(mutable_mesh) / max_num_cells;
 
		for (unsigned i=0; i<num_cells_seeded; i++) 
		{
			// Pick a node randomly avoiding boundary nodes
			unsigned num_node;
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


		// Add in all the cell writers 
		potts_population.AddCellWriter<CellIdWriter>();
		potts_population.AddCellWriter<CellAspectRatioWriter>();
		potts_population.AddCellWriter<CellMajorAxisAngleWriter>();

        potts_population.SetNumSweepsPerTimestep(1);

		// Set up Potts simulation
		OnLatticeSimulation<3> potts_simulator(potts_population);
		// potts_simulator.SetOutputDirectory("TestPottsOn2DSquare_Potts_withoutCurvature");
		potts_simulator.SetDt(0.1);
		potts_simulator.SetSamplingTimestepMultiple(1);
		
        // // Reduce temperature to avoid randomly removing cells before they grow (parameters will need proper fitting)
		potts_population.SetTemperature(1e-9);

        /*
         * Simulate cell seeding from randomly placed cell markers
        //  */

		MAKE_PTR(ArbitraryVolumeConstraintPottsUpdateRule<3>, p_volume_constraint_update_rule);
		PRINT_VARIABLE(target_volume);
		target_volume = 698e-4 *scale/100;
		p_volume_constraint_update_rule->SetMatureCellTargetVolume(target_volume);
		p_volume_constraint_update_rule->SetDeformationEnergyParameter(4000000);//4000000*1000);
		potts_simulator.AddUpdateRule(p_volume_constraint_update_rule);

		MAKE_PTR(ArbitraryAdhesionPottsUpdateRule<3>, p_adhesion_update_rule);
		p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02 / 1000);
		p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16 / 1000);
		potts_simulator.AddUpdateRule(p_adhesion_update_rule);

		MAKE_PTR(ArbitrarySurfaceAreaConstraintPottsUpdateRule<3>, p_area_constraint_update_rule);
        p_area_constraint_update_rule->SetSurfaceAreaEnergyParameter(0.016e1); //(0.016e20); //(0.016e18); // lambda_p
		double TargetCellPerimeter = 221e-4* scale/10000;// endothelial cell perimeter 2221e-4
		p_area_constraint_update_rule->SetTargetSurfaceArea(TargetCellPerimeter); // p_area_constraint_update_rule->SetTargetSurfaceArea (150e-6); // (150e-6); // P (is 150e-12 in paper but value was found wrong and corrected for 150e-6)
		potts_simulator.AddUpdateRule(p_area_constraint_update_rule);


		MAKE_PTR(PottsCellPropertiesModifier<3>, p_modifier);
        potts_simulator.AddSimulationModifier(p_modifier);


		/*
		* ADDITIONAL UPDATE RULES 
		*/

        // MAKE_PTR(ArbitraryCurvatureConstraintPottsUpdateRule<3>, p_Curvature_constraint_update_rule);
        // p_Curvature_constraint_update_rule->SetCurvatureEnergyParameter(0.016e18); //(0.016e18); // lambda_p
		// p_Curvature_constraint_update_rule->SetTargetCurvature(20993);// 7993) // p_area_constraint_update_rule->SetTargetSurfaceArea (150e-6); // (150e-6); // P (is 150e-12 in paper but value was found wrong and corrected for 150e-6)
		// potts_simulator.AddUpdateRule(p_Curvature_constraint_update_rule);

		/*
		*  Add mechanotaxis to Potts simulation for timesteps to follow
		*/

        // MAKE_PTR(MechanotaxisPottsUpdateRule<3>, p_mechanotaxis_update_rule);
        // p_mechanotaxis_update_rule->SetTractionCorrelationParameter(0.016e25);//1e-8);
        // potts_simulator.AddUpdateRule(p_mechanotaxis_update_rule);



// // Need to have a look at what this is doing, to asserts need to be DIM ==2 
//         MAKE_PTR(AspectRatioConstraintPottsUpdateRule<3>, p_aspect_ratio_update_rule);
//         p_aspect_ratio_update_rule->SetTargetAspectRatio(4);
//         p_aspect_ratio_update_rule->SetAspectRatioEnergyParameter(1e-7);
//         potts_simulator.AddUpdateRule(p_aspect_ratio_update_rule);


    

        std::string traction_file = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/CollpasedInitalState/results/Extracted/surface-tractions.xtr";


		unsigned Confluence_time = 10;
        potts_simulator.SetEndTime(Confluence_time);
        potts_simulator.SetOutputDirectory("TestPottsOn2dSquare_Potts_Setup");

		// Traction needs to be read in before the update, because the traction is reapplied at each lattive in the update step 
		p_potts_mesh->TractionDataLoader(traction_file);// "projects/VascularRemodelling/test/data/straight_vessel-surface-tractions.xtr");
		p_potts_mesh->UpdatePottsNodeLocationFromDelaunay();
		p_potts_mesh->GetRadius(Radius);

		// TRACE("solve");
        // potts_simulator.Solve();

		/*
		*  Add mechanotaxis to Potts simulation for timesteps to follow
		*/




        MAKE_PTR(MechanotaxisPottsUpdateRule<3>, p_mechanotaxis_update_rule);
        p_mechanotaxis_update_rule->SetTractionCorrelationParameter(0.01e20);//1e-8);
        potts_simulator.AddUpdateRule(p_mechanotaxis_update_rule);
		


		// potts_simulator.SetStartTime(Confluence_time)
		potts_simulator.SetEndTime(Confluence_time);
        // potts_simulator.SetOutputDirectory("TestPottsOn2dSquare_Potts_Migration");
		potts_simulator.Solve();


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




		// /*
		//  * Setup growth problem
		//  */
        // // Create representation of cells required for the mesh based simulation
        // TRACE("NewCells");
        // std::vector<CellPtr> vessel_cells;
        // CellsGenerator<FixedG1GenerationalCellCycleModel, 3> vessel_cells_generator;
        // vessel_cells_generator.GenerateBasicRandom(vessel_cells, mutable_mesh->GetNumNodes(), p_diff_type);

        // // Create a cell population linking mechanics mesh and cells
        // MeshBasedCellPopulation<2,3> vessel_cell_population(*mutable_mesh, vessel_cells);
        // vessel_cell_population.SetWriteVtkAsPoints(true);
        // vessel_cell_population.SetOutputMeshInVtk(true);

		//  // Set up cell-based simulation
        // TRACE("NewSimulation");
        // OffLatticeSimulation<2,3> vessel_simulator(vessel_cell_population);
        // vessel_simulator.SetOutputDirectory("TestPottsOnVesselGrowth_Growth");
        // vessel_simulator.SetDt(0.002);
        // vessel_simulator.SetSamplingTimestepMultiple(500);
        // vessel_simulator.SetUpdateCellPopulationRule(false); // No remeshing.

		 // This loads tractions into the lattice site objects so it can be used for mechanotaxis.
        // @todo the force law for normal stress (a few lines above) should probably be reimplemented reading from file too.
        // TractionDataLoader potts_traction_loader("projects/VascularRemodelling/test/data/straight_vessel-surface-tractions.xtr");
        // potts_traction_loader.UpdateLatticeSiteData(p_potts_mesh);

		 /*       
        -----------------------------
        Tractionforce 
        ----------------------------
        */

        // // std::string working_directory = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/DeformingCylindricalBasicHetroWallTesting/";
        // std::string working_directory = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/FluidFlowInPipe/";
		//  std::string output_directory = "AppropriateInitalConditions";//Shrinking"; // + Parameters + "/";

        


        //  // Create an Applied Force modifier to couple to Flow
        // std::string traction_file = working_directory + "results/Extracted/surface-tractions.xtr";
        // boost::shared_ptr<AppliedForceModifier<2, 3> > p_force_modifier(new AppliedForceModifier<2, 3>());

        // p_force_modifier->SetResetTractionsOnCells(true, traction_file);
		// p_force_modifier->SetupVessel(vessel_cell_population, output_directory);
        // vessel_simulator.AddSimulationModifier(p_force_modifier);

