
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
#include "FixedRegionBoundaryCondition.hpp"
#include "AbstractForce.hpp"
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


#include "projects/VascularRemodelling/src/MembraneForces/MembraneShearForce.hpp"
#include "projects/VascularRemodelling/src/MembraneForces/MembraneStiffnessForce.hpp"
#include "projects/VascularRemodelling/src/MembraneForces/MembraneSurfaceForce.hpp"




// Needed for NodesOnlyMesh.
#include "PetscSetupAndFinalize.hpp"

#include <iostream>
#include <stdio.h>


class RadialForce : public AbstractForce<2,3>
{
private:

    double mStrength;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<2,3> >(*this);
        archive & mStrength;
    }

public:
    RadialForce(double strength=1.0)
        : AbstractForce<2,3>(),
          mStrength(strength)
    {
        assert(mStrength > 0.0);
    }

    void AddForceContribution(AbstractCellPopulation<2,3>& rCellPopulation)
    {
    	// Helper variables
    	MeshBasedCellPopulation<2,3>* p_cell_population = static_cast<MeshBasedCellPopulation<2,3>*>(&rCellPopulation);


    	// Calculate midpoint
    	c_vector<double,3> centroid = zero_vector<double>(3);
   		for (AbstractCellPopulation<2,3>::Iterator cell_iter = rCellPopulation.Begin();
					 cell_iter != rCellPopulation.End();
					 ++cell_iter)
		{
			unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
			centroid += rCellPopulation.GetNode(node_index)->rGetLocation();
		}
   		centroid /= rCellPopulation.GetNumRealCells();

        for (AbstractCellPopulation<2,3>::Iterator cell_iter = rCellPopulation.Begin();
                     cell_iter != rCellPopulation.End();
                     ++cell_iter)
		{
        	unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        	Node<3>* p_node = rCellPopulation.GetNode(node_index);
        	c_vector<double,3> cell_location = p_node->rGetLocation() - centroid;
        	cell_location(2) = 0.0;
        	c_vector<double, 3> force = zero_vector<double>(3);

			//Calculate cell normal (average of element normals)
			c_vector<double,3> normal = zero_vector<double>(3);

			std::set<unsigned>&  containing_elements = p_node->rGetContainingElementIndices();
			assert(containing_elements.size()>0);
			for (std::set<unsigned>::iterator iter = containing_elements.begin();
				 iter != containing_elements.end();
				 ++iter)
			{

				normal += p_cell_population->rGetMesh().GetElement(*iter)->CalculateNormal();
			}
			normal /= norm_2(normal);

			double cell_area = rCellPopulation.GetVolumeOfCell(*cell_iter);
			// force = mStrength * cell_area * normal; // cell_location / norm_2(cell_location);

			force = -mStrength * normal; // cell_location / norm_2(cell_location);
			cell_iter->GetCellData()->SetItem("area", cell_area);

        	rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(force);

            cell_iter->GetCellData()->SetItem("applied_force_x", force[0]);
            cell_iter->GetCellData()->SetItem("applied_force_y", force[1]);
            cell_iter->GetCellData()->SetItem("applied_force_z", force[2]);
            cell_iter->GetCellData()->SetItem("applied_force_mag", norm_2(force));

            cell_iter->GetCellData()->SetItem("norm_x", normal[0]);
			cell_iter->GetCellData()->SetItem("norm_y", normal[1]);
			cell_iter->GetCellData()->SetItem("norm_z", normal[2]);

//
//
//            cell_iter->GetCellData()->SetItem("radius", norm_2(cell_location));
		}
    }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<Strength>" << mStrength << "</Strength>\n";
        AbstractForce<2,3>::OutputForceParameters(rParamsFile);
    }
};


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
double scale = 1e-4; 
            double N_D = 50;
            double N_Z = 50;
            Honeycomb3DMeshGenerator generator(N_D, N_Z, 1,1);
            MutableMesh<2, 3>* mutable_mesh = generator.GetMesh();
            mutable_mesh->Translate( unit_vector<double>(3, 2));
          mutable_mesh->Scale(5e-4 *scale, 30e-4 *scale, 1.0); // so distances are back in original scal
        //    MutableMesh<2,3> mutable_mesh = * mutable_mesh_pointer;




// Random cylinder
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


// // scale everything up 

// double scale = 1e1;
// PRINT_VARIABLE(scale);

// 			double NumberOfCells =4;
// 			unsigned N_D = 10;
// 			unsigned N_Z = 20;//N_D*1.5;
// 			double Radius = 5e-4 *scale;
// 			double Length = 30e-4 *scale;//2*96.6e-4 * scale; //12e-3;
// 			double trans = -Length/2;
			
			

// 			// MutableMesh<2, 3>* p_mesh = p_mesh_base;
// 			Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
// 			// Honeycomb3DCylinderMeshGenerator generator(N_D[N_D_index], N_Z, 1.5e-3, Length);
// 			MutableMesh<2, 3>* mutable_mesh = generator.GetMesh();
// 			mutable_mesh->Translate(trans * unit_vector<double>(3, 2));



// Birfucation
		//   unsigned N = 581;
        //   std::stringstream out;
        //     out << N;
        //     std::string mesh_size = out.str();
        //     std::string mesh_file = "projects/VascularRemodelling/test/data/bifurcation_cut/config.vtu";
        //     std::string output_directory = "BifucationPotts/" ;

        //     // This data file is in mm??
        //     VtkMeshReader<2,3> mesh_reader(mesh_file);
        //     MutableMesh<2,3> mutable_mesh;
        //     mutable_mesh.ConstructFromMeshReader(mesh_reader);
             // so distances are in m
		// 	mutable_mesh.Scale(1.0 *scale , 1.0 *scale , 1.0 * scale ); // so distances are back in original scal

				double Radius = 5e-4 *scale;
				double Length = 30e-4 *scale;//2*96.6e-4 * scale; //12e-3;




        /*
		 * Setup Potts simulation
		 */
		
        // PottsMeshFromMutableMeshGeneratorJess<3> potts_generator(mutable_mesh);
        PottsMeshFromMutableMeshGeneratorJess<2> potts_generator(* mutable_mesh);
		PottsArbitrarySurfaceIn3DMesh<2>* p_potts_mesh = potts_generator.GetMesh();
	
		

		// unsigned max_num_cells =6;  
		// unsigned num_cells_seeded = 6;
		// double CellScale =1e2;


		// double LongAxis = 10/2* 1e-6 * CellScale ;//*scale ;
		// double ShortAxis = 5/2 *1e-6 * CellScale ;//*scale ;
		// double CellArea = M_PI * LongAxis  * ShortAxis ;

		// double AxisRatio = (LongAxis -ShortAxis)* (LongAxis -ShortAxis) /(LongAxis +ShortAxis)*(LongAxis +ShortAxis); 
		// double CellPerimeter = M_PI * (LongAxis +ShortAxis)* (3 *AxisRatio * 1/(sqrt(-3*AxisRatio+4) +10) +1) ;   

		// // PRINT_3_VARIABLES(CellArea ,CellPerimeter );

		// for (unsigned i=0; i<num_cells_seeded; i++) 
		// {
		// 	// Pick a node randomly avoiding boundary nodes
		// 	unsigned num_node;
		// 	do
		// 	{
		// 		num_node = RandomNumberGenerator::Instance()->randMod(p_potts_mesh->GetNumNodes());  // Selects a random node 

		// 	}
		// 	while(p_potts_mesh->GetNode(num_node)->IsBoundaryNode()); // XXX not sure what this does, Im not sure it does anything 
		// 	PRINT_VARIABLE(i);

		// 	std::vector<Node<3>*> element_nodes;
		// 	element_nodes.push_back(p_potts_mesh->GetNode(num_node));
		// 	p_potts_mesh->AddElement(new PottsElement<3>(i, element_nodes));
		// 	std::cout << "Seeded at "<< num_node << std::endl;
		// }

        // // Randomly place the cell markers for seeding
		// std::vector<CellPtr> potts_cells;
		// MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type); // Set the cell type for the seeds 
		// CellsGenerator<FixedG1GenerationalCellCycleModel, 3> potts_cells_generator; // Set the cell cylce for the seeds 
		// potts_cells_generator.GenerateBasicRandom(potts_cells, p_potts_mesh->GetNumElements(), p_diff_type); // Assigning the cell type and the cell cycle 

        // TRACE("Generate Potts cell population");
		// // Create cell population linking potts mesh and cells
		// PottsBasedCellPopulation<3> potts_population(*p_potts_mesh, potts_cells);


		// // Add in all the cell writers 
		// potts_population.AddCellWriter<CellIdWriter>();

        // potts_population.SetNumSweepsPerTimestep(1);

		// // Set up Potts simulation
		// OnLatticeSimulation<3> potts_simulator(potts_population);
		// potts_simulator.SetDt(0.1);
		// potts_simulator.SetSamplingTimestepMultiple(2);
		
        // // // Reduce temperature to avoid randomly removing cells before they grow (parameters will need proper fitting)
		// potts_population.SetTemperature(1e-9);

        /*
         * Simulate cell seeding from randomly placed cell markers
         */

		// MAKE_PTR(ArbitraryVolumeConstraintPottsUpdateRule<3>, p_volume_constraint_update_rule);
		// p_volume_constraint_update_rule->SetMatureCellTargetVolume(CellArea);
		// p_volume_constraint_update_rule->SetDeformationEnergyParameter(4000000000);
		// potts_simulator.AddUpdateRule(p_volume_constraint_update_rule);

		// MAKE_PTR(ArbitraryAdhesionPottsUpdateRule<3>, p_adhesion_update_rule);
		// p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02 / 1000);
		// p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16 / 1000);
		// potts_simulator.AddUpdateRule(p_adhesion_update_rule);

		// MAKE_PTR(ArbitrarySurfaceAreaConstraintPottsUpdateRule<3>, p_area_constraint_update_rule);
        // p_area_constraint_update_rule->SetSurfaceAreaEnergyParameter(0.016e2);
		// p_area_constraint_update_rule->SetTargetSurfaceArea(CellPerimeter/20); 
		// potts_simulator.AddUpdateRule(p_area_constraint_update_rule);


		// /*
		// * ADDITIONAL UPDATE RULES 
		// */

		// 			// MAKE_PTR(PottsCellPropertiesModifier<3>, p_modifier);
		// 			// potts_simulator.AddSimulationModifier(p_modifier);



		// 			// MAKE_PTR(ArbitraryCurvatureConstraintPottsUpdateRule<3>, p_Curvature_constraint_update_rule);
		// 			// p_Curvature_constraint_update_rule->SetCurvatureEnergyParameter(0.016e18); //(0.016e18); // lambda_p
		// 			// p_Curvature_constraint_update_rule->SetTargetCurvature(20993);// 7993) // p_area_constraint_update_rule->SetTargetSurfaceArea (150e-6); // (150e-6); // P (is 150e-12 in paper but value was found wrong and corrected for 150e-6)
		// 			// potts_simulator.AddUpdateRule(p_Curvature_constraint_update_rule);

		// 			/*
		// 			*  Add mechanotaxis to Potts simulation for timesteps to follow
		// 			*/

		// 			// MAKE_PTR(MechanotaxisPottsUpdateRule<3>, p_mechanotaxis_update_rule);
		// 			// p_mechanotaxis_update_rule->SetTractionCorrelationParameter(0.016e25);//1e-8);
		// 			// potts_simulator.AddUpdateRule(p_mechanotaxis_update_rule);



		// 	// // Need to have a look at what this is doing, to asserts need to be DIM ==2 
		// 	//         MAKE_PTR(AspectRatioConstraintPottsUpdateRule<3>, p_aspect_ratio_update_rule);
		// 	//         p_aspect_ratio_update_rule->SetTargetAspectRatio(4);
		// 	//         p_aspect_ratio_update_rule->SetAspectRatioEnergyParameter(1e-7);
		// 	//         potts_simulator.AddUpdateRule(p_aspect_ratio_update_rule);


        
	
	
	
	
		
		// std::string traction_file = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/FlowFolder/results/Extracted/surface-tractions.xtr";
        // potts_simulator.SetOutputDirectory("TestPottsOn2dSquare_Potts_Setup");

		// // Traction needs to be read in before the update, because the traction is reapplied at each lattive in the update step 
		// p_potts_mesh->TractionDataLoader(traction_file);
		// p_potts_mesh->UpdatePottsNodeLocationFromDelaunay();
		// p_potts_mesh->GetRadius(Radius);
		// /*
		// 	*  Add mechanotaxis to Potts simulation for timesteps to follow
		// 	*/
		// 	// MAKE_PTR(MechanotaxisPottsUpdateRule<3>, p_mechanotaxis_update_rule);
		// 	// p_mechanotaxis_update_rule->SetTractionCorrelationParameter(0.01e20);//1e-8);
		// 	// potts_simulator.AddUpdateRule(p_mechanotaxis_update_rule);

		// 	// potts_simulator.SetStartTime(Confluence_time)
		// // potts_simulator.SetEndTime(Confluence_time);
	
		// // potts_simulator.Solve();



	
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



		// // for(unsigned iter_num=0; iter_num<num_iterations; iter_num++)
		// {
		// 	 std::cout << "Iteration  " << iter_num << std::endl;
		//     /*
		//      * Do one timestep of cell movement
		//      */

		//     std::cout << "Potts timestep" << std::endl;

		// 	// Destroy the simulation time that has been left from the setup step or the previous iteration 
        //     SimulationTime::Instance()->Destroy();
        //     SimulationTime::Instance()->SetStartTime(potts_current_time);
        //     if (potts_current_time > 0)
        //     {
        //         double dt = simulator.GetDt();
        //         unsigned num_time_steps = (unsigned) (potts_timestep/dt+0.5);
        //         SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(potts_current_time + potts_timestep, num_time_steps);
        //     }

        //     potts_simulator.SetEndTime(potts_current_time + potts_timestep);
        //     potts_simulator.Solve();
        //     potts_current_time += potts_timestep;

        //     /*
        //      * Do one timestep of vessel mechanics
        //      */
        //     std::cout << "Mechanics timestep" << std::endl;

		// 	SimulationTime::Instance()->Destroy();
        //     SimulationTime::Instance()->SetStartTime(growth_current_time);
        //     if (growth_current_time > 0)
        //     {
        //         double dt = simulator.GetDt();
        //         unsigned num_time_steps = (unsigned) (growth_timestep/dt+0.5);
        //         SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(growth_current_time + growth_timestep, num_time_steps);
        //      }
             
        //     simulator.SetEndTime(growth_current_time + growth_timestep);

		// 	// there is a problem here with the solve???
        //     simulator.Solve();

        //     // // // Update location of nodes in potts mesh from recently computed deformation
        //     p_potts_mesh->UpdatePottsNodeLocationFromDelaunay();

        //     growth_current_time += growth_timestep;
		// }

// Now need a step where the mesh deforms and we go again 

	}
};

#endif /*TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_*/


