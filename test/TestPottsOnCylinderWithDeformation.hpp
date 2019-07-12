
#ifndef TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_
#define TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <cstdio>
#include <ctime>


// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellsGenerator.hpp"
#include "OnLatticeSimulation.hpp"
#include "OffLatticeSimulation.hpp"
// #include "FixedG1GenerationalCellCycleModel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "CellIdWriter.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "ArbitraryVolumeConstraintPottsUpdateRule.hpp"
#include "ArbitraryAdhesionPottsUpdateRule.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PottsMeshGenerator.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "Warnings.hpp"
#include "LogFile.hpp"
#include "SmartPointers.hpp"
// #include "MechanotaxisPottsUpdateRule.hpp"
#include "AspectRatioConstraintPottsUpdateRule.hpp"
#include "PottsMeshFromMutableMeshGenerator.hpp"
#include "VtkMeshReader.hpp"
#include "RandomNumberGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "AbstractForce.hpp"
#include "ChemotaxisPottsUpdateRule.hpp"
#include "Debug.hpp"

#include "UblasCustomFunctions.hpp"

#include "PottsMeshFromMutableMeshGeneratorJess.hpp"

#include "PottsArbitrarySurfaceIn3DMesh.hpp"



#include "Honeycomb3DCylinderMeshGenerator.hpp"

#include "/Users/jcrawshaw/Documents/Chaste/projects/Jess/src/DeformableModel/MembraneSurfaceForceCylinder.hpp"
#include "/Users/jcrawshaw/Documents/Chaste/projects/Jess/src/DeformableModel/MembraneStiffnessForceCylinder.hpp"
#include "/Users/jcrawshaw/Documents/Chaste/projects/Jess/src/DeformableModel/MembraneShearForceCylinder.hpp"


// #include "TractionDataLoader.hpp"


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

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(RadialForce)
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(RadialForce)

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
        VtkMeshReader<2,3> mesh_reader("/Users/jcrawshaw/Documents/ChasteWorkingDirectory/CollpasedInitalState/SetUpData/config.vtu");
        MutableMesh<2,3> mutable_mesh;
        mutable_mesh.ConstructFromMeshReader(mesh_reader);
        double scaling = 1e-3;  // so distances are in m
		mutable_mesh.Scale(scaling,scaling,scaling);
		TRACE("A");

        /*
		 * Setup Potts simulation
		 */
		
        PottsMeshFromMutableMeshGeneratorJess<3> potts_generator(mutable_mesh);
		PottsArbitrarySurfaceIn3DMesh<3>* p_potts_mesh = potts_generator.GetMesh();
		TS_ASSERT_EQUALS(p_potts_mesh->GetNumNodes(), mutable_mesh.GetNumNodes());// XXX Why equate? they are not equal, and it doesnt break everything
        // PRINT_2_VARIABLES(p_potts_mesh->GetNumNodes(), mutable_mesh.GetNumNodes());
		

		unsigned max_num_cells = 10;  
		unsigned num_cells_seeded = 5;
		double target_volume = ComputeMeshTotalArea(mutable_mesh) / max_num_cells;
 
		for (unsigned i=0; i<num_cells_seeded; i++)  // XXX how does this give cell locations and lattice groups?? 
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

		TS_ASSERT_EQUALS(p_potts_mesh->GetNumElements(), num_cells_seeded);
        PRINT_2_VARIABLES(p_potts_mesh->GetNumElements(),num_cells_seeded);
		
	

        // Randomly place the cell markers for seeding
		std::vector<CellPtr> potts_cells;
		MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type); // Set the cell type for the seeds 
		CellsGenerator<FixedG1GenerationalCellCycleModel, 3> potts_cells_generator; // Set the cell cylce for the seeds 
		potts_cells_generator.GenerateBasicRandom(potts_cells, p_potts_mesh->GetNumElements(), p_diff_type); // Assigning the cell type and the cell cycle 

        TRACE("GenerateBasicRandom");
		// Create cell population linking potts mesh and cells
		PottsBasedCellPopulation<3> potts_population(*p_potts_mesh, potts_cells);
		potts_population.AddCellWriter<CellIdWriter>();
        potts_population.SetNumSweepsPerTimestep(1);
		// Set up Potts simulation
		OnLatticeSimulation<3> potts_simulator(potts_population);
		potts_simulator.SetOutputDirectory("TestPottsOnVesselGrowth_Potts");
		potts_simulator.SetDt(0.1);
		potts_simulator.SetSamplingTimestepMultiple(1);


         /*
		 * Setup the Mesh simulator 
		 */
		
        // Create representation of cells required for the simulator

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mutable_mesh.GetNumNodes(), p_diff_type);

        // Create a cell population. Links the Mesh and the cells
        MeshBasedCellPopulation<2,3> cell_population(mutable_mesh, cells);

        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        // cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.CalculateRestLengths();
		// cell_population.SetOutputCellIdData(true);
        //cell_population.SetDampingConstantNormal(0.1);//e5);
        //cell_population.SetDampingConstantMutant(0.1);//e5);


        // Set up cell-based simulation
        TRACE("NewSimulation");
        OffLatticeSimulation<2,3> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsOnVesselGrowth_Growth");
        simulator.SetDt(0.002);
        simulator.SetSamplingTimestepMultiple(500);
        simulator.SetUpdateCellPopulationRule(false); // No remeshing.
        



        // Add a force law for normal stress, i.e. blood pressure only
        double pressure = 0.00002;//1.0666e4; // to match 80mmhg
        MAKE_PTR_ARGS(RadialForce, p_radial_force, (pressure));
        simulator.AddForce(p_radial_force);

		/*
        -----------------------------
        Shearing Force 
        ----------------------------
        */


		double AreaDilationModulus = 1e-10; // membrane_constant[membrane_constant_index];
		double ElasticShearModulus = 4.6e-08; // membrane_constant[membrane_constant_index2];
		boost::shared_ptr<MembraneShearForceCylinder> p_shear_force(new MembraneShearForceCylinder());
		p_shear_force->SetupMembraneConfiguration(cell_population);
		p_shear_force->SetAreaDilationModulus(1e-8);
		p_shear_force->SetElasticShearModulus(5e-4);
		// p_shear_force->SetNc(N_D);
		simulator.AddForce(p_shear_force);


        //   //  -----------------------------
        //   //  Surface Area Force
        //   //  ----------------------------

            boost::shared_ptr<MembraneSurfaceForceCylinder> p_surface_force(new MembraneSurfaceForceCylinder());
            p_surface_force->SetupInitialAreas(cell_population);
            p_surface_force->SetMembraneStiffness(1e-13);
            simulator.AddForce(p_surface_force);

        //  //  -----------------------------
        //  //  Bending Force 
        //  //  ----------------------------
        
            // double Stiff_constant = 1e-15;
            // boost::shared_ptr<MembraneStiffnessForceCylinder> p_membrane_force(new MembraneStiffnessForceCylinder());
            // p_membrane_force->SetupInitialMembrane(mutable_mesh);
            // p_membrane_force->SetMembraneStiffness(Stiff_constant);
            // simulator.AddForce(p_membrane_force);




        
        //boundary conditions 
        boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Create_c_vector(0, 0, -0.006), Create_c_vector(0, 0, 1), 10));
        simulator.AddCellPopulationBoundaryCondition(p_condition_1);
        boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Create_c_vector(0, 0, 0.006), Create_c_vector(0, 0, -1), 10));
        simulator.AddCellPopulationBoundaryCondition(p_condition_2);



		// This loads tractions into the lattice site objects so it can be used for mechanotaxis.
        // @todo the force law for normal stress (a few lines above) should probably be reimplemented reading from file too.

        // TractionDataLoader potts_traction_loader("projects/VascularRemodelling/test/data/straight_vessel-surface-tractions.xtr");
        // potts_traction_loader.UpdateLatticeSiteData(p_potts_mesh);

		/*
		 * We need timestepping for two different simulations
		 */
        unsigned potts_timestep = 10;
        unsigned growth_timestep = 10;
        unsigned potts_current_time = 0;
        unsigned growth_current_time = 0;

        // Reduce temperature to avoid randomly removing cells before they grow (parameters will need proper fitting)
		potts_population.SetTemperature(1e-9);

        /*
         * Simulate cell seeding from randomly placed cell markers
         */
		MAKE_PTR(ArbitraryVolumeConstraintPottsUpdateRule<3>, p_volume_constraint_update_rule);
		PRINT_VARIABLE(target_volume);
		p_volume_constraint_update_rule->SetMatureCellTargetVolume(target_volume);
		p_volume_constraint_update_rule->SetDeformationEnergyParameter(4000000);
		potts_simulator.AddUpdateRule(p_volume_constraint_update_rule);

		MAKE_PTR(ArbitraryAdhesionPottsUpdateRule<3>, p_adhesion_update_rule);
		p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02 / 1000);
		p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16 / 1000);
		potts_simulator.AddUpdateRule(p_adhesion_update_rule);


// Need to have a look at what this is doing, to asserts need to be DIM ==2 
        // MAKE_PTR(AspectRatioConstraintPottsUpdateRule<3>, p_aspect_ratio_update_rule);
        // p_aspect_ratio_update_rule->SetTargetAspectRatio(4);
        // p_aspect_ratio_update_rule->SetAspectRatioEnergyParameter(1e-7);
        // potts_simulator.AddUpdateRule(p_aspect_ratio_update_rule);



		// Set up the inital seeding 
		unsigned Confluence_time = 1;
        potts_simulator.SetEndTime(potts_current_time + Confluence_time);
        potts_simulator.SetOutputDirectory("TestPottsOnVesselGrowth_Potts_Setup");
		potts_simulator.Solve();


		// Change the potts outdirectory now that the Potts population should have achieved confluence 
        potts_simulator.SetOutputDirectory("TestPottsOnVesselGrowth_Potts");

        // Add mechanotaxis to Potts simulation for timesteps to follow
        // MAKE_PTR(MechanotaxisPottsUpdateRule<3>, p_mechanotaxis_update_rule);
        // p_mechanotaxis_update_rule->SetTractionCorrelationParameter(1e-8);
        // potts_simulator.AddUpdateRule(p_mechanotaxis_update_rule);

		/*
		 * Now iterate between Cell Movement and Growth
		 */

	
		unsigned num_iterations = 10;
		for(unsigned iter_num=0; iter_num<num_iterations; iter_num++)
		{
			 std::cout << "Iteration  " << iter_num << std::endl;
		    /*
		     * Do one timestep of cell movement
		     */

		    std::cout << "Potts timestep" << std::endl;

			// Destroy the simulation time that has been left from the setup step or the previous iteration 
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(potts_current_time);
            if (potts_current_time > 0)
            {
                double dt = simulator.GetDt();
                unsigned num_time_steps = (unsigned) (potts_timestep/dt+0.5);
                SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(potts_current_time + potts_timestep, num_time_steps);
            }

            potts_simulator.SetEndTime(potts_current_time + potts_timestep);
            potts_simulator.Solve();
            potts_current_time += potts_timestep;

            /*
             * Do one timestep of vessel mechanics
             */
            std::cout << "Mechanics timestep" << std::endl;

			SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(growth_current_time);
            if (growth_current_time > 0)
            {
                double dt = simulator.GetDt();
                unsigned num_time_steps = (unsigned) (growth_timestep/dt+0.5);
                SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(growth_current_time + growth_timestep, num_time_steps);
            }
			growth_timestep =20;
            simulator.SetEndTime(growth_current_time + growth_timestep);
            simulator.Solve();

            // // Update location of nodes in potts mesh from recently computed deformation
            p_potts_mesh->UpdatePottsNodeLocationFromDelaunay();

            growth_current_time += growth_timestep;
		}









	}
};

#endif /*TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_*/
