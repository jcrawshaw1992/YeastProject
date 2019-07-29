/*

Copyright (c) 2005-2013, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
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

#include "PottsArbitraryMeshFromMutableMeshGenerator.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"
// #include "TractionDataLoader.hpp"


// Needed for NodesOnlyMesh.
#include "PetscSetupAndFinalize.hpp"


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
			force = mStrength * cell_area * normal; // cell_location / norm_2(cell_location);
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

class TestPottsVesselModel : public AbstractCellBasedTestSuite
{
public:

    void xTestPottsVesselMonolayer() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(100, 2, 10, 100, 2, 10, 1, 1, 1, false, true, true, false);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;

        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        // Make this pointer first as if we move it after creating the cell population the label numbers aren't tracked
        MAKE_PTR(CellLabel, p_label);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        //cell_population.SetOutputCellIdData(true);

        //cell_population.SetOutputCellMutationStates(true); // So outputs the labelled cells

        //cell_population.SetNumSweepsPerTimestep(5);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            if (RandomNumberGenerator::Instance()->ranf() < 0.5)
            {
                (*cell_iter)->AddCellProperty(p_label);
            }
        }

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsVesselMonolayer");
        simulator.SetDt(0.1);
        simulator.SetEndTime(10);


        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(100);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.2);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);

        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02);
        p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // MAKE_PTR(AspectRatioConstraintPottsUpdateRule<2>, p_aspect_ratio_update_rule);
        // p_aspect_ratio_update_rule->SetTargetAspectRatio(10);
        // p_aspect_ratio_update_rule->SetAspectRatioEnergyParameter(0.2);
        // simulator.AddUpdateRule(p_aspect_ratio_update_rule);

//        MAKE_PTR(ChemotaxisPottsUpdateRule<2>, p_chemotaxis_update_rule);
//        simulator.AddUpdateRule(p_chemotaxis_update_rule);

        // Run simulation
        simulator.Solve();
    }


    void xTestGeneratePottsMeshFromMutableMeshIn2D()
    {
    	TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
    	MutableMesh<2,2> mutable_mesh;
    	mutable_mesh.ConstructFromMeshReader(mesh_reader);

    	PottsMeshFromMutableMeshGenerator<2,2> potts_generator(mutable_mesh);

    	PottsMesh<2>* p_potts_mesh = potts_generator.GetMesh();

    	TS_ASSERT_EQUALS(p_potts_mesh->GetNumNodes(), mutable_mesh.GetNumNodes());

        // Make two triangular elements out of these nodes
        std::vector<Node<2>*> element_nodes;
        element_nodes.push_back(p_potts_mesh->GetNode(0));
//        element_nodes.push_back(p_potts_mesh->GetNode(1));
//        element_nodes.push_back(p_potts_mesh->GetNode(2));

        p_potts_mesh->AddElement(new PottsElement<2>(0, element_nodes));

        TS_ASSERT_EQUALS(p_potts_mesh->GetNumElements(), 1u);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_potts_mesh->GetNumElements(), p_diff_type);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_potts_mesh, cells);
        //cell_population.SetOutputCellIdData(true);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsOnTetMesh");
        simulator.SetDt(0.1);
        simulator.SetEndTime(10);


        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(100);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.2);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);

        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02);
        p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16);
        simulator.AddUpdateRule(p_adhesion_update_rule);

//        MAKE_PTR(AspectRatioConstraintPottsUpdateRule<2>, p_aspect_ratio_update_rule);
//        p_aspect_ratio_update_rule->SetTargetAspectRatio(10);
//        p_aspect_ratio_update_rule->SetDeformationEnergyParameter(0.2);
//        simulator.AddUpdateRule(p_aspect_ratio_update_rule);

//        MAKE_PTR(ChemotaxisPottsUpdateRule<2>, p_chemotaxis_update_rule);
//        simulator.AddUpdateRule(p_chemotaxis_update_rule);

        // Run simulation
        simulator.Solve();
    }

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

    void noTestPottsOnRealVessel()
    {
    	VtkMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/bifurcation_cut/config.vtu");
		MutableMesh<2,3> mutable_mesh;
		mutable_mesh.ConstructFromMeshReader(mesh_reader);
        //double scaling = 1e-3;  // so distances are in m
        //mutable_mesh.Scale(scaling,scaling,scaling);

        PottsArbitraryMeshFromMutableMeshGenerator<3> potts_generator(mutable_mesh);
        PottsArbitrarySurfaceIn3DMesh<3>* p_potts_mesh = potts_generator.GetMesh();

        // Original implementation where each Potts lattice is assumed a unit cube around each node in the mesh
    	// PottsMeshFromMutableMeshGenerator<2,3> potts_generator(mutable_mesh);
    	// PottsMesh<3>* p_potts_mesh = potts_generator.GetMesh();

        // TS_ASSERT_EQUALS(p_potts_mesh->GetNumNodes(), mutable_mesh.GetNumNodes());
    	PRINT_VARIABLE(mutable_mesh.GetNumNodes());
        PRINT_VARIABLE(p_potts_mesh->GetNumNodes());

        unsigned num_cells_seeded = 10;

        // Volume of a 2D mesh in 3D is area.
    	double target_volume = 1.2*ComputeMeshTotalArea(mutable_mesh) / num_cells_seeded;

    	// In the original implementation it was assumed that each node in the Potts element contributes a unit volume
    	// unsigned target_volume = p_potts_mesh->GetNumNodes()/max_num_cells;
    	for (unsigned i=0; i<num_cells_seeded; i++)
    	{
    		std::vector<Node<3>*> element_nodes;
    		element_nodes.push_back(p_potts_mesh->GetNode(RandomNumberGenerator::Instance()->randMod(p_potts_mesh->GetNumNodes())));
    		p_potts_mesh->AddElement(new PottsElement<3>(i, element_nodes));
    	}
        TS_ASSERT_EQUALS(p_potts_mesh->GetNumElements(), num_cells_seeded);

        // Randomly place the cell seeds
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_potts_mesh->GetNumElements(), p_diff_type);

        // Create cell population linking potts mesh and cells
        PottsBasedCellPopulation<3> cell_population(*p_potts_mesh, cells);
        //cell_population.SetOutputCellIdData(true);
		cell_population.AddCellWriter<CellIdWriter>();
        cell_population.SetNumSweepsPerTimestep(1);

        // Set up cell-based simulation
        OnLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsOnVessel");
        simulator.SetDt(0.01);

        /*
         * Simulate cell growth from randomly placed seeds
//          */
        MAKE_PTR(VolumeConstraintPottsUpdateRule<3>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(target_volume);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);

        MAKE_PTR(AdhesionPottsUpdateRule<3>, p_adhesion_update_rule);
        simulator.AddUpdateRule(p_adhesion_update_rule);
// //
        simulator.SetEndTime(1.0);
        simulator.Solve();

        /*
         * Simulate elongation
         */
//        MAKE_PTR(AspectRatioConstraintPottsUpdateRule<3>, p_aspect_ratio_update_rule);
//        simulator.AddUpdateRule(p_aspect_ratio_update_rule);

        // Add Killer and simulate Chemotaxis

        //Create a plane killer for inlet
        c_vector<double,3> point;
        point(0)=11.2770231292084;
        point(1)=36.0824073418958;
        point(2)=42.620856463622;
        c_vector<double,3> normal;
        normal(0)=0.119547495670065;
        normal(1)=-0.586263370820288;
        normal(2)=-0.801251306590791;
        MAKE_PTR_ARGS(PlaneBasedCellKiller<3>, p_cell_killer1, (&cell_population,point-1.0*normal,normal));
        simulator.AddCellKiller(p_cell_killer1);


        // Add Killer and simulate Chemotaxis

        //Create a plane killer for outlets

        point(0)=12.5554613754726;
        point(1)=46.8175974715283;
        point(2)=47.2958315257554;

        normal(0)=0.0212323101972157;
        normal(1)=0.996217721361657;
        normal(2)=0.0842581785269404;
        MAKE_PTR_ARGS(PlaneBasedCellKiller<3>, p_cell_killer, (&cell_population,point,normal));
        simulator.AddCellKiller(p_cell_killer);


        simulator.SetEndTime(0.5);
        simulator.Solve();

        /*
         * Simulate motion
         */
		MAKE_PTR(ChemotaxisPottsUpdateRule<3>, p_chemotaxis_update_rule);
//		p_chemotaxis_update_rule.SetShearStressFile("projects/VascularRemodelling/test/data/surface-tangentialprojectiontraction.xtr")
		simulator.AddUpdateRule(p_chemotaxis_update_rule);
//
		simulator.SetEndTime(2);
        simulator.Solve();

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }

    void TestPottsOnVesselWithDeformation()
	{
    	// For this mesh to be read in properly, one needs the patch in #2491 (until it gets committed of course)
		VtkMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/straight_vessel/config.vtu");
		MutableMesh<2,3> mutable_mesh;
		mutable_mesh.ConstructFromMeshReader(mesh_reader);
		double scaling = 1e-3;  // so distances are in m
		mutable_mesh.Scale(scaling,scaling,scaling);

		/*
		 * Setup Potts problem
		 */
		PottsArbitraryMeshFromMutableMeshGenerator<3> potts_generator(mutable_mesh);
		PottsArbitrarySurfaceIn3DMesh<3>* p_potts_mesh = potts_generator.GetMesh();
		TS_ASSERT_EQUALS(p_potts_mesh->GetNumNodes(), mutable_mesh.GetNumNodes());// XXX Why equate?

		unsigned max_num_cells = 10;  
		unsigned num_cells_seeded = 5;
		double target_volume = ComputeMeshTotalArea(mutable_mesh) / max_num_cells;
 
		for (unsigned i=0; i<num_cells_seeded; i++)  // XXX how does this give cell locations and lattice groups?? 
		{
			// Pick a node randomly avoiding boundary nodes
			unsigned num_node;
			do
			{
				num_node = RandomNumberGenerator::Instance()->randMod(p_potts_mesh->GetNumNodes());
			}
			while(p_potts_mesh->GetNode(num_node)->IsBoundaryNode());

			std::vector<Node<3>*> element_nodes;
			element_nodes.push_back(p_potts_mesh->GetNode(num_node));
			p_potts_mesh->AddElement(new PottsElement<3>(i, element_nodes));
			std::cout << "Seeded at "<< num_node << std::endl;
		}

		// TS_ASSERT_EQUALS(p_potts_mesh->GetNumElements(), num_cells_seeded);
        PRINT_2_VARIABLES(p_potts_mesh->GetNumElements(),num_cells_seeded);

        // Randomly place the cell markers for seeding
		std::vector<CellPtr> potts_cells;
		MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		CellsGenerator<FixedG1GenerationalCellCycleModel, 3> potts_cells_generator;
		potts_cells_generator.GenerateBasicRandom(potts_cells, p_potts_mesh->GetNumElements(), p_diff_type);

        TRACE("GenerateBasicRandom");
		// Create cell population linking potts mesh and cells
		PottsBasedCellPopulation<3> potts_cell_population(*p_potts_mesh, potts_cells);
		potts_cell_population.AddCellWriter<CellIdWriter>();
        cell_population.SetNumSweepsPerTimestep(1);
		// Set up Potts simulation
		OnLatticeSimulation<3> potts_simulator(potts_cell_population);
		potts_simulator.SetOutputDirectory("TestPottsOnVesselGrowth_Potts");
		potts_simulator.SetDt(0.1);
		potts_simulator.SetSamplingTimestepMultiple(10);

		/*
		 * Setup growth problem
		 */
        // Create representation of cells required for the mesh based simulation
        TRACE("NewCells");
        std::vector<CellPtr> vessel_cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> vessel_cells_generator;
        vessel_cells_generator.GenerateBasicRandom(vessel_cells, mutable_mesh.GetNumNodes(), p_diff_type);

        // Create a cell population linking mechanics mesh and cells
        MeshBasedCellPopulation<2,3> vessel_cell_population(mutable_mesh, vessel_cells);
        vessel_cell_population.SetWriteVtkAsPoints(true);
        vessel_cell_population.SetOutputMeshInVtk(true);
//        vessel_cell_population.SetOutputCellProliferativeTypes(true);
//        vessel_cell_population.SetOutputCellMutationStates(true);
        vessel_cell_population.CalculateRestLengths();
        //cell_population.SetDampingConstantNormal(0.1);//e5);
        //cell_population.SetDampingConstantMutant(0.1);//e5);


        // Set up cell-based simulation
        TRACE("NewSimulation");
        OffLatticeSimulation<2,3> vessel_simulator(vessel_cell_population);
        vessel_simulator.SetOutputDirectory("TestPottsOnVesselGrowth_Growth");
        vessel_simulator.SetDt(0.002);
        vessel_simulator.SetSamplingTimestepMultiple(500);
        vessel_simulator.SetUpdateCellPopulationRule(false); // No remeshing.

        // Add a force law for the springs
        //boost::shared_ptr<LinearSpringWithRestLengthDependentSpringConstantsForce<2,3> > p_linear_force(new LinearSpringWithRestLengthDependentSpringConstantsForce<2,3>());
        boost::shared_ptr<GeneralisedLinearSpringForce<2,3> > p_linear_force(new GeneralisedLinearSpringForce<2,3>());
        p_linear_force->SetMeinekeSpringStiffness(50.0);
        vessel_simulator.AddForce(p_linear_force);

        // Add a force law for normal stress, i.e. blood pressure only
        double pressure = 1.0666e4; // to match 80mmhg
        MAKE_PTR_ARGS(RadialForce, p_radial_force, (pressure));
        vessel_simulator.AddForce(p_radial_force);

        //Create a plane boundary to represent the inlet and pass them to the simulation
        c_vector<double,3> point;
        point(0)=-0.0288726;
        point(1)=-0.260683;
        point(2)=-0.313772;
        c_vector<double,3> normal;
        normal(0)=-0.232271;
        normal(1)= 0.95232;
        normal(2)=-0.197833;
        // boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_inlet_condition(new FixedRegionBoundaryCondition<2,3>(&vessel_cell_population,point,normal));
        // vessel_simulator.AddCellPopulationBoundaryCondition(p_inlet_condition);

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
		potts_cell_population.SetTemperature(1e-9);

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

        // MAKE_PTR(AspectRatioConstraintPottsUpdateRule<3>, p_aspect_ratio_update_rule);
        // p_aspect_ratio_update_rule->SetTargetAspectRatio(4);
        // p_aspect_ratio_update_rule->SetAspectRatioEnergyParameter(1e-7);
        // potts_simulator.AddUpdateRule(p_aspect_ratio_update_rule);

        potts_simulator.SetEndTime(potts_current_time + potts_timestep);
        potts_simulator.SetOutputDirectory("TestPottsOnVesselGrowth_Potts_Setup");
		potts_simulator.Solve();
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
		    /*
		     * Do one timestep of cell movement
		     */
		    std::cout << "Potts timestep" << std::endl;

            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(potts_current_time);
            if (potts_current_time > 0)
            {
                double dt = vessel_simulator.GetDt();
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
                double dt = vessel_simulator.GetDt();
                unsigned num_time_steps = (unsigned) (growth_timestep/dt+0.5);
                SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(growth_current_time + growth_timestep, num_time_steps);
            }

            vessel_simulator.SetEndTime(growth_current_time + growth_timestep);
            vessel_simulator.Solve();

            // Update location of nodes in potts mesh from recently computed deformation
            p_potts_mesh->UpdatePottsNodeLocationFromDelaunay();

            growth_current_time += growth_timestep;
		}
	}
};

#endif /*TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_*/
