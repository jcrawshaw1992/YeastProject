
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
#include "OffLatticeSimulation.hpp"

// Include the relevant cell writers
#include "CellAspectRatioWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellMajorAxisAngleWriter.hpp"

#include "PottsBasedCellPopulation.hpp"

#include "MeshBasedCellPopulation.hpp"

#include "ArbitraryAdhesionPottsUpdateRule.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
// New rule to make the cells wrap
#include "AppliedForceModifier.hpp"
#include "ArbitraryCurvatureConstraintPottsUpdateRule.hpp"

#include "PottsCellPropertiesModifier.hpp"

// #include "TractionDataLoader.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "PottsMeshGenerator.hpp"
// #include "WildTypeCellMutationState.hpp"
#include "AspectRatioConstraintPottsUpdateRule.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "LogFile.hpp"
#include "MechanotaxisPottsUpdateRule.hpp"
#include "PottsMeshFromMutableMeshGenerator.hpp"
#include "SmartPointers.hpp"
#include "Warnings.hpp"
// #include "VtkMeshReader.hpp"
#include "RandomNumberGenerator.hpp"
#include "AbstractForce.hpp"
#include "ArbitrarySurfaceAreaConstraintPottsUpdateRule.hpp"
#include "ArbitraryVolumeConstraintPottsUpdateRule.hpp"
#include "Debug.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "Honeycomb3DMeshGenerator.hpp"

#include "PottsArbitrarySurfaceIn3DMesh.hpp"
#include "PottsMeshFromMutableMeshGeneratorJess.hpp"
#include "UblasCustomFunctions.hpp"

// Needed for NodesOnlyMesh.
#include "PetscSetupAndFinalize.hpp"

#include <iostream>
#include <stdio.h>


class TestSetUpOnBirfication : public AbstractCellBasedTestSuite
{
public:

    void TestSetupPottsOnVessel()
    {
        // Birfucation
        unsigned N = 581;
        std::stringstream out;
        out << N;
        std::string mesh_size = out.str();
        std::string mesh_file = "projects/VascularRemodelling/test/data/bifurcation_cut/Scalled/config.vtu";
        std::string output_directory = "BifucationPottsTest/";

        // This data file is in mm??
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> mutable_mesh;
        mutable_mesh.ConstructFromMeshReader(mesh_reader);
        double scale = 1e-3; // so distances are in m
        mutable_mesh.Scale(1.0 * scale, 1.0 * scale, 1.0 * scale); // so distances are back in original scal

        double Radius = 5e-4 * scale;
        double Length = 30e-4 * scale; //2*96.6e-4 * scale; //12e-3;

        /*
		 * Setup Potts simulation
		 */

        PottsMeshFromMutableMeshGeneratorJess<3> potts_generator(mutable_mesh);
        PottsArbitrarySurfaceIn3DMesh<3>* p_potts_mesh = potts_generator.GetMesh();

        unsigned max_num_cells = 10;
        unsigned num_cells_seeded = 10;

        double CellScale =1e-3;
        double AthesticScalling = 5;
		double LongAxis = 10.0/2.0* 3.0*AthesticScalling *CellScale;//* 1e-6 * CellScale;//*scale ;
		double ShortAxis = 5.0/2.0*AthesticScalling*  CellScale;//*1e-6 * CellScale ;//*scale ;
		double CellArea = M_PI * LongAxis  * ShortAxis ;

		double AxisRatio = pow(LongAxis -ShortAxis,2)/ pow(LongAxis +ShortAxis,2);//  (LongAxis -ShortAxis)* (LongAxis -ShortAxis) /(LongAxis +ShortAxis)*(LongAxis +ShortAxis); 
		double CellPerimeter = M_PI * (LongAxis +ShortAxis)* (3 *AxisRatio * 1/(sqrt(-3*AxisRatio+4) +10) +1) ;   
        CellPerimeter/=15;//

        PRINT_2_VARIABLES(CellArea, CellPerimeter);

        for (unsigned i = 0; i < num_cells_seeded; i++)
        {
            // Pick a node randomly avoiding boundary nodes
            unsigned num_node;
            do
            {
                num_node = RandomNumberGenerator::Instance()->randMod(p_potts_mesh->GetNumNodes()); // Selects a random node

            } while (p_potts_mesh->GetNode(num_node)->IsBoundaryNode()); // XXX not sure what this does, Im not sure it does anything
            PRINT_VARIABLE(i);

            std::vector<Node<3>*> element_nodes;
            element_nodes.push_back(p_potts_mesh->GetNode(num_node));
            p_potts_mesh->AddElement(new PottsElement<3>(i, element_nodes));
            std::cout << "Seeded at " << num_node << std::endl;
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
        potts_population.SetNumSweepsPerTimestep(1);
        double Confluence_time =30;
        // Set up Potts simulation
        OnLatticeSimulation<3> potts_simulator(potts_population);


        // potts_simulator.SetDt(0.1);
        // potts_simulator.SetOutputDirectory("TestPottsOn2dSquare_Potts_Setup");
        // potts_simulator.SetSamplingTimestepMultiple(5);
        // potts_simulator.SetEndTime(Confluence_time);

        // // // Reduce temperature to avoid randomly removing cells before they grow (parameters will need proper fitting)
        // potts_population.SetTemperature(1e-9);

        /*
         * Simulate cell seeding from randomly placed cell markers
         */

        // MAKE_PTR(ArbitraryVolumeConstraintPottsUpdateRule<3>, p_volume_constraint_update_rule);
        // p_volume_constraint_update_rule->SetMatureCellTargetVolume(CellArea);
        // p_volume_constraint_update_rule->SetDeformationEnergyParameter(4e7);
        // potts_simulator.AddUpdateRule(p_volume_constraint_update_rule);

        // MAKE_PTR(ArbitraryAdhesionPottsUpdateRule<3>, p_adhesion_update_rule);
        // p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02 / 1000);
        // p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16 / 1000);
        // potts_simulator.AddUpdateRule(p_adhesion_update_rule);


        // MAKE_PTR(ArbitrarySurfaceAreaConstraintPottsUpdateRule<3>, p_area_constraint_update_rule);
        // p_area_constraint_update_rule->SetSurfaceAreaEnergyParameter(3e4);//(0.016e2);
        // p_area_constraint_update_rule->SetTargetSurfaceArea(CellPerimeter);
        // potts_simulator.AddUpdateRule(p_area_constraint_update_rule);

        /*
          Add mechanotaxis to Potts simulation for timesteps to follow
        */

        // std::string traction_file = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/Tester/results/Extracted/surface-tractions.xtr";
        // p_potts_mesh->TractionDataLoader(traction_file);
        // MAKE_PTR(MechanotaxisPottsUpdateRule<3>, p_mechanotaxis_update_rule);
        // p_mechanotaxis_update_rule->SetTractionCorrelationParameter(0.016e25); //1e-8);
        // // potts_simulator.AddUpdateRule(p_mechanotaxis_update_rule);




        // p_potts_mesh->UpdatePottsNodeLocationFromDelaunay();
        // p_potts_mesh->GetRadius(Radius);

        CellBasedSimulationArchiver<3, OnLatticeSimulation<3> >::Save(&potts_simulator);


        // void CellBasedSimulationArchiver<ELEMENT_DIM, SIM, SPACE_DIM>::Save(SIM* pSim)
        //  CellBasedSimulationArchiver<2, OnLatticeSimulation<2> >::Save(p_simulator1);

        // potts_simulator.Solve();

        
    }
};

#endif /*TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_*/





