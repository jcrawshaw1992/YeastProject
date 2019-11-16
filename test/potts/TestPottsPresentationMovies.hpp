
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
 #include "MechanotaxisPottsWrappedUpdateRule.hpp"
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
#include "ArbitraryWrappedAdhesionPottsUpdateRule.hpp"

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

#include "MathsFunctions.hpp"
#include "MechanotaxisPottsWrappedUpdateRule.hpp"

class TestShearMigration : public AbstractCellBasedTestSuite
{
public:
     void TestWrappedPottsMigration()
    {
         MathsFunctions M;;

          // Regular square --hexagonal lattie
            double scale = 1e-3; 
             double N_D = 35;
             double N_Z = 50;
			 double Width = 35; //5e-3
			 double Length = 50; //30e-3

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

      // Randomly select three nodes. 

        double num_cells_seeded = 3;
        unsigned num_node;
        std::vector<unsigned> NodesAllocated;
        // PRINT_VARIABLE(p_potts_mesh->GetNode(num_node)->IsBoundaryNode());
            bool OnEdge = 1;
            
            for (unsigned i=0; i<num_cells_seeded*2; i+=2) 
            {
                OnEdge = 1;
                double yPos =Length;
                while(OnEdge)
                {
                    num_node = RandomNumberGenerator::Instance()->ranf()*(p_potts_mesh->GetNumNodes()-60);  // Selects a random node 
                    yPos = mutable_mesh->GetNode(num_node)->rGetLocation()[1];
                
                
                    if (yPos <Length/4 && p_potts_mesh->GetNode(num_node)->IsBoundaryNode()==0 && p_potts_mesh->GetNode(num_node+1)->IsBoundaryNode() ==0 && M.IsNumberInVector(NodesAllocated, num_node)==0 && M.IsNumberInVector(NodesAllocated, num_node+1)==0)
                    {
                        OnEdge =0;
                    }
                }
                    
                NodesAllocated.push_back(num_node); NodesAllocated.push_back(num_node+1);
            
            	std::vector<Node<3>*> element_nodes1; std::vector<Node<3>*> element_nodes2;
            
                element_nodes1.push_back(p_potts_mesh->GetNode(num_node));
            	element_nodes2.push_back(p_potts_mesh->GetNode(num_node+1));
                
            	p_potts_mesh->AddElement(new PottsElement<3>(i, element_nodes1));
            	p_potts_mesh->AddElement(new PottsElement<3>(i+1, element_nodes2));

            	ElementPairing[i] = i+1;
            	ElementPairing[i+1] = i;
            	std::cout << "Seeded at "<< num_node << std::endl;
            }
      
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
        OnLatticeSimulation<3> potts_simulator(potts_population);

        double LongAxis = 10/2 ;//* 1e-6 * CellScale ;//*scale ;
		double ShortAxis = 5.0/2.0 ;//*1e-6 * CellScale ;//*scale ;
		double CellArea = M_PI * LongAxis  * ShortAxis ;

        MAKE_PTR(ArbitraryVolumeOnSurfacePottsUpdateRule<3>, p_volume_constraint_update_rule);
		p_volume_constraint_update_rule->SetMatureCellTargetVolume(CellArea);
		p_volume_constraint_update_rule->SetDeformationEnergyParameter(4000);//(4000); //r(4000000000);
		potts_simulator.AddUpdateRule(p_volume_constraint_update_rule);

        double AxisRatio = pow(LongAxis -ShortAxis,2)/ pow(LongAxis +ShortAxis,2);//  
        double CellPerimeter = M_PI * (LongAxis +ShortAxis)* (3 *AxisRatio * 1/(sqrt(-3*AxisRatio+4) +10) +1) ;  
 

        MAKE_PTR(ArbitraryPerimeterOnSurfacePottsUpdateRule<3>, p_area_constraint_update_rule);
        p_area_constraint_update_rule->SetSurfaceAreaEnergyParameter(3e4);
		p_area_constraint_update_rule->SetTargetSurfaceArea(CellPerimeter); 
		potts_simulator.AddUpdateRule(p_area_constraint_update_rule);

        // MAKE_PTR(ArbitraryWrappedAdhesionPottsUpdateRule<3>, p_adhesion_update_rule);
		// p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02 / 1000);
		// p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16 / 1000);
		// potts_simulator.AddUpdateRule(p_adhesion_update_rule);

        // Set up Potts simulation
        

        p_potts_mesh->UpdatePottsNodeLocationFromDelaunay();
        potts_population.SetTemperature(1e-9);
        double WallShearStress = -6;//
		p_potts_mesh->SetConstantWallShearStress(WallShearStress);

        /*
        *  Add mechanotaxis to Potts simulation for timesteps to follow
         */

         MAKE_PTR(MechanotaxisPottsWrappedUpdateRule<3>, p_mechanotaxis_update_rule);
		 p_mechanotaxis_update_rule->SetTractionCorrelationParameter(10e3);// that is stupidly big ---  Works for ---  update rule 
        potts_simulator.AddUpdateRule(p_mechanotaxis_update_rule);

        MAKE_PTR(PottsCellPropertiesModifier<3>, p_modifier);
        p_modifier->SetMeshDimensions(N_D, N_Z ,  Width, Length);
        potts_simulator.AddSimulationModifier(p_modifier);


        /////////////////////////////

        potts_simulator.SetOutputDirectory("PottsPresentationMovies/ShearMigration");
        potts_simulator.SetSamplingTimestepMultiple(1);
        potts_simulator.SetEndTime(10);

        potts_simulator.Solve();
        TRACE("Simulation Complete")
    }



    void NOTestPottsOnCylinderWithDeformation()
    {
         MathsFunctions M;
        double N_D = 10;
        double N_Z = 10;
        double Width = 10;

        double Length = 10;
        double NumNodes = N_D*N_Z;

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

      // Randomly select three nodes. 

        double num_cells_seeded = 3;
        unsigned num_node;
        std::vector<unsigned> NodesAllocated;
        // PRINT_VARIABLE(p_potts_mesh->GetNode(num_node)->IsBoundaryNode());
            // bool OnEdge = 1;
            // for (unsigned i=0; i<num_cells_seeded*2; i+=2) 
            // {
            //     OnEdge = 1;
            //     while(OnEdge)
            //     {
            //         num_node = RandomNumberGenerator::Instance()->ranf()*(p_potts_mesh->GetNumNodes()-60);  // Selects a random node 
                
                
            //         if (p_potts_mesh->GetNode(num_node)->IsBoundaryNode()==0 && p_potts_mesh->GetNode(num_node+1)->IsBoundaryNode() ==0 && M.IsNumberInVector(NodesAllocated, num_node)==0 && M.IsNumberInVector(NodesAllocated, num_node+1)==0)
            //         {
            //             OnEdge =0;
            //         }
            //     }
                    
            //     NodesAllocated.push_back(num_node); NodesAllocated.push_back(num_node+1);
            
            // 	std::vector<Node<3>*> element_nodes1; std::vector<Node<3>*> element_nodes2;
            
            //     element_nodes1.push_back(p_potts_mesh->GetNode(num_node));
            // 	element_nodes2.push_back(p_potts_mesh->GetNode(num_node+1));
                
            // 	p_potts_mesh->AddElement(new PottsElement<3>(i, element_nodes1));
            // 	p_potts_mesh->AddElement(new PottsElement<3>(i+1, element_nodes2));

            // 	ElementPairing[i] = i+1;
            // 	ElementPairing[i+1] = i;
            // 	std::cout << "Seeded at "<< num_node << std::endl;
            // }
      
        // // Randomly place the cell markers for seeding
        unsigned i;

        i =0;
        num_node = 24;//63-18;  // Selects a random node 
        std::vector<Node<3>*> element_nodes1; std::vector<Node<3>*> element_nodes2;

        element_nodes1.push_back(p_potts_mesh->GetNode(num_node));
        element_nodes2.push_back(p_potts_mesh->GetNode(num_node+1));

        p_potts_mesh->AddElement(new PottsElement<3>(i, element_nodes1));
        p_potts_mesh->AddElement(new PottsElement<3>(i+1, element_nodes2));

        ElementPairing[i] = i+1;
        ElementPairing[i+1] = i;

        i =2;
        num_node = 47;  // Selects a random node 
        std::vector<Node<3>*> Aelement_nodes1; std::vector<Node<3>*> Aelement_nodes2;

        Aelement_nodes1.push_back(p_potts_mesh->GetNode(num_node));
        Aelement_nodes2.push_back(p_potts_mesh->GetNode(num_node+1));

        p_potts_mesh->AddElement(new PottsElement<3>(i, Aelement_nodes1));
        p_potts_mesh->AddElement(new PottsElement<3>(i+1, Aelement_nodes2));

        ElementPairing[i] = i+1;
        ElementPairing[i+1] = i;

         i =4;
        num_node = 85;//201+45-2*18;  // Selects a random node 
        std::vector<Node<3>*> Belement_nodes1; std::vector<Node<3>*> Belement_nodes2;

        Belement_nodes1.push_back(p_potts_mesh->GetNode(num_node));
        Belement_nodes2.push_back(p_potts_mesh->GetNode(num_node+1));

        p_potts_mesh->AddElement(new PottsElement<3>(i, Belement_nodes1));
        p_potts_mesh->AddElement(new PottsElement<3>(i+1, Belement_nodes2));

        ElementPairing[i] = i+1;
        ElementPairing[i+1] = i;
		


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
        OnLatticeSimulation<3> potts_simulator(potts_population);

        MAKE_PTR(ArbitraryVolumeOnSurfacePottsUpdateRule<3>, p_volume_constraint_update_rule);
		p_volume_constraint_update_rule->SetMatureCellTargetVolume(15);
		p_volume_constraint_update_rule->SetDeformationEnergyParameter(400);//(4000); //r(4000000000);
		potts_simulator.AddUpdateRule(p_volume_constraint_update_rule);


        MAKE_PTR(ArbitraryPerimeterOnSurfacePottsUpdateRule<3>, p_area_constraint_update_rule);
        p_area_constraint_update_rule->SetSurfaceAreaEnergyParameter(0.03e4);
		p_area_constraint_update_rule->SetTargetSurfaceArea(10); 
		potts_simulator.AddUpdateRule(p_area_constraint_update_rule);

        MAKE_PTR(ArbitraryWrappedAdhesionPottsUpdateRule<3>, p_adhesion_update_rule);
		p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02 / 1000);
		p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16 / 1000);
		potts_simulator.AddUpdateRule(p_adhesion_update_rule);

        // Set up Potts simulation
        

        p_potts_mesh->UpdatePottsNodeLocationFromDelaunay();
    
        // MAKE_PTR(PottsCellPropertiesModifier<3>, p_modifier);
        // potts_simulator.AddSimulationModifier(p_modifier);


        /////////////////////////////

        potts_simulator.SetOutputDirectory("PottsPresentationMovies/SimulationWithUpdateRules");
        potts_simulator.SetSamplingTimestepMultiple(1);
        potts_simulator.SetEndTime(10);

        potts_simulator.Solve();
        TRACE("Simulation Complete")
    }
};

#endif /*TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_*/
