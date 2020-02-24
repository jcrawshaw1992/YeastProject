
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
#include "OnLatticeSimulation.hpp"

// Include the relevant cell writers
#include "CellIdWriter.hpp"
#include "PottsCellPropertiesModifier.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "LogFile.hpp"
#include "SmartPointers.hpp"
#include "Warnings.hpp"

#include "PottsMeshFromMutableMeshGenerator.hpp"

// #include "VtkMeshReader.hpp"

#include "ArbitraryWrappedAdhesionPottsUpdateRule.hpp"
#include "Debug.hpp"
#include "PottsMeshFromMutableMeshGeneratorJess.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"

#include "MutableMesh.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"

#include "PetscSetupAndFinalize.hpp"

#include <iostream>
#include <stdio.h>

// Adapted update rules

#include "ArbitraryPerimeterOnSurfacePottsUpdateRule.hpp"
#include "ArbitraryVolumeOnSurfacePottsUpdateRule.hpp"
#include "PeriodicRectangleMeshGenerator.hpp"
#include "WrappedPottsBasedCellPopulation.hpp"

// #include "Honeycomb3DMeshGenerator.hpp"
#include "MathsFunctions.hpp"

#include "CellAreaWriter.hpp"
#include "CellCenterWriter.hpp"
#include "CellPerimeterWriter.hpp"

 #include "MechanotaxisPottsWrappedUpdateRule.hpp"

class TestShearMigration : public AbstractCellBasedTestSuite
{
public:


 void TestSingleParameterPair()
    {
        MathsFunctions M;

        // Regular hexagonal lattie
        double N_D = 30;
        double N_Z = 90;
        double Width = 15;
        double Length = 45;

        PeriodicRectangleMeshGenerator generator(N_D, N_Z, Width, Length);
        MutableMesh<2, 3>* mutable_mesh = generator.GetMesh();
        std::vector<unsigned> BoundaryVector = generator.GetBoundaryVector();

        /*
		 * Setup Potts simulation
		*/
        for (double j = -1; j > -4; j--)
        {
            for (int i = 100; i < 500; i++)
            {
                    PRINT_2_VARIABLES(i, j)       
                    PottsMeshFromMutableMeshGeneratorJess<3> potts_generator(*mutable_mesh);
                    PottsArbitrarySurfaceIn3DMesh<3>* p_potts_mesh = potts_generator.GetMesh();
                    p_potts_mesh->SetBoundaries(BoundaryVector);
                    p_potts_mesh->SetMeshSize( N_D,  N_Z,  Width,  Length);
                    unsigned num_node;
                    double num_cells_seeded = 1;
                    std::map<unsigned, unsigned> ElementPairing;
                    for (unsigned i = 0; i < num_cells_seeded * 2; i += 2)
                    {
                        bool OnEdge = 1;
                        double yPos = Length;
                        
                        while (OnEdge || yPos > 15|| yPos < 8)
                        {
                            num_node = RandomNumberGenerator::Instance()->ranf() * (p_potts_mesh->GetNumNodes()); // Selects a random node

                            yPos = mutable_mesh->GetNode(num_node)->rGetLocation()[1];
                            // PRINT_VARIABLE(yPos)
                            if (p_potts_mesh->GetNode(num_node)->IsBoundaryNode() == 0)
                            {
                                OnEdge = 0;
                            }
                        }
                    
                        std::vector<Node<3>*> element_nodes;
                        std::vector<Node<3>*> element_nodes2;
                        element_nodes.push_back(p_potts_mesh->GetNode(num_node));
                        element_nodes2.push_back(p_potts_mesh->GetNode(num_node + 1));

                        p_potts_mesh->AddElement(new PottsElement<3>(i, element_nodes));
                        p_potts_mesh->AddElement(new PottsElement<3>(i + 1, element_nodes2));

                        // ElementPairing.push_back(std::make_pair(i , i+1) );
                        ElementPairing[i] = i + 1;
                        ElementPairing[i + 1] = i;
                    }

                    // Randomly select three nodes.

                    std::vector<CellPtr> potts_cells;
                    MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type); // Set the cell type for the seeds
                    CellsGenerator<FixedG1GenerationalCellCycleModel, 3> potts_cells_generator; // Set the cell cylce for the seeds
                    potts_cells_generator.GenerateBasicRandom(potts_cells, p_potts_mesh->GetNumElements(), p_diff_type); // Assigning the cell type and the cell cycle
                    WrappedPottsBasedCellPopulation<3> potts_population(*p_potts_mesh, potts_cells, ElementPairing, BoundaryVector);
            
                    // // Add in all the cell writers
                    potts_population.AddCellWriter<CellAreaWriter>();
                    potts_population.AddCellWriter<CellIdWriter>();
                    potts_population.AddCellWriter<CellCenterWriter>();
                    potts_population.AddCellWriter<CellPerimeterWriter>();

                    potts_population.SetNumSweepsPerTimestep(1);
                    OnLatticeSimulation<3> potts_simulator(potts_population);

                    double LongAxis = 10 / 2;
                    double ShortAxis = 5.0 / 2.0;
                    double CellArea = M_PI * LongAxis * ShortAxis;

                    MAKE_PTR(ArbitraryVolumeOnSurfacePottsUpdateRule<3>, p_volume_constraint_update_rule);
                    p_volume_constraint_update_rule->SetMatureCellTargetVolume(CellArea);
                    p_volume_constraint_update_rule->SetDeformationEnergyParameter(1); //0.1);
                    potts_simulator.AddUpdateRule(p_volume_constraint_update_rule);

                    double AxisRatio = pow(LongAxis - ShortAxis, 2) / pow(LongAxis + ShortAxis, 2); //
                    double CellPerimeter = M_PI * (LongAxis + ShortAxis) * (3 * AxisRatio * 1 / (sqrt(-3 * AxisRatio + 4) + 10) + 1);

                    MAKE_PTR(ArbitraryPerimeterOnSurfacePottsUpdateRule<3>, p_area_constraint_update_rule);
                    p_area_constraint_update_rule->SetSurfaceAreaEnergyParameter(0.1); //0.1);
                    p_area_constraint_update_rule->SetTargetSurfaceArea(CellPerimeter);

                    potts_simulator.AddUpdateRule(p_area_constraint_update_rule);

                    double WallShearStress = -6;//
                    p_potts_mesh->SetConstantWallShearStress(WallShearStress);

                    MAKE_PTR(ArbitraryWrappedAdhesionPottsUpdateRule<3>, p_adhesion_update_rule);
                    p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02 / 1000);
                    p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16 / 1000);
                    potts_simulator.AddUpdateRule(p_adhesion_update_rule);

                    // Set up Potts simulation

                    p_potts_mesh->UpdatePottsNodeLocationFromDelaunay();
                    potts_population.SetTemperature(0.1);

                    MAKE_PTR(PottsCellPropertiesModifier<3>, p_modifier);
                    p_modifier->SetMeshDimensions(N_D, N_Z, Width, Length);
                    potts_simulator.AddSimulationModifier(p_modifier);

                    potts_simulator.SetSamplingTimestepMultiple(100);
                    potts_simulator.SetEndTime(10);

                    std::stringstream out;
                    out << "_" << i;
                    std::string Iteration = out.str();

                    std::stringstream out2;
                    out2 << "SSParamter_" << j;
                    std::string Parameters = out2.str();
                    // potts_simulator.SetOutputDirectory("PottsMetrics/Periodic/FindingBug/Trial" + Iteration);
                   
                    potts_simulator.SetOutputDirectory("PottsMetrics/Periodic/ShearStressTest/" + Parameters + "/Trial" + Iteration);
                    // potts_simulator.SetOutputDirectory("PottsMetrics/Periodic/");
                    potts_simulator.Solve();


                    /*
                    *  Add mechanotaxis to Potts simulation for timesteps to follow
                    */

                    MAKE_PTR(MechanotaxisPottsWrappedUpdateRule<3>, p_mechanotaxis_update_rule);
                    p_mechanotaxis_update_rule->SetTractionCorrelationParameter(pow(10,j));
                    potts_simulator.AddUpdateRule(p_mechanotaxis_update_rule);
                    potts_simulator.SetSamplingTimestepMultiple(10);
                    potts_simulator.SetEndTime(30);
                    potts_simulator.Solve();
                    


                    SimulationTime::Instance()->Destroy();
                    SimulationTime::Instance()->SetStartTime(0.0);
            }
        }
    }

};

#endif /*TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_*/
