
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
// #include "CellMajorAxisAngleWriter.hpp"
// #include "PottsBasedCellPopulation.hpp"
#include "ArbitraryWrappedAdhesionPottsUpdateRule.hpp"
#include "MeshBasedCellPopulation.hpp"
// New rule to make the cells wrap
// #include "ArbitraryCurvatureConstraintPottsUpdateRule.hpp"

#include "PottsCellPropertiesModifier.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "LogFile.hpp"
#include "PottsMeshGenerator.hpp"
#include "SmartPointers.hpp"
#include "Warnings.hpp"

#include "AspectRatioConstraintPottsUpdateRule.hpp"
#include "PottsMeshFromMutableMeshGenerator.hpp"

#include "VtkMeshReader.hpp"

#include "Debug.hpp"
#include "RandomNumberGenerator.hpp"

#include "UblasCustomFunctions.hpp"

#include "ArbitraryWrappedAdhesionPottsUpdateRule.hpp"
#include "PottsMeshFromMutableMeshGeneratorJess.hpp"

#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "MutableMesh.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"

#include "Honeycomb3DMeshGenerator.hpp"

// Needed for NodesOnlyMesh.
#include "PetscSetupAndFinalize.hpp"

#include <iostream>
#include <stdio.h>

// Adapted update rules

#include "ArbitraryPerimeterOnSurfacePottsUpdateRule.hpp"
#include "ArbitraryVolumeOnSurfacePottsUpdateRule.hpp"

#include "WrappedPottsBasedCellPopulation.hpp"

// #include "PeriodicRectangleMeshGenerator.hpp"

#include "MathsFunctions.hpp"
#include "MechanotaxisPottsWrappedUpdateRule.hpp"

#include "CellAreaWriter.hpp"
#include "CellCenterWriter.hpp"
#include "CellPerimeterWriter.hpp"

class TestShearMigration : public AbstractCellBasedTestSuite
{
public:
  
     void TestNonperiodicPottsSweep()
    {
        MathsFunctions M;


        // Regular hexagonal lattie

        double N_D = 20;
        double N_Z = 20;
        double Width = 20;
        double Length = 20;

        Honeycomb3DMeshGenerator generator(N_D, N_Z, Width, Length);
        MutableMesh<2, 3>* mutable_mesh = generator.GetMesh();
        std::vector<unsigned> BoundaryVector = generator.GetBoundaryVector();

        /*
		 * Setup Potts simulation
		*/

        for (double k =-10; k<3; k++)
        {
          for (double j =-10; j<3; j++)
          {
             for(int i =0; i<100; i++)
             {
                    PottsMeshFromMutableMeshGeneratorJess<3> potts_generator(*mutable_mesh);
                    PottsArbitrarySurfaceIn3DMesh<3>* p_potts_mesh = potts_generator.GetMesh();
                    p_potts_mesh->SetBoundaries(BoundaryVector);

                    std::map<unsigned, unsigned> ElementPairing;

                    // Randomly select three nodes.

                    double num_cells_seeded = 1;
                    unsigned num_node;
                    // PRINT_VARIABLE(p_potts_mesh->GetNode(num_node)->IsBoundaryNode());
                    bool OnEdge = 1;
                    // double yPos = Length;
                    for (unsigned i = 0; i < num_cells_seeded; ++i) // i+=2)
                    {
                        OnEdge = 1;
                        while (OnEdge)
                        {
                            num_node = RandomNumberGenerator::Instance()->ranf() * (p_potts_mesh->GetNumNodes() - 60); // Selects a random node
                            // yPos = mutable_mesh->GetNode(num_node)->rGetLocation()[1];
                            //yPos < Length / 4 &&
                            if ( p_potts_mesh->GetNode(num_node)->IsBoundaryNode() == 0)
                            {
                                OnEdge = 0;
                            }
                        }
                        std::vector<Node<3>*> element_nodes1;
                        element_nodes1.push_back(p_potts_mesh->GetNode(num_node));
                        p_potts_mesh->AddElement(new PottsElement<3>(i, element_nodes1));
                        // std::cout << "Seeded at " << num_node << std::endl;
                    }

                    std::vector<CellPtr> potts_cells;
                    MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type); // Set the cell type for the seeds
                    CellsGenerator<FixedG1GenerationalCellCycleModel, 3> potts_cells_generator; // Set the cell cylce for the seeds
                    potts_cells_generator.GenerateBasicRandom(potts_cells, p_potts_mesh->GetNumElements(), p_diff_type); // Assigning the cell type and the cell cycle

                    WrappedPottsBasedCellPopulation<3> potts_population(*p_potts_mesh, potts_cells, BoundaryVector);

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
                    p_volume_constraint_update_rule->SetDeformationEnergyParameter(pow(10,i) );//0.1);
                    potts_simulator.AddUpdateRule(p_volume_constraint_update_rule);

                    double AxisRatio = pow(LongAxis - ShortAxis, 2) / pow(LongAxis + ShortAxis, 2); //
                    double CellPerimeter = M_PI * (LongAxis + ShortAxis) * (3 * AxisRatio * 1 / (sqrt(-3 * AxisRatio + 4) + 10) + 1); // M_PI*( LongAxis  + ShortAxis );//
               
                    MAKE_PTR(ArbitraryPerimeterOnSurfacePottsUpdateRule<3>, p_area_constraint_update_rule);
                    p_area_constraint_update_rule->SetSurfaceAreaEnergyParameter(pow(10,j) );//0.1);
                    p_area_constraint_update_rule->SetTargetSurfaceArea(CellPerimeter);
                    potts_simulator.AddUpdateRule(p_area_constraint_update_rule);

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

                    potts_simulator.SetSamplingTimestepMultiple(10);
                    potts_simulator.SetEndTime(11);

                    std::stringstream out;
                    out << "_" << i ;
                    std::string Iteration = out.str();

                    std::stringstream out2;
                    out2 << "Area_" << j << "Perimeter_" <<k;
                    std::string Parameters = out2.str();

                    potts_simulator.SetOutputDirectory("PottsMetrics/AreaAndPerimeterCourse/"+Parameters+"/Trial"+Iteration );
                    potts_simulator.Solve();
                    SimulationTime::Instance()->Destroy();
                    SimulationTime::Instance()->SetStartTime(0.0);
                }
            }
        }
    }

};

#endif /*TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_*/
