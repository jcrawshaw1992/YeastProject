/*

Copyright (c) 2005-2015, University of Oxford.
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

#ifndef TESTCPMADHESIVECOUPLING_HPP_
#define TESTCPMADHESIVECOUPLING_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "GTPasePDEMultipleMeshesModifier.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "SmartPointers.hpp"
#include "PottsMeshGenerator.hpp"
#include "OnLatticeSimulation.hpp"
#include "OffLatticeSimulation.hpp"

#include "VolumeConstraintPottsUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"

#include "TransitCellProliferativeType.hpp"
#include "CellLabelWriter.hpp"
#include "RhoBasedContractionUpdateRule.hpp"
#include "BarbedEndsUpdateRule.hpp"
#include "PottsArbitraryMeshFromMutableMeshGenerator.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "ArbitraryVolumeConstraintPottsUpdateRule.hpp"
#include "ArbitrarySurfaceAreaConstraintPottsUpdateRule.hpp"
#include "ArbitraryAdhesionPottsUpdateRule.hpp"

#include "VolumeTrackingModifier.hpp"
#include "ContactInhibitionCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "NagaiHondaForce.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "Warnings.hpp"
//#include "CellOrientationWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "PottsBasedCellPopulation.hpp"

#include "MeshBasedCellPopulation.hpp"


#include "CellMutationStatesCountWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellLabelWriter.hpp"
#include "PottsArbitraryMeshFromMutableMeshGenerator.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "ArbitraryVolumeConstraintPottsUpdateRule.hpp"
#include "ArbitrarySurfaceAreaConstraintPottsUpdateRule.hpp"
#include "ArbitraryAdhesionPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "LogFile.hpp"
#include "PottsMeshGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "RandomNumberGenerator.hpp"
#include "AbstractForce.hpp"
#include "CellIdWriter.hpp"
#include <math.h>
#include <iostream>
#include <fstream>
#include "AspectRatioConstraintPottsUpdateRule.hpp"

#include "RacDependentAdhesionPottsUpdateRule.hpp"
#include "BoundaryNodeWriter.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestGTPaseAdhesivityCoupling : public AbstractCellBasedTestSuite
{
public:


    void TestRacDependentUpdateRule() throw (Exception)
    {
     EXIT_IF_PARALLEL;

     // Generate a monolayer of 50x50um
     MutableMesh<2,2> original_delaunay_mesh;
     original_delaunay_mesh.ConstructRectangularMesh(50,50);
     original_delaunay_mesh.Scale(1e-6, 1e-6);

     // Writes out the generated mesh. Useful for determining the element indices used for cell creation
     VtkMeshWriter<2,2> mesh_writer("GTPaseCPMMonolayer_generated_mesh", "generated_mesh", false);
     mesh_writer.WriteFilesUsingMesh(original_delaunay_mesh);

     // Generate a potts mesh from an arbitrary mesh of triangles
     PottsArbitraryMeshFromMutableMeshGenerator<2> potts_generator(original_delaunay_mesh);
     PottsArbitrarySurfaceIn3DMesh<2>* p_mesh = potts_generator.GetMesh();

     // Create a cell from a number of triangles in the original mesh
     std::vector<Node<2>*> element_nodes_cell_0;
     element_nodes_cell_0.push_back(p_mesh->GetNode(3418));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3419));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3420));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3421));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3422));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3423));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3424));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3425));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3426));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3427));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3428));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3429));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3430));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3431));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3432));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3519));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3520));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3521));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3522));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3523));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3524));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3525));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3526));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3527));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3528));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3529));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3530));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3531));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3532));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3619));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3620));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3621));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3622));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3623));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3624));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3625));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3626));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3627));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3628));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3629));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3630));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3631));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3632));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3719));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3720));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3721));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3722));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3723));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3724));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3725));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3726));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3727));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3728));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3729));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3731));
     element_nodes_cell_0.push_back(p_mesh->GetNode(3732));
     p_mesh->AddElement(new PottsElement<2>(0, element_nodes_cell_0));

     // Create cell 1
     std::vector<Node<2>*> element_nodes_cell_1;
      element_nodes_cell_1.push_back(p_mesh->GetNode(2055));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2056));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2057));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2058));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2059));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2060));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2061));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2062));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2063));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2064));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2065));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2066));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2067));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2068));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2069));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2155));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2156));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2157));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2158));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2159));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2160));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2161));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2162));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2163));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2164));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2165));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2166));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2167));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2168));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2169));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2255));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2256));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2257));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2258));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2259));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2260));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2261));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2262));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2263));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2264));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2265));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2266));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2267));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2268));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2269));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2355));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2356));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2357));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2358));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2359));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2360));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2361));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2362));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2363));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2364));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2365));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2366));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2367));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2368));
      element_nodes_cell_1.push_back(p_mesh->GetNode(2369));
     p_mesh->AddElement(new PottsElement<2>(1,element_nodes_cell_1));

     // Create cells
     std::vector<CellPtr> cells;
     MAKE_PTR(DifferentiatedCellProliferativeType, p_transit_type);
     CellsGenerator<UniformCellCycleModel, 2> cells_generator;
     cells_generator.GenerateBasicRandom (cells, p_mesh->GetNumElements(), p_transit_type);

     // Make this pointer first as if we move it after creating the cell population the label numbers aren't tracked
    MAKE_PTR(CellLabel, p_label);

    // Create cell population
    PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
    cell_population.SetTemperature (0.008e9); //(0.008e9);

    //cell_population.AddCellWriter<CellIdWriter>();
    //cell_population.AddCellWriter<CellOrientationWriter>();

     // Set up cell-based simulation
     OnLatticeSimulation<2> simulator(cell_population);
     simulator.SetOutputDirectory("GTPaseAdhesiveCoupling");
     simulator.SetDt(0.1);
     simulator.SetEndTime(3.0);

     MAKE_PTR(ArbitraryVolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
     p_volume_constraint_update_rule->SetMatureCellTargetVolume (300e-12); //(300e-12); // A
     p_volume_constraint_update_rule->SetDeformationEnergyParameter (4e30); // (4e30); // lambda_a (is 4e18 in paper but value was found wrong and corrected for 4e30)
     simulator.AddUpdateRule(p_volume_constraint_update_rule);

     MAKE_PTR(ArbitrarySurfaceAreaConstraintPottsUpdateRule<2>, p_area_constraint_update_rule);
     p_area_constraint_update_rule->SetTargetSurfaceArea (150e-6); // (150e-6); // P (is 150e-12 in paper but value was found wrong and corrected for 150e-6)
     p_area_constraint_update_rule->SetSurfaceAreaEnergyParameter(0.016e18); //(0.016e18); // lambda_p
     simulator.AddUpdateRule(p_area_constraint_update_rule);

     MAKE_PTR(ArbitraryAdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
     p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.75e6 / 200e-9); // Coupling energy per boundary site per unit length J_CM/Delta_x
	 p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.5e6 / 200e-9); // Coupling energy per boundary site per unit length J_CM/Delta_x
     simulator.AddUpdateRule(p_adhesion_update_rule);

     std::cout << "LET THE CELLS REACH TARGET PERIMETER/AREA" << std::endl;
     simulator.Solve();

     simulator.RemoveAllUpdateRules();
     simulator.AddUpdateRule(p_volume_constraint_update_rule);
     simulator.AddUpdateRule(p_area_constraint_update_rule);

     MAKE_PTR(GTPasePDEMultipleMeshesModifier<2>, p_pde_modifier);
        simulator.AddSimulationModifier(p_pde_modifier);

     /// CAMBIAR NOMBRES A LOS PARAMETERS
     MAKE_PTR(RacDependentAdhesionPottsUpdateRule<2>, p_rac_based_adhesion_update_rule);
     p_rac_based_adhesion_update_rule->SetLowLowRacAdhesionEnergyParameter(0.5*0.75e6 / 200e-9); //LowLow
     p_rac_based_adhesion_update_rule->SetHighLowRacAdhesionEnergyParameter(0.1*0.75e6 / 200e-9); //Low/High
     p_rac_based_adhesion_update_rule->SetHighHighRacAdhesionEnergyParameter(0.01*0.75e6 / 200e-9); //HighHigh
     p_rac_based_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.75e6*0.75e6 / 200e-9); //CAMBIAR
     simulator.AddUpdateRule(p_rac_based_adhesion_update_rule);

     simulator.SetEndTime(20.0);
     std::cout << "RUN WITH GTPASE UPDATE RULES" << std::endl;
     simulator.Solve();
    }

    /*
     * EMPTYLINE
     *
     *To compile: scons b=GccOpt_ndebug_4 ts=projects/VascularRemodelling/test/TestGTPaseAdhesivityCoupling.hpp
     *
     *
     * To visualize the results, we need to use Paraview.
     * See UserTutorials/VisualizingWithParaview for more information.
     *
     * Load the file {{{/tmp/$USER/testoutput/PottsCellSorting3D/results_from_time_0/results.pvd}}}, and click apply
     *
     * Add box "Glyphs" to represent lattice sites. You will need to adjust the size so they don't overlap.
     *
     * Select the "Display" tab and select "color by" cell label. (you can also "color by" cell index to see individual cells)
     *
     * Add a "Threshold" filter, filter by cell type and make the lower threshold 0 or greater (unocupied lattice sites are labelled with -1). This will allow you to view only the cells.
     *
     * Click play to see the evolution of the simulation.
     *
     * You should see that the cells sort into ones of the same type.
     *
     *Visualize the entire mesh tmp/mobernabeu/testoutput/GTPaseCPMMonolayer_generated_mesh  generated_mesh.vtu
     *
     * EMPTYLINE
     */
};

#endif /* TESTCPMADHESIVECOUPLING_HPP_ */
