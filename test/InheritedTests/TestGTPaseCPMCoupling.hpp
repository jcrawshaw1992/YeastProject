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

#ifndef TESTCPMFEMCOUPLING_HPP_
#define TESTCPMFEMCOUPLING_HPP_

#include <cxxtest/TestSuite.h>
#include "GTPasePDEMultipleMeshesModifier.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "UniformCellCycleModel.hpp"
#include "PottsMeshGenerator.hpp"
#include "OnLatticeSimulation.hpp"
#include "PottsBasedCellPopulation.hpp"

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

class TestGTPaseCPMCoupling : public AbstractCellBasedTestSuite
{
public:

    void TestMonolayer() throw(Exception)
    {
        /** The next line is needed because we cannot currently run Potts simulations in parallel. */
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
        element_nodes_cell_0.push_back(p_mesh->GetNode(2228));
        element_nodes_cell_0.push_back(p_mesh->GetNode(2229));
        element_nodes_cell_0.push_back(p_mesh->GetNode(2230));
        element_nodes_cell_0.push_back(p_mesh->GetNode(2231));
        element_nodes_cell_0.push_back(p_mesh->GetNode(2328));
        element_nodes_cell_0.push_back(p_mesh->GetNode(2329));
        element_nodes_cell_0.push_back(p_mesh->GetNode(2330));
        element_nodes_cell_0.push_back(p_mesh->GetNode(2331));
        p_mesh->AddElement(new PottsElement<2>(0, element_nodes_cell_0));
/*
        // Create cell 1
        std::vector<Node<2>*> element_nodes_cell_1;
           element_nodes_cell_1.push_back(p_mesh->GetNode(1872));
           element_nodes_cell_1.push_back(p_mesh->GetNode(1873));
           element_nodes_cell_1.push_back(p_mesh->GetNode(1874));
           element_nodes_cell_1.push_back(p_mesh->GetNode(1875));
           element_nodes_cell_1.push_back(p_mesh->GetNode(1773));
           element_nodes_cell_1.push_back(p_mesh->GetNode(1772));
           element_nodes_cell_1.push_back(p_mesh->GetNode(1774));
           element_nodes_cell_1.push_back(p_mesh->GetNode(1775));
           p_mesh->AddElement(new PottsElement<2>(1, element_nodes_cell_1));


        // Create cell 2
        std::vector<Node<2>*> element_nodes_cell_1;
        element_nodes_cell_1.push_back(p_mesh->GetNode(4538));
        element_nodes_cell_1.push_back(p_mesh->GetNode(4539));
        element_nodes_cell_1.push_back(p_mesh->GetNode(4541));
        element_nodes_cell_1.push_back(p_mesh->GetNode(4540));
        element_nodes_cell_1.push_back(p_mesh->GetNode(4440));
        element_nodes_cell_1.push_back(p_mesh->GetNode(4438));
        element_nodes_cell_1.push_back(p_mesh->GetNode(4441));
        element_nodes_cell_1.push_back(p_mesh->GetNode(4439));
        p_mesh->AddElement(new PottsElement<2>(1, element_nodes_cell_1));

        // Create cell 3
        std::vector<Node<2>*> element_nodes_cell_3;
        element_nodes_cell_3.push_back(p_mesh->GetNode(558));
        element_nodes_cell_3.push_back(p_mesh->GetNode(536));
        element_nodes_cell_3.push_back(p_mesh->GetNode(537));
        element_nodes_cell_3.push_back(p_mesh->GetNode(539));
        element_nodes_cell_3.push_back(p_mesh->GetNode(436));
        element_nodes_cell_3.push_back(p_mesh->GetNode(437));
        element_nodes_cell_3.push_back(p_mesh->GetNode(438));
        element_nodes_cell_3.push_back(p_mesh->GetNode(439));
        p_mesh->AddElement(new PottsElement<2>(3, element_nodes_cell_3));

        // Create cell 4
        std::vector<Node<2>*> element_nodes_cell_4;
        element_nodes_cell_4.push_back(p_mesh->GetNode(3010));
        element_nodes_cell_4.push_back(p_mesh->GetNode(3011));
        element_nodes_cell_4.push_back(p_mesh->GetNode(3012));
        element_nodes_cell_4.push_back(p_mesh->GetNode(3013));
        element_nodes_cell_4.push_back(p_mesh->GetNode(2910));
        element_nodes_cell_4.push_back(p_mesh->GetNode(2911));
        element_nodes_cell_4.push_back(p_mesh->GetNode(2913));
        element_nodes_cell_4.push_back(p_mesh->GetNode(2912));
        p_mesh->AddElement(new PottsElement<2>(4, element_nodes_cell_4));

        // Create cell 5
        std::vector<Node<2>*> element_nodes_cell_5;
        element_nodes_cell_5.push_back(p_mesh->GetNode(720));
        element_nodes_cell_5.push_back(p_mesh->GetNode(722));
        element_nodes_cell_5.push_back(p_mesh->GetNode(723));
        element_nodes_cell_5.push_back(p_mesh->GetNode(721));
        element_nodes_cell_5.push_back(p_mesh->GetNode(620));
        element_nodes_cell_5.push_back(p_mesh->GetNode(621));
        element_nodes_cell_5.push_back(p_mesh->GetNode(622));
        element_nodes_cell_5.push_back(p_mesh->GetNode(623));
        p_mesh->AddElement(new PottsElement<2>(5, element_nodes_cell_5));

        // Create cell 6
        std::vector<Node<2>*> element_nodes_cell_6;
        element_nodes_cell_6.push_back(p_mesh->GetNode(4027));
        element_nodes_cell_6.push_back(p_mesh->GetNode(4026));
        element_nodes_cell_6.push_back(p_mesh->GetNode(4028));
        element_nodes_cell_6.push_back(p_mesh->GetNode(4029));
        element_nodes_cell_6.push_back(p_mesh->GetNode(3926));
        element_nodes_cell_6.push_back(p_mesh->GetNode(3927));
        element_nodes_cell_6.push_back(p_mesh->GetNode(3929));
        element_nodes_cell_6.push_back(p_mesh->GetNode(3928));
        p_mesh->AddElement(new PottsElement<2>(6, element_nodes_cell_6));
*/
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom (cells, p_mesh->GetNumElements(), p_transit_type); //cells

        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetTemperature(0.008e9);

        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("GTPaseCPMMonolayer");
        // Do 40 time steps during the first second to let the cells reach target area. OK cause next three update rules don't specify temporal dynamics
        simulator.SetEndTime(1.0);
        simulator.SetDt(1.0/40);

        MAKE_PTR(ArbitraryVolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume (300e-12); // A
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(2e31); // lambda_a (maps delta_a=0.5um^2 to delta_H=5e6)
        simulator.AddUpdateRule(p_volume_constraint_update_rule);

        MAKE_PTR(ArbitrarySurfaceAreaConstraintPottsUpdateRule<2>, p_area_constraint_update_rule);
        // Ellipse a=16.9, b~=a/3~=5.65, A~=300, P~=75. Because of the regular grid discretisation we increase it to 120um
        p_area_constraint_update_rule->SetTargetSurfaceArea (120e-6); // P (is 150e-12 in paper but value was found wrong and corrected for 150e-6)
        p_area_constraint_update_rule->SetSurfaceAreaEnergyParameter(5e18); // lambda_p (maps delta_p=1um to delta_H=5e6)
        simulator.AddUpdateRule(p_area_constraint_update_rule);

        MAKE_PTR(ArbitraryAdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.75e6 / 200e-9); // Coupling energy per boundary site per unit length J_CM/Delta_x
		p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.75e6 / 200e-9); // Coupling energy per boundary site per unit length J_CM/Delta_x
        simulator.AddUpdateRule(p_adhesion_update_rule);

        std::cout << "LET THE CELLS REACH TARGET PERIMETER/AREA" << std::endl;
        simulator.Solve();

        MAKE_PTR(GTPasePDEMultipleMeshesModifier<2>, p_pde_modifier);
        simulator.AddSimulationModifier(p_pde_modifier);
        /// \todo MIGUEL
        // cell_population.AddUpdateLocationCallback(&p_pde_modifier->UpdateLocationCallback); // or (p_pde_modifier);

        MAKE_PTR(RhoBasedContractionUpdateRule<2>, p_rho_based_contraction_update_rule);
        p_rho_based_contraction_update_rule->SetMatureCellTargetRhoContractionForces(1.25); // P_th  Rho contraction threshold
        p_rho_based_contraction_update_rule->SetDeformationEnergyParameter(0.0025e9); // Effect of Rho on contraction
        simulator.AddUpdateRule(p_rho_based_contraction_update_rule);

        MAKE_PTR(BarbedEndsUpdateRule<2>, p_barbed_ends_update_rule);
        p_barbed_ends_update_rule->SetDeformationEnergyParameter(1.0);
        simulator.AddUpdateRule(p_barbed_ends_update_rule);

        simulator.SetEndTime(125.0);
        simulator.SetDt(1.0);
        std::cout << "RUN WITH GTPASE UPDATE RULES" << std::endl;
        simulator.Solve();
    }





    /*
     * EMPTYLINE
     *
     *To compile: scons ts=projects/VascularRemodelling/test/TestGTPaseCPMCoupling.hpp
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

#endif /* TESTCPMFEMCOUPLING_HPP_ */
