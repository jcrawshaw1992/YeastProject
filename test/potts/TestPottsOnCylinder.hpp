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
#include "UblasCustomFunctions.hpp"

#include <cmath>
#include <cstdio>
#include <ctime>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellsGenerator.hpp"
#include "OnLatticeSimulation.hpp"
#include "OffLatticeSimulation.hpp"
// #include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "ArbitraryVolumeConstraintPottsUpdateRule.hpp"
#include "ArbitraryAdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PottsMeshGenerator.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "Warnings.hpp"
#include "LogFile.hpp"
#include "SmartPointers.hpp"
#include "ChemotaxisPottsUpdateRule.hpp"
#include "AspectRatioConstraintPottsUpdateRule.hpp"
#include "PottsMeshFromMutableMeshGenerator.hpp"
#include "VtkMeshReader.hpp"
#include "RandomNumberGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "AbstractForce.hpp"
#include "Debug.hpp"

#include "PottsArbitraryMeshFromMutableMeshGenerator.hpp"

#include "FixedG1GenerationalCellCycleModel.hpp"
// #include "PottsArbitrarySurfaceIn3DMesh.hpp"
// #include "PottsArbitrarySurfaceIn3DMesh.hpp"

#include "PottsArbitraryMeshFromMutableMeshGenerator.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"

#include "VtkMeshReader.hpp"
#include "MutableMesh.hpp"


// Needed for NodesOnlyMesh.
#include "PetscSetupAndFinalize.hpp"



class TestPottsOnSurface : public AbstractCellBasedTestSuite
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




    void TestPottsOn2DIn3DDisk() throw(Exception)
     {
    // 	TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/disk_in_3d");
    // 	MutableMesh<2,3> original_delaunay_mesh;
    // 	original_delaunay_mesh.ConstructFromMeshReader(mesh_reader);
    VtkMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/SimpleTestCylinder/config.vtu");
    	// VtkMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/bifurcation_cut/config.vtu");
		MutableMesh<2,3> mutable_mesh;
		mutable_mesh.ConstructFromMeshReader(mesh_reader);


    PottsArbitraryMeshFromMutableMeshGenerator<3> potts_generator(mutable_mesh);
    PottsArbitrarySurfaceIn3DMesh<3>* p_potts_mesh = potts_generator.GetMesh();

        // Make two triangular elements out of these nodes
        // std::vector<Node<3>*> element_nodes;
        // element_nodes.push_back(p_potts_mesh->GetNode(100));

        // p_potts_mesh->AddElement(new PottsElement<3>(0, element_nodes));



        unsigned num_cells_seeded = 2; //50;

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
      

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_potts_mesh->GetNumElements(), p_diff_type);

        // Create cell population
        PottsBasedCellPopulation<3> cell_population(*p_potts_mesh, cells);
        //cell_population.SetTemperature(1e-9);

        // Set up cell-based simulation
        OnLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("PottsOn3dDisc");
        simulator.SetDt(0.1);
        simulator.SetEndTime(1);


        // Create update rules and pass to the simulation
        MAKE_PTR(ArbitraryVolumeConstraintPottsUpdateRule<3>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(2.5);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(400);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        // simulator.AddUpdateRule(p_volume_constraint_update_rule);


      //   MAKE_PTR(ArbitraryAdhesionPottsUpdateRule<3>, p_adhesion_update_rule);
      //   p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02);
      //   p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16 * 10);
      //   // simulator.AddPottsUpdateRule(p_adhesion_update_rule);
      //   simulator.AddUpdateRule(p_adhesion_update_rule);

      //  MAKE_PTR(AspectRatioConstraintPottsUpdateRule<3>, p_aspect_ratio_update_rule);
      //  p_aspect_ratio_update_rule->SetTargetAspectRatio(1);
      //  p_aspect_ratio_update_rule->SetAspectRatioEnergyParameter(1);
      // //  simulator.AddPottsUpdateRule(p_aspect_ratio_update_rule);
      //   simulator.AddUpdateRule(p_aspect_ratio_update_rule);

      //  MAKE_PTR(ChemotaxisPottsUpdateRule<2>, p_chemotaxis_update_rule);
      //  simulator.AddPottsUpdateRule(p_chemotaxis_update_rule);

        // Run simulation
        simulator.Solve();

        // // original_delaunay_mesh.Scale(1.5,3.0,1.0);
        // p_potts_mesh->UpdatePottsNodeLocationFromDelaunay();

        // simulator.SetEndTime(2);
        // simulator.Solve();

    }
};

#endif /*TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_*/
