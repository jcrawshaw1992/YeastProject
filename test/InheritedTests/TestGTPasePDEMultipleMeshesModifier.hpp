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

#ifndef TESTGTPASEPDEMULTIPLEMESHESMODIFIER_HPP_
#define TESTGTPASEPDEMULTIPLEMESHESMODIFIER_HPP_

#include <cxxtest/TestSuite.h>
#include "PetscSetupAndFinalize.hpp"
#include "GTPasePDEMultipleMeshesModifier.hpp"
#include "GTPasePDESystemParameters.hpp"

class TestGTPasePDEMultipleMeshesModifier : public CxxTest::TestSuite
{

private:

    static const unsigned DIM = 2;

    TetrahedralSubsetMesh<DIM,DIM>* mpBiggerSubsetMesh;
    TetrahedralSubsetMesh<DIM,DIM>* mpSmallerSubsetMesh;

    GTPasePDEMultipleMeshesModifier<DIM> mGTPaseModifier;

    void setUp()
    {
        TetrahedralMesh<DIM,DIM> original_mesh;
        original_mesh.ConstructRectangularMesh(5,5);

        // 0,1,2,11,13
        std::vector<unsigned> element_mask;
        element_mask.push_back(0);
        element_mask.push_back(1);
        element_mask.push_back(2);
        element_mask.push_back(11);
        element_mask.push_back(13);

        // 2,11
        std::vector<unsigned> element_mask_subset;
        element_mask_subset.push_back(2);
        element_mask_subset.push_back(11);

        mpBiggerSubsetMesh = new TetrahedralSubsetMesh<DIM,DIM>(original_mesh, element_mask);
        mpSmallerSubsetMesh = new TetrahedralSubsetMesh<DIM,DIM>(original_mesh, element_mask_subset);
    }

    void tearDown()
    {
        delete mpBiggerSubsetMesh;
        delete mpSmallerSubsetMesh;
    }

    Vec CreateNodeWiseSolution()
    {
        std::vector<double> data;

        // Node 0 is the one shared by both triangles
        for (unsigned var_index=0; var_index<GTPASE_PROBLEM_DIM; ++var_index)
        {
            data.push_back(1.0);
        }

        // Nodes 1 and 2 belong to the first element
        for (unsigned var_index=0; var_index<2*GTPASE_PROBLEM_DIM; ++var_index)
        {
            data.push_back(4.0);
        }

        // Nodes 3 and 4 belong to the second element
        for (unsigned var_index=0; var_index<2*GTPASE_PROBLEM_DIM; ++var_index)
        {
            data.push_back(7.0);
        }

        return PetscTools::CreateVec(data);
    }

public:

    void TestComputeGlobalToLocalMap()
    {

        const std::map<unsigned, unsigned>& global_local_node_map = mGTPaseModifier.ComputeGlobalToLocalMap(mpBiggerSubsetMesh);
        const std::vector<unsigned>& global_node_indices = mpBiggerSubsetMesh->GetOriginalNodeIndices();

        // Test that global_local_node_map is the inverse of global_node_indices
        for (unsigned node_local_index = 0; node_local_index < mpBiggerSubsetMesh->GetNumNodes(); ++node_local_index)
        {
            unsigned global_node_index = global_node_indices[node_local_index];
            TS_ASSERT_EQUALS(global_local_node_map.at(global_node_index), node_local_index);
        }
    }

    void TestIntersectionAndDifference()
    {
        std::vector<unsigned> common_nodes;
        std::vector<unsigned> new_nodes;
        mGTPaseModifier.GetNodeSetIntersectionAndDifference(mpBiggerSubsetMesh, mpSmallerSubsetMesh, common_nodes, new_nodes);

        TS_ASSERT_EQUALS(common_nodes.size(), 5u);
        TS_ASSERT_EQUALS(common_nodes[0], 1u);
        TS_ASSERT_EQUALS(common_nodes[1], 6u);
        TS_ASSERT_EQUALS(common_nodes[2], 7u);
        TS_ASSERT_EQUALS(common_nodes[3], 8u);
        TS_ASSERT_EQUALS(common_nodes[4], 13u);

        TS_ASSERT_EQUALS(new_nodes.size(), 1u);
        TS_ASSERT_EQUALS(new_nodes[0], 0u);
    }

    void TestAverageNodeWiseSolutionAtElements()
    {
        Vec node_wise_solution = CreateNodeWiseSolution();
        Vec element_wise_solution = mGTPaseModifier.AverageNodeWiseSolutionAtElements(node_wise_solution, mpSmallerSubsetMesh);

        unsigned num_elements = PetscVecTools::GetSize(element_wise_solution)/GTPASE_PROBLEM_DIM;
        TS_ASSERT_EQUALS(num_elements, mpSmallerSubsetMesh->GetNumElements());

        for (unsigned var_index=0; var_index<GTPASE_PROBLEM_DIM; ++var_index)
        {
            // Node values of 1, 4, and 4 average to 3
            TS_ASSERT_EQUALS(3.0, PetscVecTools::GetElement(element_wise_solution, var_index));

            // Node values of 1, 7, and 7 average to 5
            TS_ASSERT_EQUALS(5.0, PetscVecTools::GetElement(element_wise_solution, GTPASE_PROBLEM_DIM+var_index));
        }

        PetscTools::Destroy(node_wise_solution);
        PetscTools::Destroy(element_wise_solution);
    }

    void TestInterpolationBetweenIdenticalMeshes()
    {
        Vec original_solution = CreateNodeWiseSolution();
        Vec interpolated_solution = mGTPaseModifier.InterpolateSolutionBetweenMeshes(mpSmallerSubsetMesh, mpSmallerSubsetMesh, original_solution);

        for (unsigned solution_vec_index=0; solution_vec_index<PetscVecTools::GetSize(original_solution); ++solution_vec_index)
        {
            TS_ASSERT_EQUALS(PetscVecTools::GetElement(original_solution, solution_vec_index),
                             PetscVecTools::GetElement(interpolated_solution, solution_vec_index));
        }

        PetscTools::Destroy(original_solution);
        PetscTools::Destroy(interpolated_solution);
    }

    void TestInterpolationBetweenDifferentMeshes()
    {
        // Test that interpolating from a smaller to a bigger mesh populates the new entries with default values

        Vec smaller_mesh_solution = CreateNodeWiseSolution();
        Vec bigger_mesh_solution = mGTPaseModifier.InterpolateSolutionBetweenMeshes(mpBiggerSubsetMesh, mpSmallerSubsetMesh, smaller_mesh_solution);

        std::vector<unsigned> node_map(5);
        node_map[0] = 2; // local node 0 in smaller mesh is local node 2 in bigger mesh
        node_map[1] = 1;
        node_map[2] = 4;
        node_map[3] = 5;
        node_map[4] = 0;

        for (unsigned smaller_mesh_node_index=0; smaller_mesh_node_index < node_map.size(); ++smaller_mesh_node_index)
        {
            unsigned bigger_mesh_node_index = node_map[smaller_mesh_node_index];

            for (unsigned var_index=0; var_index<GTPASE_PROBLEM_DIM; ++var_index)
            {
                TS_ASSERT_EQUALS(PetscVecTools::GetElement(smaller_mesh_solution, smaller_mesh_node_index*GTPASE_PROBLEM_DIM + var_index),
                                 PetscVecTools::GetElement(bigger_mesh_solution, bigger_mesh_node_index*GTPASE_PROBLEM_DIM + var_index));
            }
        }

        unsigned bigger_mesh_node_index = 3;
        for (unsigned var_index=0; var_index<GTPASE_PROBLEM_DIM; ++var_index)
        {
            TS_ASSERT_EQUALS(PetscVecTools::GetElement(bigger_mesh_solution, bigger_mesh_node_index*GTPASE_PROBLEM_DIM + var_index),
                            0.5 * (PetscVecTools::GetElement(bigger_mesh_solution, 0*GTPASE_PROBLEM_DIM + var_index) +
                                   PetscVecTools::GetElement(bigger_mesh_solution, 1*GTPASE_PROBLEM_DIM + var_index)));
        }


        PetscTools::Destroy(smaller_mesh_solution);
        PetscTools::Destroy(bigger_mesh_solution);
    }

    // TODO: test that interpolating from a smaller to a bigger mesh populates the new entries while conserving mass
    void TestInterpolationConservesFromSmallToBig()
    {
        TS_FAIL("Not implemented yet");
    }

    // TODO: test that interpolating from a bigger to a smaller mesh populates the new entries while conserving mass
    void TestInterpolationConservesFromBigToSmall()
    {
        TS_FAIL("Not implemented yet");
    }

};

#endif /* TESTGTPASEPDEMULTIPLEMESHESMODIFIER_HPP_ */
