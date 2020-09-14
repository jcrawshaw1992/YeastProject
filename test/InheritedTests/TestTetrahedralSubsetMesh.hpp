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

#ifndef TESTTETRAHEDRALSUBSETMESH_HPP_
#define TESTTETRAHEDRALSUBSETMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "TetrahedralSubsetMesh.hpp"
#include "VtkMeshWriter.hpp"

class TestTetrahedralSubsetMesh : public CxxTest::TestSuite
{

    TetrahedralMesh<2,2> mOriginalMesh;
    TetrahedralSubsetMesh<2,2>* mpSubsetMesh;

	void setUp()
	{
        mOriginalMesh.ConstructRectangularMesh(5,5);

        // Writes out the generated mesh. Useful for determining the element indices used for cell creation
        VtkMeshWriter<2,2> mesh_writer("TestTetrahedralSubsetMesh", "mOriginalMesh");
        mesh_writer.WriteFilesUsingMesh(mOriginalMesh);

        std::vector<unsigned> element_mask;
        element_mask.push_back(0);
        element_mask.push_back(1);
        element_mask.push_back(2);
        element_mask.push_back(11);
        element_mask.push_back(13);

        mpSubsetMesh = new TetrahedralSubsetMesh<2,2>(mOriginalMesh, element_mask);
	}

	void tearDown()
	{
		delete mpSubsetMesh;
	}

public:

	void TestNumberOfNodesElements()
	{
        TS_ASSERT_EQUALS(mpSubsetMesh->GetNumElements(), 5u);
        TS_ASSERT_EQUALS(mpSubsetMesh->GetNumNodes(), 6u);
	}

	void TestNodeLocations()
	{
		const c_vector<double, 2>& original_mesh_node_0_location = mpSubsetMesh->GetNode(0)->rGetLocation();
		const c_vector<double, 2>& subset_mesh_node_1_location = mOriginalMesh.GetNode(6)->rGetLocation();

        TS_ASSERT_EQUALS(original_mesh_node_0_location[0], subset_mesh_node_1_location[0]);
        TS_ASSERT_EQUALS(original_mesh_node_0_location[1], subset_mesh_node_1_location[1]);
	}

	void TestWorkingOutBoundaryNodesElements()
	{
        TS_ASSERT_EQUALS(mpSubsetMesh->GetNumBoundaryElements(), 5u);
        TS_ASSERT_EQUALS(mpSubsetMesh->GetNumBoundaryNodes(), 5u);
	}

	void TestGetOriginalElementIndices()
	{
		const std::vector<unsigned>& original_indices = mpSubsetMesh->GetOriginalElementIndices();

		TS_ASSERT_EQUALS(original_indices[0], 0u);
		TS_ASSERT_EQUALS(original_indices[1], 1u);
		TS_ASSERT_EQUALS(original_indices[2], 2u);
		TS_ASSERT_EQUALS(original_indices[3], 11u);
		TS_ASSERT_EQUALS(original_indices[4], 13u);
	}

	void TestGetOriginalNodeIndices()
	{
		const std::vector<unsigned>& original_indices = mpSubsetMesh->GetOriginalNodeIndices();

		TS_ASSERT_EQUALS(original_indices[0], 6u);
		TS_ASSERT_EQUALS(original_indices[1], 1u);
		TS_ASSERT_EQUALS(original_indices[2], 7u);
		TS_ASSERT_EQUALS(original_indices[3], 0u);
		TS_ASSERT_EQUALS(original_indices[4], 8u);
		TS_ASSERT_EQUALS(original_indices[5], 13u);
	}

	void TestGetNeighbouringNodes()
	{
		const std::set<unsigned>& neighbours = mpSubsetMesh->GetNeighbouringNodes(0);

		TS_ASSERT_EQUALS(neighbours.size(), 4u);

		std::set<unsigned>::const_iterator neighbour_iter = neighbours.begin();
		TS_ASSERT_EQUALS(*(neighbour_iter++), 1u);
		TS_ASSERT_EQUALS(*(neighbour_iter++), 2u);
		TS_ASSERT_EQUALS(*(neighbour_iter++), 3u);
		TS_ASSERT_EQUALS(*(neighbour_iter++), 5u);
	}

	void TestBoundaryNormals()
	{
	    TS_ASSERT_THROWS_NOTHING(mpSubsetMesh->CheckOutwardNormals());
	}

};

#endif /* TESTTETRAHEDRALSUBSETMESH_HPP_ */

