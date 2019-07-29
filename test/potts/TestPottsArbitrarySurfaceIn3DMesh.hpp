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
#ifndef TESTPOTTSARBITRARYSURFACEIN3DMESH_HPP_
#define TESTPOTTSARBITRARYSURFACEIN3DMESH_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "PottsArbitraryMeshFromMutableMeshGenerator.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"
#include "RandomNumberGenerator.hpp"
#include "VtkMeshWriter.hpp"
#include "Debug.hpp"

class TestPottsArbitrarySurfaceIn3DMesh : public AbstractCellBasedTestSuite
{
private:

    double mSite3 Volume;
    double mSite5Volume;
    double mSite3SurfaceArea;
    double mSite5SurfaceArea;
    double mContactAreaBetween3And5;
    double mElement0SurfaceArea;

    void setUp()
    {
        mSite3Volume = 0.25;
        mSite5Volume = 0.25;
        mSite3SurfaceArea = 1 + 2*sqrt(0.5);
        mSite5SurfaceArea = mSite3SurfaceArea;
        mContactAreaBetween3And5 = 1.0;
        mElement0SurfaceArea = 4*sqrt(0.5);
    }

public:

    void TestSimpleGenerateArbitraryPottsMesh() throw(Exception)
    {
        TrianglesMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/simple_2d_in_3d");
        MutableMesh<2,3> original_delaunay_mesh;
        original_delaunay_mesh.ConstructFromMeshReader(mesh_reader);

        PottsArbitraryMeshFromMutableMeshGenerator<3> potts_generator(original_delaunay_mesh);
        PottsArbitrarySurfaceIn3DMesh<3>* p_potts_mesh = potts_generator.GetMesh();

        TS_ASSERT_EQUALS(p_potts_mesh->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(p_potts_mesh->GetNumAllNodes(), 8u);
        TS_ASSERT_EQUALS(p_potts_mesh->GetNumElements(), 0u);

        TS_ASSERT(p_potts_mesh->GetNode(0)->IsBoundaryNode());

        // Test Neighbourhoods
        std::set<unsigned> moore_neighbouring_sites = p_potts_mesh->GetMooreNeighbouringNodeIndices(4);
        std::set<unsigned> von_neumann_neighbouring_sites = p_potts_mesh->GetVonNeumannNeighbouringNodeIndices(4);

        TS_ASSERT_EQUALS(moore_neighbouring_sites.size(), 2u);
        TS_ASSERT_EQUALS(von_neumann_neighbouring_sites.size(), 2u);

        std::set<unsigned> expected_neighbouring_sites;
        expected_neighbouring_sites.insert(1);
        expected_neighbouring_sites.insert(7);

        TS_ASSERT_EQUALS(moore_neighbouring_sites, expected_neighbouring_sites);
        TS_ASSERT_EQUALS(von_neumann_neighbouring_sites, expected_neighbouring_sites);
   }

    void TestVolumeAndSurfaceAreaCalculations() throw(Exception)
    {
        TrianglesMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/simple_2d_in_3d");
        MutableMesh<2,3> original_delaunay_mesh;
        original_delaunay_mesh.ConstructFromMeshReader(mesh_reader);

        PottsArbitraryMeshFromMutableMeshGenerator<3> potts_generator(original_delaunay_mesh);
        PottsArbitrarySurfaceIn3DMesh<3>* p_potts_mesh = potts_generator.GetMesh();

        TS_ASSERT_DELTA(p_potts_mesh->GetVolumeOfLatticeSite(3), mSite3Volume, 1e-6);
        TS_ASSERT_DELTA(p_potts_mesh->GetVolumeOfLatticeSite(5), mSite5Volume, 1e-6)

        TS_ASSERT_DELTA(p_potts_mesh->GetSurfaceAreaOfLatticeSite(3), mSite3SurfaceArea, 1e-6);
        TS_ASSERT_DELTA(p_potts_mesh->GetSurfaceAreaOfLatticeSite(5), mSite5SurfaceArea, 1e-6);

        TS_ASSERT_DELTA(p_potts_mesh->GetContactAreaBetweenLatticeSite(3,5), mContactAreaBetween3And5, 1e-6);

        // Not neighbours but sharing a node
        TS_ASSERT_DELTA(p_potts_mesh->GetContactAreaBetweenLatticeSite(3,6), 0., 1e-6);

        // Not neighbours and sharing no node
        TS_ASSERT_DELTA(p_potts_mesh->GetContactAreaBetweenLatticeSite(4,2), 0., 1e-6);

        // Add a PottsElement made of three nodes
        std::vector<Node<3>*> element_nodes;
        element_nodes.push_back(p_potts_mesh->GetNode(3));
        element_nodes.push_back(p_potts_mesh->GetNode(5));
        p_potts_mesh->AddElement(new PottsElement<3>(0, element_nodes));

        TS_ASSERT_DELTA(p_potts_mesh->GetVolumeOfElement(0), mSite3Volume+mSite5Volume, 1e-6);
        TS_ASSERT_DELTA(p_potts_mesh->GetSurfaceAreaOfElement(0), mElement0SurfaceArea, 1e-6);
    }

    void TestUpdatePottsNodeLocationFromDelaunay()
    {
        TrianglesMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/simple_2d_in_3d");
        MutableMesh<2,3> original_delaunay_mesh;
        original_delaunay_mesh.ConstructFromMeshReader(mesh_reader);

        PottsArbitraryMeshFromMutableMeshGenerator<3> potts_generator(original_delaunay_mesh);
        PottsArbitrarySurfaceIn3DMesh<3>* p_potts_mesh = potts_generator.GetMesh();

        // Add a PottsElement made of three nodes
        std::vector<Node<3>*> element_nodes;
        element_nodes.push_back(p_potts_mesh->GetNode(3));
        element_nodes.push_back(p_potts_mesh->GetNode(5));
        p_potts_mesh->AddElement(new PottsElement<3>(0, element_nodes));

        // Use copy constructor rather than taking a reference, otherwise the location will be updated by Scale()
        const c_vector<double, 3> node_0_original_location(p_potts_mesh->GetNode(0)->rGetLocation());

        original_delaunay_mesh.Scale(2.0,2.0,2.0);
        p_potts_mesh->UpdatePottsNodeLocationFromDelaunay();

        const c_vector<double, 3>& node_0_updated_location = p_potts_mesh->GetNode(0)->rGetLocation();

        TS_ASSERT_EQUALS(2*node_0_original_location[0], node_0_updated_location[0]);
        TS_ASSERT_EQUALS(2*node_0_original_location[1], node_0_updated_location[1]);
        TS_ASSERT_EQUALS(2*node_0_original_location[2], node_0_updated_location[2]);

        TS_ASSERT_DELTA(p_potts_mesh->GetVolumeOfLatticeSite(3), 4*mSite3Volume, 1e-6);
        TS_ASSERT_DELTA(p_potts_mesh->GetVolumeOfLatticeSite(5), 4*mSite5Volume, 1e-6)

        TS_ASSERT_DELTA(p_potts_mesh->GetSurfaceAreaOfLatticeSite(3), 2*mSite3SurfaceArea, 1e-6);
        TS_ASSERT_DELTA(p_potts_mesh->GetSurfaceAreaOfLatticeSite(5), 2*mSite5SurfaceArea, 1e-6)

        TS_ASSERT_DELTA(p_potts_mesh->GetVolumeOfElement(0), 4.0*(mSite3Volume+mSite5Volume), 1e-6);
        TS_ASSERT_DELTA(p_potts_mesh->GetSurfaceAreaOfElement(0), 2.0*mElement0SurfaceArea, 1e-6);

        TS_ASSERT_DELTA(p_potts_mesh->GetContactAreaBetweenLatticeSite(3,5), 2.0*mContactAreaBetween3And5, 1e-6);

        // Useful for debugging
        // VtkMeshWriter<2,3> mesh_writer("TestPottsArbitrarySurfaceIn3DMesh<3>", "delaunay", false);
        // mesh_writer.WriteFilesUsingMesh(original_delaunay_mesh);
    }
};

#endif /*TESTPOTTSARBITRARYSURFACEIN3DMESH_HPP_*/
