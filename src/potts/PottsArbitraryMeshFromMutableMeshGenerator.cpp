/*

Copyright (c) 2005-2014, University of Oxford.
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

#include "PottsArbitraryMeshFromMutableMeshGenerator.hpp"
#include "Debug.hpp"

template<unsigned SPACE_DIM>
PottsArbitraryMeshFromMutableMeshGenerator<SPACE_DIM>::PottsArbitraryMeshFromMutableMeshGenerator(MutableMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
{
	assert(SPACE_DIM==2 || SPACE_DIM==3);

    std::vector<Node<SPACE_DIM>*> nodes;
    std::vector<std::set<unsigned> > node_neighbours;

    // Loop over all elements in the triangular mesh to create lattice sites
    for(typename MutableMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator iter = rMesh.GetElementIteratorBegin();
        iter != rMesh.GetElementIteratorEnd();
        ++iter)
    {
        // We consider a lattice site boundary if any of the nodes in the original triangle was a boundary node.
        bool is_element_boundary(false);

        // Neighbours are defined as all surrounding triangles sharing an edge
        std::set<unsigned> neighbouring_lattice_sites;

        // Loop over the edges of the triangle and assign as neighbouring all the triangles that share an edge with iter
        for (unsigned local_node_index=0; local_node_index<ELEMENT_DIM+1; ++local_node_index)
        {
            unsigned node_a_local_index = local_node_index;
            unsigned node_b_local_index = (local_node_index+1) % (ELEMENT_DIM+1);

            std::set<unsigned> node_a_elements = iter->GetNode(node_a_local_index)->rGetContainingElementIndices();
            std::set<unsigned> node_b_elements = iter->GetNode(node_b_local_index)->rGetContainingElementIndices();

            // Compute intersection of mesh_node->rGetContainingElementIndices()
            std::vector<unsigned> element_set_intersection(std::min(node_a_elements.size(), node_b_elements.size()));
            std::vector<unsigned>::iterator element_set_intersection_end;

            element_set_intersection_end = std::set_intersection (node_a_elements.begin(), node_a_elements.end(),
                                                                  node_b_elements.begin(), node_b_elements.end(),
                                                                  element_set_intersection.begin());
            element_set_intersection.resize(element_set_intersection_end - element_set_intersection.begin());

            neighbouring_lattice_sites.insert(element_set_intersection.begin(), element_set_intersection.end());

            is_element_boundary |= iter->GetNode(node_a_local_index)->IsBoundaryNode();
            // p_potts_mesh->GetNumNodes()
        }
        // At this point neighbouring_lattice_sites will contain the current element, remove it.
        neighbouring_lattice_sites.erase(iter->GetIndex());
        // PRINT_VECTOR(iter->CalculateCentroid());
        // PRINT_VARIABLE(is_element_boundary);

        // A node in a Potts mesh is the centroid of a lattice site.
        nodes.push_back(new Node<SPACE_DIM>(iter->GetIndex(), iter->CalculateCentroid(), is_element_boundary));


        // PRINT_VECTOR(nodes);
        node_neighbours.push_back(neighbouring_lattice_sites);
    }

    // We create the Potts mesh without any elements (cells). This will be defined later by the simulation.
    std::vector<PottsElement<SPACE_DIM>*> elements;
    // TRACE("E");
    mpMesh = new PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>(nodes, elements, node_neighbours, node_neighbours, &rMesh);
    // PRINT_VARIABLE(mpMesh->GetNumNodes());
}

template<unsigned SPACE_DIM>
PottsArbitraryMeshFromMutableMeshGenerator<SPACE_DIM>::~PottsArbitraryMeshFromMutableMeshGenerator()
{
    delete mpMesh;
}

template<unsigned SPACE_DIM>
PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>* PottsArbitraryMeshFromMutableMeshGenerator<SPACE_DIM>::GetMesh()
{
    return mpMesh;
}

template class PottsArbitraryMeshFromMutableMeshGenerator<2>;
template class PottsArbitraryMeshFromMutableMeshGenerator<3>;

//#include "SerializationExportWrapperForCpp.hpp"
//EXPORT_TEMPLATE_CLASS1(PottsArbitraryMeshFromMutableMeshGenerator, 2)
//EXPORT_TEMPLATE_CLASS1(PottsArbitraryMeshFromMutableMeshGenerator, 3)
