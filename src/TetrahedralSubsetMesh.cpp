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

#include "TetrahedralSubsetMesh.hpp"
#include <map>
#include <boost/bimap.hpp>

using boost::bimap;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TetrahedralSubsetMesh<ELEMENT_DIM, SPACE_DIM>::TetrahedralSubsetMesh(const TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rOriginalMesh, const std::vector<unsigned>& rElementSubset):
	mOriginalElementIndices(rElementSubset)
{
    unsigned node_index(0);
    unsigned element_index(0);
    unsigned boundary_element_index(0);

    // Auxiliary data structure used to keep track of newly create nodes. We need to avoid creating multiple
    // instances of the same node when the original one is visited several times via multiple mesh elements. It
    // maps node objects in the original mesh and node objects in the subset mesh
    typedef bimap<Node<SPACE_DIM>*, Node<SPACE_DIM>*> NodeMapType;
    typedef typename NodeMapType::value_type NodePair;
    NodeMapType original_vs_subset_nodes_map;

    // Create new nodes and elements
    for (std::vector<unsigned>::const_iterator original_element_index = rElementSubset.begin(); original_element_index != rElementSubset.end(); ++original_element_index)
    {
        std::vector<Node<SPACE_DIM>*> nodes_in_element(ELEMENT_DIM+1);
        for (unsigned node_local_index=0; node_local_index<ELEMENT_DIM+1; ++node_local_index)
        {
            // These are the node in the original mesh and its index
            Node<SPACE_DIM>* node_in_original_mesh = rOriginalMesh.GetElement(*original_element_index)->GetNode(node_local_index);

            // If this node in the original mesh was already visited via a different element, we will retrieve it from an auxiliary
            // data structure. Otherwise create and store it and update the auxiliary data structure.
            Node<SPACE_DIM>* node_in_element;
            try
            {
                node_in_element = original_vs_subset_nodes_map.left.at(node_in_original_mesh);
            }
            catch(const std::out_of_range& exception)
            {
                node_in_element = new Node<SPACE_DIM>(node_index, node_in_original_mesh->rGetLocation());
                this->mNodes.push_back(node_in_element);

                original_vs_subset_nodes_map.insert(NodePair(node_in_original_mesh, node_in_element));
                mOriginalNodeIndices.push_back(node_in_original_mesh->GetIndex());

                node_index++;
            }
            nodes_in_element[node_local_index] = node_in_element;
        }
        this->mElements.push_back(new Element<ELEMENT_DIM, SPACE_DIM>(element_index, nodes_in_element));
        element_index++;
    }

    // Work out which nodes are boundary nodes
    typedef typename NodeMapType::left_const_iterator NodeMapLeftIterator;
    for (NodeMapLeftIterator index_node_pair = original_vs_subset_nodes_map.left.begin();
         index_node_pair != original_vs_subset_nodes_map.left.end(); ++index_node_pair)
    {
        Node<SPACE_DIM>* original_node = index_node_pair->first;
        Node<SPACE_DIM>* subset_node = index_node_pair->second;

        // The condition to be a boundary node in the new mesh is having been a boundary node in the original mesh or being
        // a node in the new mesh contained in fewer elements that its original counterpart.
        bool is_subset_node_boundary = original_node->IsBoundaryNode() ||
        							   (subset_node->GetNumContainingElements() < original_node->GetNumContainingElements());

        if (is_subset_node_boundary)
        {
            subset_node->SetAsBoundaryNode();
            this->mBoundaryNodes.push_back(subset_node);
        }
    }

    // Create boundary elements
    for (typename TetrahedralSubsetMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator subset_element = this->GetElementIteratorBegin();
         subset_element != this->GetElementIteratorEnd(); ++subset_element)
    {
    	// The following two nested loops will visit all the edges of a triangle just once.
    	///\todo This operation needs to be generalised to arbitrary ELEMENT_DIM, i.e. current code doesn't work for tets
    	assert(ELEMENT_DIM==2);
    	for (unsigned node_a_index = 0; node_a_index != subset_element->GetNumNodes(); ++node_a_index)
    	{
        	for (unsigned node_b_index = node_a_index+1; node_b_index != subset_element->GetNumNodes(); ++node_b_index)
        	{
        		Node<SPACE_DIM>* subset_node_a = subset_element->GetNode(node_a_index);
        		Node<SPACE_DIM>* subset_node_b = subset_element->GetNode(node_b_index);

        		bool both_boundary_in_subset = subset_node_a->IsBoundaryNode() && subset_node_b->IsBoundaryNode();
        		unsigned cardinality_subset_containing_elements_intersection = GetIntersectionContainintElements(subset_node_a, subset_node_b).size();

        		Node<SPACE_DIM>* original_node_a = original_vs_subset_nodes_map.right.at(subset_node_a);
        		Node<SPACE_DIM>* original_node_b = original_vs_subset_nodes_map.right.at(subset_node_b);

        		bool both_boundary_in_original = original_node_a->IsBoundaryNode() && original_node_b->IsBoundaryNode();
        		unsigned cardinality_original_containing_elements_intersection = GetIntersectionContainintElements(original_node_a, original_node_b).size();

        		// The condition for a pair of nodes to define a boundary edge is being boundary nodes in the original mesh and sharing a single element
        		// (i.e. having already defined a boundary element in the original mesh) or being boundary nodes in the new mesh and sharing a single
        		// element (i.e. being at the boundary of the new mesh without having been in the original mesh
        		if ((both_boundary_in_original && cardinality_original_containing_elements_intersection==1) ||
        			(both_boundary_in_subset && cardinality_subset_containing_elements_intersection==1))
        		{
        		    this->mBoundaryElements.push_back(CreateEdgeWithCorrectNormal(boundary_element_index, subset_node_a, subset_node_b, *subset_element));
        		    boundary_element_index++;
        		}
        	}
    	}
    }

    this->RefreshMesh();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TetrahedralSubsetMesh<ELEMENT_DIM, SPACE_DIM>::TetrahedralSubsetMesh()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>* TetrahedralSubsetMesh<ELEMENT_DIM, SPACE_DIM>::CreateEdgeWithCorrectNormal(unsigned boundaryElementIndex, Node<SPACE_DIM>* pEdgeNodeA, Node<SPACE_DIM>* pEdgeNodeB, Element<ELEMENT_DIM, SPACE_DIM>& rContainingElement) const
{
    Node<SPACE_DIM>* p_opposite_node = NULL;

    for (unsigned node_index = 0; node_index != rContainingElement.GetNumNodes(); ++node_index)
    {
        p_opposite_node = rContainingElement.GetNode(node_index);
        if (p_opposite_node != pEdgeNodeA && p_opposite_node != pEdgeNodeB)
        {
            break;
        }
    }
    assert(p_opposite_node);

    std::vector<Node<SPACE_DIM>*> nodes;
    nodes.push_back(pEdgeNodeA);
    nodes.push_back(pEdgeNodeB);

    BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>* p_boundary_elem = new BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>(boundaryElementIndex, nodes);

    c_vector<double, SPACE_DIM> into_mesh = p_opposite_node->rGetLocation() - p_boundary_elem->CalculateCentroid();
    c_vector<double, SPACE_DIM> normal = p_boundary_elem->CalculateNormal();

    // If the edge normal is pointing inwards, swap the node ordering
    if (inner_prod(into_mesh, normal) > 0.0)
    {
        Node<SPACE_DIM>* temp0 = p_boundary_elem->GetNode(0);
        Node<SPACE_DIM>* temp1 = p_boundary_elem->GetNode(1);
        p_boundary_elem->UpdateNode(1, temp0);
        p_boundary_elem->UpdateNode(0, temp1);
    }

    return p_boundary_elem;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> TetrahedralSubsetMesh<ELEMENT_DIM, SPACE_DIM>::GetIntersectionContainintElements(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB) const
{
	std::set<unsigned>& node_a_elements = pNodeA->rGetContainingElementIndices();
	std::set<unsigned>& node_b_elements = pNodeB->rGetContainingElementIndices();

	std::set<unsigned> element_intersecion;
	std::set_intersection(node_a_elements.begin(), node_a_elements.end(),
						  node_b_elements.begin(), node_b_elements.end(),
						  std::inserter(element_intersecion, element_intersecion.begin()));

	return element_intersecion;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const std::vector<unsigned>& TetrahedralSubsetMesh<ELEMENT_DIM, SPACE_DIM>::GetOriginalElementIndices() const
{
	return mOriginalElementIndices;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const std::vector<unsigned>& TetrahedralSubsetMesh<ELEMENT_DIM, SPACE_DIM>::GetOriginalNodeIndices() const
{
	return mOriginalNodeIndices;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> TetrahedralSubsetMesh<ELEMENT_DIM, SPACE_DIM>::GetNeighbouringNodes(unsigned nodeIndex) const
{
	std::set<unsigned> neighbouring_nodes;

	const std::set<unsigned>& containing_elements = this->GetNode(nodeIndex)->rGetContainingElementIndices();

	for (std::set<unsigned>::const_iterator element_index = containing_elements.begin(); element_index != containing_elements.end(); ++element_index)
	{
		Element<ELEMENT_DIM, SPACE_DIM>* element = this->GetElement(*element_index);

		for (unsigned node_local_index = 0; node_local_index<element->GetNumNodes(); ++node_local_index)
		{
			unsigned node_global_index = element->GetNode(node_local_index)->GetIndex();

			if (node_global_index != nodeIndex)
			{
				neighbouring_nodes.insert(node_global_index);
			}
		}
	}

	return neighbouring_nodes;
}


// Explicit instantiation
template class TetrahedralSubsetMesh<2, 2>;

// Serialisation
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS2(TetrahedralSubsetMesh, 2, 2)
