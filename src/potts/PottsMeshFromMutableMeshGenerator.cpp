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

#include "PottsMeshFromMutableMeshGenerator.hpp"

#include <boost/scoped_array.hpp>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PottsMeshFromMutableMeshGenerator<ELEMENT_DIM,SPACE_DIM>::PottsMeshFromMutableMeshGenerator(MutableMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
{
    std::vector<Node<SPACE_DIM>*> nodes;
    std::vector<PottsElement<SPACE_DIM>*>  elements;
    std::vector<std::set<unsigned> > node_neighbours;

    for(typename MutableMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator iter = rMesh.GetNodeIteratorBegin();
    	iter != rMesh.GetNodeIteratorEnd();
    	++iter)
    {
    	nodes.push_back(new Node<SPACE_DIM>(iter->GetIndex(),iter->rGetLocation(),iter->IsBoundaryNode()));

    	std::set<unsigned> neighbouring_nodes;

    	for(std::set<unsigned>::const_iterator elem_iter = iter->rGetContainingElementIndices().begin();
    		elem_iter != iter->rGetContainingElementIndices().end();
    		++elem_iter)
    	{
    		Element<ELEMENT_DIM, SPACE_DIM>* p_containing_element = rMesh.GetElement(*elem_iter);

    		for (unsigned node_index=0; node_index<ELEMENT_DIM+1; ++node_index)
    		{
    			neighbouring_nodes.insert(p_containing_element->GetNodeGlobalIndex(node_index));
    		}
    	}

    	node_neighbours.push_back(neighbouring_nodes);
    }

    mpMesh = new PottsMesh<SPACE_DIM>(nodes, elements, node_neighbours, node_neighbours);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PottsMeshFromMutableMeshGenerator<ELEMENT_DIM,SPACE_DIM>::~PottsMeshFromMutableMeshGenerator()
{
    delete mpMesh;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PottsMesh<SPACE_DIM>* PottsMeshFromMutableMeshGenerator<ELEMENT_DIM,SPACE_DIM>::GetMesh()
{
    return mpMesh;
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class PottsMeshFromMutableMeshGenerator<1,1>;
template class PottsMeshFromMutableMeshGenerator<2,2>;
template class PottsMeshFromMutableMeshGenerator<2,3>;
template class PottsMeshFromMutableMeshGenerator<3,3>;
