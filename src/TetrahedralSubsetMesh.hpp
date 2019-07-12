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

#ifndef TETRAHEDRALSUBSETMESH_HPP_
#define TETRAHEDRALSUBSETMESH_HPP_

#include "TetrahedralMesh.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class TetrahedralSubsetMesh: public TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>
{

private:

	/// Vector containing the indices in the original mesh for all the elements in this mesh.
	std::vector<unsigned> mOriginalElementIndices;

	/// Vector containing the indices in the original mesh for all the elements in this mesh.
	std::vector<unsigned> mOriginalNodeIndices;

	/**
	 * Helper method that takes two nodes in a mesh and returns the set of mesh elements containing both nodes.
	 *
	 * @param pNodeA pointer to first node
	 * @param pNodeB pointer to second node
	 * @return set of indices of mesh elements containing both nodes
	 */
	std::set<unsigned> GetIntersectionContainintElements(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB) const;

	/**
	 * Helper method that, given two nodes, creates a boundary edge taking care that the edge normal is pointing outwards wrt the mesh
	 *
	 * @param boundaryElementIndex the index of the new boundary element, i.e. the edge
	 * @param pEdgeNodeA pointer to first node defining the edge
	 * @param pEdgeNodeB pointer to second node defining the edge
	 * @return the new boundary element defining an edge between pEdgeNodeA and pEdgeNodeB
	 */
	BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>* CreateEdgeWithCorrectNormal(unsigned boundaryElementIndex, Node<SPACE_DIM>* pEdgeNodeA, Node<SPACE_DIM>* pEdgeNodeB, Element<ELEMENT_DIM, SPACE_DIM>& rContainingElement) const;

public:
	/**
	 * Constructor. Creates a mesh made of a subset of elements of another mesh.
	 *
	 * @param rOriginalMesh the original mesh
	 * @param rElementSubset the subset to elements defining the new mesh
	 */
	TetrahedralSubsetMesh(const TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rOriginalMesh, const std::vector<unsigned>& rElementSubset);

	/**
	 * Constructor required for archiving.
	 */
	TetrahedralSubsetMesh();

	/**
	 * Return the indices in the original mesh for all the elements in this mesh (mOriginalElementIndices)
	 *
	 * @return const reference to mOriginalElementIndices
	 */
	const std::vector<unsigned>& GetOriginalElementIndices() const;

	/**
	 * Return the indices in the original mesh for all the nodes in this mesh (mOriginalNodeIndices)
	 *
	 * @return const reference to mOriginalNodeIndices
	 */
	const std::vector<unsigned>& GetOriginalNodeIndices() const;

	/**
	 * For a given node index return the indices of its neighbouring nodes.
	 *
	 * @return set of neighbouring nodes
	 */
	std::set<unsigned> GetNeighbouringNodes(unsigned nodeIndex) const;

};

#endif /* TETRAHEDRALSUBSETMESH_HPP_ */
