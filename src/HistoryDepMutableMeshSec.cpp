/*

Copyright (c) 2005-2018, University of Oxford.
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

#include <map>
#include <cstring>

#include "HistoryDepMutableMesh.hpp"
#include "OutputFileHandler.hpp"

//Jonathan Shewchuk's triangle and Hang Si's tetgen
#define REAL double
#define VOID void
#include "triangle.h"
#include "tetgen.h"
#undef REAL
#undef VOID
#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>::HistoryDepMutableMesh()
                                              :MutableMesh<ELEMENT_DIM, SPACE_DIM>()
{
     TRACE("Constructor 1 HistoryDepMutableMesh.cpp")
    this->mMeshChangesDuringSimulation = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>::HistoryDepMutableMesh(std::vector<Node<SPACE_DIM> *> nodes)
                                              :MutableMesh<ELEMENT_DIM, SPACE_DIM>(nodes)
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>::~HistoryDepMutableMesh()
{
    //  TRACE("HistoryDep Mesh  Destructor")
    //  this->Clear();
    //  TRACE("Cleared hist dep mesh")
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>::DeleteMesh()
{
    this->mNodes.clear();
    this->mElements.clear();
    this->mBoundaryNodes.clear();
    this->mNodePermutation.clear();
    this->mMeshChangesDuringSimulation = true;
    
    PRINT_2_VARIABLES(this->mElements.size(), this->mNodes.size());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>::AssignNewMesh(HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>* New_Mesh)
{

    // Size of mesh is about to change
    // this->mpDistributedVectorFactory.clear();

    // this->mDeletedBoundaryElementIndices.clear();
    // this->mDeletedNodeIndices.clear();

    // TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Clear();
    // this->mpDistributedVectorFactory = new DistributedVectorFactory(New_Mesh->GetNumNodes());

    this->mpDistributedVectorFactory =  New_Mesh->mpDistributedVectorFactory;


    // for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = New_Mesh->GetNodeIteratorBegin(); 
    //      node_iter != New_Mesh->GetNodeIteratorEnd();
    //       ++node_iter)
    // {
    //     // unsigned index =  node_iter->GetIndex(); // assume that the ordering matches
    //     // AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::AddCellUsingLocationIndex(index,*iter);
    //     Node<SPACE_DIM>* pNewNode = &*(node_iter);//New_Mesh->mNodes[i];
    //     this->mNodes.push_back(pNewNode);
        
    // }



    for (int i=0; i< New_Mesh->mNodes.size(); ++i)
    {
        Node<SPACE_DIM>* pNewNode = New_Mesh->mNodes[i];
        this->mNodes.push_back(pNewNode);
    }
    for (int i=0; i< New_Mesh->mElements.size(); ++i)
    {  
         Element<ELEMENT_DIM,SPACE_DIM>* pNewElement = New_Mesh->mElements[i];
         this->mElements.push_back(pNewElement);

    }
     for (int i=0; i< New_Mesh->mBoundaryNodes.size(); ++i)
    {
        Node<SPACE_DIM>* pNewNode = New_Mesh->mBoundaryNodes[i];
        pNewNode->SetAsBoundaryNode();
        this->mBoundaryNodes.push_back(pNewNode);
    }
     New_Mesh->mElements.clear();
     New_Mesh->mNodes.clear();
    New_Mesh->mBoundaryNodes.clear();


    this->RefreshMesh();

    // TRACE("get rid of the new mesh that is redundant now")
    // New_Mesh->Clear();
    // New_Mesh->~HistoryDepMutableMesh();

}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>::AddANewNodeBehindBoundary()
{

    // Node<3>* p_node_3 = new Node<3>(3, false, 1.0, 1.0, 0.0);
    // unsigned new_node_index = node_based_cell_population.AddNode(p_node2);
    // unsigned new_node_index = MutableMesh<2,2>::AddNode(new Node<2>(0, location));

    c_vector<double,SPACE_DIM> point;
    //  1e-3,10e-3);
        point[2] = 12e-3;
        point[1] = 0;
        point[0] = 1e-3;
        Node<SPACE_DIM>* p_node = new Node<SPACE_DIM>(this->GetNumNodes(), point);

        unsigned new_index = this->AddNode(p_node);

    // I need to get the boudary nodes -- if possible grab a pair, but im not sure I have pair info... 
    // Add an extra node that is off the edge of the cylinder
    // this->mBoundaryNodes
    // this->mElements


    // Add the element connecting the new node

}





// Explicit instantiation
template class HistoryDepMutableMesh<2,2>;
template class HistoryDepMutableMesh<2,3>;
template class HistoryDepMutableMesh<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(HistoryDepMutableMesh)
