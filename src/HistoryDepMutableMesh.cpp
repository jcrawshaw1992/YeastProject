/*


*/

#include <cstring>
#include <map>

#include "HistoryDepMutableMesh.hpp"
#include "OutputFileHandler.hpp"

//Jonathan Shewchuk's triangle and Hang Si's tetgen
#define REAL double
#define VOID void
#include "tetgen.h"
#include "triangle.h"
#undef REAL
#undef VOID
#include "Debug.hpp"
#include "VtkMeshWriter.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>::HistoryDepMutableMesh()
        : MutableMesh<ELEMENT_DIM, SPACE_DIM>()
{
    this->mMeshChangesDuringSimulation = true;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>::HistoryDepMutableMesh(std::vector<Node<SPACE_DIM>*> nodes)
        : MutableMesh<ELEMENT_DIM, SPACE_DIM>(nodes)
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>::~HistoryDepMutableMesh()
{
    //  this->Clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>::DeleteMesh()
{
    this->mNodes.clear();
    this->mElements.clear();
    this->mBoundaryNodes.clear();
    this->mNodePermutation.clear();
    this->mMeshChangesDuringSimulation = true;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>::AssignNewMesh(HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>* New_Mesh)
{

    // this->mNodes = New_Mesh->mNodes;
    // this->mElements  = New_Mesh->mElements;


    this->mpDistributedVectorFactory = New_Mesh->mpDistributedVectorFactory;
    for (int i = 0; i < New_Mesh->mNodes.size(); ++i)
    {
        Node<SPACE_DIM>* pNewNode = New_Mesh->mNodes[i];
        this->mNodes.push_back(pNewNode);
    }
    for (int i = 0; i < New_Mesh->mElements.size(); ++i)
    {
        Element<ELEMENT_DIM, SPACE_DIM>* pNewElement = New_Mesh->mElements[i];
        this->mElements.push_back(pNewElement);
    }
    for (int i = 0; i < New_Mesh->mBoundaryNodes.size(); ++i)
    {
        Node<SPACE_DIM>* pNewNode = New_Mesh->mBoundaryNodes[i];
        // pNewNode->SetAsBoundaryNode();
        this->mBoundaryNodes.push_back(pNewNode);
    }
    New_Mesh->mElements.clear();
    New_Mesh->mNodes.clear();
    New_Mesh->mBoundaryNodes.clear();

    this->RefreshMesh();
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>::CreateNewMesh(HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>* New_Mesh, std::map<unsigned, c_vector<double, SPACE_DIM> > InitalPositionOfRemeshedNodes)
{
    // this->mpDistributedVectorFactory = New_Mesh->mpDistributedVectorFactory;
    for (unsigned i = 0; i < New_Mesh->mNodes.size(); ++i)
    {
        PRINT_VECTOR(InitalPositionOfRemeshedNodes[i])
        Node<SPACE_DIM>* pNewNode = new Node<SPACE_DIM>(i, InitalPositionOfRemeshedNodes[i], false); // never on boundary
        this->mNodes.push_back(pNewNode);
    }
    for (unsigned i = 0; i < New_Mesh->mElements.size(); ++i)
    {
        Element<ELEMENT_DIM, SPACE_DIM>* pNewElement = New_Mesh->mElements[i];
        this->mElements.push_back(pNewElement);
    }
}


 



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>::AddANewNodeBehindBoundary()
{

    // THis function lets me see what happens if we add one extra node, and then later redo all the cells
    c_vector<double, SPACE_DIM> point;
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
template class HistoryDepMutableMesh<2, 2>;
template class HistoryDepMutableMesh<2, 3>;
template class HistoryDepMutableMesh<3, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(HistoryDepMutableMesh)

// TRACE("get rid of the new mesh that is redundant now")
// New_Mesh->Clear();
// New_Mesh->~HistoryDepMutableMesh();

// Node<3>* p_node_3 = new Node<3>(3, false, 1.0, 1.0, 0.0);
// unsigned new_node_index = node_based_cell_population.AddNode(p_node2);
// unsigned new_node_index = MutableMesh<2,2>::AddNode(new Node<2>(0, location));
// for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = New_Mesh->GetNodeIteratorBegin();
    //      node_iter != New_Mesh->GetNodeIteratorEnd();
    //       ++node_iter)
    // {
    //     // unsigned index =  node_iter->GetIndex(); // assume that the ordering matches
    //     // AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::AddCellUsingLocationIndex(index,*iter);
    //     Node<SPACE_DIM>* pNewNode = &*(node_iter);//New_Mesh->mNodes[i];
    //     this->mNodes.push_back(pNewNode);

    // }


    // Size of mesh is about to change
    // this->mpDistributedVectorFactory.clear();

    // this->mDeletedBoundaryElementIndices.clear();
    // this->mDeletedNodeIndices.clear();

    // TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::Clear();
    // this->mpDistributedVectorFactory = new DistributedVectorFactory(New_Mesh->GetNumNodes());
