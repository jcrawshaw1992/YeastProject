

#include "ElementAnglesWriter.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementAnglesWriter<ELEMENT_DIM, SPACE_DIM>::ElementAnglesWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("ElementAnglesWriter.dat")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ElementAnglesWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{

    assert(SPACE_DIM==2 || SPACE_DIM==3); // LCOV_EXCL_LINE
    HistoryDepMeshBasedCellPopulation< ELEMENT_DIM,   SPACE_DIM>* p_cell_population = static_cast<HistoryDepMeshBasedCellPopulation<  ELEMENT_DIM,   SPACE_DIM>*>(pCellPopulation);

    // Iterate over all springs and add force contributions
    for (typename MeshBasedCellPopulation< ELEMENT_DIM,   SPACE_DIM>::SpringIterator spring_iterator = p_cell_population->SpringsBegin();
         spring_iterator != p_cell_population->SpringsEnd();
         ++spring_iterator)
    {
        Node<SPACE_DIM>* pNode1 = spring_iterator.GetNodeA();
        Node<SPACE_DIM>* pNode3 = spring_iterator.GetNodeB();

        std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> edge = std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*>(pNode1, pNode3);

        std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > nonUnitNormals;
        std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> otherNodes;

        bool boundary_edge_found = CalculateElementNormals(p_cell_population->rGetMesh(),
                                                           edge, nonUnitNormals, otherNodes);

        if (boundary_edge_found)
        {
            continue;
        }

        c_vector<double, SPACE_DIM> normal_1 = nonUnitNormals.first;
        normal_1 /= norm_2(normal_1);

        c_vector<double, SPACE_DIM> normal_2 = nonUnitNormals.second;
        normal_2 /= norm_2(normal_2);


        
        double Angle = acos(inner_prod(normal_1, normal_2));

        *this->mpOutStream << Angle << " ";
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ElementAnglesWriter<ELEMENT_DIM, SPACE_DIM>::CalculateElementNormals(MutableMesh<ELEMENT_DIM, SPACE_DIM>& rMesh, std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> edge,
                                                     std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> >& nonUnitNormals,
                                                     std::pair<Node<SPACE_DIM>*,  Node<SPACE_DIM>*>& otherNodes)
{
    Node<SPACE_DIM>* pNode1 = edge.first;
    Node<SPACE_DIM>* pNode3 = edge.second;

    /*
     *  Find common triangles
     */
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_node1 = pNode1->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_node3 = pNode3->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_node1.begin(),
                          elements_containing_node1.end(),
                          elements_containing_node3.begin(),
                          elements_containing_node3.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    switch (shared_elements.size())
    {
    case 1:
        // We found a boundary edge and we finish here
        return true;
    case 2:
        break;
    default:
        NEVER_REACHED;
    }

    std::set<unsigned>::iterator set_iter = shared_elements.begin();
    Element<ELEMENT_DIM,SPACE_DIM>* pElement1 = rMesh.GetElement(*set_iter);
    ++set_iter;
    Element<ELEMENT_DIM,SPACE_DIM>* pElement2 = rMesh.GetElement(*set_iter);

    // Find additional nodes
    Node<SPACE_DIM>* pNode2 = NULL;
    Node<SPACE_DIM>* pNode4 = NULL;
    for (unsigned local_index = 0; local_index < SPACE_DIM; ++local_index)
    {
        unsigned index_for_node2 = pElement1->GetNodeGlobalIndex(local_index);
        unsigned index_for_node4 = pElement2->GetNodeGlobalIndex(local_index);

        if ((index_for_node2 != pNode1->GetIndex()) && (index_for_node2 != pNode3->GetIndex()))
        {
            pNode2 = pElement1->GetNode(local_index);
        }

        if ((index_for_node4 != pNode1->GetIndex()) && (index_for_node4 != pNode3->GetIndex()))
        {
            pNode4 = pElement2->GetNode(local_index);
        }
    }
    assert(pNode2 != NULL);
    assert(pNode4 != NULL);

    // Calculate the force acting on each node
    c_vector<double, SPACE_DIM> vector_A = pNode1->rGetLocation() - pNode3->rGetLocation();
    c_vector<double,SPACE_DIM> vector_B = pNode2->rGetLocation() - pNode3->rGetLocation();
    c_vector<double, SPACE_DIM> normal_1 = VectorProduct(vector_A,vector_B);

    vector_A = pNode4->rGetLocation() - pNode3->rGetLocation();
    vector_B = pNode1->rGetLocation() - pNode3->rGetLocation();
    c_vector<double, SPACE_DIM> normal_2 = VectorProduct(vector_A,vector_B);

    nonUnitNormals = std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> >(normal_1, normal_2);
    otherNodes = std::pair<Node<SPACE_DIM>*,  Node<SPACE_DIM>*>(pNode2, pNode4);

    return false;
}






template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ElementAnglesWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("ElementAnglesWriter cannot be used with a CaBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ElementAnglesWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("ElementAnglesWriter cannot be used with a NodeBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ElementAnglesWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("ElementAnglesWriter cannot be used with a PottsBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ElementAnglesWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("ElementAnglesWriter cannot be used with a VertexBasedCellPopulation");
}

// Explicit instantiation
template class ElementAnglesWriter<1,1>;
template class ElementAnglesWriter<1,2>;
template class ElementAnglesWriter<2,2>;
template class ElementAnglesWriter<1,3>;
template class ElementAnglesWriter<2,3>;
template class ElementAnglesWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ElementAnglesWriter)
