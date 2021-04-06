#include "PottsArbitrarySurfaceIn3DMesh.hpp"
#include <cxxtest/TestSuite.h>
#include "MathsFunctions.hpp"

template <unsigned SPACE_DIM>
PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::PottsArbitrarySurfaceIn3DMesh(std::vector<Node<SPACE_DIM>*> nodes,
                                                                        std::vector<PottsElement<SPACE_DIM>*> pottsElements,
                                                                        std::vector<std::set<unsigned> > vonNeumannNeighbouringNodeIndices,
                                                                        std::vector<std::set<unsigned> > mooreNeighbouringNodeIndices,
                                                                        MutableMesh<2, SPACE_DIM>* pDelaunayMesh) : PottsMesh<SPACE_DIM>(nodes, pottsElements, vonNeumannNeighbouringNodeIndices, mooreNeighbouringNodeIndices),
                                                                                                                    mpDelaunayMesh(pDelaunayMesh)
{
    assert(SPACE_DIM == 2 || SPACE_DIM == 3);

    mMeshElementMidPoints.clear();
    mLatticeVolume.clear();
    mPeriodic = 1;
    // TRACE("Potts mesh is constructed");
}

template <unsigned SPACE_DIM>
PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::~PottsArbitrarySurfaceIn3DMesh()
{
}

template <unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::SetBoundaries(std::vector<unsigned> BoundaryVector)
{
    mBoundaryVector = BoundaryVector;
}

template <unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::SetMeshSize(double N_C, double N_D, double width, double length)
{
    mN_C = N_C;
    mN_D = N_D;
    mWidth = width;
    mLength = length;
}

template <unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::UpdatePottsNodeLocationFromDelaunay()
{

    // This code was edited by Jess.
    // This function updates the node locations in the potts mesh based on where the nodes have been moved to in the Delaunay mesh
    typename MutableMesh<2, SPACE_DIM>::NodeIterator node_iter = mpDelaunayMesh->GetNodeIteratorBegin();
    typename PottsMesh<SPACE_DIM>::NodeIterator potts_node_iter = this->GetNodeIteratorBegin();

    for (;
         node_iter != mpDelaunayMesh->GetNodeIteratorEnd();
         ++node_iter, ++potts_node_iter)
    {

        // PRINT_2_VARIABLES(potts_node_iter->GetIndex(), node_iter->GetIndex());
        assert(potts_node_iter->GetIndex() == node_iter->GetIndex());
        potts_node_iter->rGetModifiableLocation() = node_iter->rGetLocation();
    }

    // TRACE("Get Edges");
    GetEdges();
    EdgeLenghts();

    mAngleRange = 0;
    // MapCylinderToPlane();
    // TRACE("Second");

    // TRACE("Third");
    CalculateLatticeVolumes();

    // TRACE("Fifth");
    // CalculateTraction();
}

template <unsigned SPACE_DIM>
inline double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::ComputeTriangleArea(const c_vector<double, SPACE_DIM>& vertexA, const c_vector<double, SPACE_DIM>& vertexB, const c_vector<double, SPACE_DIM>& vertexC) const
{
    const c_vector<double, SPACE_DIM> AC = vertexC - vertexA;
    const c_vector<double, SPACE_DIM> AB = vertexB - vertexA;
    return 0.5 * norm_2(VectorProduct(AC, AB));
}

template <unsigned SPACE_DIM>
bool PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::DoNodesShareElement(unsigned indexNodeA, unsigned indexNodeB)
{
    std::set<unsigned> node_a_elements = this->GetNode(indexNodeA)->rGetContainingElementIndices();
    std::set<unsigned> node_b_elements = this->GetNode(indexNodeB)->rGetContainingElementIndices();

    bool is_a_medium = node_a_elements.empty();
    bool is_b_medium = node_b_elements.empty();

    if (is_a_medium xor is_b_medium)
    {
        return false;
    }

    if (is_a_medium and is_b_medium)
    {
        return true;
    }
    else
    {
        // PRINT_2_VARIABLES(node_a_elements.size(), node_b_elements.size());

        assert(node_a_elements.size() == 1);
        assert(node_b_elements.size() == 1);

        return *node_a_elements.begin() == *node_b_elements.begin();
    }
}

template <unsigned SPACE_DIM>
bool PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::IsLatticeInElement(unsigned Latticeindex, unsigned pottsElementIndex)
{
    std::set<unsigned> Lattice_Element = this->GetNode(Latticeindex)->rGetContainingElementIndices();

    bool is_a_medium = Lattice_Element.empty();

    if (is_a_medium)
    {
        // Is the lattice part of the medium
        return false;
    }
    else
    {
        assert(Lattice_Element.size() == 1); // double check this lattice site is only in one element
        return *Lattice_Element.begin() == pottsElementIndex;
    }
}

template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetContactAreaBetweenLatticeSite(unsigned index_A, unsigned index_B)
{
    // TRACE("GetContactAreaBetweenLatticeSite")
    assert(index_A != index_B);

    // In this piece of code I need to get the contact area between two lattices sites
    // to do this I will find the Neighbouring element pairs that these lattice sites have in common
    //
    std::vector<std::pair<unsigned, unsigned> > ElementPairsToLattiveA = mMapNodesToAssociateElementPairs[index_A];
    std::vector<std::pair<unsigned, unsigned> > ElementPairsToLattiveB = mMapNodesToAssociateElementPairs[index_B];

    std::vector<std::pair<unsigned, unsigned> > CommonElements;

    std::set_intersection(ElementPairsToLattiveA.begin(), ElementPairsToLattiveA.end(),
                          ElementPairsToLattiveB.begin(), ElementPairsToLattiveB.end(),
                          back_inserter(CommonElements));

    if (CommonElements.size() == 1)
    {
        // TRACE("One common elemnt pair");
        return mDistanceBetweenElements[CommonElements[0]];
    }
    else
    {
        return 0;
    }
}

template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetCurvatureOfElement(unsigned pottsElementIndex)
{
    assert(SPACE_DIM == 3);
    MathsFunctions M;
    // TRACE(" UPDATED ------  GetVolumeOfElement");
    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);

    double potts_element_ave_curvature = 0.0;

    // An element is made of a number of lattice sites, which are centered around a number of nodes in the Delaunay mesh

    // learn how to iterate ove the container. need the caller in the map below to be a unsigned

    // Im going to adapt this so that I am taking the curvature of the elemen
    std::vector<std::pair<unsigned, unsigned> > VectorOfElementPairs;
    for (unsigned node_index = 0; node_index < p_potts_element->GetNumNodes(); ++node_index)
    {

        unsigned lattice_site_index = p_potts_element->GetNodeGlobalIndex(node_index);

        for (std::set<unsigned>::const_iterator neighbour_lattice_index = this->mMooreNeighbouringNodeIndices[lattice_site_index].begin();
             neighbour_lattice_index != this->mMooreNeighbouringNodeIndices[lattice_site_index].end();
             ++neighbour_lattice_index)
        {
            if (DoNodesShareElement(lattice_site_index, *neighbour_lattice_index))
            {
                // In Neighbourhood and are in the same CPM cell -- need to remove the pair to length
                // check this pair hasnt already been coundted
                std::pair<unsigned, unsigned> NodePair = M.Create_pair(lattice_site_index, *neighbour_lattice_index);
                bool AlreadyCounted = M.IsPairInVector(VectorOfElementPairs, NodePair);
                if (AlreadyCounted == 0)
                {
                    // potts_element_ave_curvature += mAngleBetweenNodes[NodePair];
                    potts_element_ave_curvature += mInstantaneousCurvature[NodePair];
                    VectorOfElementPairs.push_back(NodePair);
                }
            }
        }
    }
    potts_element_ave_curvature /= p_potts_element->GetNumNodes();
    // PRINT_VARIABLE(potts_element_ave_curvature);
    return potts_element_ave_curvature;
}

template <unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetEdges()
{
    MathsFunctions M;

    /*
        Loop over the Mesh nodes and get the conneting edges, label them in a map of connecting edges 
        e.g edge 1 -> pair <Node A, Node B> 
        also have something that gives me if this edge is a boundary 
        -- It will be obviosuy a boudary becuase only one Mesh element will be associated with it, or will be label as such 
    */

    std::vector<std::pair<unsigned, unsigned> > ListOfNodePairs;
    unsigned counter = 0;

    for (typename MutableMesh<2, SPACE_DIM>::NodeIterator node_iter = mpDelaunayMesh->GetNodeIteratorBegin();
         node_iter != mpDelaunayMesh->GetNodeIteratorEnd();
         ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();
        std::set<unsigned> neighbour_lattice_index = this->mMooreNeighbouringNodeIndices[node_index];
        for (std::set<unsigned>::iterator iter = neighbour_lattice_index.begin(); iter != neighbour_lattice_index.end(); ++iter)
        {
            if (node_index != *iter)
            {
                std::pair<unsigned, unsigned> NodePair = std::make_pair(std::min(node_index, *iter), std::max(node_index, *iter));
                if (M.IsPairInVector(ListOfNodePairs, NodePair) == 0)
                {
                    ListOfNodePairs.push_back(NodePair);
                    mEdges[counter] = NodePair;

                    // Get the mesh elements associated with each node
                    std::set<unsigned>& containing_elements1 = node_iter->rGetContainingElementIndices();
                    std::set<unsigned>& containing_elements2 = mpDelaunayMesh->GetNode(*iter)->rGetContainingElementIndices();

                    // Find the Mesh elements that are unique to this edge (there will only be one or two and they will be the intersection)
                    std::set<unsigned> SharedMeshElements;
                    set_intersection(containing_elements2.begin(), containing_elements2.end(), containing_elements1.begin(), containing_elements1.end(),
                                     std::inserter(SharedMeshElements, SharedMeshElements.begin()));

                    std::vector<unsigned> VectorOfContainingElements(SharedMeshElements.begin(), SharedMeshElements.end());

                    mGetElementsFromEdge[counter] = VectorOfContainingElements;
                    // mGetElementsFromNodesPair[NodePair] = VectorOfContainingElements;

                    // Determine if either node is a boundary
                    // unsigned Boundary1 = mpDelaunayMesh->GetNode(node_index)->IsBoundaryNode();
                    // unsigned Boundary2 = mpDelaunayMesh->GetNode(*iter)->IsBoundaryNode();

                    std::vector<unsigned>::iterator Bound_iter = mBoundaryVector.begin();
                    std::vector<unsigned>::iterator Bound_iter2 = mBoundaryVector.begin();

                    std::advance(Bound_iter, node_index);
                    std::advance(Bound_iter2, *iter);
                    unsigned Boundary1 = *Bound_iter;
                    unsigned Boundary2 = *Bound_iter2;
                    mIsEdgeBoundary[counter] = 1000;
                    mIsEdgeBoundaryMap[NodePair] = 1000;

                    // 0 means its in the centeral region and everything is easy
                    // 1 means its on one of the zipping edges
                    // 2 means its on the top or bottom edge
                    // 3 means the edge is at one of the corners and will be a pain

                    if (Boundary1 == 1 && Boundary2 == 1 && VectorOfContainingElements.size() == 2) // THis should be along the edge
                    {
                        // This is a wrapping edge
                        mIsEdgeBoundary[counter] = 1;
                        mIsEdgeBoundaryMap[NodePair] = 1;
                    }
                    else if (Boundary1 == 0 || Boundary2 == 0) // THis should be in the central region
                    {
                        TS_ASSERT(VectorOfContainingElements.size() == 2);
                        mIsEdgeBoundary[counter] = 0;
                        mIsEdgeBoundaryMap[NodePair] = 0;
                    }
                    else if (Boundary1 == 2 && Boundary2 == 2) // this edge is along the top or bottom, should have  VectorOfContainingElements.size() == 1
                    {
                        // TS_ASSERT(VectorOfContainingElements.size() == 1);
                        mIsEdgeBoundary[counter] = 2;
                        mIsEdgeBoundaryMap[NodePair] = 2;
                    }
                    else if ((Boundary1 == 2 && Boundary2 == 3) || (Boundary1 == 3 && Boundary2 == 2)) // this edge is along the top or bottom, should have  VectorOfContainingElements.size() == 1
                    {
                        // TS_ASSERT(VectorOfContainingElements.size() == 1);
                        mIsEdgeBoundary[counter] = 2;
                        mIsEdgeBoundaryMap[NodePair] = 2;
                    }
                    else if ((Boundary1 == 1 && Boundary2 == 3) || (Boundary1 == 3 && Boundary2 == 1)) // One is on the top bottom edge, and the other is along the zip edge
                    {
                        mIsEdgeBoundary[counter] = 4;
                        mIsEdgeBoundaryMap[NodePair] = 4;
                    }
                    else if ((Boundary1 == 1 && Boundary2 == 2 && VectorOfContainingElements.size() == 2) || (Boundary1 == 2 && Boundary2 == 1 && VectorOfContainingElements.size() == 2)) // THis is for the other two nodes on the elements at the edge, this edge has two elements associated with it, so its fine to count as a zero edge
                    {
                        mIsEdgeBoundary[counter] = 0;
                        mIsEdgeBoundaryMap[NodePair] = 0;
                    }
                    else if ((Boundary1 == 1 && Boundary2 == 2 && VectorOfContainingElements.size() == 1) || (Boundary1 == 2 && Boundary2 == 1 && VectorOfContainingElements.size() == 1)) // THis is for the other two nodes on the elements at the edge, this edge has two elements associated with it, so its fine to count as a zero edge
                    {
                        mIsEdgeBoundary[counter] = 3;
                        mIsEdgeBoundaryMap[NodePair] = 3;
                    }
                    else if (Boundary1 == 3 && Boundary2 == 3) // THis is for the other two nodes on the elements at the edge, this edge has two elements associated with it, so its fine to count as a zero edge
                    {
                        mIsEdgeBoundary[counter] = 5;
                        mIsEdgeBoundaryMap[NodePair] = 5;
                    }

                    // PRINT_3_VARIABLES(node_index, *iter, mIsEdgeBoundary[counter] );
                    if (VectorOfContainingElements.size() == 1)
                    {
                        // This edge is on the top or bottom, so I must mark to add extra to the perimeter here

                        mGetNodePairOnBoundary[NodePair] = 1;
                    }
                    else
                    {
                        mGetNodePairOnBoundary[NodePair] = 0;
                    }
                    counter = counter + 1;
                }
            }
        }
    }
}

template <unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::EdgeLenghts()
{
    assert(SPACE_DIM == 3);
    double storedEdgeLongXEdge;
    double storedEdgeShortEdge;
    MathsFunctions M;
    // Calculate the Midpoint for each Mesh element
    for (typename MutableMesh<2, SPACE_DIM>::ElementIterator elem_iter = mpDelaunayMesh->GetElementIteratorBegin(); elem_iter != mpDelaunayMesh->GetElementIteratorEnd(); ++elem_iter)
    {
        unsigned Element_index = elem_iter->GetIndex();
        c_vector<double, 3> CenterPoint = Create_c_vector(0, 0, 0);
        for (int i = 0; i < 3; i++)
        {
            Node<SPACE_DIM>* pNode = mpDelaunayMesh->GetNode(mpDelaunayMesh->GetElement(Element_index)->GetNodeGlobalIndex(i));
            CenterPoint += pNode->rGetLocation();
        }
        mMeshElementMidPoints[Element_index] = CenterPoint / 3;
    }

    double NumberOfEdge = mEdges.size();

    for (unsigned i = 0; i < NumberOfEdge; i++)
    {
        if (mIsEdgeBoundary[i] == 0)
        {
            TS_ASSERT(mGetElementsFromEdge[i].size() == 2);
            Node<SPACE_DIM>* pNode0 = mpDelaunayMesh->GetNode(mEdges[i].first);
            Node<SPACE_DIM>* pNode1 = mpDelaunayMesh->GetNode(mEdges[i].second);

            unsigned Element1 = mGetElementsFromEdge[i][0];
            unsigned Element2 = mGetElementsFromEdge[i][0];
            // Node 0 will be origin, make poisition relative
            c_vector<double, 3> PositionVector = pNode1->rGetLocation() - pNode0->rGetLocation();

            c_vector<double, 3> MidPoint1 = mMeshElementMidPoints[Element1] - pNode0->rGetLocation();
            c_vector<double, 3> MidPoint2 = mMeshElementMidPoints[Element2] - pNode0->rGetLocation();

            c_vector<double, 3> b1 = inner_prod(PositionVector, MidPoint1) / inner_prod(PositionVector, PositionVector) * PositionVector - MidPoint1;
            c_vector<double, 3> b2 = inner_prod(PositionVector, MidPoint2) / inner_prod(PositionVector, PositionVector) * PositionVector - MidPoint2;

            mLengthOfLatticeEdge[mEdges[i]] = norm_2(b1) + norm_2(b2);
            if (PositionVector[1] == 0)
            {
                storedEdgeLongXEdge = norm_2(b1) + norm_2(b2);
            }
            else
            {
                storedEdgeShortEdge = norm_2(b1) + norm_2(b2);
            }
        }
        else if (mIsEdgeBoundary[i] == 2) // THis both of the nodes  are on the top or bottom???
        {
            // PRINT_VARIABLE(mGetElementsFromEdge[i].size() ); // This need to be 2? so the top or the bottom
            // TS_ASSERT(mGetElementsFromEdge[i].size() == 1);

            Node<SPACE_DIM>* pNode0 = mpDelaunayMesh->GetNode(mEdges[i].first);
            Node<SPACE_DIM>* pNode1 = mpDelaunayMesh->GetNode(mEdges[i].second);

            unsigned Element1 = mGetElementsFromEdge[i][0];
            // Node 0 will be origin, make poisition relative
            c_vector<double, 3> PositionVector = pNode1->rGetLocation() - pNode0->rGetLocation();
            mBoundaryEdgeUnitLength = norm_2(PositionVector);
            // TS_ASSERT(mBoundaryEdgeUnitLength ==1); // This assertion not valid for the non periodic mesh
            c_vector<double, 3> MidPoint1 = mMeshElementMidPoints[Element1] - pNode0->rGetLocation();
            c_vector<double, 3> b1 = inner_prod(PositionVector, MidPoint1) / inner_prod(PositionVector, PositionVector) * PositionVector - MidPoint1;
            mLengthOfLatticeEdge[mEdges[i]] = norm_2(b1) + 0.5 * norm_2(PositionVector);
        }
    }
    for (unsigned i = 0; i < NumberOfEdge; i++)
    {
        if (mIsEdgeBoundary[i] == 1) // This will be the edges along the zip
        {
            TS_ASSERT(mGetElementsFromEdge[i].size() == 2); // Need the edge here to be associated with 2 Mesh elements
            Node<SPACE_DIM>* pNode0 = mpDelaunayMesh->GetNode(mEdges[i].first);
            Node<SPACE_DIM>* pNode1 = mpDelaunayMesh->GetNode(mEdges[i].second);
            c_vector<double, 3> PositionVector = pNode1->rGetLocation() - pNode0->rGetLocation();

            if (PositionVector[1] == 0) // No change in y, so just connecting bottom elements or top elements
            {
                mLengthOfLatticeEdge[mEdges[i]] = storedEdgeLongXEdge;
            }
            else
            {
                mLengthOfLatticeEdge[mEdges[i]] = storedEdgeShortEdge;
            }
        }
        else if (mIsEdgeBoundary[i] == 3)
        {
            // TS_ASSERT(mGetElementsFromEdge[i].size() == 1); // Pretty sure these ones should only have 1 element associated
            // mLengthOfLatticeEdge[mEdges[i]] = storedEdgeLongXEdge / 2 + 0.5 * mBoundaryEdgeUnitLength;
            mLengthOfLatticeEdge[mEdges[i]] = storedEdgeShortEdge / 2 + 0.5 * mBoundaryEdgeUnitLength;
            //  PRINT_VARIABLE(storedEdgeShortEdge / 2 + 0.5 * mBoundaryEdgeUnitLength);
            // PRINT_VARIABLE(storedEdgeLongXEdge / 2 + 0.5 * mBoundaryEdgeUnitLength);
        }
        else if (mIsEdgeBoundary[i] == 4)
        {
            mLengthOfLatticeEdge[mEdges[i]] = storedEdgeShortEdge;
        }
        else if (mIsEdgeBoundary[i] == 5)
        {
            mLengthOfLatticeEdge[mEdges[i]] = storedEdgeLongXEdge / 2 + 0.5 * mBoundaryEdgeUnitLength;
        }
        // PRINT_VARIABLE(mLengthOfLatticeEdge[mEdges[i]] )
    }

    //   for (unsigned i = 0; i <NumberOfEdge; i++)
    // {
    //     PRINT_VARIABLE(mLengthOfLatticeEdge[mEdges[i]]);
    // }
}

// Determine graphs connectivity
template <unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::FindElementNeighbours()
{
    TRACE("FindElementNeighbours")
    MathsFunctions M;

    double j = 0;
    assert(SPACE_DIM == 3);
    // Iterate over all of the elements
    for (typename MutableMesh<2, SPACE_DIM>::ElementIterator elem_iter = mpDelaunayMesh->GetElementIteratorBegin();
         elem_iter != mpDelaunayMesh->GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned Element_index = elem_iter->GetIndex();

        // For each element we are getting its 3 nodes and then from these we can find the elements that are contained with these nodes
        // We can then select for the neighbouroing mesh elements that are sharing 2 nodes

        std::set<unsigned> NodesInElement1;

        Node<SPACE_DIM>* pNode0 = mpDelaunayMesh->GetNode(mpDelaunayMesh->GetElement(Element_index)->GetNodeGlobalIndex(0));
        std::set<unsigned>& containing_elements0 = pNode0->rGetContainingElementIndices();
        Node<SPACE_DIM>* pNode1 = mpDelaunayMesh->GetNode(mpDelaunayMesh->GetElement(Element_index)->GetNodeGlobalIndex(1));
        std::set<unsigned>& containing_elements1 = pNode1->rGetContainingElementIndices();
        Node<SPACE_DIM>* pNode2 = mpDelaunayMesh->GetNode(mpDelaunayMesh->GetElement(Element_index)->GetNodeGlobalIndex(2));
        std::set<unsigned>& containing_elements2 = pNode2->rGetContainingElementIndices();

        NodesInElement1.insert(pNode0->GetIndex());
        NodesInElement1.insert(pNode1->GetIndex());
        NodesInElement1.insert(pNode2->GetIndex());

        std::set<unsigned> intersect;
        set_intersection(containing_elements0.begin(), containing_elements0.end(), containing_elements1.begin(), containing_elements1.end(),
                         std::inserter(intersect, intersect.begin()));

        std::set<unsigned> intersect1;
        set_intersection(containing_elements2.begin(), containing_elements2.end(), containing_elements0.begin(), containing_elements0.end(),
                         std::inserter(intersect1, intersect1.begin()));

        std::set<unsigned> intersect2;
        set_intersection(containing_elements2.begin(), containing_elements2.end(), containing_elements1.begin(), containing_elements1.end(),
                         std::inserter(intersect2, intersect2.begin()));

        std::vector<unsigned> VectorOfContainingElements;
        for (std::set<unsigned>::iterator iter = intersect.begin(); iter != intersect.end(); ++iter)
        {
            if (*iter != Element_index)
            {
                VectorOfContainingElements.push_back(*iter);
            }
        }
        for (std::set<unsigned>::iterator iter = intersect1.begin(); iter != intersect1.end(); ++iter)
        {
            if (*iter != Element_index)
            {
                VectorOfContainingElements.push_back(*iter);
            }
        }
        for (std::set<unsigned>::iterator iter = intersect2.begin(); iter != intersect2.end(); ++iter)
        {
            if (*iter != Element_index)
            {
                VectorOfContainingElements.push_back(*iter);
            }
        }

        // Now need to loop over this set and get a the nodes and b the elements i want
        for (std::vector<unsigned>::iterator iter = VectorOfContainingElements.begin(); iter != VectorOfContainingElements.end(); ++iter)
        {
            std::pair<unsigned, unsigned> ElementPair = std::make_pair(std::min(*iter, Element_index), std::max(*iter, Element_index));
            if (M.IsPairInVector(mMeshElementPairs, ElementPair) == 0)
            {
                mMeshElementPairs.push_back(ElementPair);
                std::set<unsigned> NodesInElement2;
                Node<SPACE_DIM>* pNode;
                for (int i = 0; i < 3; i++)
                {
                    pNode = mpDelaunayMesh->GetNode(mpDelaunayMesh->GetElement(*iter)->GetNodeGlobalIndex(i));
                    NodesInElement2.insert(pNode->GetIndex());
                }

                std::set<unsigned> SharedLatticeSites;
                set_intersection(NodesInElement2.begin(), NodesInElement2.end(), NodesInElement1.begin(), NodesInElement1.end(),
                                 std::inserter(SharedLatticeSites, SharedLatticeSites.begin()));

                //

                std::vector<unsigned> IntersectionVector;
                std::set<unsigned>::iterator Set_iter = SharedLatticeSites.begin(); // get iterator to 1st element
                IntersectionVector.push_back(*Set_iter);
                std::advance(Set_iter, 1); // advance by 1
                IntersectionVector.push_back(*Set_iter);
                sort(IntersectionVector.begin(), IntersectionVector.end());

                std::pair<unsigned, unsigned> AssociatedNodes = M.Create_pair(IntersectionVector[0], IntersectionVector[1]);

                assert(SharedLatticeSites.size() == 2);

                mMapElementPairsToNodes[ElementPair] = IntersectionVector;
                mMapNodesPairsToElementPairs[AssociatedNodes] = ElementPair;
                // mMeshElementNeighbours[Element_index].push_back(Neighbour_Element_index);

                // Also save the paired elements for each of the nodes
                mMapNodesToAssociateElementPairs[IntersectionVector[0]].push_back(ElementPair);
                mMapNodesToAssociateElementPairs[IntersectionVector[1]].push_back(ElementPair);
            }
        }
    }

    // Have now saved all of the mesh element neighbours into member vectors and maps

    // ------------------------------

    // Check that the neighbour element pairs associated with each node is the same as the number of elements contained at the node

    for (typename MutableMesh<2, SPACE_DIM>::NodeIterator node_iter = mpDelaunayMesh->GetNodeIteratorBegin();
         node_iter != mpDelaunayMesh->GetNodeIteratorEnd();
         ++node_iter)
    {

        unsigned node_index = node_iter->GetIndex();
        std::vector<std::pair<unsigned, unsigned> > VectorOfElementPairs = mMapNodesToAssociateElementPairs[node_index];

        std::vector<unsigned> ElementsFromPairs;

        for (std::vector<std::pair<unsigned, unsigned> >::iterator it = VectorOfElementPairs.begin(); it != VectorOfElementPairs.end(); ++it)
        {

            ElementsFromPairs.push_back(it->first);
            ElementsFromPairs.push_back(it->second);
            // PRINT_PAIR(*it);
        }

        std::set<unsigned>& containing_elements = node_iter->rGetContainingElementIndices();
        // Vector of the contained elements
        std::vector<unsigned> ContainingElements(containing_elements.begin(), containing_elements.end());
        sort(ContainingElements.begin(), ContainingElements.end());
        sort(ElementsFromPairs.begin(), ElementsFromPairs.end());
        ElementsFromPairs.erase(unique(ElementsFromPairs.begin(), ElementsFromPairs.end()), ElementsFromPairs.end());
        // Check if the contained elements and the Assocaited element pairs are the same and throw an error if they are not
        // The Asscoiated vector pairs will be empty of there is only one element
        if (M.AreVectorsSame(ContainingElements, ElementsFromPairs) != 0 && ElementsFromPairs.size() != 0)
        {
            EXCEPTION("Elements in element pairs are not the same as the elements contained at this node");
        }
    }
    // TRACE("Have neighbours");
}

template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetPerimeterOfElement(unsigned pottsElementIndex) // perimiter
{
    assert(SPACE_DIM == 3);
    /*  1) Loop over the lattice sites and see which ones have neighbours that are in the Sister element 
        2) Should only be counted if they are around the center.  --- only one seam     
    */
    MathsFunctions M;
    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);
    double potts_element_surface = 0.0;

    //  Loop over each of the lattice sites in the cell

    for (unsigned local_lattice_index = 0; local_lattice_index < p_potts_element->GetNumNodes(); ++local_lattice_index)
    {
        unsigned lattice_site_index = p_potts_element->GetNodeGlobalIndex(local_lattice_index);

        // iterate over all neighbours for this element to check if they are in the same element
        for (std::set<unsigned>::const_iterator neighbour_lattice_index = this->mMooreNeighbouringNodeIndices[lattice_site_index].begin();
             neighbour_lattice_index != this->mMooreNeighbouringNodeIndices[lattice_site_index].end();
             ++neighbour_lattice_index)
        {
            std::pair<unsigned, unsigned> NodePair = M.Create_pair(lattice_site_index, *neighbour_lattice_index);
            if (DoNodesShareElement(lattice_site_index, *neighbour_lattice_index) == 0)
            {
                potts_element_surface += mLengthOfLatticeEdge[NodePair];
                TS_ASSERT(mLengthOfLatticeEdge[NodePair] != 0)
            }
            else if (mIsEdgeBoundaryMap[NodePair] == 2)
            {
                // Node pair is on top or bottom and needs to be counted and edges
                potts_element_surface += mBoundaryEdgeUnitLength / 2; // Divide by two because each length is going to get counted twice in iterating around each lattice's neighbourhood
                double BoundL = mBoundaryEdgeUnitLength / 2;
            }
        }
    }

    return potts_element_surface;
}

template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetPerimeterOfCoupledElements(unsigned pottsElementIndex, unsigned pottsSisterElementIndex, double Center) // perimiter
{

    /*  
        1) Loop over the lattice sites and see which ones have neighbours that are in the Sister element 
        2) Should only be counted if they are around the center.  --- only one seam     
    */

    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);
    PottsElement<SPACE_DIM>* p_potts_Sister_element = this->GetElement(pottsSisterElementIndex);

    double potts_element_surface = 0.0;

    unsigned EdgesAddedA = 0;
    unsigned EdgesAddedB = 0;

    MathsFunctions M;
    //  Loop over each of the lattice sites in the cell
    std::vector<std::pair<unsigned, unsigned> > EdgeLattices;
    for (unsigned local_lattice_index = 0; local_lattice_index < p_potts_element->GetNumNodes(); ++local_lattice_index)
    {
        unsigned lattice_site_index = p_potts_element->GetNodeGlobalIndex(local_lattice_index);

        // iterate over all neighbours for this element to check if they are in the same element
        for (std::set<unsigned>::const_iterator neighbour_lattice_index = this->mMooreNeighbouringNodeIndices[lattice_site_index].begin();
             neighbour_lattice_index != this->mMooreNeighbouringNodeIndices[lattice_site_index].end();
             ++neighbour_lattice_index)
        {
            std::pair<unsigned, unsigned> NodePair = M.Create_pair(lattice_site_index, *neighbour_lattice_index);
            if (DoNodesShareElement(lattice_site_index, *neighbour_lattice_index) == 0)
            {
                potts_element_surface += mLengthOfLatticeEdge[NodePair];
                if (DoNodesShareElement(p_potts_Sister_element->GetNodeGlobalIndex(0), *neighbour_lattice_index) == 1)
                {
                    // If the neighbour is in the sister elememt, then I want to save it in this Edge lattice vector -- then I can Loop through this vector, and take the edges of the
                    // ones I dont needs
                    EdgeLattices.push_back(NodePair);
                }

                TS_ASSERT(mLengthOfLatticeEdge[NodePair] != 0)
            }
            else if (mIsEdgeBoundaryMap[NodePair] == 2)
            {
                // Node pair is on top or bottom and needs to be counted
                potts_element_surface += mBoundaryEdgeUnitLength / 2; // Divide by two because each length is going to get counted twice in iterating around each lattice's neighbourhood
                double BoundL = mBoundaryEdgeUnitLength / 2;
            }
        }
    }
    // EdgeLattices.erase(unique(EdgeLattices.begin(), EdgeLattices.end()), EdgeLattices.end());
    // Now do the same for the sister

    for (unsigned local_sister_lattice_index = 0; local_sister_lattice_index < p_potts_Sister_element->GetNumNodes(); ++local_sister_lattice_index)
    {

        unsigned lattice_site_index = p_potts_Sister_element->GetNodeGlobalIndex(local_sister_lattice_index);

        // iterate over all neighbours for this element to check if they are in the same element
        for (std::set<unsigned>::const_iterator neighbour_lattice_index = this->mMooreNeighbouringNodeIndices[lattice_site_index].begin();
             neighbour_lattice_index != this->mMooreNeighbouringNodeIndices[lattice_site_index].end();
             ++neighbour_lattice_index)
        {
            std::pair<unsigned, unsigned> NodePair = M.Create_pair(lattice_site_index, *neighbour_lattice_index);
            if (DoNodesShareElement(lattice_site_index, *neighbour_lattice_index) == 0)
            {
                potts_element_surface += mLengthOfLatticeEdge[NodePair];
                TS_ASSERT(mLengthOfLatticeEdge[NodePair] != 0)
            }
            else if (mIsEdgeBoundaryMap[NodePair] == 2)
            {
                // Node pair is on top or bottom and needs to be counted
                potts_element_surface += mBoundaryEdgeUnitLength / 2; // Divide by two because each length is going to get counted twice in iterating around each lattice's neighbourhood
                double BoundL = mBoundaryEdgeUnitLength / 2;

            }
        }
    }

    /*
        Now I need to lop over everything and get the central region out.  How to do this --- I have the central part. I should loop over this and look for the lattices in Element 1 within 
        Within a region close to the centeral element will be defined as  UnitStep = Width/(N_D-1);
    */

    double UnitStep = mWidth / (2*mN_C - 2);
    for (std::vector<std::pair<unsigned, unsigned> >::iterator NodePair = EdgeLattices.begin(); NodePair != EdgeLattices.end(); ++NodePair)
    {
        
        Node<SPACE_DIM>* p_node1 = mpDelaunayMesh->GetNode(NodePair->first);
        Node<SPACE_DIM>* p_node2 = mpDelaunayMesh->GetNode(NodePair->second);
        c_vector<double, SPACE_DIM> LatticeLocation = p_node1->rGetLocation();
        c_vector<double, SPACE_DIM> NeighbourLocation = p_node2->rGetLocation();

        
        // double Distance1 =std::abs(LatticeLocation[0] - Center);
        // double Distance2 =std::abs(NeighbourLocation[0] - Center);

        if (std::abs(LatticeLocation[0] - Center) < 5*UnitStep || std::abs(NeighbourLocation[0] - Center) <  5*UnitStep) // We also need to consider the zip edge here!!!!!
        {   //Might want to put this back to 2* UnitStep 

            potts_element_surface -= (2 * mLengthOfLatticeEdge[*NodePair]);

            // TRACE("At midline - Must remove from perimeter")
            // M.PRINT_PAIR(*NodePair);
            // PRINT_3_VARIABLES(Distance1, Distance2, UnitStep)
            // TRACE("...")
            if (mIsEdgeBoundaryMap[*NodePair] == 2 || mIsEdgeBoundaryMap[*NodePair] == 5) // Is on top or bottom?
            {
                // This pair is on the edge, i need to take out the perimeter dividing them, but leave the edge length along the top or bottom of the mesh
                potts_element_surface += mBoundaryEdgeUnitLength; // Divide by two because each length is going to get counted twice in iterating around each lattice's neighbourhood
            }
        }
        // else
        // {
        //     TRACE("At Zip")
        //     M.PRINT_PAIR(*NodePair);
        //     PRINT_3_VARIABLES(Distance1, Distance2, UnitStep)
        //     TRACE("...")
        // }
    }
    return potts_element_surface;
}

template <unsigned SPACE_DIM>
bool PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::IsNodeNearCenter(unsigned LatticeIndex, double Center) // perimiter
{

    // Need to write a function giving ditance to center
    Node<SPACE_DIM>* pNode = mpDelaunayMesh->GetNode(LatticeIndex);
    double DistanceToCenter = std::abs(Center - pNode->rGetLocation()[0]);

    if (DistanceToCenter > 4.5 * mLatticeSpaceing && DistanceToCenter != NAN)
    {
        return false;
        // node is not nearu cemter line this should happen if the cell has wrapped around on itself.
    }
    else
    {
        // Node is near center line
        return true;
    }
}

// I am going to need this funciton to split the cell in half
template <unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetCellCentre(std::vector<unsigned> Edges, unsigned pottsElementIndex)
{
    assert(SPACE_DIM == 3);
    // MathsFunctions M;

    c_vector<double, SPACE_DIM> Mean_Location2 = zero_vector<double>(SPACE_DIM);
    std::vector<c_vector<double, SPACE_DIM> > Locations;

    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);

    for (std::vector<unsigned>::iterator it = Edges.begin(); it != Edges.end(); ++it)
    {
        // std::vector<double> LocationOfLatticeX = p_potts_element->GetNodeLocation(*it);

        c_vector<double, SPACE_DIM> LocationOfLatticeX = p_potts_element->GetNodeLocation(*it);
        LocationOfLatticeX[2] = 0;
        if (IsVectorInVector(Locations, LocationOfLatticeX) == 0)
        {
            Locations.push_back(LocationOfLatticeX);
        }
        Mean_Location2 += p_potts_element->GetNodeLocation(*it);
    }
    // Locations.erase(unique(Locations.begin(), Locations.end()), Locations.end());
    c_vector<double, 3> Mean_Location = zero_vector<double>(SPACE_DIM);
    for (typename std::vector<c_vector<double, SPACE_DIM> >::iterator it = Locations.begin(); it != Locations.end(); ++it)
    {

        Mean_Location += *it;
        // p_potts_element->GetNodeLocation(node_index);
    }
    Mean_Location /= Locations.size();
    Mean_Location2 /= Locations.size();
    return Mean_Location2;
}

template <unsigned SPACE_DIM>
std::map<unsigned, c_vector<double, 2> > PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::MapCellTo2DPlane(unsigned pottsElementIndex)
{

    assert(SPACE_DIM == 3);

    std::map<unsigned, c_vector<double, 2> > MappedCellToPlane;

    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);
    double potts_element_volume = 0.0;

    for (unsigned node_index = 0; node_index < p_potts_element->GetNumNodes(); ++node_index)
    {
        c_vector<double, 3> Location = p_potts_element->GetNodeLocation(node_index);
        unsigned lattice_site_index = p_potts_element->GetNodeGlobalIndex(node_index);
        double X_Pos = Location[0];
        double Y_Pos = Location[1];
        double Z_Pos = Location[2];

        double Theta = atan(Y_Pos / X_Pos);

        if (X_Pos < 0 && Y_Pos > 0)
        {
            Theta = Theta + M_PI;
        }
        else if (X_Pos < 0 && Y_Pos < 0)
        {
            Theta = Theta - M_PI;
        }
        else if (Y_Pos == 0)
        {
            if (X_Pos > 0)
            {
                Theta = 0;
            }
            else if (X_Pos < 0)
            {
                Theta = M_PI;
            }
        }
        double C = mRadius * Theta;
        MappedCellToPlane[lattice_site_index] = Create_c_vector(C, Z_Pos);
    }
    return MappedCellToPlane;
}

template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::CellLength(unsigned pottsElementIndex)
{

    assert(SPACE_DIM == 3);

    double X_max;
    double X_min;

    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);
    double potts_element_volume = 0.0;

    for (unsigned node_index = 0; node_index < p_potts_element->GetNumNodes(); ++node_index)
    {
        c_vector<double, 3> Location = p_potts_element->GetNodeLocation(node_index);
        double X_Pos = Location[0];
        double Y_Pos = Location[1];

        double Theta = atan(Y_Pos / X_Pos);

        if (X_Pos < 0 && Y_Pos > 0)
        {
            Theta = Theta + M_PI;
        }
        else if (X_Pos < 0 && Y_Pos < 0)
        {
            Theta = Theta - M_PI;
        }
        else if (Y_Pos == 0)
        {
            if (X_Pos > 0)
            {
                Theta = 0;
            }
            else if (X_Pos < 0)
            {
                Theta = M_PI;
            }
        }
        double C = mRadius * Theta;
        if (node_index == 0)
        {
            X_max = C;
            X_min = C;
        }

        if (C > X_max)
        {
            X_max = C;
        }
        else if (C < X_min)
        {
            X_min = C;
        }
    }

    double length = X_max - X_min;

    return length;
}

template <unsigned SPACE_DIM>
Node<SPACE_DIM>* PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetPottsLatticeSite(unsigned latticeSiteIndex)
{
    Node<SPACE_DIM>* p_node = mpDelaunayMesh->GetNode(latticeSiteIndex);
    return p_node;
}

template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetVolumeOfLatticeSite(unsigned NodeIndex)
{

    return mLatticeVolume[NodeIndex];
}

template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetNumberOfLatticesInElement(unsigned pottsElementIndex)
{
    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);
    return p_potts_element->GetNumNodes();
}

template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetVolumeOfElement(unsigned pottsElementIndex)
{
    // TRACE(" UPDATED ------  GetVolumeOfElement");
    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);

    double potts_element_volume = 0.0;

    // An element is made of a number of lattice sites, which are centered around a number of nodes in the Delaunay mesh

    // learn how to iterate ove the container. need the caller in the map below to be a unsigned
    for (unsigned node_index = 0; node_index < p_potts_element->GetNumNodes(); ++node_index)
    {

        unsigned latticeSiteIndex = p_potts_element->GetNodeGlobalIndex(node_index);
        potts_element_volume += mLatticeVolume[p_potts_element->GetNodeGlobalIndex(node_index)];
    }

    return potts_element_volume;
}

template <unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::CalculateLatticeVolumes()
{
    //  TRACE("Jess is the best");
    MathsFunctions M;

    assert(SPACE_DIM == 3);
    double AverageArea = 0;

    // Have iterated over all of the nodes and saved the area of the lattice site,
    for (typename MutableMesh<2, SPACE_DIM>::NodeIterator node_iter = mpDelaunayMesh->GetNodeIteratorBegin();
         node_iter != mpDelaunayMesh->GetNodeIteratorEnd();
         ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();
        std::set<unsigned>& containing_elements = node_iter->rGetContainingElementIndices();
        assert(containing_elements.size() > 0);

        std::vector<unsigned>::iterator Bound_iter = mBoundaryVector.begin();
        std::advance(Bound_iter, node_index);

        double LatticeArea = 0;
        c_vector<double, 3> Normal = Create_c_vector(0, 0, 0);
        if (*Bound_iter == 1)
        {
            std::vector<unsigned> containing_elements_Vector(containing_elements.size());
            std::copy(containing_elements.begin(), containing_elements.end(), containing_elements_Vector.begin());

            std::vector<unsigned> SideA;
            std::vector<unsigned> SideB;
            // I have all the neigbours, need one that is on the same side of the mesh, to do this, I only need to find the ones with a boundary condition of 0
            std::set<unsigned> neighbour_lattice_index = this->mMooreNeighbouringNodeIndices[node_index];
            // std::vector<unsigned> NodesAssociatedWith_node_iter;
            for (std::set<unsigned>::iterator iter = neighbour_lattice_index.begin(); iter != neighbour_lattice_index.end(); ++iter)
            {
                // Find the element that is common to both
                Node<SPACE_DIM>* pNeighbourNode = mpDelaunayMesh->GetNode(*iter);
                std::set<unsigned>& NeighbourELements = pNeighbourNode->rGetContainingElementIndices();

                std::vector<unsigned> NeighbourELementsVector(NeighbourELements.size());
                std::copy(NeighbourELements.begin(), NeighbourELements.end(), NeighbourELementsVector.begin());

                c_vector<double, 3> Direction = pNeighbourNode->rGetLocation() - node_iter->rGetLocation();

                if (norm_2(Direction) < 4.5)
                {
                    // We found a neighbour on this side of the mesh iterate over the common elements
                    std::set_intersection(containing_elements_Vector.begin(), containing_elements_Vector.end(),
                                          NeighbourELementsVector.begin(), NeighbourELementsVector.end(),
                                          back_inserter(SideA));
                }
                else
                {
                    std::set_intersection(containing_elements_Vector.begin(), containing_elements_Vector.end(),
                                          NeighbourELementsVector.begin(), NeighbourELementsVector.end(),
                                          back_inserter(SideB));
                }
            }
            sort(SideB.begin(), SideB.end());
            SideB.erase(unique(SideB.begin(), SideB.end()), SideB.end());
            sort(SideA.begin(), SideA.end());
            SideA.erase(unique(SideA.begin(), SideA.end()), SideA.end());

            SideA.erase(std::set_difference(SideA.begin(), SideA.end(), SideB.begin(), SideB.end(), SideA.begin()), SideA.end());

            double counter = 0;
            for (std::vector<unsigned>::iterator iter = SideA.begin(); iter != SideA.end(); ++iter)
            {
                counter += 1;
                unsigned elem_index = *iter;
                // PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(elem_index);

                c_vector<c_vector<double, 3>, 2> Vectors;
                std::vector<unsigned> LocalNodes;
                double j = 0;

                // Loops over the nodes in the element to get the two position vectors describing the element form the inital node
                for (int i = 0; i < 3; i++)
                {
                    Node<SPACE_DIM>* pNode = mpDelaunayMesh->GetNode(mpDelaunayMesh->GetElement(*iter)->GetNodeGlobalIndex(i));
                    unsigned node_index_i = pNode->GetIndex();
                    if (node_index_i != node_index)
                    {
                        LocalNodes.push_back(node_index_i);
                        Vectors[j] = pNode->rGetLocation() - node_iter->rGetLocation();
                        j += 1;
                    }
                }
                LatticeArea += 0.5 * norm_2(VectorProduct(Vectors[0], Vectors[1]));
                Normal += VectorProduct(Vectors[0], Vectors[1]);
                if (counter == 2)
                {
                    break;
                }
            }

            Normal /= norm_2(Normal);
        }
        else if (*Bound_iter == 3)
        {
            LatticeArea = (mWidth / mN_C * mLength / mN_D) / 2; //0.5;
        }
        else
        {
            // Loop over the elements for this node, to find the perimiter of the lines connecting to the midpoint
            for (std::set<unsigned>::iterator iter = containing_elements.begin();
                 iter != containing_elements.end();
                 ++iter)
            {
                unsigned elem_index = *iter;
                c_vector<c_vector<double, 3>, 2> Vectors;
                std::vector<unsigned> LocalNodes;
                double j = 0;

                // Loops over the nodes in the element to get the two position vectors describing the element form the inital node
                for (int i = 0; i < 3; i++)
                {
                    Node<SPACE_DIM>* pNode = mpDelaunayMesh->GetNode(mpDelaunayMesh->GetElement(*iter)->GetNodeGlobalIndex(i));
                    unsigned node_index_i = pNode->GetIndex();
                    if (node_index_i != node_index)
                    {
                        LocalNodes.push_back(node_index_i);
                        Vectors[j] = pNode->rGetLocation() - node_iter->rGetLocation();
                        j += 1;
                    }
                }

                c_vector<double, 3> MidPoint = mMeshElementMidPoints[elem_index] - node_iter->rGetLocation();

                c_vector<double, 3> a1 = inner_prod(Vectors[0], MidPoint) / inner_prod(Vectors[0], Vectors[0]) * Vectors[0];
                c_vector<double, 3> a2 = inner_prod(Vectors[1], MidPoint) / inner_prod(Vectors[1], Vectors[1]) * Vectors[1];

                double Area1 = 0.5 * norm_2(VectorProduct(a1, MidPoint));
                double Area2 = 0.5 * norm_2(VectorProduct(a2, MidPoint));
                LatticeArea += Area1 + Area2;

                c_vector<double, 3> MiniElementNormal = VectorProduct(a1, a2);
                MiniElementNormal *= M.MaintainOutwardsPointingNormal(MiniElementNormal, node_iter->rGetLocation());
                MiniElementNormal /= norm_2(MiniElementNormal);
                Normal += MiniElementNormal;
            }

            Normal /= norm_2(Normal);
        }
        mLatticeVolume[node_index] = LatticeArea;
        AverageArea += LatticeArea;
        // PRINT_VARIABLE(LatticeArea);
        mLatticeNormal[node_index] = Normal;
    }
    mAverageArea = AverageArea / mpDelaunayMesh->GetNumNodes();
}

template <unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetRadius(double Radius)
{
    mRadius = Radius;
}

template <unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::MapCylinderToPlane()
{
    assert(SPACE_DIM == 3);
    double MaxC = 0;
    double MinC = 0;
    double MinTheta = 0;
    double MaxTheta = 0;
    // double Radius =1.5e-3;

    for (typename MutableMesh<2, SPACE_DIM>::NodeIterator node_iter = mpDelaunayMesh->GetNodeIteratorBegin();
         node_iter != mpDelaunayMesh->GetNodeIteratorEnd();
         ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();
        c_vector<double, 3> Location = node_iter->rGetLocation();

        double X_Pos = node_iter->rGetLocation()[0];
        double Y_Pos = node_iter->rGetLocation()[1];
        double Z_Pos = node_iter->rGetLocation()[2];

        // Mapped -> Circumferencial lenght

        double Theta = atan(Y_Pos / X_Pos);

        if (X_Pos < 0 && Y_Pos > 0)
        {
            Theta = Theta + M_PI;
        }
        else if (X_Pos < 0 && Y_Pos < 0)
        {
            Theta = Theta - M_PI;
        }
        else if (Y_Pos == 0)
        {
            if (X_Pos > 0)
            {
                Theta = 0;
            }
            else if (X_Pos < 0)
            {
                Theta = M_PI;
            }
        }
        double C = mRadius * Theta;
        if (C > MaxC)
        {
            MaxC = C;
        }
        else if (C < MinC)
        {
            MinC = C;
        }
        if (Theta < MinTheta)
        {
            MinTheta = Theta;
        }
        else if (Theta > MaxTheta)
        {
            MaxTheta = Theta;
        }

        mMappedLocation[node_index] = Create_c_vector(C, Z_Pos);
    }
    mMaxC = MaxC;
    mMinC = MinC;
    mLength = MaxC - MinC;
    //  PRINT_3_VARIABLES(MinTheta, MaxTheta,  mLength );
}

template <unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetCellCentre(unsigned pottsElementIndex)
{

    //  assert(SPACE_DIM == 3);
    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);
    c_vector<double, SPACE_DIM> Mean_Location = zero_vector<double>(SPACE_DIM);

    for (unsigned i = 0; i < p_potts_element->GetNumNodes(); i++)
    {
        Mean_Location += p_potts_element->GetNodeLocation(i);
    }
    Mean_Location /= p_potts_element->GetNumNodes();

    return Mean_Location;
}

// Gives aspect ratio of cell on non-periodic mesh
template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetAspectRatio(unsigned pottsElementIndex)
{

    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);
    assert(p_potts_element->GetNumNodes() != 0);

    if (p_potts_element->GetNumNodes() <= 2)
    {
        return 1.0;
    }

    double eig_max;
    double eig_min;

    // See http://stackoverflow.com/questions/7059841/estimating-aspect-ratio-of-a-convex-hull for how to do it.

    // Calculate entries of covariance matrix (var_x,cov_xy;cov_xy,var_y)
    double mean_x = 0;
    double mean_y = 0;

    for (unsigned i = 0; i < p_potts_element->GetNumNodes(); i++)
    {

        Node<SPACE_DIM>* pNode = p_potts_element->GetNode(i);
        unsigned node_index = pNode->GetIndex();
        // c_vector<double, 2> MappedPosition =  mMappedLocation[node_index];

        double X_Pos = p_potts_element->GetNode(i)->rGetLocation()[0]; // mMappedLocation[node_index][0]; // ;
        double Y_Pos = p_potts_element->GetNode(i)->rGetLocation()[1]; // mMappedLocation[node_index][1]; //

        mean_x += X_Pos;
        mean_y += Y_Pos;
    }

    mean_x /= (p_potts_element->GetNumNodes());
    mean_y /= (p_potts_element->GetNumNodes());

    double variance_x = 0;
    double variance_y = 0;
    double covariance_xy = 0;

    for (unsigned i = 0; i < p_potts_element->GetNumNodes(); i++)
    {
        Node<SPACE_DIM>* pNode = p_potts_element->GetNode(i);
        unsigned node_index = pNode->GetIndex();
        // c_vector<double, 2> MappedPosition =  mMappedLocation[node_index];

        double X_Pos = p_potts_element->GetNode(i)->rGetLocation()[0]; // mMappedLocation[node_index][0]; //
        double Y_Pos = p_potts_element->GetNode(i)->rGetLocation()[1]; //mMappedLocation[node_index][1]; //

        variance_x += pow((X_Pos - mean_x), 2);
        variance_y += pow((Y_Pos - mean_y), 2);
        covariance_xy += (X_Pos - mean_x) * (Y_Pos - mean_y);
    }

    variance_x /= (p_potts_element->GetNumNodes());
    variance_y /= (p_potts_element->GetNumNodes());
    covariance_xy /= (p_potts_element->GetNumNodes());

    // Calculate max/min eigenvalues
    double trace = variance_x + variance_y;
    double det = variance_x * variance_y - covariance_xy * covariance_xy;

    eig_max = 0.5 * (trace + sqrt(trace * trace - 4 * det));
    eig_min = 0.5 * (trace - sqrt(trace * trace - 4 * det));

    if (eig_min == 0)
    {
        // TRACE("All nodes in an element lie in the same line/plane (2D/3D) so aspect ratio is infinite. This interferes with calculation of the Hamiltonian.");
        return -1;
        // EXCEPTION("All nodes in an element lie in the same line/plane (2D/3D) so aspect ratio is infinite. This interferes with calculation of the Hamiltonian.");
    }

    // TO get the major axis will be the eigenvector corresponding to to the major eigenvaule eig_max -- basic eigenvector calcuations

    // double MajorAxisX = (variance_x - eig_max) * (variance_y - eig_max) / (covariance_xy * covariance_xy);
    // double MajorAxisY = -(variance_x - eig_max) * (variance_x - eig_max) * (variance_y - eig_max) / (covariance_xy * covariance_xy * covariance_xy);

    c_vector<double, 2> MajorAxis = Create_c_vector((eig_max - variance_y) / covariance_xy, 1);

    MajorAxis /= norm_2(MajorAxis);
    if (covariance_xy == 0)
    {
        MajorAxis = Create_c_vector(0, 0);
    }
    mMajorAxis[pottsElementIndex] = MajorAxis; //Create_c_vector(MajorAxisX, MajorAxisY);
    // TS_ASSERT(!isnan(mMajorAxis[0]));
    return eig_max / eig_min;
}

// Gives aspect ratio of cell on wrapped mesh
template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetAspectRatio(unsigned pottsElementIndex, unsigned pottsSisterIndex)
{

    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);
    PottsElement<SPACE_DIM>* p_potts_sister = this->GetElement(pottsSisterIndex);
    assert(p_potts_element->GetNumNodes() != 0);
    assert(p_potts_sister->GetNumNodes() != 0);

    if (p_potts_element->GetNumNodes() + p_potts_sister->GetNumNodes() <= 2)
    {
        return 1.0;
    }

    double eig_max;
    double eig_min;

    // See http://stackoverflow.com/questions/7059841/estimating-aspect-ratio-of-a-convex-hull for how to do it.

    // Calculate entries of covariance matrix (var_x,cov_xy;cov_xy,var_y)
    double mean_x = 0;
    double mean_y = 0;

    for (unsigned i = 0; i < p_potts_element->GetNumNodes(); i++)
    {

        Node<SPACE_DIM>* pNode = p_potts_element->GetNode(i);
        unsigned node_index = pNode->GetIndex();
        // c_vector<double, 2> MappedPosition =  mMappedLocation[node_index];

        double X_Pos = p_potts_element->GetNode(i)->rGetLocation()[0]; // mMappedLocation[node_index][0]; // ;
        double Y_Pos = p_potts_element->GetNode(i)->rGetLocation()[1]; // mMappedLocation[node_index][1]; //

        mean_x += X_Pos;
        mean_y += Y_Pos;
    }

    for (unsigned i = 0; i < p_potts_sister->GetNumNodes(); i++)
    {

        Node<SPACE_DIM>* pNode = p_potts_sister->GetNode(i);
        unsigned node_index = pNode->GetIndex();
        // c_vector<double, 2> MappedPosition =  mMappedLocation[node_index];

        double X_Pos = p_potts_sister->GetNode(i)->rGetLocation()[0]; // mMappedLocation[node_index][0]; // ;
        double Y_Pos = p_potts_sister->GetNode(i)->rGetLocation()[1]; // mMappedLocation[node_index][1]; //

        mean_x += X_Pos;
        mean_y += Y_Pos;
    }

    mean_x /= (p_potts_element->GetNumNodes() + p_potts_sister->GetNumNodes());
    mean_y /= (p_potts_element->GetNumNodes() + p_potts_sister->GetNumNodes());

    double variance_x = 0;
    double variance_y = 0;
    double covariance_xy = 0;

    for (unsigned i = 0; i < p_potts_element->GetNumNodes(); i++)
    {
        Node<SPACE_DIM>* pNode = p_potts_element->GetNode(i);
        unsigned node_index = pNode->GetIndex();
        // c_vector<double, 2> MappedPosition =  mMappedLocation[node_index];

        double X_Pos = p_potts_element->GetNode(i)->rGetLocation()[0]; // mMappedLocation[node_index][0]; //
        double Y_Pos = p_potts_element->GetNode(i)->rGetLocation()[1]; //mMappedLocation[node_index][1]; //

        variance_x += pow((X_Pos - mean_x), 2);
        variance_y += pow((Y_Pos - mean_y), 2);
        covariance_xy += (X_Pos - mean_x) * (Y_Pos - mean_y);
    }

    for (unsigned i = 0; i < p_potts_sister->GetNumNodes(); i++)
    {
        Node<SPACE_DIM>* pNode = p_potts_sister->GetNode(i);
        unsigned node_index = pNode->GetIndex();
        // c_vector<double, 2> MappedPosition =  mMappedLocation[node_index];

        double X_Pos = p_potts_sister->GetNode(i)->rGetLocation()[0]; // mMappedLocation[node_index][0]; //
        double Y_Pos = p_potts_sister->GetNode(i)->rGetLocation()[1]; //mMappedLocation[node_index][1]; //

        variance_x += pow((X_Pos - mean_x), 2);
        variance_y += pow((Y_Pos - mean_y), 2);
        covariance_xy += (X_Pos - mean_x) * (Y_Pos - mean_y);
    }
    // PRINT_VARIABLE(covariance_xy)

    variance_x /= (p_potts_element->GetNumNodes() + p_potts_sister->GetNumNodes());
    variance_y /= (p_potts_element->GetNumNodes() + p_potts_sister->GetNumNodes());
    covariance_xy /= (p_potts_element->GetNumNodes() + p_potts_sister->GetNumNodes());

    // Calculate max/min eigenvalues
    double trace = variance_x + variance_y;
    double det = variance_x * variance_y - covariance_xy * covariance_xy;

    eig_max = 0.5 * (trace + sqrt(trace * trace - 4 * det));
    eig_min = 0.5 * (trace - sqrt(trace * trace - 4 * det));

    if (eig_min == 0)
    {
        // TRACE("All nodes in an element lie in the same line/plane (2D/3D) so aspect ratio is infinite. This interferes with calculation of the Hamiltonian.");
        return -1;
        // EXCEPTION("All nodes in an element lie in the same line/plane (2D/3D) so aspect ratio is infinite. This interferes with calculation of the Hamiltonian.");
    }
    // double MajorAxis_1 = sqrt( (variance_y+ eig_max)*(variance_y+ eig_max) / (covariance_xy*covariance_xy + (variance_y+ eig_max)*(variance_y+ eig_max)));
    // double MajorAxis_2 = sqrt( 1- MajorAxis_1 * MajorAxis_1 );
    c_vector<double, 2> MajorAxis = Create_c_vector((eig_max - variance_y) / covariance_xy, 1);
    // c_vector<double, 2> MajorAxis = Create_c_vector(MajorAxis_1, MajorAxis_2);
    // PRINT_VECTOR(MajorAxis)
    if (covariance_xy == 0)
    {
        // The major axis lies on the  or y axis  https://cookierobotics.com/007/
        // Can check which one
        if (variance_x >= variance_y)
        {
            //    theta =0;
            MajorAxis = Create_c_vector(1, 0);
        }
        else if (variance_x <= variance_y)
        {
            //    theta = M_PI/2;
            MajorAxis = Create_c_vector(0, 1);
        }
    }
    MajorAxis /= norm_2(MajorAxis);

    // doubletheta= asin( ( (variance_y+ eig_max)*(variance_y+ eig_max) / (covariance_xy*covariance_xy + (variance_y+ eig_max)*(variance_y+ eig_max)))/ (1- MajorAxis_1 * MajorAxis_1));

    double Theta = atan(MajorAxis[0] / MajorAxis[1]);
    double Theta_check = atan((eig_max - variance_x) / covariance_xy);

    double Theta_check5 = Theta_check;
    if (!(std::abs(std::abs(Theta) - std::abs(Theta_check)) < 0.15 || std::abs(std::abs(Theta) - std::abs(Theta_check - M_PI / 2)) < 0.15 || std::abs(std::abs(Theta) - std::abs(Theta_check - M_PI)) < 0.15 || std::abs(std::abs(Theta) - std::abs(Theta_check + M_PI / 2)) < 0.15 || std::abs(std::abs(Theta) - std::abs(Theta_check + M_PI)) < 0.15 || covariance_xy == 0))
    {
        // TRACE("Did not Passed")
        PRINT_2_VARIABLES(Theta, Theta_check)
        TRACE("Angle between major axis is incorrect  -- Should be the same both ways of caluclating it")
        //
        // EXCEPTION("Angle between major axis is incorrect");
    }
    assert(!isnan(norm_2(MajorAxis)));

    mMajorAxis[pottsElementIndex] = MajorAxis;
    mMajorAxis[pottsSisterIndex] = MajorAxis;

    return eig_max / eig_min;
}

template <unsigned SPACE_DIM>
c_vector<double, 2> PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetMajorAxisVector(unsigned pottsElementIndex)
{
    // PRINT_VECTOR(mMajorAxis[pottsElementIndex]);
    // TS_ASS  RT(!isnan(mMajorAxis[pottsElementIndex][0]));
    // assert(!isnan(mMajorAxis[pottsElementIndex][0]));
    return mMajorAxis[pottsElementIndex];
}

template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetMajorAxisAngle(unsigned pottsElementIndex)
{

    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);
    assert(p_potts_element->GetNumNodes() != 0);

    if (p_potts_element->GetNumNodes() <= 2)
    {
        return 0;
    }

    c_vector<double, 2> MajorAxis = mMajorAxis[pottsElementIndex];
    // now I need to find the angle between the MajorAxis and the x axis

    double theta = atan(MajorAxis[1] / MajorAxis[0]);
    if (MajorAxis[0] == 0)
    {
        theta = M_PI / 2;
    }
    return theta;
}

template <unsigned SPACE_DIM>
bool PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::IsVectorInVector(std::vector<c_vector<double, SPACE_DIM> > Vector, c_vector<double, SPACE_DIM> Location)
{
    assert(SPACE_DIM == 3);

    bool IsInVector = 0;
    for (typename std::vector<c_vector<double, SPACE_DIM> >::iterator it = Vector.begin(); it != Vector.end(); ++it)
    {
        c_vector<double, SPACE_DIM> AVector = *it;
        if ((AVector[0] == Location[0]) && (AVector[1] == Location[1]))
        {
            //  IsInVector = 1;
            return 1;
        }
    }
    return IsInVector;
}

/**************************************
          Fluids Stuff
  ***************************************/

template <unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::SetConstantWallShearStress(double WallShearStress)
{
    assert(SPACE_DIM == 3);
    for (typename MutableMesh<2, SPACE_DIM>::NodeIterator iter = mpDelaunayMesh->GetNodeIteratorBegin();
         iter != mpDelaunayMesh->GetNodeIteratorEnd();
         ++iter)
    {
        unsigned node_index = iter->GetIndex();
        c_vector<double, 3> location = iter->rGetLocation();
        // PRINT_2_VARIABLES(location[1], location[1] * 1/49+1);
        mForceOnLattice[node_index] = Create_c_vector(0, 0, 1);
        mTractionOnLattice[node_index] = Create_c_vector(0, location[1] * WallShearStress + 20, 0);

        // Need to get cell iter
        // cell_iter->GetCellData()->SetItem("WallShearMag", norm_2(mTractionOnLattice[node_index]) );
    }

}

template <unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::TractionDataLoader(const std::string& tractionFilename)
{
    TRACE("Load traction file");
    FILE* traction_file = fopen((char*)tractionFilename.c_str(), "r");
    assert(traction_file != NULL);
    hemelb::io::writers::xdr::XdrFileReader reader(traction_file);

    // File format described in http://pauli.chem.ucl.ac.uk/trac/hemelb/wiki/ExtractionFiles

    unsigned hemelb_magic_number, extraction_magic_number, extraction_version_number;
    reader.readUnsignedInt(hemelb_magic_number);
    reader.readUnsignedInt(extraction_magic_number);
    reader.readUnsignedInt(extraction_version_number);

    assert(hemelb_magic_number == 0x686c6221);
    assert(extraction_magic_number == 0x78747204);
    assert(extraction_version_number == 4);

    double voxel_size;
    reader.readDouble(voxel_size);

    c_vector<double, 3> origin;
    reader.readDouble(origin[0]);
    reader.readDouble(origin[1]);
    reader.readDouble(origin[2]);

    unsigned long long number_fluid_sites;
    reader.readUnsignedLong(number_fluid_sites);

    unsigned field_count;
    reader.readUnsignedInt(field_count);
    assert(field_count == 2); // Traction and tangetial component of traction

    unsigned header_length;
    reader.readUnsignedInt(header_length);

    // Traction field header
    std::string field_name;
    unsigned number_floats;
    double traction_offset;

    reader.readString(field_name, header_length);
    reader.readUnsignedInt(number_floats);
    assert(number_floats == 3);
    reader.readDouble(traction_offset);

    // Tangential traction field header
    double tanget_traction_offset;

    reader.readString(field_name, header_length);
    reader.readUnsignedInt(number_floats);
    assert(number_floats == 3);
    reader.readDouble(tanget_traction_offset);

    // Data section (we are reading a single timestep)
    unsigned long long timestep_num;
    reader.readUnsignedLong(timestep_num);

    mAppliedPosition.clear();
    mAppliedTractions.clear();

    for (unsigned fluid_site_index = 0; fluid_site_index < number_fluid_sites; fluid_site_index++)
    {
        {
            c_vector<unsigned, 3> coords;
            reader.readUnsignedInt(coords[0]);
            reader.readUnsignedInt(coords[1]);
            reader.readUnsignedInt(coords[2]);

            mAppliedPosition.push_back(origin + voxel_size * coords);
        }

        {
            c_vector<float, 3> traction;
            reader.readFloat(traction[0]);
            traction[0] += traction_offset;
            reader.readFloat(traction[1]);
            traction[1] += traction_offset;
            reader.readFloat(traction[2]);
            traction[2] += traction_offset;

            assert(fabs(traction[0]) < 1e10);
            assert(fabs(traction[1]) < 1e10);
            assert(fabs(traction[2]) < 1e10);

            mAppliedTractions.push_back(traction);
        }

        {
            c_vector<float, 3> tangent_traction;
            reader.readFloat(tangent_traction[0]);
            tangent_traction[0] += tanget_traction_offset;
            reader.readFloat(tangent_traction[1]);
            tangent_traction[1] += tanget_traction_offset;
            reader.readFloat(tangent_traction[2]);
            tangent_traction[2] += tanget_traction_offset;

            assert(fabs(tangent_traction[0]) < 1e10);
            assert(fabs(tangent_traction[1]) < 1e10);
            assert(fabs(tangent_traction[2]) < 1e10);

            mAppliedTangentTractions.push_back(tangent_traction);
        }
    }

    assert(mAppliedPosition.size() == number_fluid_sites);
    assert(mAppliedTractions.size() == number_fluid_sites);
    assert(mAppliedTangentTractions.size() == number_fluid_sites);
}

template <unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::CalculateTraction()
{
    assert(SPACE_DIM == 3);
    TRACE("Jess Is here");

    for (typename MutableMesh<2, SPACE_DIM>::NodeIterator iter = mpDelaunayMesh->GetNodeIteratorBegin();
         iter != mpDelaunayMesh->GetNodeIteratorEnd();
         ++iter)
    {
        unsigned node_index = iter->GetIndex();
        // Ask james about loadin in the potts mesh here, I think it is already floting around somewhere in here, I just dont know where to get it

        if (!iter->IsBoundaryNode())
        {
            c_vector<double, 3> location = iter->rGetLocation();

            TRACE("Jess needs to remove the division of 4.8");
            unsigned nearest_fluid_site = UNSIGNED_UNSET;
            double distance_to_fluid_site = DBL_MAX;
            for (unsigned fluid_site_index = 0; fluid_site_index < mAppliedPosition.size(); fluid_site_index++)
            {

                double distance = norm_2(location - (mAppliedPosition[fluid_site_index]) / 4.8);
                if (distance < distance_to_fluid_site)
                {
                    distance_to_fluid_site = distance;
                    nearest_fluid_site = fluid_site_index;
                }
            }
            assert(nearest_fluid_site != UNSIGNED_UNSET);

            // Calculate the approximate area of the voronoi region around the cell by including a third of the area
            // of all surrounding triangles. Useful for turning stresses into forces.
            // double lattice_cell_area = pPottsMesh->GetVolumeOfLatticeSite(iter->GetIndex());

            c_vector<double, 3> force = mAppliedTractions[nearest_fluid_site] * GetVolumeOfLatticeSite(node_index); //*lattice_cell_area;
            c_vector<double, 3> shear_stress = mAppliedTangentTractions[nearest_fluid_site]; //* GetVolumeOfLatticeSite(node_index);
            // // c_vector<double,3> shear_stress = 1000*Create_c_vector(shear_stress2[0]* shear_stress2[0], shear_stress2[1]* shear_stress2[1], shear_stress2[2]* shear_stress2[2]  );
            // PRINT_VECTOR(shear_stress );

            // shear_stress =  Create_c_vector(0,30,0);
            assert(fabs(force[0]) < 1e10);
            assert(fabs(force[1]) < 1e10);
            assert(fabs(force[2]) < 1e10);
            assert(fabs(shear_stress[0]) < 1e10);
            assert(fabs(shear_stress[1]) < 1e10);
            assert(fabs(shear_stress[2]) < 1e10);

            // Now I need to collect these into their own maps
            mForceOnLattice[iter->GetIndex()] = force;
            mTractionOnLattice[iter->GetIndex()] = shear_stress;
            PRINT_VECTOR(shear_stress);
            PRINT_VECTOR(force);

            double normForce = norm_2(force);
            double normShear = norm_2(shear_stress);
            // PRINT_2_VARIABLES(normForce,normShear);
        }
    }
}

template <unsigned SPACE_DIM>
unsigned PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetSizeOfElement(unsigned pottsElementIndex)
{
    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);
    return p_potts_element->GetNumNodes();
}

template <unsigned SPACE_DIM>
c_vector<double, 3> PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetTractionOnElement(unsigned pottsElementIndex)
{
    // // Average traction on element
    // TRACE("getting traction on element");
    // PRINT_VARIABLE(pottsElementIndex);
    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);

    c_vector<double, 3> potts_element_traction = Create_c_vector(0, 0, 0);

    // An element is made of a number of lattice sites, which are centered around a number of nodes in the Delaunay mesh
    for (unsigned node_index = 0; node_index < p_potts_element->GetNumNodes(); ++node_index)
    {

        unsigned latticeSiteIndex = p_potts_element->GetNodeGlobalIndex(node_index);
        potts_element_traction += mTractionOnLattice[latticeSiteIndex];
    }
    potts_element_traction /= p_potts_element->GetNumNodes();

    return potts_element_traction;
}

template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetSumShearStressOnElement(unsigned pottsElementIndex)
{
    // // Average |shear stress| acting on the potts element
    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);

    double T_mag = 0;

    // An element is made of a number of lattice sites, which are centered around a number of nodes in the Delaunay mesh
    for (unsigned node_index = 0; node_index < p_potts_element->GetNumNodes(); ++node_index)
    {

        unsigned latticeSiteIndex = p_potts_element->GetNodeGlobalIndex(node_index);
        T_mag += norm_2(mTractionOnLattice[latticeSiteIndex]);
    }
    // potts_element_traction /= p_potts_element->GetNumNodes();

    return T_mag;
}

template <unsigned SPACE_DIM>
c_vector<double, 3> PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetTractionOnLattice(unsigned latticeSiteIndex)
{
    // Average traction on element

    c_vector<double, 3> potts_lattice_traction = mTractionOnLattice[latticeSiteIndex];

    return potts_lattice_traction;
}

template <unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::UnwrapCylinder(double radius)
{

    TRACE("Trying to unwrap cylinder")
    // Loop over the nodes and redust them
    assert(SPACE_DIM == 3);
    ChastePoint<SPACE_DIM> NewLocation; // Needs to be a chaste point im setting
    // double Radius =1.5e-3;

    for (typename MutableMesh<2, SPACE_DIM>::NodeIterator node_iter = mpDelaunayMesh->GetNodeIteratorBegin();
         node_iter != mpDelaunayMesh->GetNodeIteratorEnd();
         ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();
        c_vector<double, 3> Location = node_iter->rGetLocation();

        double X_Pos = node_iter->rGetLocation()[0];
        double Y_Pos = node_iter->rGetLocation()[1];
        double Z_Pos = node_iter->rGetLocation()[2];

        // Mapped -> Circumferencial lenght

        double Theta = atan(Y_Pos / X_Pos);

        double C = radius * Theta;

        PRINT_2_VARIABLES(C, Z_Pos)
        //  node_iter->rGetLocation()[0] =C;
        mMappedLocation[node_index] = Create_c_vector(C, Z_Pos);
        NewLocation.rGetLocation()[0] = C; //NewPositions[node_index][0];
        NewLocation.rGetLocation()[1] = Z_Pos; //NewPositions[node_index][1];
        NewLocation.rGetLocation()[2] = 0; // NewPositions[node_index][2];
        mpDelaunayMesh->SetNode(node_index, NewLocation, false);
    }

    typename MutableMesh<2, SPACE_DIM>::NodeIterator node_iter = mpDelaunayMesh->GetNodeIteratorBegin();
    typename PottsMesh<SPACE_DIM>::NodeIterator potts_node_iter = this->GetNodeIteratorBegin();

    for (;
         node_iter != mpDelaunayMesh->GetNodeIteratorEnd();
         ++node_iter, ++potts_node_iter)
    {

        // PRINT_2_VARIABLES(potts_node_iter->GetIndex(), node_iter->GetIndex());
        assert(potts_node_iter->GetIndex() == node_iter->GetIndex());
        potts_node_iter->rGetModifiableLocation() = node_iter->rGetLocation();
    }
}

template <unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::CorrectedPeriodicEdges()
{

    // Need to loop over the edges and the areas and correct the lenghts, just make them the average

    double AverageEdgeLength = 0;
    for (std::vector<std::pair<unsigned, unsigned> >::iterator it = mMeshElementPairs.begin(); it != mMeshElementPairs.end(); ++it)
    {

        // Get the two nodes associated with each element
        std::vector<unsigned> CommonNodes = mMapElementPairsToNodes[*it];
        if (mPerimeterBetweenLatticeSites[CommonNodes] > mAverageEdgeLength)
        {
            mPerimeterBetweenLatticeSites[CommonNodes] = 0.5 * mAverageEdgeLength;

            // I need to do this properly
        }
    }

    for (typename MutableMesh<2, SPACE_DIM>::NodeIterator node_iter = mpDelaunayMesh->GetNodeIteratorBegin();
         node_iter != mpDelaunayMesh->GetNodeIteratorEnd();
         ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();

        if (mLatticeVolume[node_index] > mAverageArea)
        {
            mLatticeVolume[node_index] = 0.9 * mAverageArea;
        }
    }
}

template <unsigned SPACE_DIM>
MutableMesh<2, SPACE_DIM>* PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetDelaunayMesh()
{
    return mpDelaunayMesh;
}

template <unsigned SPACE_DIM>
std::vector<std::pair<unsigned, unsigned> > PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::RemoveInternalEdges(std::vector<std::pair<unsigned, unsigned> > ElementPairsAssociatedWithPottsCell)
{

    sort(ElementPairsAssociatedWithPottsCell.begin(), ElementPairsAssociatedWithPottsCell.end());
    std::vector<std::pair<unsigned, unsigned> > DuplicatesRemoved = ElementPairsAssociatedWithPottsCell;
    DuplicatesRemoved.erase(unique(DuplicatesRemoved.begin(), DuplicatesRemoved.end()), DuplicatesRemoved.end());

    std::vector<std::pair<unsigned, unsigned> > InternalPairs;
    std::set_difference(ElementPairsAssociatedWithPottsCell.begin(), ElementPairsAssociatedWithPottsCell.end(), DuplicatesRemoved.begin(), DuplicatesRemoved.end(),
                        std::inserter(InternalPairs, InternalPairs.begin()));

    // I now have the original vector, a vector with the duplicates removed, and a vector with only the repeating elements
    // I can use the set_difference to find the difference between the InternalPairs and the DuplicatesRemoved vectors
    std::vector<std::pair<unsigned, unsigned> > ExternalPairs;
    std::set_difference(DuplicatesRemoved.begin(), DuplicatesRemoved.end(), InternalPairs.begin(), InternalPairs.end(),
                        std::inserter(ExternalPairs, ExternalPairs.begin()));

    return ExternalPairs;
}

template class PottsArbitrarySurfaceIn3DMesh<2>;
template class PottsArbitrarySurfaceIn3DMesh<3>;

// #include "SerializationExportWrapperForCpp.hpp"
// EXPORT_TEMPLATE_CLASS1(PottsArbitrarySurfaceIn3DMesh, 2)
// EXPORT_TEMPLATE_CLASS1(PottsArbitrarySurfaceIn3DMesh, 3)

//  Code graveyard
// ---------------

//     // Figure out which elements are the nighbours to this one -- should be one or two neighbours
//     double a= 0;
//     c_vector<std::vector<unsigned> , 6 > NeighbourELement;
//     for (std::vector<unsigned>::iterator  Neighbour_iter = Elements.begin();
//      Neighbour_iter != Elements.end();
//     ++ Neighbour_iter)

//     {
//         c_vector<double , 3> NeighbourElementNodeIncides;
//         std::vector<unsigned> CommonNodes;

//         // TRACE("Iterating over the nodes in each neighbouring mesh element ");
//         Element<2,SPACE_DIM>* Neighbour_element = mpDelaunayMesh->GetElement(*  Neighbour_iter);
//         unsigned N_Element_index = Neighbour_element->GetIndex();

//         for (int i=0; i<3; i++)
//         {
//             Node<SPACE_DIM>* pNode = mpDelaunayMesh->GetNode(Neighbour_element->GetNodeGlobalIndex(i));
//             NeighbourElementNodeIncides[i]  = pNode->GetIndex();
//         }

//         // See if there are any nodes from element iter intersecting with the original nodes
//          std::set_intersection(NeighbourElementNodeIncides.begin(), NeighbourElementNodeIncides.end(),
//                           NodeIncides.begin(), NodeIncides.end(),
//                          std::back_inserter(CommonNodes));

//         if(CommonNodes.size()>1 & CommonNodes.size()!=3) // IF the two elements share 2 nodes, they are neighbours
//         {
//             TRACE("Common nodes greater than 1");
//            NeighbourELementMap[element_index].push_back(N_Element_index );
//            NeighbourELementMap[N_Element_index ].push_back(element_index);

//         }
//         // Elements.erase(Neighbour_iter);
//         // PRINT_VECTOR(NeighbourElementNodeIncides);
//         // PRINT_VECTOR(NodeIncides);
//         // PRINT_VECTOR(CommonNodes);
//     }
//     TRACE("Done looping");

// }

// for (std::set<unsigned>::iterator iter = containing_elements.begin();
//     iter != containing_elements.end();
//     ++iter)

// {
//         Element<2,SPACE_DIM>* Neighbour_element = mpDelaunayMesh->GetElement(* iter);
//         unsigned N_Element_index = Neighbour_element->GetIndex();
//        PRINT_VARIABLE(N_Element_index);
//         sort(NeighbourELementMap[N_Element_index].begin(), NeighbourELementMap[N_Element_index].end());
//         NeighbourELementMap[N_Element_index].erase(unique(NeighbourELementMap[N_Element_index].begin(), NeighbourELementMap[N_Element_index].end()), NeighbourELementMap[N_Element_index].end());
//         PRINT_VECTOR(NeighbourELementMap[N_Element_index]);

//     }

// I dont know why this didnt work

// double PottsCellPerimeter1 = 0;
// std::vector<std::pair<unsigned, unsigned> > ElementPairsAssociatedWithPottsCell;
// // loop over the lattice sites in this CPM cell
// for (unsigned local_lattice_index = 0; local_lattice_index < p_potts_element->GetNumNodes(); ++local_lattice_index)
// {
//     unsigned lattice_site_index = p_potts_element->GetNodeGlobalIndex(local_lattice_index);

//     // Collect the element pairs associate with this node/lattice site.
//     std::vector<std::pair<unsigned, unsigned> > ElementPairsAssociatedWithLatticeSite = mMapNodesToAssociateElementPairs[lattice_site_index];

//     // Loop over a vector here XXXX

//     for(std::vector<std::pair<unsigned, unsigned> >::iterator it = ElementPairsAssociatedWithLatticeSite.begin(); it != ElementPairsAssociatedWithLatticeSite.end(); ++it)
//         {
//               /* std::cout << *it; ... */
//               ElementPairsAssociatedWithPottsCell.push_back(*it);
//         }
//     // for (int j = 0; j < ElementPairsAssociatedWithLatticeSite.size(); j++)
//     // {
//     //     ElementPairsAssociatedWithPottsCell.push_back(ElementPairsAssociatedWithLatticeSite[j]);
//     // }
// }

// ElementPairsAssociatedWithPottsCell = RemoveInternalEdges(ElementPairsAssociatedWithPottsCell);
// // PRINT_VARIABLE(ElementPairsAssociatedWithPottsCell.size());

// // loop over the relevant element pairs and get the edge lenght between then
// // - I think i need to be looping over a map here -- this might be causing a size problem
// double i =0;
//  for(std::vector<std::pair<unsigned, unsigned> >::iterator it = ElementPairsAssociatedWithPottsCell.begin(); it != ElementPairsAssociatedWithPottsCell.end(); ++it)
//     {
//         PottsCellPerimeter1 += mDistanceBetweenElements[*it];
//         PRINT_VARIABLE(i);
//         i+=1;
//     }

// // return PottsCellPerimeter;
// // PRINT_VARIABLE(PottsCellPerimeter);

// // -------------------------

// template <unsigned SPACE_DIM>
// double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetSurfaceAreaOfElement(unsigned pottsElementIndex) // perimiter
// {
//     /*  1) Have the CPM cell the lattice site belongs to
//         2) Loop over the lattice sites in thi CPM cell and collect all of the associated mesh element pairs
//            2a) -  remove anything that appears more than once
//         3) Loop over all of the associated mesh element pairs and add the perimeter contribution from each.

//         */
//     //
//     PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);
//     //   std::map< unsigned, std::map< unsigned, bool > > mEdgeLatticeSitesList;

//     double potts_element_surface = 0.0;
//     for (unsigned local_lattice_index = 0; local_lattice_index < p_potts_element->GetNumNodes(); ++local_lattice_index)
//     {
//         unsigned lattice_site_index = p_potts_element->GetNodeGlobalIndex(local_lattice_index);
//         // potts_element_surface += mLatticePerimeters[lattice_site_index];

//         // Removing the internal elements
//         double counter =0;
//         for (std::set<unsigned>::const_iterator neighbour_lattice_index = this->mMooreNeighbouringNodeIndices[lattice_site_index].begin();
//              neighbour_lattice_index != this->mMooreNeighbouringNodeIndices[lattice_site_index].end();
//              ++neighbour_lattice_index)
//         {

//             // // if both the lattices are in the element and they are both edges, then this element has wrapped around on itself
//             // if (DoNodesShareElement(lattice_site_index, *neighbour_lattice_index))
//             // {
//             //     counter +=1;
//             //     if (mEdgeLatticeSitesList[pottsElementIndex][lattice_site_index] == 1)
//             //     {
//             //         mEdgeLatticeSitesList[pottsElementIndex][lattice_site_index] = 0;
//             //     }
//             // }

//             // If lattice sites dont share an element, then it is the edge
//             if (DoNodesShareElement(lattice_site_index, *neighbour_lattice_index) ==0)
//             {
//                 // save these two indices
//                 // TRACE("Lattice sites on the edge");
//                 // // PerimeterPairs.push_back(MathsFunctions.Create_pair(lattice_site_index, *neighbour_lattice_index) );
//                 potts_element_surface  += mDistanceBetweenElements[ mMapNodesPairsToElementPairs[MathsFunctions.Create_pair(lattice_site_index, *neighbour_lattice_index)  ]];
//                 // save the lattice site as an edge lattive of the Potts cell
//                 // mEdgeLatticeSites[pottsElementIndex].push_back(lattice_site_index);
//                 // Save this as an edge lattise site pair (one in the cell/element, the other on the outside)

//             }
//             // else if (IsNumberInVector(mEdgeLatticeSites[pottsElementIndex], lattice_site_index) && IsNumberInVector(mEdgeLatticeSites[pottsElementIndex], *neighbour_lattice_index))
//             // {
//             //     TRACE("Lattice sites share element, and are both marked as edges");
//             //     // if the nodes share the same potts element but one or more of them is an edge element I need to add it too
//             //    // this happens where the cell has wrapped around the vessel and attached to itself

//             //     // PerimeterPairs.push_back(MathsFunctions.Create_pair(lattice_site_index, *neighbour_lattice_index) );
//             //     potts_element_surface  += 2* mDistanceBetweenElements[ mMapNodesPairsToElementPairs[MathsFunctions.Create_pair(lattice_site_index, *neighbour_lattice_index)  ]];
//             // }else
//             // {
//             //     TRACE("Both on the inside ");
//             // }
//         }

//         // if (counter <6)
//         //     {
//         //         // at least one of the neighbours to this lattice site if not in the element, so it is an edge
//         //         mEdgeLatticeSitesList[pottsElementIndex][lattice_site_index] == 1;
//         //     }
//         //      if (counter == 6 && mEdgeLatticeSitesList[pottsElementIndex][lattice_site_index] == 1 )
//         //     {
//         //         if( neighbourEdges.size()!= 0)
//         //         // this lattice site has all of its neighbours, but is labeled as an edge, either needs updating, or cell has wrapped around on itself
//         //         mEdgeLatticeSitesList[pottsElementIndex][lattice_site_index] == 1;
//         //     }

//     // mEdgeLatticeSites[pottsElementIndex].erase(unique(mEdgeLatticeSites[pottsElementIndex].begin(), mEdgeLatticeSites[pottsElementIndex].end()), mEdgeLatticeSites[pottsElementIndex].end());
//     }
//     return potts_element_surface;
// }

/*  1) loop over the lattice sites of the Potts element
        2) check if the lattice site neighbours are in the cell. 
        3) For the neighbour not also in the cell, save it as a edge lattice

   **** 4) I need to think of a way to label the lattice sites of two edges of the element that have wrapped around to connect
    */

// template <unsigned SPACE_DIM>
// double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetSurfaceAreaOfElement(unsigned pottsElementIndex) // perimiter
// {
//     /*  1) Have the CPM cell the lattice site belongs to
//         2) Loop over the lattice sites in thi CPM cell and collect all of the associated mesh element pairs
//            2a) -  remove anything that appears more than once
//         3) Loop over all of the associated mesh element pairs and add the perimeter contribution from each.

//         */
//     //
//     PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);

//     double potts_element_surface = 0.0;

// std::vector<std::pair<unsigned, unsigned>> PerimeterPairs;
//     for (unsigned local_lattice_index = 0; local_lattice_index < p_potts_element->GetNumNodes(); ++local_lattice_index)
//     {
//         unsigned lattice_site_index = p_potts_element->GetNodeGlobalIndex(local_lattice_index);
//         // potts_element_surface += mLatticePerimeters[lattice_site_index];

//         // Removing the internal elements
//         for (std::set<unsigned>::const_iterator neighbour_lattice_index = this->mMooreNeighbouringNodeIndices[lattice_site_index].begin();
//              neighbour_lattice_index != this->mMooreNeighbouringNodeIndices[lattice_site_index].end();
//              ++neighbour_lattice_index)
//         {
//             // if (DoNodesShareElement(lattice_site_index, *neighbour_lattice_index))
//             // {
//             //     double distanceFromNeighbour = mDistanceBetweenElements[ mMapNodesPairsToElementPairs[MathsFunctions.Create_pair(lattice_site_index, *neighbour_lattice_index) ]];

//             //     // In Neighbourhood and are in the same CPM cell -- need to remove the pair to length
//             //     potts_element_surface -= distanceFromNeighbour;//mDistanceBetweenElements[std::make_pair(std::min(lattice_site_index, *neighbour_lattice_index), std::max(lattice_site_index, *neighbour_lattice_index))]; //GetContactAreaBetweenLatticeSite(lattice_site_index, *neighbour_lattice_index);
//             // }
//             if (DoNodesShareElement(lattice_site_index, *neighbour_lattice_index) ==0)
//             {
//                 // save these two indices
//                 // PerimeterPairs.push_back(MathsFunctions.Create_pair(lattice_site_index, *neighbour_lattice_index) );
//                 potts_element_surface2  += mDistanceBetweenElements[ mMapNodesPairsToElementPairs[MathsFunctions.Create_pair(lattice_site_index, *neighbour_lattice_index)  ]];
//                 // save the lattice site as an edge lattive of the Potts cell
//                 mEdgeLatticeSites[pottsElementIndex].push_back(lattice_site_index);

//             } else if (IsNumberInVector(mEdgeLatticeSites[pottsElementIndex], lattice_site_index))
//             {
//                 // if the nodes share the same potts element but one or more of them is an edge element I need to add it too
//               // this happens where the cell has wrapped around the vessel and attached to itself

//                 // PerimeterPairs.push_back(MathsFunctions.Create_pair(lattice_site_index, *neighbour_lattice_index) );
//                 potts_element_surface2  += mDistanceBetweenElements[ mMapNodesPairsToElementPairs[MathsFunctions.Create_pair(lattice_site_index, *neighbour_lattice_index)  ]];
//             } else if (IsNumberInVector(mEdgeLatticeSites[pottsElementIndex], *neighbour_lattice_index))
//             {
//                 // PerimeterPairs.push_back(MathsFunctions.Create_pair(lattice_site_index, *neighbour_lattice_index) );
//                 potts_element_surface2  += mDistanceBetweenElements[ mMapNodesPairsToElementPairs[MathsFunctions.Create_pair(lattice_site_index, *neighbour_lattice_index)  ]];
//             }
//         }
//     }
//     double potts_element_surface2 = 0;

//     for(std::vector<std::pair<unsigned, unsigned> >::iterator it = PerimeterPairs.begin(); it != PerimeterPairs.end(); ++it)
//         {
//               /* std::cout << *it; ... */

//               potts_element_surface2  += mDistanceBetweenElements[ mMapNodesPairsToElementPairs[*it ]];
//         }

//     // PRINT_2_VARIABLES(potts_element_surface2 , potts_element_surface2);
//     return potts_element_surface2;

// }

// template <unsigned SPACE_DIM>
// double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetSurfaceAreaOfElement(unsigned pottsElementIndex) // perimiter
// {

//     // For now just use the basic perimeter calculation

//     /*  1) Have the CPM cell the lattice site belongs to
//         2) Loop over the lattice sites in thi CPM cell and collect all of the associated mesh element pairs
//            2a) -  remove anything that appears more than once
//         3) Loop over all of the associated mesh element pairs and add the perimeter contribution from each.

//         */
//     //
//     PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);
//     c_vector<double, SPACE_DIM> CenterOfCell = this->GetCentroidOfElement(pottsElementIndex);
//     // PRINT_VECTOR(CenterOfCell);
//     // std::map<unsigned, double> Relative_theta = GetRelativeLattiePositionsWithinCell(pottsElementIndex);

//     double edgecounter = 0;
//     double potts_element_surface = 0.0;
//     std::vector<unsigned> ListOfEdges;

//     //  Loop over each of the lattice sites in the cell
//     for (unsigned local_lattice_index = 0; local_lattice_index < p_potts_element->GetNumNodes(); ++local_lattice_index)
//     {

//         unsigned lattice_site_index = p_potts_element->GetNodeGlobalIndex(local_lattice_index);

//         // potts_element_surface += mLatticePerimeters[lattice_site_index];

//         // Removing the internal elements
//         double counter = 0;
//         // iterate over all neighbours for this element to check if they are in the same element

//         for (std::set<unsigned>::const_iterator neighbour_lattice_index = this->mMooreNeighbouringNodeIndices[lattice_site_index].begin();
//              neighbour_lattice_index != this->mMooreNeighbouringNodeIndices[lattice_site_index].end();
//              ++neighbour_lattice_index)
//         {

//             if (DoNodesShareElement(lattice_site_index, *neighbour_lattice_index) == 0)
//             {

//                 potts_element_surface += mDistanceBetweenElements[mMapNodesPairsToElementPairs[MathsFunctions.Create_pair(lattice_site_index, *neighbour_lattice_index)]];
//                 counter += 1;
//                 //  } else if ( (Relative_theta[*neighbour_lattice_index] < -0.001 && Relative_theta[lattice_site_index]  > 0.001 )|| (Relative_theta[*neighbour_lattice_index] > 0.001 && Relative_theta[lattice_site_index]  < -0.001 )   )
//             }
//         }

//     }
//     return potts_element_surface;

//     // for (unsigned local_lattice_index = 0; local_lattice_index < p_potts_element->GetNumNodes(); ++local_lattice_index)
//     // {

//     //     unsigned lattice_site_index = p_potts_element->GetNodeGlobalIndex(local_lattice_index);

//     //     // potts_element_surface += mLatticePerimeters[lattice_site_index];

//     //     // Removing the internal elements
//     //     double counter = 0;
//     //     // iterate over all neighbours for this element to check if they are in the same element

//     //     for (std::set<unsigned>::const_iterator neighbour_lattice_index = this->mMooreNeighbouringNodeIndices[lattice_site_index].begin();
//     //          neighbour_lattice_index != this->mMooreNeighbouringNodeIndices[lattice_site_index].end();
//     //          ++neighbour_lattice_index)
//     //     {

//     //         if (DoNodesShareElement(lattice_site_index, *neighbour_lattice_index) == 0)
//     //         {

//     //             potts_element_surface += mDistanceBetweenElements[mMapNodesPairsToElementPairs[MathsFunctions.Create_pair(lattice_site_index, *neighbour_lattice_index)]];
//     //             counter += 1;
//     //             //  } else if ( (Relative_theta[*neighbour_lattice_index] < -0.001 && Relative_theta[lattice_site_index]  > 0.001 )|| (Relative_theta[*neighbour_lattice_index] > 0.001 && Relative_theta[lattice_site_index]  < -0.001 )   )
//     //         }

//     //         else if (DoNodesShareElement(lattice_site_index, *neighbour_lattice_index) == 1)
//     //         {
//     //              TRACE("A");
//     //             if ((Relative_theta[*neighbour_lattice_index] != Relative_theta[lattice_site_index]) && Relative_theta[lattice_site_index]!=2 &&Relative_theta[*neighbour_lattice_index]!=2)
//     //             // if (Relative_theta[*neighbour_lattice_index] != Relative_theta[lattice_site_index])
//     //             {
//     //                  TRACE("B");
//     //                 // PRINT_2_VARIABLES(Relative_theta[*neighbour_lattice_index],Relative_theta[lattice_site_index]);
//     //                 // they share an element, but they are on different ends of the cell, meaning it has wrapped around
//     //                 edgecounter += 1;
//     //             }
//     //         }
//     //     }

//     //     // if (counter > 0)
//     //     // {

//     //     //     ListOfEdges.push_back(lattice_site_index);

//     //     // }
//     // }
//     // double CentroidRadialLoction;

//     // if (p_potts_element->GetNumNodes() != ListOfEdges.size())
//     // {
//     //     TRACE("c");
//     //     c_vector<double, 3> Centroid = GetCellCentre(ListOfEdges, pottsElementIndex);
//     //     TRACE("d");
//     //     // PRINT_VECTOR(Centroid);
//     //     CentroidRadialLoction = norm_2(Create_c_vector(Centroid[0], Centroid[1]));
//     //     TRACE("e");
//     // }
//     // else
//     // {
//     //     CentroidRadialLoction = 10;
//     // }
//     // if (CentroidRadialLoction < 0.0009)
//     // {
//     //     PRINT_2_VARIABLES(CentroidRadialLoction, pottsElementIndex);
// // }
//     // bool wrapped = TestRandomWalkAlongEdges(ListOfEdges, pottsElementIndex);

//     // if (edgecounter > 0  )
//     // {
//     //     TRACE("Wrapped accroding to tests");
//     //     PRINT_2_VARIABLES(pottsElementIndex, mQuaters );
//     // }
//     // // // else if ( edgecounter > 0 )
//     // // {
//     // //      TRACE("False Positive? ");
//     // //     PRINT_2_VARIABLES(CentroidRadialLoction, pottsElementIndex);

//     // // }
//     // if (  wrapped  ==1 )
//     // {
//     //      TRACE("Random walk shows wrapping ");
//     //     PRINT_VARIABLE( pottsElementIndex);

//     // }

//     // bool wrapped = TestRandomWalk(pottsElementIndex);

//     // else if (wrapped )
//     // {
//     //     TRACE("False negative ");
//     // }else if (edgecounter > 0)
//     // {
//     //         TRACE("False positive");
//     // }
//     // else
//     // {
//     //     TRACE("Not Wrapped");
//     // }

//     // double Length = CellLength(pottsElementIndex);
//     //  double Diff1 = 2*M_PI - (2*M_PI/10);
//     //         // double Diff2 = 2*M_PI - (M_PI/10)- mAngleRange;
//     //         double Diff3 = 2*M_PI -  mAngleRange;
//     //         if (mAngleRange <= Diff1  &&  edgecounter > 0 )//Diff1 ==0 ||Diff2 ==0 || Diff3 ==0)
//     //         {
//     //             // PRINT_3_VARIABLES( mAngleRange, Diff1,  Diff3 );

//     //             // TRACE("----------  Both measures say wrapped");
//     //             bool wrapped = TestRandomWalk(pottsElementIndex);
//     //             if (wrapped ==1)
//     //             {
//     //                 TRACE("A & B - RIGHT ");
//     //             }
//     //             else {
//     //                 TRACE(" A & B - WRONG");
//     //             }

//     //         }else if  (mAngleRange <= Diff1 ) // Method B
//     //         {
//     //             // TRACE(" ///////   Only Angle range");
//     //              bool wrapped = TestRandomWalk(pottsElementIndex);
//     //             if (wrapped ==1)
//     //             {
//     //                 TRACE("B - RIGHT ");
//     //             }
//     //             else {
//     //                 TRACE("B - WRONG");
//     //             }
//     //         } else if ( edgecounter > 0 ) // Method A
//     //         {
//     //             // PRINT_3_VARIABLES( mAngleRange, Diff1,  Diff3 );
//     //             // TRACE(" ***********   only left-right examination");
//     //              bool wrapped = TestRandomWalk(pottsElementIndex);
//     //             if (wrapped ==1)
//     //             {
//     //                 TRACE("A - RIGHT ");
//     //             }
//     //             else {
//     //                 TRACE("A - WRONG");
//     //             }
//     //             //  PRINT_VARIABLE(pottsElementIndex);
//     //         }
//     //         else
//     //         {
//     //             TRACE(" Not wrapped ");
//     //         }

//     // / need someway to determine if the cell wraps right around

// }
