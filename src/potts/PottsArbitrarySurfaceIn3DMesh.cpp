#include "PottsArbitrarySurfaceIn3DMesh.hpp"

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
    TRACE("Potts mesh is constructed");
}

template <unsigned SPACE_DIM>
PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::~PottsArbitrarySurfaceIn3DMesh()
{
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
        assert(node_a_elements.size() == 1);
        assert(node_b_elements.size() == 1);
        return *node_a_elements.begin() == *node_b_elements.begin();
    }
}

template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetContactAreaBetweenLatticeSite(unsigned index_A, unsigned index_B)
{
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

    // switch (CommonElements.size())
    // {
    //     case 3:
    //         /// TODO code up for the ELEMENT_DIM=3 case. Use ComputeTriangleArea helper method
    //         NEVER_REACHED;
    //     case 2:
    //     {
    //         TRACE("More than one commone element pair, this is a problem ");
    //         NEVER_REACHED;
    //         // const c_vector<double, SPACE_DIM>& node_a_location = node_set_intersecion[0]->rGetLocation();
    //         // const c_vector<double, SPACE_DIM>& node_b_location = node_set_intersecion[1]->rGetLocation();
    //         // // TRACE("Two nodes shared between thesments");
    //         // return norm_2(node_a_location - node_b_location);
    //     }
    //     case 1:
    //     TRACE("One common elemnt pair");
    //     return mDistanceBetweenElements[ CommonElements[0]];

    //     case 0:
    //     // TRACE("No nodes shared between the elements --- this is bad??");
    //         return 0;
    //     default:
    //         NEVER_REACHED;
    // }

    // Now find the intersection

    // mDistanceBetweenElements

    // // TRACE("Need to check this");

    // /// TODO make this a template parameter and code up the ELEMENT_DIM=3 case
    // const unsigned ELEMENT_DIM = 2;

    // // Get the elements for node a and b
    // Element<2,SPACE_DIM>* mesh_element_a = mpDelaunayMesh->GetElement(index_A);
    // Element<2,SPACE_DIM>* mesh_element_b = mpDelaunayMesh->GetElement(index_B);

    // std::vector<Node<SPACE_DIM>*> node_set_a(ELEMENT_DIM+1);
    // std::vector<Node<SPACE_DIM>*> node_set_b(ELEMENT_DIM+1);

    // // Get the pointers  to all of the nodes in each of the elements -- vectors of nodes
    // for (unsigned local_node_index=0; local_node_index<ELEMENT_DIM+1; ++local_node_index)
    // {
    //     node_set_a[local_node_index] = mesh_element_a->GetNode(local_node_index);
    //     node_set_b[local_node_index] = mesh_element_b->GetNode(local_node_index);
    // }

    // // Order the nodes in each element
    // std::sort(node_set_a.begin(), node_set_a.end());
    // std::sort(node_set_b.begin(), node_set_b.end());

    // // PRINT_VECTOR(node_set_a);
    // // PRINT_VECTOR(node_set_b);

    // // Find the common nodes
    // std::vector<Node<SPACE_DIM>*> node_set_intersecion(ELEMENT_DIM);
    // typename std::vector<Node<SPACE_DIM>*>::iterator node_set_intersection_end;
    // node_set_intersection_end = std::set_intersection(node_set_a.begin(), node_set_a.end(),
    //                                                   node_set_b.begin(), node_set_b.end(),
    //                                                   node_set_intersecion.begin());

    // node_set_intersecion.resize(node_set_intersection_end - node_set_intersecion.begin());
    // // PRINT_VECTOR(node_set_intersecion);

    // //node_set_intersecion has returned the common nodes betweem the two elements -- what

    // switch (node_set_intersecion.size())
    // {
    //     case 3:
    //         /// TODO code up for the ELEMENT_DIM=3 case. Use ComputeTriangleArea helper method
    //         NEVER_REACHED;
    //     case 2:
    //     {
    //         const c_vector<double, SPACE_DIM>& node_a_location = node_set_intersecion[0]->rGetLocation();
    //         const c_vector<double, SPACE_DIM>& node_b_location = node_set_intersecion[1]->rGetLocation();
    //         // TRACE("Two nodes shared between thesments");
    //         return norm_2(node_a_location - node_b_location);
    //     }
    //     case 1:
    //     // TRACE("One nodes shared between the elements -- which is defult becuase there neighbours ");
    //     case 0:
    //     // TRACE("No nodes shared between the elements --- this is bad??");
    //         return 0.;
    //     default:
    //         NEVER_REACHED;
    // }
}

// Determine graphs connectivity
template <unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::FindElementNeighbours()
{

    TRACE("Set up neighbouring element pairs");
    double j = 0;
    assert(SPACE_DIM == 3);

    for (typename MutableMesh<2, SPACE_DIM>::ElementIterator elem_iter = mpDelaunayMesh->GetElementIteratorBegin();
         elem_iter != mpDelaunayMesh->GetElementIteratorEnd();
         ++elem_iter)
    {

        unsigned Element_index = elem_iter->GetIndex();

        // Store the node indices from this element in a vector
        std::vector<unsigned> NodeIncides;

        Node<SPACE_DIM>* pNode;

        // loop over the indices from this element and save them in the vector NodeIncides
        for (int i = 0; i < 3; i++)
        {

            pNode = mpDelaunayMesh->GetNode(mpDelaunayMesh->GetElement(Element_index)->GetNodeGlobalIndex(i));
            NodeIncides.push_back(pNode->GetIndex());
        }

        // Loop over elements to find the elements sharing nodes, indicating they are neighbours.
        for (typename MutableMesh<2, SPACE_DIM>::ElementIterator N_elem_iter = mpDelaunayMesh->GetElementIteratorBegin();
             N_elem_iter != mpDelaunayMesh->GetElementIteratorEnd();
             ++N_elem_iter)
        {
            unsigned Neighbour_Element_index = N_elem_iter->GetIndex();
            // std::vector<unsigned > TotalNodes;

            // Store the node indices of the neighbour element in a vector
            std::vector<unsigned> NeighbourNodeIncides;

            // loop over the indices from the neighbour element and save them in the vector NeighbourNodeIncides
            for (int k = 0; k < 3; k++)
            {
                // TRACE("Looping over neighbor indices");
                Node<SPACE_DIM>* pNeighbourNode = mpDelaunayMesh->GetNode(mpDelaunayMesh->GetElement(Neighbour_Element_index)->GetNodeGlobalIndex(k));
                NeighbourNodeIncides.push_back(pNeighbourNode->GetIndex());
            }

            // Determine if any nodes are in both elements, If the share 2 nodes, they are neighbours, if they share three nodes, it is the same element
            std::vector<unsigned> IntersectionVector = Intersection(NeighbourNodeIncides, NodeIncides);

            if (IntersectionVector.size() == 2)
            {
                //  TRACE("Share two nodes, so save this as neigbouring elements");
                std::pair<unsigned, unsigned> ElementPair = std::make_pair(std::min(Neighbour_Element_index, Element_index), std::max(Neighbour_Element_index, Element_index));

                if (IsPairInVector(mMeshElementPairs, ElementPair) == 0)
                {
                    sort(IntersectionVector.begin(), IntersectionVector.end());
                    std::pair<unsigned, unsigned> AssociatedNodes = Create_pair(IntersectionVector[0], IntersectionVector[1]);
                    mMeshElementPairs.push_back(ElementPair);
                    mMapElementPairsToNodes[ElementPair] = IntersectionVector;

                    mMapNodesPairsToElementPairs[AssociatedNodes] = ElementPair;
                    // mMeshElementNeighbours[Element_index].push_back(Neighbour_Element_index);

                    // Also save the paired elements for each of the nodes
                    mMapNodesToAssociateElementPairs[IntersectionVector[0]].push_back(ElementPair);
                    mMapNodesToAssociateElementPairs[IntersectionVector[1]].push_back(ElementPair);
                }
            }
        }

    } // Have now saved all of the mesh element neighbours into member vectors and maps

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
        if (AreVectorsSame(ContainingElements, ElementsFromPairs) != 0 && ElementsFromPairs.size() != 0)
        {
            EXCEPTION("Elements in element pairs are not the same as the elements contained at this node");
        }
    }
}

template <unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::CalculateEdgeLenghts()
{
    // TRACE("Recalculating the edge lengths");

    assert(SPACE_DIM == 3);

    // Calculate the Midpoint for each element
    for (typename MutableMesh<2, SPACE_DIM>::ElementIterator elem_iter = mpDelaunayMesh->GetElementIteratorBegin();
         elem_iter != mpDelaunayMesh->GetElementIteratorEnd();
         ++elem_iter)
    {
        // // Calculate the center point for each node
        unsigned Element_index = elem_iter->GetIndex();
        c_vector<double, 3> CenterPoint = Create_c_vector(0, 0, 0);
        c_vector<double, 3> NodeIncides;

        for (int i = 0; i < 3; i++)
        {
            Node<SPACE_DIM>* pNode = mpDelaunayMesh->GetNode(mpDelaunayMesh->GetElement(Element_index)->GetNodeGlobalIndex(i));
            CenterPoint += pNode->rGetLocation();
            NodeIncides[i] = pNode->GetIndex();
        }
        mMeshElementMidPoints[Element_index] = CenterPoint / 3;
    }

    //---------------

    // For each Mesh neighbour element pair calculate the lenghts between the center points, remebering this is in 3D and the elements probably dont lie in a 2D plane
    //   std::vector< std::pair <unsigned,unsigned> > mMeshElementPairs;

    for (std::vector<std::pair<unsigned, unsigned> >::iterator it = mMeshElementPairs.begin(); it != mMeshElementPairs.end(); ++it)
    {
        //    PRINT_2_VARIABLES(it->first ,it->second);  /* std::cout << *i; ... */
        // Get each of the elements in the elementpair
        unsigned Element1 = it->first;
        unsigned Element2 = it->second;

        // Get the two nodes associated with each element
        std::vector<unsigned> CommonNodes = mMapElementPairsToNodes[*it];
        Node<SPACE_DIM>* pNode0 = mpDelaunayMesh->GetNode(CommonNodes[0]);
        Node<SPACE_DIM>* pNode1 = mpDelaunayMesh->GetNode(CommonNodes[1]);

        // CommonNode[0] will be origin, make poisition relative
        c_vector<double, 3> PositionVector = pNode1->rGetLocation() - pNode0->rGetLocation();

        c_vector<double, 3> MidPoint1 = mMeshElementMidPoints[Element1] - pNode0->rGetLocation();
        c_vector<double, 3> MidPoint2 = mMeshElementMidPoints[Element2] - pNode0->rGetLocation();

        c_vector<double, 3> b1 = inner_prod(PositionVector, MidPoint1) / inner_prod(PositionVector, PositionVector) * PositionVector - MidPoint1;
        c_vector<double, 3> b2 = inner_prod(PositionVector, MidPoint2) / inner_prod(PositionVector, PositionVector) * PositionVector - MidPoint2;
        // PRINT_2_VARIABLES(norm_2(b1), norm_2(b2) );
        mDistanceBetweenElements[*it] = norm_2(b1) + norm_2(b2);

        mPerimeterBetweenLatticeSites[CommonNodes] = norm_2(b1) + norm_2(b2);
    }

    //---------------

    // Have iterated over all of the nodes and saved the area of the lattice site, proabaly not necessary
    for (typename MutableMesh<2, SPACE_DIM>::NodeIterator node_iter = mpDelaunayMesh->GetNodeIteratorBegin();
         node_iter != mpDelaunayMesh->GetNodeIteratorEnd();
         ++node_iter)
    {

        double Perimeter = 0;
        unsigned node_index = node_iter->GetIndex();
        std::set<unsigned>& containing_elements = node_iter->rGetContainingElementIndices();
        assert(containing_elements.size() > 0);
        // Loop over the elements for this node, to find the perimiter of the lines connecting to the midpoint
        for (std::set<unsigned>::iterator iter = containing_elements.begin();
             iter != containing_elements.end();
             ++iter)
        {
            unsigned elem_index = *iter;
            c_vector<c_vector<double, 3>, 2> Vectors;
            std::vector<unsigned> LocalNodes;
            double j = 0;

            // Loops over the nodes in the element and saves them as neighbours
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

            c_vector<double, 3> b1 = inner_prod(Vectors[0], MidPoint) / inner_prod(Vectors[0], Vectors[0]) * Vectors[0] - MidPoint;
            c_vector<double, 3> b2 = inner_prod(Vectors[1], MidPoint) / inner_prod(Vectors[1], Vectors[1]) * Vectors[1] - MidPoint;
            Perimeter += norm_2(b1) + norm_2(b2);
        }
        // PRINT_VARIABLE(Perimeter);
        mLatticePerimeters[node_index] = Perimeter;
        // PRINT_VARIABLE(Perimeter);
    }
}

template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetSurfaceAreaOfElement(unsigned pottsElementIndex) // perimiter
{
    /*  1) Have the CPM cell the lattice site belongs to
        2) Loop over the lattice sites in thi CPM cell and collect all of the associated mesh element pairs
           2a) -  remove anything that appears more than once 
        3) Loop over all of the associated mesh element pairs and add the perimeter contribution from each. 

        */
    //
    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);
    c_vector<double, SPACE_DIM> CenterOfCell = this->GetCentroidOfElement(pottsElementIndex);
    // PRINT_VECTOR(CenterOfCell);
    std::map<unsigned, double> Relative_theta = GetRelativeLattiePositionsWithinCell(pottsElementIndex);

    double edgecounter = 0;
    double potts_element_surface = 0.0;
    std::vector<unsigned> ListOfEdges;

    for (unsigned local_lattice_index = 0; local_lattice_index < p_potts_element->GetNumNodes(); ++local_lattice_index)
    {

        unsigned lattice_site_index = p_potts_element->GetNodeGlobalIndex(local_lattice_index);

        // potts_element_surface += mLatticePerimeters[lattice_site_index];

        // Removing the internal elements
        double counter = 0;
        // iterate over all neighbours for this element to check if they are in the same element

        for (std::set<unsigned>::const_iterator neighbour_lattice_index = this->mMooreNeighbouringNodeIndices[lattice_site_index].begin();
             neighbour_lattice_index != this->mMooreNeighbouringNodeIndices[lattice_site_index].end();
             ++neighbour_lattice_index)
        {

            if (DoNodesShareElement(lattice_site_index, *neighbour_lattice_index) == 0)
            {

                potts_element_surface += mDistanceBetweenElements[mMapNodesPairsToElementPairs[Create_pair(lattice_site_index, *neighbour_lattice_index)]];
                counter += 1;
                //  } else if ( (Relative_theta[*neighbour_lattice_index] < -0.001 && Relative_theta[lattice_site_index]  > 0.001 )|| (Relative_theta[*neighbour_lattice_index] > 0.001 && Relative_theta[lattice_site_index]  < -0.001 )   )
            }

            else if (DoNodesShareElement(lattice_site_index, *neighbour_lattice_index) == 1)
            {
                 TRACE("A");
                if ((Relative_theta[*neighbour_lattice_index] != Relative_theta[lattice_site_index]) && Relative_theta[lattice_site_index]!=2 &&Relative_theta[*neighbour_lattice_index]!=2)
                // if (Relative_theta[*neighbour_lattice_index] != Relative_theta[lattice_site_index])
                {
                     TRACE("B");
                    // PRINT_2_VARIABLES(Relative_theta[*neighbour_lattice_index],Relative_theta[lattice_site_index]);
                    // they share an element, but they are on different ends of the cell, meaning it has wrapped around
                    edgecounter += 1;
                }
            }
        }

        if (counter > 0)
        {
      
            ListOfEdges.push_back(lattice_site_index);
      
        }
    }
    // double CentroidRadialLoction;

    // if (p_potts_element->GetNumNodes() != ListOfEdges.size())
    // {
    //     TRACE("c");
    //     c_vector<double, 3> Centroid = GetCellCentre(ListOfEdges, pottsElementIndex);
    //     TRACE("d");
    //     // PRINT_VECTOR(Centroid);
    //     CentroidRadialLoction = norm_2(Create_c_vector(Centroid[0], Centroid[1]));
    //     TRACE("e");
    // }
    // else
    // {
    //     CentroidRadialLoction = 10;
    // }
    // if (CentroidRadialLoction < 0.0009)
    // {
    //     PRINT_2_VARIABLES(CentroidRadialLoction, pottsElementIndex);
// }
    // bool wrapped = TestRandomWalkAlongEdges(ListOfEdges, pottsElementIndex);

    if (edgecounter > 0  ) 
    {
        TRACE("Wrapped accroding to tests");
        PRINT_2_VARIABLES(pottsElementIndex, mQuaters );
    }
    // // else if ( edgecounter > 0 )
    // // {
    // //      TRACE("False Positive? ");
    // //     PRINT_2_VARIABLES(CentroidRadialLoction, pottsElementIndex);

    // // }
    // if (  wrapped  ==1 )
    // {
    //      TRACE("Random walk shows wrapping ");
    //     PRINT_VARIABLE( pottsElementIndex);

    // }

    

    // bool wrapped = TestRandomWalk(pottsElementIndex);

    // else if (wrapped )
    // {
    //     TRACE("False negative ");
    // }else if (edgecounter > 0)
    // {
    //         TRACE("False positive");
    // }
    // else
    // {
    //     TRACE("Not Wrapped");
    // }

    // double Length = CellLength(pottsElementIndex);
    //  double Diff1 = 2*M_PI - (2*M_PI/10);
    //         // double Diff2 = 2*M_PI - (M_PI/10)- mAngleRange;
    //         double Diff3 = 2*M_PI -  mAngleRange;
    //         if (mAngleRange <= Diff1  &&  edgecounter > 0 )//Diff1 ==0 ||Diff2 ==0 || Diff3 ==0)
    //         {
    //             // PRINT_3_VARIABLES( mAngleRange, Diff1,  Diff3 );

    //             // TRACE("----------  Both measures say wrapped");
    //             bool wrapped = TestRandomWalk(pottsElementIndex);
    //             if (wrapped ==1)
    //             {
    //                 TRACE("A & B - RIGHT ");
    //             }
    //             else {
    //                 TRACE(" A & B - WRONG");
    //             }

    //         }else if  (mAngleRange <= Diff1 ) // Method B
    //         {
    //             // TRACE(" ///////   Only Angle range");
    //              bool wrapped = TestRandomWalk(pottsElementIndex);
    //             if (wrapped ==1)
    //             {
    //                 TRACE("B - RIGHT ");
    //             }
    //             else {
    //                 TRACE("B - WRONG");
    //             }
    //         } else if ( edgecounter > 0 ) // Method A
    //         {
    //             // PRINT_3_VARIABLES( mAngleRange, Diff1,  Diff3 );
    //             // TRACE(" ***********   only left-right examination");
    //              bool wrapped = TestRandomWalk(pottsElementIndex);
    //             if (wrapped ==1)
    //             {
    //                 TRACE("A - RIGHT ");
    //             }
    //             else {
    //                 TRACE("A - WRONG");
    //             }
    //             //  PRINT_VARIABLE(pottsElementIndex);
    //         }
    //         else
    //         {
    //             TRACE(" Not wrapped ");
    //         }

    // / need someway to determine if the cell wraps right around
    return potts_element_surface;
}

template <unsigned SPACE_DIM>
bool PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::TestRandomWalkAlongEdges(std::vector<unsigned> ListOfEdges, unsigned pottsElementIndex)
{
    bool wrapped = 0;
    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);
    std::map<unsigned, c_vector<double, 2> > CellOnPlane = MapCellTo2DPlane(pottsElementIndex);

    // testing random walk from lattice index 0
    if (p_potts_element->GetNumNodes() == ListOfEdges.size())
    {
        TRACE(" Element is only made up of edges ");
        return 0; // all of the lattice sites are edges
    }
    else if (p_potts_element->GetNumNodes() <= 1)
    {
        TRACE(" Element has less than two lattices ");
        return 0;
    }

    // TRACE(" More then eddges");
 
    unsigned lattice_site_index = ListOfEdges[0]; // this is a global index

    std::set<unsigned> NeighbourSet = this->mMooreNeighbouringNodeIndices[lattice_site_index];
    std::vector<unsigned> Neighbours(NeighbourSet.begin(), NeighbourSet.end());

    std::vector<unsigned> EdgeNeighbours = Intersection(ListOfEdges, Neighbours);
    // PRINT_VECTOR(ListOfEdges);
    // PRINT_VECTOR(Neighbours);

    if (AreVectorsSame(ListOfEdges, EdgeNeighbours) == 0)
    { 
        TRACE(" Neighbour edgers and all edges are the same ");
        return 0; // all of the edges are within the neighborhood, therefore it cant wrap around.
    }

    unsigned endpoint = lattice_site_index + 10;

    std::vector<unsigned> listOfLeftNodes;
    listOfLeftNodes.push_back(lattice_site_index);

    double x = 0;
    while (endpoint != lattice_site_index)
    {
        std::vector<unsigned> NewlistOfLeftNodes;

        for (std::vector<unsigned>::iterator it = listOfLeftNodes.begin() + x; it != listOfLeftNodes.end(); ++it)
        {
            c_vector<double, 3> CurrentPosition = CellOnPlane[*it];

            NeighbourSet = this->mMooreNeighbouringNodeIndices[*it];
            std::vector<unsigned> NextNeighbours(NeighbourSet.begin(), NeighbourSet.end());
            
            EdgeNeighbours = Intersection(ListOfEdges, NextNeighbours);
          
            for (std::vector<unsigned>::iterator neighbour_lattice_index  = EdgeNeighbours.begin(); neighbour_lattice_index  != EdgeNeighbours.end(); ++neighbour_lattice_index )
            {
                c_vector<double, 3> direction = CellOnPlane[*neighbour_lattice_index] - CurrentPosition; 
                if (direction[0] > 0)
                    {
                        NewlistOfLeftNodes.push_back(*neighbour_lattice_index);
                    }
                    else if (CellOnPlane[*neighbour_lattice_index][0] < 0 && CurrentPosition[0] > 0)
                    {
                        NewlistOfLeftNodes.push_back(*neighbour_lattice_index);
                        // Has looped around here
                    }
                
            }


    
            // for (std::set<unsigned>::const_iterator neighbour_lattice_index = this->mMooreNeighbouringNodeIndices[*it].begin();
            //      neighbour_lattice_index != this->mMooreNeighbouringNodeIndices[*it].end();
            //      ++neighbour_lattice_index)
            // {
            //     // TRACE("Here is a neighbours");                              // check I havent already been using this lattice to step on
            //     if (DoNodesShareElement(*it, *neighbour_lattice_index) == 1 && IsNumberInVector(listOfLeftNodes, *neighbour_lattice_index) == 0)
            //     {
            //         // TRACE("Neighbour in element");
            //         c_vector<double, 3> direction = CellOnPlane[*neighbour_lattice_index] - CurrentPosition;
            //         // if the next lattie is in the positive direction, then take it
            //         if (direction[0] > 0)
            //         {
            //             NewlistOfLeftNodes.push_back(*neighbour_lattice_index);

            //         }
            //         else if (CellOnPlane[*neighbour_lattice_index][0] < 0 && CurrentPosition[0] > 0)
            //         {
            //             NewlistOfLeftNodes.push_back(*neighbour_lattice_index);
            //         }
            //     }
            // }
        }
    
        if (NewlistOfLeftNodes.size() == 0)
        {
          return 0;
        }

        NewlistOfLeftNodes.erase(unique(NewlistOfLeftNodes.begin(), NewlistOfLeftNodes.end()), NewlistOfLeftNodes.end());
        bool check = IsNumberInVector(NewlistOfLeftNodes, lattice_site_index);
        if (check == 1)
        {
            // endpoint = lattice_site_index;
            // wrapped = 1;
            return 1;
        }
        x = listOfLeftNodes.size() - 1; // before the new nodes are added

        for (std::vector<unsigned>::iterator it = NewlistOfLeftNodes.begin(); it != NewlistOfLeftNodes.end(); ++it)
        {
            listOfLeftNodes.push_back(*it);
        }
    }

return wrapped;

}

template <unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetCellCentre(std::vector<unsigned> Edges, unsigned pottsElementIndex)
{
    assert(SPACE_DIM == 3);

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
    // PRINT_VECTOR(Mean_Location2 );
    // PRINT_VECTOR(Mean_Location );
    return Mean_Location2;
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

/*  1) loop over the lattice sites of the Potts element
        2) check if the lattice site neighbours are in the cell. 
        3) For the neighbour not also in the cell, save it as a edge lattice

   **** 4) I need to think of a way to label the lattice sites of two edges of the element that have wrapped around to connect
    */

//   if ( (Relative_theta[*neighbour_lattice_index] < -M_PI/8 && Relative_theta[lattice_site_index]  > M_PI/8 ))
//                     {
//                         // they share an element, but they are on different ends of the cell, meaning it has wrapped around
//                         edgecounter +=1;
//                     }else if (Relative_theta[*neighbour_lattice_index] > M_PI/8 && Relative_theta[lattice_site_index]  < -M_PI/8 )
//                     {
//                         edgecounter +=1;
//                     }
//                  }

//  TestRandomWalkAlongEdges(ListOfEdges, pottsElementIndex);


template <unsigned SPACE_DIM>
std::map<unsigned, double> PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetRelativeLattiePositionsWithinCell(unsigned pottsElementIndex)
{
assert(SPACE_DIM == 3);
     TRACE("C");
     std::map<unsigned, double> D_theta;
    std::map<unsigned, double> AbsoluteAngles;

    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);

if (p_potts_element->GetNumNodes() >10)
{
    
    std::vector<double> Angles;

    double Theta = 0;
    // now loop over everything and determine its relative angle in the circumfrencial direction
    for (unsigned node_index = 0; node_index < p_potts_element->GetNumNodes(); ++node_index)
    {
        c_vector<double, 3> Location = p_potts_element->GetNodeLocation(node_index);
        unsigned lattice_site_index = p_potts_element->GetNodeGlobalIndex(node_index); // the relative positions in the map are saved based on their global node index as I need to check they relative position between neighbours
        double X_Pos = Location[0];
        double Y_Pos = Location[1];

        Theta = atan(Y_Pos / X_Pos);

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
        // // double C =  mRadius * Theta;
        // D_theta[lattice_site_index] =  -Centroid_Theta + Theta;

        AbsoluteAngles[lattice_site_index] = Theta;
        Angles.push_back(Theta);
    }

    
// translate everything towards the orgin defined at
    bool FirstQuater = 0 ;  bool SecondQuater = 0 ;  bool ThirdQuater = 0 ;   bool FourthQuater = 0 ;  
    
    for (std::vector<double>::iterator it  = Angles.begin(); it  != Angles.end(); ++it )
            {
                // PRINT_VARIABLE(*it);
                TRACE("D");
                if ( (0 <= *it) &&  ( *it<=M_PI/2) )// IN first quater
                {   FirstQuater = 1;
                 TRACE("E");
                }
                else if ( (M_PI/2 < *it ) && ( *it<=M_PI ) ) // IN second quater
                { 
                      SecondQuater = 1;
                }
                else if ( (-M_PI < *it) && ( *it < -M_PI/2 )) // IN thid quater
                {   ThirdQuater = 1;
                }
                else if ( (-M_PI/2 <= *it)  && ( *it< 0) ) // IN foruth quater
                {   FourthQuater = 1; 
                } 
                 TRACE("F");
            }
//   PRINT_4_VARIABLES(FirstQuater, SecondQuater, ThirdQuater, FourthQuater);
 TRACE("G");
        if ( SecondQuater && ThirdQuater  && FirstQuater == 0    && FourthQuater == 0 )
        {
            for (std::vector<double>::iterator it  = Angles.begin(); it  != Angles.end(); ++it )
            { //Everything needs rotating by -M_PI
                *it = AddAngles( *it, -M_PI) ;
                
            }
            // TRACE("Added -M_PI");
  
           
        }else if ( SecondQuater && ThirdQuater  && FirstQuater   && FourthQuater == 0 )
        {
            for (std::vector<double>::iterator it  = Angles.begin(); it  != Angles.end(); ++it )
            { //Everything needs rotating by -M_PI
               *it = AddAngles( *it, -M_PI) ;
            }
        //    TRACE("Added -M_PI");
        } else if ( SecondQuater && ThirdQuater  && FirstQuater  == 0  && FourthQuater )
        {
            for (std::vector<double>::iterator it  = Angles.begin(); it  != Angles.end(); ++it )
            {//Everything needs rotating by M_PI
                *it = AddAngles( *it, M_PI) ;
            }
            // TRACE("Added M_PI");
            
        }
         TRACE("H");
         PRINT_VARIABLE(Angles.size())

    // double average = accumulate(Angles.begin(),Angles.end(), 0)/ Angles.size();
  
    double max = Angles[0];   double min = Angles[0];
    for (std::vector<double>::iterator it  = Angles.begin(); it  != Angles.end(); ++it )
    {     
        if (*it > max)
        {
            max = *it;
        }  else if (*it < min)
        {
            min = *it;
        } 
    }
     TRACE("I");
    double MiddleAngle = (max - min)/2;
    // PRINT_3_VARIABLES(MiddleAngle, max ,min);
    //   cout << "_ "  << endl;

 TRACE("J");
    double AverageDifference = (max- min) / Angles.size();
    double counter =0;

    for (unsigned node_index = 0; node_index < p_potts_element->GetNumNodes(); ++node_index)
    {
         TRACE("K");
        unsigned lattice_site_index = p_potts_element->GetNodeGlobalIndex(node_index); // the relative positions in the map are saved based on their global node index as I need to check they relative position between neighbours
        // AbsoluteAngles[lattice_site_index] = AddAngles(AbsoluteAngles[lattice_site_index], - MiddleAngle);  
        Angles.push_back(Theta);
        if (Angles[counter] > MiddleAngle + 2*AverageDifference ) // the AverageDifference is the buffer so im not picking up neighbours near theta == 0
        {
            D_theta[lattice_site_index] = 1;
        }
        else if (Angles[counter] <= MiddleAngle - 2* AverageDifference )
        { 
             TRACE("L");
            D_theta[lattice_site_index] = 0;
        }else 
        { 
            D_theta[lattice_site_index] = 2;
        }
        counter += 1;
         TRACE("M");
    }
 TRACE("N");
    mQuaters = FirstQuater+ SecondQuater+ ThirdQuater+ FourthQuater;
    return D_theta;
    
}else 
{
    for (unsigned node_index = 0; node_index < p_potts_element->GetNumNodes(); ++node_index)
    {
        unsigned lattice_site_index = p_potts_element->GetNodeGlobalIndex(node_index); 
        D_theta[lattice_site_index]  = 2;
    }
    there are size problems here with the number of lattices in the element

}
    
}







template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::AddAngles(double alpha, double beta)
{
    double eta = alpha + beta;

    if (eta > M_PI)
    {  // This anlge is now negative 
        eta = eta -2*M_PI;
    }else if(eta < -M_PI)
     {  // This anlge is now positive
        eta = eta + 2*M_PI;
    }
    return eta;
}

// template <unsigned SPACE_DIM>
// double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::MiddleAngle(std::vector<double> Vector)
// {
//     std::vector<double> Angles = Vector;
//     double initalAngle = Angles[0];

//     double max = 0;
//     double min = 0;
//     // translate everything towards the orgin defined at
//     bool FirstQuater = 0 ;  bool SecondQuater = 0 ;  bool ThirdQuater = 0 ;   bool FourthQuater = 0 ;
//      for (std::vector<double>::iterator it  = Angles.begin(); it  != Angles.end(); ++it )
//         {
//            if (0 <= *it<M_PI/2 )// IN first quater
//            {   FirstQuater = 1;
//            }else if (M_PI/2 <= *it<M_PI ) // IN second quater
//            {   SecondQuater = 1;
//            }else if (-M_PI <= *it < -M_PI/2 ) // IN thid quater
//            {   ThirdQuater = 1;
//            }else if (-M_PI/2 <= *it < 0) // IN foruth quater
//            {
//                 FourthQuater = 1; 
//            }
            
//         }


//         if ( SecondQuater && ThirdQuater  && FirstQuater == 0    && FourthQuater == 0 )
//         {
//             for (std::vector<double>::iterator it  = Angles.begin(); it  != Angles.end(); ++it )
//             {
//                 *it = AddAngles( *it, -M_PI) ;
//             }
//             //Everything needs rotating by -M_PI
//         }else if ( SecondQuater && ThirdQuater  && FirstQuater   && FourthQuater == 0 )
//         {
//             for (std::vector<double>::iterator it  = Angles.begin(); it  != Angles.end(); ++it )
//             {
//                 *it = AddAngles( *it, -M_PI) ;
//             }
//             //Everything needs rotating by -M_PI
//         } else if ( SecondQuater && ThirdQuater  && FirstQuater  == 0  && FourthQuater )
//         {
//             for (std::vector<double>::iterator it  = Angles.begin(); it  != Angles.end(); ++it )
//             {
//                 *it = AddAngles( *it, M_PI) ;
//             }
//             //Everything needs rotating by M_PI
//         }

//     // double average = accumulate(Angles.begin(),Angles.end(), 0)/ Angles.size();
//     // cout << "The average is " << average << endl;
//     for (std::vector<double>::iterator it  = Angles.begin(); it  != Angles.end(); ++it )
//     {     
//            if (*it > max)
//             {
//                 max = *it;
//             }  else if (*it < min)
//             {
//                 min = *it;
//             } 
//     }
//     double Middle = (max - min)/2;
    
//     return Middle;
// }

template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::PeriodicAngle(double eta)
{
    if (eta > M_PI)
    {  // This anlge is now negative 
        eta = eta -2*M_PI;
    }else if(eta < -M_PI)
     {  // This anlge is now positive
        eta = eta + 2*M_PI;
    }
    return eta;
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
bool PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::TestRandomWalk(unsigned pottsElementIndex)
{
    bool wrapped = 0;

    PottsElement<SPACE_DIM>* p_potts_element = this->GetElement(pottsElementIndex);

    std::map<unsigned, c_vector<double, 2> > CellOnPlane = MapCellTo2DPlane(pottsElementIndex);

    // testing random walk from lattice index 0

    if (p_potts_element->GetNumNodes() > 1)
    {

        unsigned lattice_site_index = p_potts_element->GetNodeGlobalIndex(0);
        // TRACE("B");

        unsigned endpoint = lattice_site_index + 10;

        double counter = 0;
        std::vector<std::vector<unsigned> > NextInSequence;
        std::vector<unsigned> lattice_site_index_vector;

        lattice_site_index_vector.push_back(lattice_site_index);
        // TRACE("D.7");
        NextInSequence.push_back(lattice_site_index_vector);
        std::vector<unsigned> listOfLeftNodes; //      refresh the list
        listOfLeftNodes.push_back(lattice_site_index);
        // TRACE("E");
        double x = 0;
        while (endpoint != lattice_site_index)
        {

            counter += 1;

            std::vector<unsigned> NewlistOfLeftNodes;
            double LeftNeighbours = 0 ;
            for (std::vector<unsigned>::iterator it = listOfLeftNodes.begin() + x; it != listOfLeftNodes.end(); ++it)
            {
                c_vector<double, 3> CurrentPosition = CellOnPlane[*it];
                //  PRINT_VARIABLE(*/

                for (std::set<unsigned>::const_iterator neighbour_lattice_index = this->mMooreNeighbouringNodeIndices[*it].begin();
                     neighbour_lattice_index != this->mMooreNeighbouringNodeIndices[*it].end();
                     ++neighbour_lattice_index)
                {
                    // TRACE("Here is a neighbours");                              // check I havent already been using this lattice to step on
                    if (DoNodesShareElement(*it, *neighbour_lattice_index) == 1 && IsNumberInVector(listOfLeftNodes, *neighbour_lattice_index) == 0)
                    {
                        // TRACE("Neighbour in element");
                        c_vector<double, 3> direction = CellOnPlane[*neighbour_lattice_index] - CurrentPosition;
                        // PRINT_VECTOR(direction);
                        // if the next lattie is in the positive direction, then take it
                        if (direction[0] > 0)
                        {
                            // list.push_back(*neighbour_lattice_index);

                            NewlistOfLeftNodes.push_back(*neighbour_lattice_index);
                            // NextInSequence[counter].push_back(*neighbour_lattice_index);

                            LeftNeighbours += 1;
                            // x +=1;
                        }
                        else if (CellOnPlane[*neighbour_lattice_index][0] < 0 && CurrentPosition[0] > 0)
                        {
                            // this is where the cylinder would wrap around on itself.
                            // list.push_back(*neighbour_lattice_index);
                            // NextInSequence[counter].push_back(*neighbour_lattice_index);
                            NewlistOfLeftNodes.push_back(*neighbour_lattice_index);
                            LeftNeighbours += 1;
                            // x +=1;
                            // TRACE("Have wrapped ");
                        }
                    }
                }
            }
            // PRINT_VARIABLE(NextInSequence[counter].size());
            // PRINT_VARIABLE(listOfLeftNodes.size());z
            // PRINT_VECTOR( NewlistOfLeftNodes);
            // PRINT_VECTOR( listOfLeftNodes);
            if (LeftNeighbours == 0)
            {
                //  TRACE("we have run out of things on the left to step to, thus have come to the end ");
                return 0;
                // break;
            }
            // listOfLeftNodes.clear();

            NewlistOfLeftNodes.erase(unique(NewlistOfLeftNodes.begin(), NewlistOfLeftNodes.end()), NewlistOfLeftNodes.end());
            // std::vector<unsigned> listOfLeftNodesOfLeftNodesLeft     // PRINT_VECTOR(listOfLeftNodesOfLeftN
            bool check = IsNumberInVector(NewlistOfLeftNodes, lattice_site_index);
            if (check == 1)
            {
                endpoint = lattice_site_index;
                // TRACE("Loops around");
                wrapped = 1;
                // PRINT_VARIABLE(listOfLeftNodes.size());
                return wrapped;
            }
            x = listOfLeftNodes.size() - 1; // before the new nodes are added

            for (std::vector<unsigned>::iterator it = NewlistOfLeftNodes.begin(); it != NewlistOfLeftNodes.end(); ++it)
            {
                listOfLeftNodes.push_back(*it);
            }
        }
    }
    else if (p_potts_element->GetNumNodes() <= 1)
    {
        return 0;
    }

    return wrapped;
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

    assert(SPACE_DIM == 3);

    // Have iterated over all of the nodes and saved the area of the lattice site,
    for (typename MutableMesh<2, SPACE_DIM>::NodeIterator node_iter = mpDelaunayMesh->GetNodeIteratorBegin();
         node_iter != mpDelaunayMesh->GetNodeIteratorEnd();
         ++node_iter)
    {
        // Set up the area and normal to add to;
        double LatticeArea = 0;
        c_vector<double, 3> Normal = Create_c_vector(0, 0, 0);

        unsigned node_index = node_iter->GetIndex();
        std::set<unsigned>& containing_elements = node_iter->rGetContainingElementIndices();
        assert(containing_elements.size() > 0);

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
            MiniElementNormal *= MaintainOutwardsPointingNormal(MiniElementNormal, node_iter->rGetLocation());
            MiniElementNormal /= norm_2(MiniElementNormal);
            Normal += MiniElementNormal;
        }

        Normal /= norm_2(Normal);
        mLatticeVolume[node_index] = LatticeArea;
        mLatticeNormal[node_index] = Normal;
    }
}

template <unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::CalculateCurvature()
{
    assert(SPACE_DIM == 3);

    for (typename MutableMesh<2, SPACE_DIM>::NodeIterator node_iter = mpDelaunayMesh->GetNodeIteratorBegin();
         node_iter != mpDelaunayMesh->GetNodeIteratorEnd();
         ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();
        double MeanCurvature = 0;
        unsigned NumberOfNeighbours = 0;

        // Iterate ove the more neirghourhood to get the instantaenous curvature between this node and the neighbours

        // Removing the internal elements
        for (std::set<unsigned>::const_iterator neighbour_node_index = this->mMooreNeighbouringNodeIndices[node_index].begin();
             neighbour_node_index != this->mMooreNeighbouringNodeIndices[node_index].end();
             ++neighbour_node_index) // iterating over the set
        {
            NumberOfNeighbours += 1;
            // Need to find the angle between the 2 nodes --can do this using the definition of the curvature
            // k =1/R

            Node<SPACE_DIM>* pNeighbourNode = mpDelaunayMesh->GetNode(*neighbour_node_index);

            // Find the difference between the two points

            c_vector<double, 3> Vector01 = pNeighbourNode->rGetLocation() - node_iter->rGetLocation();

            // Ang;e netweem the unit normals

            double angle = acos(inner_prod(mLatticeNormal[node_index], mLatticeNormal[*neighbour_node_index]));

            // Basic trig gives

            double R = norm_2(Vector01) / (2 * sin(angle / 2));
            double Curv = 1 / R;

            //  PRINT_2_VARIABLES( Curv , R);

            std::pair<unsigned, unsigned> NodePair = Create_pair(node_index, *neighbour_node_index);
            mInstantaneousCurvature[NodePair] = 1 / R;
            mAngleBetweenNodes[NodePair] = angle;
            MeanCurvature += 1 / R;
        }
        // if NumberOfNeighbours
        c_vector<double, 3> Location = node_iter->rGetLocation();
        if (node_iter->IsBoundaryNode() && abs(Location[2]) > 5.99)
        { // Need to consider here that not all nodes labeld as the boundaries are exactly on the boundary, there is the connectivly n]line in the
            // along the lenght of the cylinder -- need to fix
            // find the neighbours that are the furthest
            double Distance1 = 0;
            double Distance2 = 0;
            unsigned Neighbour1;
            unsigned Neighbour2;

            for (std::set<unsigned>::const_iterator neighbour_node_index = this->mMooreNeighbouringNodeIndices[node_index].begin();
                 neighbour_node_index != this->mMooreNeighbouringNodeIndices[node_index].end();
                 ++neighbour_node_index) // iterating over the set
            {
                Node<SPACE_DIM>* pNeighbourNode = mpDelaunayMesh->GetNode(*neighbour_node_index);
                double distance = norm_2(pNeighbourNode->rGetLocation() - node_iter->rGetLocation());

                if (distance > Distance1)
                {
                    distance = Distance1;
                    Neighbour1 = *neighbour_node_index;
                }
                else if (distance > Distance2)
                {
                    distance = Distance2;
                    Neighbour2 = *neighbour_node_index;
                }
            }

            // now have the two neighbours I need, go get the curvature

            std::pair<unsigned, unsigned> NodePair1 = Create_pair(node_index, Neighbour1);
            std::pair<unsigned, unsigned> NodePair2 = Create_pair(node_index, Neighbour2);
            MeanCurvature += mInstantaneousCurvature[NodePair1] + mInstantaneousCurvature[NodePair2];
            NumberOfNeighbours += 2;
        }

        MeanCurvature /= NumberOfNeighbours;
        // PRINT_VARIABLE();
        mMeanCurvature[node_index] = MeanCurvature;
        bool Bound = node_iter->IsBoundaryNode();

        // PRINT_2_VARIABLES( MeanCurvature, Bound );
    }
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

        double X_Pos = mMappedLocation[node_index][0]; //  p_potts_element->GetNode(i)->rGetLocation()[0];
        double Y_Pos = mMappedLocation[node_index][1]; //   p_potts_element->GetNode(i)->rGetLocation()[1];

        mean_x += X_Pos;
        mean_y += Y_Pos;
    }
    mean_x /= p_potts_element->GetNumNodes();
    mean_y /= p_potts_element->GetNumNodes();

    double variance_x = 0;
    double variance_y = 0;
    double covariance_xy = 0;

    for (unsigned i = 0; i < p_potts_element->GetNumNodes(); i++)
    {
        Node<SPACE_DIM>* pNode = p_potts_element->GetNode(i);
        unsigned node_index = pNode->GetIndex();
        // c_vector<double, 2> MappedPosition =  mMappedLocation[node_index];

        double X_Pos = mMappedLocation[node_index][0]; //  p_potts_element->GetNode(i)->rGetLocation()[0];
        double Y_Pos = mMappedLocation[node_index][1]; //   p_potts_element->GetNode(i)->rGetLocation()[1];

        variance_x += pow((X_Pos - mean_x), 2);
        variance_y += pow((Y_Pos - mean_y), 2);
        covariance_xy += (X_Pos - mean_x) * (Y_Pos - mean_y);
    }
    variance_x /= p_potts_element->GetNumNodes();
    variance_y /= p_potts_element->GetNumNodes();
    covariance_xy /= p_potts_element->GetNumNodes();

    // Calculate max/min eigenvalues
    double trace = variance_x + variance_y;
    double det = variance_x * variance_y - covariance_xy * covariance_xy;

    eig_max = 0.5 * (trace + sqrt(trace * trace - 4 * det));
    eig_min = 0.5 * (trace - sqrt(trace * trace - 4 * det));

    if (eig_min == 0)
    {
        TRACE("All nodes in an element lie in the same line/plane (2D/3D) so aspect ratio is infinite. This interferes with calculation of the Hamiltonian.");
        return -1;
        // EXCEPTION("All nodes in an element lie in the same line/plane (2D/3D) so aspect ratio is infinite. This interferes with calculation of the Hamiltonian.");
    }

    // TO get the major axis will be the eigenvector corresponding to to the major eigenvaule eig_max -- basic eigenvector calcuations

    double MajorAxisX = (variance_x - eig_max) * (variance_y - eig_max) / (covariance_xy * covariance_xy);
    double MajorAxisY = -(variance_x - eig_max) * (variance_x - eig_max) * (variance_y - eig_max) / (covariance_xy * covariance_xy * covariance_xy);

    mMajorAxis[pottsElementIndex] = Create_c_vector(MajorAxisX, MajorAxisY);
    return eig_max / eig_min;
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
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::TractionDataLoader(const std::string& tractionFilename)
{
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

    for (typename MutableMesh<2, SPACE_DIM>::NodeIterator iter = mpDelaunayMesh->GetNodeIteratorBegin();
         iter != mpDelaunayMesh->GetNodeIteratorEnd();
         ++iter)
    {
        unsigned node_index = iter->GetIndex();
        // Ask james about loadin in the potts mesh here, I think it is already floting around somewhere in here, I just dont know where to get it

        if (!iter->IsBoundaryNode())
        {
            c_vector<double, 3> location = iter->rGetLocation();

            unsigned nearest_fluid_site = UNSIGNED_UNSET;
            double distance_to_fluid_site = DBL_MAX;
            for (unsigned fluid_site_index = 0; fluid_site_index < mAppliedPosition.size(); fluid_site_index++)
            {
                double distance = norm_2(location - mAppliedPosition[fluid_site_index]);
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

            double normForce = norm_2(force);
            double normShear = norm_2(shear_stress);
            // PRINT_2_VARIABLES(normForce,normShear);
        }
    }
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
c_vector<double, 3> PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetTractionOnLattice(unsigned latticeSiteIndex)
{
    // Average traction on element

    c_vector<double, 3> potts_lattice_traction = mTractionOnLattice[latticeSiteIndex];

    return potts_lattice_traction;
}

template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetCurvatureOfElement(unsigned pottsElementIndex)
{
    assert(SPACE_DIM == 3);
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
                std::pair<unsigned, unsigned> NodePair = Create_pair(lattice_site_index, *neighbour_lattice_index);
                bool AlreadyCounted = IsPairInVector(VectorOfElementPairs, NodePair);
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

    FindElementNeighbours();
    mAngleRange = 0;
    MapCylinderToPlane();
    CalculateEdgeLenghts();
    CalculateLatticeVolumes();
    CalculateCurvature();
    CalculateTraction();
}

template <unsigned SPACE_DIM>
MutableMesh<2, SPACE_DIM>* PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::GetDelaunayMesh()
{
    return mpDelaunayMesh;
}

// Axciallary functions

template <unsigned SPACE_DIM>
std::pair<double, double> PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::Create_pair(double x, double y)
{
    std::pair<double, double> NodePair = std::make_pair(std::min(x, y), std::max(x, y));
    return NodePair;
}

template <unsigned SPACE_DIM>
std::pair<unsigned, unsigned> PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::Create_pair(unsigned x, unsigned y)
{
    std::pair<unsigned, unsigned> NodePair = std::make_pair(std::min(x, y), std::max(x, y));
    return NodePair;
}

template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::MaintainOutwardsPointingNormal(c_vector<double, 3> Normal, c_vector<double, 3> x1)
{

    double direction = 1;
    c_vector<double, 2> Normal2D = Create_c_vector(Normal[0], Normal[1]);

    c_vector<double, 2> x12D = Create_c_vector(x1[0], x1[1]);

    c_vector<double, 2> Extension = Normal2D + x12D;

    double absExtensin = norm_2(Extension);
    double absx1 = norm_2(x12D);
    if (absExtensin < absx1) // meaning the normal points towards the center
    {
        direction = -1;
        //  Normal = -Normal;// reverse the normal so it points out
    }
    return direction;

    // 2d Plane, so this need to be different

    //  double direction  =1;

    // c_vector<double, 3> Extension = Normal+ x1;

    //  if (Extension[2] < x1[2]) // meaning the normal points down
    //  {
    //     direction  =-1;
    //     //  Normal = -Normal;// reverse the normal so it points out
    //  }
    //  return direction ;
}

template <unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::PRINT_PAIR(std::pair<unsigned, unsigned> Pair)
{
    std::cout << "DEBUG: Pair = [" << Pair.first << ", " << Pair.second << "] " << endl;
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

template <unsigned SPACE_DIM>
void PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::TestRemoveInternalEdgesFunction()
{
    std::vector<std::pair<unsigned, unsigned> > TestVector;
    std::vector<std::pair<unsigned, unsigned> > ExpectedPairsVector;
    std::pair<unsigned, unsigned> Pair1 = std::make_pair(4, 6);
    std::pair<unsigned, unsigned> Pair2 = std::make_pair(1, 2); // unique
    ExpectedPairsVector.push_back(Pair2);
    std::pair<unsigned, unsigned> Pair3 = std::make_pair(4, 6);
    std::pair<unsigned, unsigned> Pair4 = std::make_pair(3, 4);
    std::pair<unsigned, unsigned> Pair5 = std::make_pair(1, 5);
    std::pair<unsigned, unsigned> Pair6 = std::make_pair(3, 4);
    std::pair<unsigned, unsigned> Pair7 = std::make_pair(1, 5);
    std::pair<unsigned, unsigned> Pair8 = std::make_pair(3, 8); // unique
    ExpectedPairsVector.push_back(Pair8);
    std::pair<unsigned, unsigned> Pair9 = std::make_pair(3, 5); // unique
    ExpectedPairsVector.push_back(Pair9);
    std::pair<unsigned, unsigned> Pair10 = std::make_pair(1, 7); // unique
    ExpectedPairsVector.push_back(Pair10);
    sort(ExpectedPairsVector.begin(), ExpectedPairsVector.end());

    TestVector.push_back(Pair1);
    TestVector.push_back(Pair2);
    TestVector.push_back(Pair3);
    TestVector.push_back(Pair4);
    TestVector.push_back(Pair5);
    TestVector.push_back(Pair6);
    TestVector.push_back(Pair7);
    TestVector.push_back(Pair8);
    TestVector.push_back(Pair9);
    TestVector.push_back(Pair10);

    TestVector = RemoveInternalEdges(TestVector);

    assert(AreVectorsSame(TestVector, ExpectedPairsVector) == 0);
}

template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::AreVectorsSame(std::vector<unsigned> Vector1, std::vector<unsigned> Vector2)
{
    sort(Vector1.begin(), Vector1.end());
    sort(Vector2.begin(), Vector2.end());

    std::vector<unsigned> DifferenceVector;
    std::set_difference(Vector1.begin(), Vector1.end(), Vector2.begin(), Vector2.end(),
                        std::inserter(DifferenceVector, DifferenceVector.begin()));

    double difference = DifferenceVector.size();
    return difference;
}

template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::AreVectorsSame(std::vector<std::pair<unsigned, unsigned> > Vector1, std::vector<std::pair<unsigned, unsigned> > Vector2)
{
    sort(Vector1.begin(), Vector1.end());
    sort(Vector2.begin(), Vector2.end());

    std::vector<std::pair<unsigned, unsigned> > DifferenceVector;
    std::set_difference(Vector1.begin(), Vector1.end(), Vector2.begin(), Vector2.end(),
                        std::inserter(DifferenceVector, DifferenceVector.begin()));

    double difference = DifferenceVector.size();
    return difference;
}

template <unsigned SPACE_DIM>
std::vector<unsigned> PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::Intersection(std::vector<unsigned> Vector1, std::vector<unsigned> Vector2)
{
    std::vector<unsigned> IntersectionVector;

    // sort(Vector1.begin(), Vector1.end());
    // sort(Vector2.begin(), Vector2.end());

    for (std::vector<unsigned>::iterator it = Vector1.begin(); it != Vector1.end(); ++it)
    {
        for (std::vector<unsigned>::iterator it2 = Vector2.begin(); it2 != Vector2.end(); ++it2)
        {
            if (*it == *it2)
            {
                IntersectionVector.push_back(*it);
            }
        }
    }

    return IntersectionVector;
}

template <unsigned SPACE_DIM>
std::vector<unsigned> PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::RemoveElement(std::vector<unsigned> Vector1, unsigned number)
{
    // std::vector<unsigned> NewVector;

    // for (std::vector<unsigned>::iterator it = Vector1.begin(); it != Vector1.end(); ++it)
    // {
    //     if (*it != number)
    //     {
    //     NewVector.push_back(*it);
    //     }
    // }

    // return NewVector;

    unsigned ElementToRemove;
    for (std::vector<unsigned>::iterator it = Vector1.begin(); it != Vector1.end(); ++it)
    {
        if (*it == number)
        {
            ElementToRemove = *it;
        }
    }
    Vector1.erase(Vector1.begin() + ElementToRemove);

    return Vector1;
}

template <unsigned SPACE_DIM>
std::vector<double> PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::RemoveElement(std::vector<double> Vector1, double number)
{
    unsigned ElementToRemove;
    for (std::vector<double>::iterator it = Vector1.begin(); it != Vector1.end(); ++it)
    {
        if (*it == number)
        {
            ElementToRemove = *it;
        }
    }
    Vector1.erase(Vector1.begin() + ElementToRemove);

    return Vector1;
}

template <unsigned SPACE_DIM>
bool PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::IsPairInVector(std::vector<std::pair<unsigned, unsigned> > Vector, std::pair<unsigned, unsigned> Pair)
{
    bool IsInVector = 0;

    for (std::vector<std::pair<unsigned, unsigned> >::iterator it = Vector.begin(); it != Vector.end(); ++it)
    {
        if (*it == Pair)
        {
            IsInVector = 1;
            return IsInVector;
        }
        /* std::cout << *it; ... */
    }

    // for (int i = 0; i < Vector.size(); i++)
    // {
    //     if (Vector[i] == Pair)
    //     {
    //         IsInVector = 1;
    //         return IsInVector;
    //     }
    // }
    return IsInVector;
}

template <unsigned SPACE_DIM>
bool PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::IsNumberInVector(std::vector<unsigned> Vector, unsigned number)
{
    bool IsInVector = 0;
    for (std::vector<unsigned>::iterator it = Vector.begin(); it != Vector.end(); ++it)
    {
        if (*it == number)
        {
            IsInVector = 1;
        } /* std::cout << *it; ... */
    }

    // for (int i = 0; i < Vector.size(); i++)
    // {
    //     if (Vector[i] == number)
    //     {
    //         IsInVector = 1;
    //     }
    // }
    return IsInVector;
}

template <unsigned SPACE_DIM>
bool PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::IsNumberInVector(std::vector<double> Vector, double number)
{
    bool IsInVector = 0;
    for (std::vector<double>::iterator it = Vector.begin(); it != Vector.end(); ++it)
    {
        if (*it == number)
        {
            IsInVector = 1;
        } /* std::cout << *it; ... */
    }

    // for (int i = 0; i < Vector.size(); i++)
    // {
    //     if (Vector[i] == number)
    //     {
    //         IsInVector = 1;
    //     }
    // }
    return IsInVector;
}

template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::min_value(std::vector<double> Vector)
{
    double min = Vector[0];
    for (std::vector<double>::iterator it = Vector.begin(); it != Vector.end(); ++it)
    {
        /* std::cout << *it; ... */
        if (*it < min)
        {
            min = *it;
        } //else if (min + *it < -M_PI)//AddAngles(*it,min)
        // {
        //  // Have crossed the periodic boundary 
        //     min = *it;
        // }
    }
    return min;
}

template <unsigned SPACE_DIM>
double PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>::max_value(std::vector<double> Vector)
{
    double max = Vector[0];
    for (std::vector<double>::iterator it = Vector.begin(); it != Vector.end(); ++it)
    {
        /* std::cout << *it; ... */
        if (*it > max)
        {
            max = *it;
        }
    }
    return max;
}

template class PottsArbitrarySurfaceIn3DMesh<2>;
template class PottsArbitrarySurfaceIn3DMesh<3>;

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
//                 // // PerimeterPairs.push_back(Create_pair(lattice_site_index, *neighbour_lattice_index) );
//                 potts_element_surface  += mDistanceBetweenElements[ mMapNodesPairsToElementPairs[Create_pair(lattice_site_index, *neighbour_lattice_index)  ]];
//                 // save the lattice site as an edge lattive of the Potts cell
//                 // mEdgeLatticeSites[pottsElementIndex].push_back(lattice_site_index);
//                 // Save this as an edge lattise site pair (one in the cell/element, the other on the outside)

//             }
//             // else if (IsNumberInVector(mEdgeLatticeSites[pottsElementIndex], lattice_site_index) && IsNumberInVector(mEdgeLatticeSites[pottsElementIndex], *neighbour_lattice_index))
//             // {
//             //     TRACE("Lattice sites share element, and are both marked as edges");
//             //     // if the nodes share the same potts element but one or more of them is an edge element I need to add it too
//             //    // this happens where the cell has wrapped around the vessel and attached to itself

//             //     // PerimeterPairs.push_back(Create_pair(lattice_site_index, *neighbour_lattice_index) );
//             //     potts_element_surface  += 2* mDistanceBetweenElements[ mMapNodesPairsToElementPairs[Create_pair(lattice_site_index, *neighbour_lattice_index)  ]];
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
//             //     double distanceFromNeighbour = mDistanceBetweenElements[ mMapNodesPairsToElementPairs[Create_pair(lattice_site_index, *neighbour_lattice_index) ]];

//             //     // In Neighbourhood and are in the same CPM cell -- need to remove the pair to length
//             //     potts_element_surface -= distanceFromNeighbour;//mDistanceBetweenElements[std::make_pair(std::min(lattice_site_index, *neighbour_lattice_index), std::max(lattice_site_index, *neighbour_lattice_index))]; //GetContactAreaBetweenLatticeSite(lattice_site_index, *neighbour_lattice_index);
//             // }
//             if (DoNodesShareElement(lattice_site_index, *neighbour_lattice_index) ==0)
//             {
//                 // save these two indices
//                 // PerimeterPairs.push_back(Create_pair(lattice_site_index, *neighbour_lattice_index) );
//                 potts_element_surface2  += mDistanceBetweenElements[ mMapNodesPairsToElementPairs[Create_pair(lattice_site_index, *neighbour_lattice_index)  ]];
//                 // save the lattice site as an edge lattive of the Potts cell
//                 mEdgeLatticeSites[pottsElementIndex].push_back(lattice_site_index);

//             } else if (IsNumberInVector(mEdgeLatticeSites[pottsElementIndex], lattice_site_index))
//             {
//                 // if the nodes share the same potts element but one or more of them is an edge element I need to add it too
//               // this happens where the cell has wrapped around the vessel and attached to itself

//                 // PerimeterPairs.push_back(Create_pair(lattice_site_index, *neighbour_lattice_index) );
//                 potts_element_surface2  += mDistanceBetweenElements[ mMapNodesPairsToElementPairs[Create_pair(lattice_site_index, *neighbour_lattice_index)  ]];
//             } else if (IsNumberInVector(mEdgeLatticeSites[pottsElementIndex], *neighbour_lattice_index))
//             {
//                 // PerimeterPairs.push_back(Create_pair(lattice_site_index, *neighbour_lattice_index) );
//                 potts_element_surface2  += mDistanceBetweenElements[ mMapNodesPairsToElementPairs[Create_pair(lattice_site_index, *neighbour_lattice_index)  ]];
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
