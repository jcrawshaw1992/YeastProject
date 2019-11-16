

#include "ArbitrarySurfaceAreaConstraintPottsUpdateRule.hpp"
#include "CPMGTPaseEventHandler.hpp"
#include "Debug.hpp"


template <unsigned DIM>
ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::ArbitrarySurfaceAreaConstraintPottsUpdateRule()
        : AbstractPottsUpdateRule<DIM>(),
          mSurfaceAreaEnergyParameter(0.5), // Educated guess
          mTargetSurfaceArea(16.0) // Defaults to a 4*4 cell size

{
    /// \todo Default values don't apply in 3D.
}

template <unsigned DIM>
ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::~ArbitrarySurfaceAreaConstraintPottsUpdateRule()
{
}
//********



// 



template <unsigned DIM>
double ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::DetermineAReasonableTargetPerimeter(double x, double y)
{
   // 1) Decide how many lattice sites I want in my cell, which is x wide and y long. Im thinking 5*10 (5 accross and 10 up)

    unsigned lattice_site_index = 22; 
    std::vector<std::pair<unsigned, unsigned> > ElementPairsAssociatedWithLatticeSite = mMapNodesToAssociateElementPairs[lattice_site_index];
    std::vector<std::pair<unsigned, unsigned> > ElementPairsAssociatedWithPottsCell;
   
    std::vector<unsigned> EdgeLenghts ;
    for (int j = 0; j < ElementPairsAssociatedWithLatticeSite.size(); j++)
    {
        EdgeLenghts.push_back(mDistanceBetweenElements[ElementPairsAssociatedWithLatticeSite[j] ]);
        
    }
    PRINT_VARIABLE(ElementPairsAssociatedWithLatticeSite.size());
    PRINT_VECTOR(EdgeLenghts); // These are giving zeros, dont know why 


    unsigned up = *std::min_element(EdgeLenghts.begin(), EdgeLenghts.end());
    unsigned diagonal = *std::max_element(EdgeLenghts.begin(), EdgeLenghts.end());

    // This is calculated based on a rectangular Potts cell 
    double UpwardsLength = y * up + y * 2 * diagonal;
    double HorizontalLenght = x * 2 * diagonal;

    double TargetPerimeter = 2*(UpwardsLength  + HorizontalLenght  );
    PRINT_VARIABLE(TargetPerimeter);
    return TargetPerimeter;


}




//---------------------

template <unsigned DIM>
double ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                                           unsigned targetNodeIndex,
                                                                                           PottsBasedCellPopulation<DIM>& rCellPopulation)
{
    CPMGTPaseEventHandler::BeginEvent(CPMGTPaseEventHandler::SURFACE);
  
//   TRACE("Evaluate Surface area Hamiltonian");
   PottsArbitrarySurfaceIn3DMesh<DIM>* p_static_cast_potts_mesh = static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*> (&(rCellPopulation.rGetMesh())); // Jess needs to edit to be  PottsArbitrarySurfaceIn3DMesh<DIM>
 
 
// note the other way to do this is do the static cast each time you need it 
//static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*>(&rCellPopulation.rGetMesh())->GetVolumeOfLatticeSite(targetNodeIndex);


    double delta_H = 0.0;
    std::set<unsigned> containing_elements = p_static_cast_potts_mesh->GetNode(currentNodeIndex)->rGetContainingElementIndices();
    // TRACE("B");
    std::set<unsigned> new_location_containing_elements = p_static_cast_potts_mesh->GetNode(targetNodeIndex)->rGetContainingElementIndices();
    // TRACE("C");
    bool current_node_contained = !containing_elements.empty();
    bool target_node_contained = !new_location_containing_elements.empty();

    // Every node must each be in at most one element
    assert(new_location_containing_elements.size() < 2);

    if (!current_node_contained && !target_node_contained)
    {
        EXCEPTION("At least one of the current node or target node must be in an element.");
    }

    if (current_node_contained && target_node_contained)
    {
        if (*(new_location_containing_elements.begin()) == *(containing_elements.begin()))
        {
            EXCEPTION("The current node and target node must not be in the same element.");
        }
    }

    //double surface_area_target_lattice_site = mLatticePerimeters[targetNodeIndex] ;//static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*>(&rCellPopulation.rGetMesh())->GetSurfaceAreaOfLatticeSite(targetNodeIndex);

    if (current_node_contained) // current node is in an element
    {
        unsigned current_element = (*containing_elements.begin());
        PottsElement<DIM>* pCurrentElement = p_static_cast_potts_mesh->GetElement(current_element);


        double current_surface_area = p_static_cast_potts_mesh->GetSurfaceAreaOfElement(current_element);
        // PRINT_VARIABLE(current_surface_area);
        // GetSurfaceAreaOfElement(current_element, pCurrentElement);
        
        double current_surface_area_difference = current_surface_area - mTargetSurfaceArea;
        // PRINT_3_VARIABLES(current_surface_area , mTargetSurfaceArea,  current_surface_area_difference);
        // PRINT_VARIABLE(current_surface_area_difference);

        // Mock change of spin to call GetSurfaceAreaOfElement with new configuration
        Node<DIM>* pTargetNode = p_static_cast_potts_mesh->GetNode(targetNodeIndex);
        pCurrentElement->AddNode(pTargetNode);
        if (target_node_contained)
        {
            unsigned target_element = (*new_location_containing_elements.begin());
            PottsElement<DIM>* pTargetElement = p_static_cast_potts_mesh->GetElement(target_element);
            pTargetElement->DeleteNode(pTargetElement->GetNodeLocalIndex(pTargetNode->GetIndex()));
        }

        double current_surface_area_after_switch = p_static_cast_potts_mesh->GetSurfaceAreaOfElement(current_element);//, pCurrentElement);
        
        double current_surface_area_difference_after_switch = current_surface_area_after_switch - mTargetSurfaceArea;
        // PRINT_2_VARIABLES(current_surface_area_difference_after_switch, current_surface_area_after_switch);
        // Undo mocked change
        pCurrentElement->DeleteNode(pCurrentElement->GetNodeLocalIndex(pTargetNode->GetIndex()));
        if (target_node_contained)
        {
            unsigned target_element = (*new_location_containing_elements.begin());
            PottsElement<DIM>* pTargetElement = p_static_cast_potts_mesh->GetElement(target_element);
            pTargetElement->AddNode(pTargetNode);
        }

        delta_H += mSurfaceAreaEnergyParameter * (current_surface_area_difference_after_switch * current_surface_area_difference_after_switch - current_surface_area_difference * current_surface_area_difference);
    }


    
    if (target_node_contained) // target node is in an element
    {
        unsigned target_element = (*new_location_containing_elements.begin());
        PottsElement<DIM>* pTargetElement = p_static_cast_potts_mesh->GetElement(target_element);

        double target_surface_area = p_static_cast_potts_mesh->GetSurfaceAreaOfElement(target_element);//, pTargetElement);
        // GetSurfaceAreaOfElement(target_element, pTargetElement);// p_static_cast_potts_mesh->GetSurfaceAreaOfElement(target_element);
        double target_surface_area_difference = target_surface_area - mTargetSurfaceArea;

        Node<DIM>* pTargetNode = p_static_cast_potts_mesh->GetNode(targetNodeIndex);
        pTargetElement->DeleteNode(pTargetElement->GetNodeLocalIndex(pTargetNode->GetIndex()));

        double target_surface_area_after_switch = p_static_cast_potts_mesh->GetSurfaceAreaOfElement(target_element);//, pTargetElement);
        // GetSurfaceAreaOfElement(target_element, pTargetElement); 
        double target_surface_area_difference_after_switch = target_surface_area_after_switch - mTargetSurfaceArea;

        // Add node back after calculating aspect ratio
        pTargetElement->AddNode(pTargetNode);

        delta_H += mSurfaceAreaEnergyParameter * (target_surface_area_difference_after_switch * target_surface_area_difference_after_switch - target_surface_area_difference * target_surface_area_difference);
    }

    CPMGTPaseEventHandler::EndEvent(CPMGTPaseEventHandler::SURFACE);

    return delta_H;

}

template <unsigned DIM>
double ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::GetSurfaceAreaEnergyParameter()
{
    return mSurfaceAreaEnergyParameter;
}

template <unsigned DIM>
void ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::SetSurfaceAreaEnergyParameter(double surfaceAreaEnergyParameter)
{
    mSurfaceAreaEnergyParameter = surfaceAreaEnergyParameter;
}












// template <unsigned DIM>
// void ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::PRINT_PAIR(std::pair<unsigned, unsigned> Pair)
// {
//     std::cout << "DEBUG: Pair = [" << Pair.first << ", " << Pair.second << "] " << endl;
// }




template <unsigned DIM>
double ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::GetTargetSurfaceArea() const
{
    return mTargetSurfaceArea;
}

template <unsigned DIM>
void ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::SetTargetSurfaceArea(double targetSurfaceArea)
{
    assert(targetSurfaceArea >= 0.0);
    mTargetSurfaceArea = targetSurfaceArea;
}

   
template <unsigned DIM>
void ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<SurfaceAreaEnergyParameter >" << mSurfaceAreaEnergyParameter << "</SurfaceAreaEnergyParameter >\n";
    *rParamsFile << "\t\t\t<TargetSurfaceArea>" << mTargetSurfaceArea << "</TargetSurfaceArea>\n";

    // Call method on direct parent class
    AbstractPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

//template class ArbitraryVolumeConstraintPottsUpdateRule<1>;
template class ArbitrarySurfaceAreaConstraintPottsUpdateRule<2>;
template class ArbitrarySurfaceAreaConstraintPottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
//EXPORT_TEMPLATE_CLASS_SAME_DIMS(ArbitraryVolumeConstraintPottsUpdateRule)
EXPORT_TEMPLATE_CLASS1(ArbitrarySurfaceAreaConstraintPottsUpdateRule, 2)
EXPORT_TEMPLATE_CLASS1(ArbitrarySurfaceAreaConstraintPottsUpdateRule, 3)

// if (IsNumberInVector(ElementsFromPairs, VectorOfElementPairs[i].first) ==0 )
// {
//     ElementsFromPairs.push_back(VectorOfElementPairs[i].first);
// }
// if (IsNumberInVector(ElementsFromPairs, VectorOfElementPairs[i].second) ==0 )
// {
//     ElementsFromPairs.push_back(VectorOfElementPairs[i].second);
// }

  // for (int j = 0; j < ElementPairsAssociatedWithPottsCell.size(); j++)
    //     {
    //         PRINT_PAIR(ElementPairsAssociatedWithPottsCell[j]);
    //     }

    // 

    //     // Need a nicer way to do this
    //     potts_element_surface += mLatticePerimeters[lattice_site_index]; //GetPerimeterOfLatticeSite(lattice_site_index);

    //     // Have added the perimeter from all of the nodes in the Potts cell, now need to remove all of the internal perimeters.

    //     // std::set<unsigned>& containing_elements = node_iter->rGetContainingElementIndices();

    //     // for (std::set<unsigned>::const_iterator neighbour_lattice_index = p_static_cast_potts_mesh->mMooreNeighbouringNodeIndices[lattice_site_index].begin();
    //     //          neighbour_lattice_index != p_static_cast_potts_mesh->mMooreNeighbouringNodeIndices[lattice_site_index].end();
    //     //          ++neighbour_lattice_index)
    //     //     {
    //     //         if (DoNodesShareElement(lattice_site_index, *neighbour_lattice_index))
    //     //         {
    //     //             potts_element_surface -= GetContactAreaBetweenLatticeSite(lattice_site_index,
    //     //                                                                       *neighbour_lattice_index);
    //     //         }
    //     //     }
    // }

    // Now i just need to iterate over the element paris and get the associated lengths

    // PRINT_VARIABLE(potts_element_surface);




// template <unsigned DIM>
// double ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::AreVectorsSame(std::vector< unsigned> Vector1, std::vector< unsigned> Vector2)
// {
//     sort(Vector1.begin(), Vector1.end());
//     sort(Vector2.begin(), Vector2.end());

//     std::vector<unsigned> DifferenceVector;
//     std::set_difference(Vector1.begin(), Vector1.end(), Vector2.begin(), Vector2.end(),
//                         std::inserter(DifferenceVector, DifferenceVector.begin()));


//     double difference = DifferenceVector.size();
//     return difference;
// }

// template <unsigned DIM>
// double ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::AreVectorsSame(std::vector<std::pair<unsigned, unsigned >> Vector1, std::vector<std::pair<unsigned, unsigned >> Vector2)
// {
//     sort(Vector1.begin(), Vector1.end());
//     sort(Vector2.begin(), Vector2.end());

//     std::vector<std::pair<unsigned, unsigned >> DifferenceVector;
//     std::set_difference(Vector1.begin(), Vector1.end(), Vector2.begin(), Vector2.end(),
//                         std::inserter(DifferenceVector, DifferenceVector.begin()));

//     double difference = DifferenceVector.size();
//     return difference;
// }

// template <unsigned DIM>
// std::vector<unsigned> ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::Intersection(std::vector<unsigned> Vector1, std::vector<unsigned> Vector2)
// {
//     std::vector<unsigned> IntersectionVector;
//     for (int i = 0; i < Vector1.size(); i++)
//     {
//         for (int j = 0; j < Vector2.size(); j++)
//         {
//             if (Vector1[i] == Vector2[j])
//             {
//                 IntersectionVector.push_back(Vector1[i]);
//             }
//         }
//     }
//     return IntersectionVector;
// }

// template <unsigned DIM>
// bool ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::IsPairInVector(std::vector<std::pair<unsigned, unsigned> > Vector, std::pair<unsigned, unsigned> Pair)
// {
//     bool IsInVector = 0;
//     for (int i = 0; i < Vector.size(); i++)
//     {
//         if (Vector[i] == Pair)
//         {
//             IsInVector = 1;
//         }
//     }
//     return IsInVector;
// }

// template <unsigned DIM>
// bool ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::IsNumberInVector(std::vector<unsigned> Vector, unsigned number)
// {
//     bool IsInVector = 0;
//     for (int i = 0; i < Vector.size(); i++)
//     {
//         if (Vector[i] == number)
//         {
//             IsInVector = 1;
//         }
//     }
//     return IsInVector;
// }

// template <unsigned DIM>
// bool ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::IsNumberInVector(std::vector<double> Vector, double number)
// {
//     bool IsInVector = 0;
//     for (int i = 0; i < Vector.size(); i++)
//     {
//         if (Vector[i] == number)
//         {
//             IsInVector = 1;
//         }
//     }
//     return IsInVector;
// }


// template <unsigned DIM>
// std::vector<std::pair<unsigned, unsigned> > ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::RemoveInternalEdges(std::vector<std::pair<unsigned, unsigned> > ElementPairsAssociatedWithPottsCell)
// {

//     sort(ElementPairsAssociatedWithPottsCell.begin(), ElementPairsAssociatedWithPottsCell.end());
//     std::vector<std::pair<unsigned, unsigned> > DuplicatesRemoved = ElementPairsAssociatedWithPottsCell;
//     DuplicatesRemoved.erase(unique(DuplicatesRemoved.begin(), DuplicatesRemoved.end()), DuplicatesRemoved.end());

//     std::vector<std::pair<unsigned, unsigned> > InternalPairs;
//     std::set_difference(ElementPairsAssociatedWithPottsCell.begin(), ElementPairsAssociatedWithPottsCell.end(), DuplicatesRemoved.begin(), DuplicatesRemoved.end(),
//                         std::inserter(InternalPairs, InternalPairs.begin()));

//     // I now have the original vector, a vector with the duplicates removed, and a vector with only the repeating elements
//     // I can use the set_difference to find the difference between the InternalPairs and the DuplicatesRemoved vectors
//     std::vector<std::pair<unsigned, unsigned> > ExternalPairs;
//     std::set_difference(DuplicatesRemoved.begin(), DuplicatesRemoved.end(), InternalPairs.begin(), InternalPairs.end(),
//                         std::inserter(ExternalPairs, ExternalPairs.begin()));

//     return ExternalPairs;
// }


// template <unsigned DIM>
// void ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::FindElementNeighbours(MutableMesh<2, DIM>& rMesh)
// {

//     TRACE("Set up neighbouring element pairs");
//     double j = 0;
//     assert(DIM == 3);
//     for (typename MutableMesh<2, DIM>::ElementIterator elem_iter = rMesh.GetElementIteratorBegin();
//          elem_iter != rMesh.GetElementIteratorEnd();
//          ++elem_iter)
//     {
//         unsigned Element_index = elem_iter->GetIndex();

//         // Store the node indices from this element in a vector
//         std::vector<unsigned> NodeIncides;

//         // loop over the indices from this element and save them in the vector NodeIncides
//         for (int i = 0; i < 3; i++)
//         {
//             Node<DIM>* pNode = rMesh.GetNode(rMesh.GetElement(Element_index)->GetNodeGlobalIndex(i));
//             NodeIncides.push_back(pNode->GetIndex());
//         }

//         // Loop over elements to find the elements sharing nodes, indicating they are neighbours.
//         for (typename MutableMesh<2, DIM>::ElementIterator N_elem_iter = rMesh.GetElementIteratorBegin();
//              N_elem_iter != rMesh.GetElementIteratorEnd();
//              ++N_elem_iter)
//         {
//             unsigned Neighbour_Element_index = N_elem_iter->GetIndex();
//             // std::vector<unsigned > TotalNodes;

//             // Store the node indices of the neighbour element in a vector
//             std::vector<unsigned> NeighbourNodeIncides;

//             // loop over the indices from the neighbour element and save them in the vector NeighbourNodeIncides
//             for (int k = 0; k < 3; k++)
//             {
//                 Node<DIM>* pNeighbourNode = rMesh.GetNode(rMesh.GetElement(Neighbour_Element_index)->GetNodeGlobalIndex(k));
//                 NeighbourNodeIncides.push_back(pNeighbourNode->GetIndex());
//                 // TotalNodes.push_back(pNeighbourNode->GetIndex());
//                 // TotalNodes.push_back(NodeIncides[k]);
//             }
//             // Determine if any nodes are in both elements, If the share 2 nodes, they are neighbours, if they share three nodes, it is the same element
//             std::vector<unsigned> IntersectionVector = Intersection(NeighbourNodeIncides, NodeIncides);
//             if (IntersectionVector.size() == 2)
//             {
//                 std::pair<unsigned, unsigned> ElementPair = std::make_pair(std::min(Neighbour_Element_index, Element_index), std::max(Neighbour_Element_index, Element_index));

//                 // Check I havent already counted this pair
//                 if (IsPairInVector(mMeshElementPairs, ElementPair) == 0)
//                 {
//                     mMeshElementPairs.push_back(ElementPair);
//                     mMapElementPairsToNodes[ElementPair] = IntersectionVector;
//                     // mMeshElementNeighbours[Element_index].push_back(Neighbour_Element_index);

//                     // Also save the paired elements for each of the nodes
//                     mMapNodesToAssociateElementPairs[IntersectionVector[0]].push_back(ElementPair);
//                     mMapNodesToAssociateElementPairs[IntersectionVector[1]].push_back(ElementPair);
//                 }
//             }
//         }
//     } // Have now saved all of the mesh element neighbours into member vectors and maps

    // ------------------------------

    // Check that the neighbour element pairs associated with each node is the same as the number of elements contained at the node

//     for (typename MutableMesh<2, DIM>::NodeIterator node_iter = rMesh.GetNodeIteratorBegin();
//          node_iter != rMesh.GetNodeIteratorEnd();
//          ++node_iter)
//     {
//         unsigned node_index = node_iter->GetIndex();
//         std::vector<std::pair<unsigned, unsigned> > VectorOfElementPairs = mMapNodesToAssociateElementPairs[node_index];
//         std::vector<unsigned> ElementsFromPairs;
//         std::vector<unsigned> ElementsFromPairsAlt;

//         // Makes some vector of each of the elements from the pairs, remove the doubles
//         for (int i = 0; i < VectorOfElementPairs.size(); i++)
//         {
//             ElementsFromPairs.push_back(VectorOfElementPairs[i].first);
//             ElementsFromPairs.push_back(VectorOfElementPairs[i].second);
//         }

//         std::set<unsigned>& containing_elements = node_iter->rGetContainingElementIndices();
//         // Vector of the contained elements
//         std::vector<unsigned> ContainingElements(containing_elements.begin(), containing_elements.end());

//         sort(ContainingElements.begin(), ContainingElements.end());

//         sort(ElementsFromPairs.begin(), ElementsFromPairs.end());
//         ElementsFromPairs.erase(unique(ElementsFromPairs.begin(), ElementsFromPairs.end()), ElementsFromPairs.end());

//         // Check if the contained elements and the Assocaited element pairs are the same and throw an error if they are not
//         // The Asscoiated vector pairs will be empty of there is only one element
//         if (AreVectorsSame(ContainingElements, ElementsFromPairs) != 0 && ElementsFromPairs.size() != 0)
//         {
//             EXCEPTION("Elements in elementpairs are not the same as the elements contained at this node");
//         }
//     }
// }


// template <unsigned DIM>
// double ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::GetSurfaceAreaOfElement(unsigned pottsElementIndex, PottsElement<DIM>* pCurrentElement)
// {
//     /*  1) Have the CPM cell the lattice site belongs to
//         2) Loop over the lattice sites in thi CPM cell and collect all of the associated mesh element pairs
//            2a) -  remove anything that appears more than once 
//         3) Loop over all of the associated mesh element pairs and add the perimeter contribution from each. 

//         */
//     //

//     double PottsCellPerimeter = 0;
//     std::vector<std::pair<unsigned, unsigned> > ElementPairsAssociatedWithPottsCell;
//     // loop over the lattice sites in this CPM cell
//     for (unsigned local_lattice_index = 0; local_lattice_index < pCurrentElement->GetNumNodes(); ++local_lattice_index)
//     {
//         unsigned lattice_site_index = pCurrentElement->GetNodeGlobalIndex(local_lattice_index);

//         // Collect the element pairs associate with this node/lattice site.
//         std::vector<std::pair<unsigned, unsigned> > ElementPairsAssociatedWithLatticeSite = mMapNodesToAssociateElementPairs[lattice_site_index];
//         for (int j = 0; j < ElementPairsAssociatedWithLatticeSite.size(); j++)
//         {
//             ElementPairsAssociatedWithPottsCell.push_back(ElementPairsAssociatedWithLatticeSite[j]);
//         }
//     }
    
//     ElementPairsAssociatedWithPottsCell = RemoveInternalEdges(ElementPairsAssociatedWithPottsCell);
//     // PRINT_VARIABLE(ElementPairsAssociatedWithPottsCell.size());

//     // loop over the relevant element pairs and get the edge lenght between then 
//     for(int i = 0; i<ElementPairsAssociatedWithPottsCell.size(); i++ )
//     {
//         PottsCellPerimeter += mDistanceBetweenElements[ElementPairsAssociatedWithPottsCell[i]];
//     }
//     // PRINT_VARIABLE(PottsCellPerimeter);

//     return PottsCellPerimeter;
// }



// template <unsigned DIM>
// void ArbitrarySurfaceAreaConstraintPottsUpdateRule<DIM>::CalculateEdgeLenghts(MutableMesh<2, DIM>& rMesh)
// {
//     //  TRACE("Jess is the best");

//     assert(DIM == 3);

//     // Calculate the Midpoint for each element
//     for (typename MutableMesh<2, DIM>::ElementIterator elem_iter = rMesh.GetElementIteratorBegin();
//          elem_iter != rMesh.GetElementIteratorEnd();
//          ++elem_iter)
//     {
//         // // Calculate the center point for each node
//         unsigned Element_index = elem_iter->GetIndex();
//         c_vector<double, 3> CenterPoint = Create_c_vector(0, 0, 0);
//         c_vector<double, 3> NodeIncides;

//         for (int i = 0; i < 3; i++)
//         {
//             Node<DIM>* pNode = rMesh.GetNode(rMesh.GetElement(Element_index)->GetNodeGlobalIndex(i));
//             CenterPoint += pNode->rGetLocation();
//             NodeIncides[i] = pNode->GetIndex();
//         }
//         mMeshElementMidPoints[Element_index] = CenterPoint / 3;
//     }

//     //---------------

//     // For each Mesh neighbour element pair calculate the lenghts between the center points, remebering this is in 3D and the elements probably dont lie in a 2D plane
//     for (int i = 0; i < mMeshElementPairs.size(); i++)
//     {
//         // Get each of the elements in the elementpair
//         unsigned Element1 = mMeshElementPairs[i].first;
//         unsigned Element2 = mMeshElementPairs[i].second;

//         // Get the two nodes associated with each element
//         std::vector<unsigned> CommonNodes = mMapElementPairsToNodes[mMeshElementPairs[i]];
//         Node<DIM>* pNode0 = rMesh.GetNode(CommonNodes[0]);
//         Node<DIM>* pNode1 = rMesh.GetNode(CommonNodes[1]);

//         // CommonNode[0] will be origin, make poisition relative
//         c_vector<double, 3> PositionVector = pNode1->rGetLocation() - pNode0->rGetLocation();

//         c_vector<double, 3> MidPoint1 = mMeshElementMidPoints[Element1] - pNode0->rGetLocation();
//         c_vector<double, 3> MidPoint2 = mMeshElementMidPoints[Element2] - pNode0->rGetLocation();

//         c_vector<double, 3> b1 = inner_prod(PositionVector, MidPoint1) / inner_prod(PositionVector, PositionVector) * PositionVector - MidPoint1;
//         c_vector<double, 3> b2 = inner_prod(PositionVector, MidPoint2) / inner_prod(PositionVector, PositionVector) * PositionVector - MidPoint2;
//         // PRINT_2_VARIABLES(norm_2(b1), norm_2(b2) );
//         mDistanceBetweenElements[mMeshElementPairs[i]] = norm_2(b1) + norm_2(b2);
//         // TRACE("Stuff set here ");
//     }

//     //---------------

//     // Have iterated over all of the nodes and saved the area of the lattice site, proabaly not necessary
//     for (typename MutableMesh<2, DIM>::NodeIterator node_iter = rMesh.GetNodeIteratorBegin();
//          node_iter != rMesh.GetNodeIteratorEnd();
//          ++node_iter)
//     {

//         double Perimeter = 0;
//         unsigned node_index = node_iter->GetIndex();
//         std::set<unsigned>& containing_elements = node_iter->rGetContainingElementIndices();
//         assert(containing_elements.size() > 0);
//         // Loop over the elements for this node, to find the perimiter of the lines connecting to the midpoint
//         for (std::set<unsigned>::iterator iter = containing_elements.begin();
//              iter != containing_elements.end();
//              ++iter)
//         {
//             unsigned elem_index = *iter;
//             c_vector<c_vector<double, 3>, 2> Vectors;
//             std::vector<unsigned> LocalNodes;
//             double j = 0;

//             // Loops over the nodes in the element and saves them as neighbours
//             for (int i = 0; i < 3; i++)
//             {
//                 Node<DIM>* pNode = rMesh.GetNode(rMesh.GetElement(*iter)->GetNodeGlobalIndex(i));
//                 unsigned node_index_i = pNode->GetIndex();
//                 if (node_index_i != node_index)
//                 {
//                     LocalNodes.push_back(node_index_i);
//                     Vectors[j] = pNode->rGetLocation() - node_iter->rGetLocation();
//                     j += 1;
//                 }
//             }

//             c_vector<double, 3> MidPoint = mMeshElementMidPoints[elem_index] - node_iter->rGetLocation();
//             ;

//             c_vector<double, 3> b1 = inner_prod(Vectors[0], MidPoint) / inner_prod(Vectors[0], Vectors[0]) * Vectors[0] - MidPoint;
//             c_vector<double, 3> b2 = inner_prod(Vectors[1], MidPoint) / inner_prod(Vectors[1], Vectors[1]) * Vectors[1] - MidPoint;
//             Perimeter += norm_2(b1) + norm_2(b2);
//         }
//         // PRINT_VARIABLE(Perimeter);
//         mLatticePerimeters[node_index] = Perimeter;
//     }
// }