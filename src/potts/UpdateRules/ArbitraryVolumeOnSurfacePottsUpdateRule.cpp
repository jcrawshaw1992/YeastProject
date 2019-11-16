
#include "ArbitraryVolumeOnSurfacePottsUpdateRule.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"
#include "Debug.hpp"
#include "CPMGTPaseEventHandler.hpp"
#include "Debug.hpp"

template<unsigned DIM>
ArbitraryVolumeOnSurfacePottsUpdateRule<DIM>::ArbitraryVolumeOnSurfacePottsUpdateRule()
    : AbstractWrappedPottsUpdateRule<DIM>(),
      mDeformationEnergyParameter(0.5), // Educated guess
      mMatureCellTargetVolume(16.0) // Defaults to a 4*4 cell size
{
        /// \todo Default values don't apply in 3D.
}



template<unsigned DIM>
ArbitraryVolumeOnSurfacePottsUpdateRule<DIM>::~ArbitraryVolumeOnSurfacePottsUpdateRule()
{
}

template<unsigned DIM>
double ArbitraryVolumeOnSurfacePottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                        unsigned targetNodeIndex,
                                                                        WrappedPottsBasedCellPopulation<DIM>& rCellPopulation)
{
   
     PottsArbitrarySurfaceIn3DMesh<DIM>* p_static_cast_potts_mesh = static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*> (&(rCellPopulation.rGetMesh())); // Jess needs to edit to be  PottsArbitrarySurfaceIn3DMesh<DIM>
 
    CPMGTPaseEventHandler::BeginEvent(CPMGTPaseEventHandler::VOLUME);
    // TRACE("EvaluateHamiltonianContribution");

    double delta_H = 0.0;
    // Element is the cell number/ identifier in the Potts model 
    
   
    std::set<unsigned> containing_elements = p_static_cast_potts_mesh->GetNode(currentNodeIndex)->rGetContainingElementIndices();

    std::set<unsigned> new_location_containing_elements = p_static_cast_potts_mesh->GetNode(targetNodeIndex)->rGetContainingElementIndices();
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

    // Returns the volume of the lattice site -> This has been updated by jess to work on a Node located CPM 
    double volume_target_lattice_site =  p_static_cast_potts_mesh->GetVolumeOfLatticeSite(targetNodeIndex);

    // Need to consider the energy change in the volume change in both cells that will be considered in this change 

    if (current_node_contained) // current node is in an element or is in an element
    {
        unsigned current_element = (*containing_elements.begin());
     
        // Can I get the sister from the map in in the potts population 
        double current_volume;
        // If this is a wrapped mesh, then I need to consider the sister volume as well 
        if(rCellPopulation.IsPottsSimulationPeriodic())
        {
            unsigned Sister = rCellPopulation.GetSister(current_element);
            current_volume = p_static_cast_potts_mesh->GetVolumeOfElement(current_element) +p_static_cast_potts_mesh->GetVolumeOfElement(Sister); //getareaofelement -> basically the sum of all the enclosed lattice sites.
        }else
        {
            current_volume = p_static_cast_potts_mesh->GetVolumeOfElement(current_element) ; //getareaofelement -> basically the sum of all the enclosed lattice sites.    
        }
        double current_volume_difference = current_volume - mMatureCellTargetVolume; // Difference from target 
        double target_volume = current_volume +volume_target_lattice_site;

        // Have checked this, it checks out mathematically -- Jess
        delta_H += mDeformationEnergyParameter*(   (current_volume_difference + volume_target_lattice_site)*(current_volume_difference + volume_target_lattice_site) - current_volume_difference*current_volume_difference);

    }
    if (target_node_contained) // target node is in an element
    {
        unsigned target_element = (*new_location_containing_elements.begin());
        double target_volume;
         if(rCellPopulation.IsPottsSimulationPeriodic())
        {
          unsigned Sister = rCellPopulation.GetSister(target_element);
           target_volume =  p_static_cast_potts_mesh->GetVolumeOfElement(target_element) +p_static_cast_potts_mesh->GetVolumeOfElement(Sister);
        }else
        {
          target_volume =  p_static_cast_potts_mesh->GetVolumeOfElement(target_element); 
        }
        double target_volume_difference = target_volume - mMatureCellTargetVolume;

        delta_H += mDeformationEnergyParameter*((target_volume_difference - volume_target_lattice_site)*(target_volume_difference - volume_target_lattice_site) - target_volume_difference*target_volume_difference);
    }

    CPMGTPaseEventHandler::EndEvent(CPMGTPaseEventHandler::VOLUME);
    return delta_H;
}

template<unsigned DIM>
double ArbitraryVolumeOnSurfacePottsUpdateRule<DIM>::GetDeformationEnergyParameter()
{
    return mDeformationEnergyParameter;
    // TRACE("GetDeformationEnergyParameter");
}

template<unsigned DIM>
void ArbitraryVolumeOnSurfacePottsUpdateRule<DIM>::SetDeformationEnergyParameter(double deformationEnergyParameter)
{
    mDeformationEnergyParameter = deformationEnergyParameter;
    // TRACE("SetDeformationEnergyParameter");
}

template<unsigned DIM>
double ArbitraryVolumeOnSurfacePottsUpdateRule<DIM>::GetMatureCellTargetVolume() const
{
    return mMatureCellTargetVolume;
    // TRACE("GetMatureCellTargetVolume");
}

template<unsigned DIM>
void ArbitraryVolumeOnSurfacePottsUpdateRule<DIM>::SetMatureCellTargetVolume(double matureCellTargetVolume)
{
    assert(matureCellTargetVolume >= 0.0);
    mMatureCellTargetVolume = matureCellTargetVolume;
    // TRACE("SetMatureCellTargetVolume");
}

template<unsigned DIM>
void ArbitraryVolumeOnSurfacePottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    
    *rParamsFile << "\t\t\t<DeformationEnergyParameter>" << mDeformationEnergyParameter << "</DeformationEnergyParameter>\n";
    *rParamsFile << "\t\t\t<MatureCellTargetVolume>" << mMatureCellTargetVolume << "</MatureCellTargetVolume>\n";

    // Call method on direct parent class
    AbstractWrappedPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

//template class ArbitraryVolumeOnSurfacePottsUpdateRule<1>;
template class ArbitraryVolumeOnSurfacePottsUpdateRule<2>;
template class ArbitraryVolumeOnSurfacePottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
//EXPORT_TEMPLATE_CLASS_SAME_DIMS(ArbitraryVolumeOnSurfacePottsUpdateRule)
EXPORT_TEMPLATE_CLASS1(ArbitraryVolumeOnSurfacePottsUpdateRule, 2)
EXPORT_TEMPLATE_CLASS1(ArbitraryVolumeOnSurfacePottsUpdateRule, 3)



// template <unsigned DIM>
// void ArbitraryVolumeOnSurfacePottsUpdateRule<DIM>::CalculateLatticeVolumes(MutableMesh<2, DIM>& rMesh)
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

//     // Have iterated over all of the nodes and saved the area of the lattice site,
//     for (typename MutableMesh<2, DIM>::NodeIterator node_iter = rMesh.GetNodeIteratorBegin();
//          node_iter != rMesh.GetNodeIteratorEnd();
//          ++node_iter)
//     {

//         double LatticeArea = 0;
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

//             // Loops over the nodes in the element to get the two position vectors describing the element form the inital node
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
      
//             c_vector<double, 3> a1 = inner_prod(Vectors[0], MidPoint) / inner_prod(Vectors[0], Vectors[0]) * Vectors[0];
//             c_vector<double, 3> a2 = inner_prod(Vectors[1], MidPoint) / inner_prod(Vectors[1], Vectors[1]) * Vectors[1];
//             double Area1 = 0.5*norm_2(VectorProduct(a1, MidPoint));
//             double Area2 = 0.5*norm_2(VectorProduct(a2, MidPoint));
//             LatticeArea += Area1+Area2;
//         }
//         mLatticeVolume[node_index] = LatticeArea ;
        
//     }



// }




// template<unsigned DIM>
// double PottsArbitrarySurfaceIn3DMesh<DIM>::GetVolumeOfElement(unsigned pottsElementIndex, PottsElement<DIM>* p_potts_element)
// {
//     // TRACE(" UPDATED ------  GetVolumeOfElement");
//     // PottsElement<DIM>* p_potts_element = this->GetElement(pottsElementIndex);

//     double potts_element_volume = 0.0;

//     // An element is made of a number of lattice sites, which are centered around a number of nodes in the Delaunay mesh
//     for (unsigned node_index = 0; node_index < p_potts_element->GetNumNodes(); ++node_index)
//     {
//         potts_element_volume += mLatticeVolume(p_potts_element->GetNodeGlobalIndex(node_index));
//     }

//     return potts_element_volume;
// }


