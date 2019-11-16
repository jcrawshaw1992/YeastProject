

#include "ArbitraryCurvatureConstraintPottsUpdateRule.hpp"
#include "CPMGTPaseEventHandler.hpp"
#include "Debug.hpp"


template <unsigned DIM>
ArbitraryCurvatureConstraintPottsUpdateRule<DIM>::ArbitraryCurvatureConstraintPottsUpdateRule()
        : AbstractPottsUpdateRule<DIM>(),
          mCurvatureEnergyParameter(0.5), // Educated guess
          mTargetCurvature(16.0) // Defaults to a 4*4 cell size

{
    /// \todo Default values don't apply in 3D.
}

template <unsigned DIM>
ArbitraryCurvatureConstraintPottsUpdateRule<DIM>::~ArbitraryCurvatureConstraintPottsUpdateRule()
{
}
//********



//---------------------

template <unsigned DIM>
double ArbitraryCurvatureConstraintPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
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

    
    if (current_node_contained) // current node is in an element
    {
        unsigned current_element = (*containing_elements.begin());
        PottsElement<DIM>* pCurrentElement = p_static_cast_potts_mesh->GetElement(current_element);


        double current_surface_area = p_static_cast_potts_mesh->GetCurvatureOfElement(current_element);
        // PRINT_VARIABLE(current_surface_area);
        // GetCurvatureOfElement(current_element, pCurrentElement);
        
        double current_surface_area_difference = current_surface_area - mTargetCurvature;
        // PRINT_VARIABLE(current_surface_area_difference);

        // Mock change of spin to call GetCurvatureOfElement with new configuration
        Node<DIM>* pTargetNode = p_static_cast_potts_mesh->GetNode(targetNodeIndex);
        pCurrentElement->AddNode(pTargetNode);
        if (target_node_contained)
        {
            unsigned target_element = (*new_location_containing_elements.begin());
            PottsElement<DIM>* pTargetElement = p_static_cast_potts_mesh->GetElement(target_element);
            pTargetElement->DeleteNode(pTargetElement->GetNodeLocalIndex(pTargetNode->GetIndex()));
        }

        double current_surface_area_after_switch = p_static_cast_potts_mesh->GetCurvatureOfElement(current_element);//, pCurrentElement);
        
        double current_surface_area_difference_after_switch = current_surface_area_after_switch - mTargetCurvature;
        // PRINT_2_VARIABLES(current_surface_area_difference_after_switch, current_surface_area_after_switch);
        // Undo mocked change
        pCurrentElement->DeleteNode(pCurrentElement->GetNodeLocalIndex(pTargetNode->GetIndex()));
        if (target_node_contained)
        {
            unsigned target_element = (*new_location_containing_elements.begin());
            PottsElement<DIM>* pTargetElement = p_static_cast_potts_mesh->GetElement(target_element);
            pTargetElement->AddNode(pTargetNode);
        }

        delta_H += mCurvatureEnergyParameter * (current_surface_area_difference_after_switch * current_surface_area_difference_after_switch - current_surface_area_difference * current_surface_area_difference);
    }
    if (target_node_contained) // target node is in an element
    {
        unsigned target_element = (*new_location_containing_elements.begin());
        PottsElement<DIM>* pTargetElement = p_static_cast_potts_mesh->GetElement(target_element);

        double target_surface_area = p_static_cast_potts_mesh->GetCurvatureOfElement(target_element);//, pTargetElement);
        // GetCurvatureOfElement(target_element, pTargetElement);// p_static_cast_potts_mesh->GetCurvatureOfElement(target_element);
        double target_surface_area_difference = target_surface_area - mTargetCurvature;

        Node<DIM>* pTargetNode = p_static_cast_potts_mesh->GetNode(targetNodeIndex);
        pTargetElement->DeleteNode(pTargetElement->GetNodeLocalIndex(pTargetNode->GetIndex()));

        double target_surface_area_after_switch = p_static_cast_potts_mesh->GetCurvatureOfElement(target_element);//, pTargetElement);
        // GetCurvatureOfElement(target_element, pTargetElement); 
        double target_surface_area_difference_after_switch = target_surface_area_after_switch - mTargetCurvature;

        // Add node back after calculating aspect ratio
        pTargetElement->AddNode(pTargetNode);

        delta_H += mCurvatureEnergyParameter * (target_surface_area_difference_after_switch * target_surface_area_difference_after_switch - target_surface_area_difference * target_surface_area_difference);
    }

    CPMGTPaseEventHandler::EndEvent(CPMGTPaseEventHandler::SURFACE);

    return delta_H;

}

template <unsigned DIM>
double ArbitraryCurvatureConstraintPottsUpdateRule<DIM>::GetCurvatureEnergyParameter()
{
    return mCurvatureEnergyParameter;
}

template <unsigned DIM>
void ArbitraryCurvatureConstraintPottsUpdateRule<DIM>::SetCurvatureEnergyParameter(double CurvatureEnergyParameter)
{
    mCurvatureEnergyParameter = CurvatureEnergyParameter;
}



template <unsigned DIM>
double ArbitraryCurvatureConstraintPottsUpdateRule<DIM>::GetTargetCurvature() const
{
    return mTargetCurvature;
}

template <unsigned DIM>
void ArbitraryCurvatureConstraintPottsUpdateRule<DIM>::SetTargetCurvature(double targetCurvature)
{
    assert(targetCurvature >= 0.0);
    mTargetCurvature = targetCurvature;
}

   
template <unsigned DIM>
void ArbitraryCurvatureConstraintPottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CurvatureEnergyParameter >" << mCurvatureEnergyParameter << "</CurvatureEnergyParameter >\n";
    *rParamsFile << "\t\t\t<TargetCurvature>" << mTargetCurvature << "</TargetCurvature>\n";

    // Call method on direct parent class
    AbstractPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

//template class ArbitraryVolumeConstraintPottsUpdateRule<1>;
template class ArbitraryCurvatureConstraintPottsUpdateRule<2>;
template class ArbitraryCurvatureConstraintPottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
//EXPORT_TEMPLATE_CLASS_SAME_DIMS(ArbitraryVolumeConstraintPottsUpdateRule)
EXPORT_TEMPLATE_CLASS1(ArbitraryCurvatureConstraintPottsUpdateRule, 2)
EXPORT_TEMPLATE_CLASS1(ArbitraryCurvatureConstraintPottsUpdateRule, 3)
