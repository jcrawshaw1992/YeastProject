
#include "AspectRatioConstraintPottsUpdateRule.hpp"
#include "Debug.hpp"

template<unsigned DIM>
AspectRatioConstraintPottsUpdateRule<DIM>::AspectRatioConstraintPottsUpdateRule()
    : AbstractPottsUpdateRule<DIM>(),
      mAspectRatioEnergyParameter(0.5), // @todo Made up defaults
      mTargetAspectRatio(40.0)
{
}

template<unsigned DIM>
AspectRatioConstraintPottsUpdateRule<DIM>::~AspectRatioConstraintPottsUpdateRule()
{
}

template<unsigned DIM>
double AspectRatioConstraintPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                        unsigned targetNodeIndex,
                                                                        PottsBasedCellPopulation<DIM>& rCellPopulation)
{
	assert(DIM==2 || DIM==3);
    TRACE("EvaluateHamiltonianContribution");
    DIM==2;

    double delta_H = 0.0;

    std::set<unsigned> containing_elements = rCellPopulation.GetNode(currentNodeIndex)->rGetContainingElementIndices();
    std::set<unsigned> new_location_containing_elements = rCellPopulation.GetNode(targetNodeIndex)->rGetContainingElementIndices();
    TRACE("Get to here");
    bool current_node_contained = !containing_elements.empty();
    bool target_node_contained = !new_location_containing_elements.empty();

    // Every node must each be in at most one element
    assert(new_location_containing_elements.size() < 2);
    TRACE("A");
    if (!current_node_contained && !target_node_contained)
    {
        TRACE("!current_node_contained && !target_node_contained");
        EXCEPTION("At least one of the current node or target node must be in an element.");
    }
TRACE("B");
    if (current_node_contained && target_node_contained)
    {   
        TRACE("current_node_contained && target_node_contained");
        if (*(new_location_containing_elements.begin()) == *(containing_elements.begin()))
        {
            EXCEPTION("The current node and target node must not be in the same element.");
        }
    }



TRACE("c");
    if (current_node_contained) // current node is in an element
    {
        unsigned current_element = (*containing_elements.begin());
        TRACE("here is the problem");
        PottsElement<DIM>* pCurrentElement = rCellPopulation.rGetMesh().GetElement(current_element);
        TRACE("Get to here");
        double current_aspect_ratio = pCurrentElement->GetAspectRatio();
        double current_aspect_ratio_difference = current_aspect_ratio - mTargetAspectRatio;

        // Add node to element here
        Node<DIM>* pTargetNode = rCellPopulation.rGetMesh().GetNode(targetNodeIndex);
		pCurrentElement->AddNode(pTargetNode);

        double current_aspect_ratio_after_switch =  pCurrentElement->GetAspectRatio();
        double current_aspect_ratio_difference_after_switch= current_aspect_ratio_after_switch - mTargetAspectRatio;

        // Remove node after calculating aspect ratio
        pCurrentElement->DeleteNode(pCurrentElement->GetNodeLocalIndex(pTargetNode->GetIndex()));

        delta_H += mAspectRatioEnergyParameter*(current_aspect_ratio_difference_after_switch*current_aspect_ratio_difference_after_switch - current_aspect_ratio_difference*current_aspect_ratio_difference);
    }
    if (target_node_contained) // target node is in an element
    {
        TRACE("D");
        unsigned target_element = (*new_location_containing_elements.begin());
        TRACE("E");
        PottsElement<DIM>* pTargetElement = rCellPopulation.rGetMesh().GetElement(target_element);
        TRACE("F");


        double target_aspect_ratio = pTargetElement->GetAspectRatio();
        
        double target_aspect_ratio_difference = target_aspect_ratio - mTargetAspectRatio;
        double target_aspect_ratio_after_switch;

        // Remove node from element here
        if (pTargetElement->GetNumNodes() > 1)
        {
        	Node<DIM>* pTargetNode = rCellPopulation.rGetMesh().GetNode(targetNodeIndex);
        	pTargetElement->DeleteNode(pTargetElement->GetNodeLocalIndex(pTargetNode->GetIndex()));

        	target_aspect_ratio_after_switch = pTargetElement->GetAspectRatio();

        	// Add node back after calculating aspect ratio
        	pTargetElement->AddNode(pTargetNode);
        }
        else
        {
        	assert(target_aspect_ratio ==1.0);
        	target_aspect_ratio_after_switch = 1.0;
        }
        double target_aspect_ratio_difference_after_switch= target_aspect_ratio_after_switch - mTargetAspectRatio;



        delta_H += mAspectRatioEnergyParameter*(target_aspect_ratio_difference_after_switch*target_aspect_ratio_difference_after_switch - target_aspect_ratio_difference*target_aspect_ratio_difference);
    }

    return delta_H;
}

template<unsigned DIM>
double AspectRatioConstraintPottsUpdateRule<DIM>::GetAspectRatioEnergyParameter()
{
    return mAspectRatioEnergyParameter;
}

template<unsigned DIM>
void AspectRatioConstraintPottsUpdateRule<DIM>::SetAspectRatioEnergyParameter(double aspectRatioEnergyParameter)
{
    mAspectRatioEnergyParameter = aspectRatioEnergyParameter;
}

template<unsigned DIM>
double AspectRatioConstraintPottsUpdateRule<DIM>::GetTargetAspectRatio() const
{
    return mTargetAspectRatio;
}

template<unsigned DIM>
void AspectRatioConstraintPottsUpdateRule<DIM>::SetTargetAspectRatio(double targetAspectRatio)
{
    assert(targetAspectRatio >= 0.0);
    mTargetAspectRatio = targetAspectRatio;
}

template<unsigned DIM>
void AspectRatioConstraintPottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<AspectRatioEnergyParameter>" << mAspectRatioEnergyParameter << "</AspectRatioEnergyParameter>\n";
    *rParamsFile << "\t\t\t<TargetAspectRatio>" << mTargetAspectRatio << "</TargetAspectRatio>\n";

    // Call method on direct parent class
    AbstractPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class AspectRatioConstraintPottsUpdateRule<1>;
template class AspectRatioConstraintPottsUpdateRule<2>;
template class AspectRatioConstraintPottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AspectRatioConstraintPottsUpdateRule)
