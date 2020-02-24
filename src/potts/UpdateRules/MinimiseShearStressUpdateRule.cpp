/*

Shear Stress update rule 

*/

#include "MinimiseShearStressUpdateRule.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"

template <unsigned DIM>
MinimiseShearStressUpdateRule<DIM>::MinimiseShearStressUpdateRule()
        : AbstractWrappedPottsUpdateRule<DIM>(),
          mShearMinimisationCorrelationParameter(0.01) // @todo Made up defaults

{
    // TRACE("Shear constructor")
}

template <unsigned DIM>
MinimiseShearStressUpdateRule<DIM>::~MinimiseShearStressUpdateRule()
{
}

template <unsigned DIM>
double MinimiseShearStressUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                                unsigned targetNodeIndex,
                                                                                WrappedPottsBasedCellPopulation<DIM>& rCellPopulation)
{
    double delta_H = 0;
    PottsArbitrarySurfaceIn3DMesh<DIM>* p_static_cast_potts_mesh = static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*>(&(rCellPopulation.rGetMesh())); // Jess needs to edit to be  PottsArbitrarySurfaceIn3DMesh<DIM>

    std::set<unsigned> containing_elements = p_static_cast_potts_mesh->GetNode(currentNodeIndex)->rGetContainingElementIndices(); 
    std::set<unsigned> new_location_containing_elements = p_static_cast_potts_mesh->GetNode(targetNodeIndex)->rGetContainingElementIndices();
    double AverageSS_0;

    bool current_node_contained = !containing_elements.empty();
    bool target_node_contained = !new_location_containing_elements.empty();
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
 
 
 
    //
    
    
    if (current_node_contained) // current node is in an element or is in an element
    {
        unsigned current_element = (*containing_elements.begin());
        PottsElement<DIM>* pCurrentElement = p_static_cast_potts_mesh->GetElement(current_element);

        double SS_0 = norm_2(p_static_cast_potts_mesh->GetTractionOnElement(current_element));
        double NumberOfLattices = pCurrentElement->GetNumNodes();

        double ShearStressOnTarget = norm_2(p_static_cast_potts_mesh->GetTractionOnLattice(targetNodeIndex));

        if(rCellPopulation.IsPottsSimulationPeriodic())
        {
             unsigned Sister_Element = rCellPopulation.GetSister(current_element);
             PottsElement<DIM>* pSisterElement = p_static_cast_potts_mesh->GetElement(Sister_Element);
             SS_0 += norm_2(p_static_cast_potts_mesh->GetTractionOnElement(Sister_Element));
             NumberOfLattices += pSisterElement->GetNumNodes();
        }

        AverageSS_0  = SS_0/NumberOfLattices;
        
        double AverageSS_1 = (SS_0 + ShearStressOnTarget)/(NumberOfLattices +1);

        delta_H += mShearMinimisationCorrelationParameter * (AverageSS_1 - AverageSS_0);
    }
     if (target_node_contained) // current node is in an element or is in an element
    {

        unsigned target_element = (*new_location_containing_elements.begin());
        PottsElement<DIM>* pTargetElement = p_static_cast_potts_mesh->GetElement(target_element);

        double SS_0 = norm_2(p_static_cast_potts_mesh->GetTractionOnElement(target_element));
        double NumberOfLattices = pTargetElement->GetNumNodes();

        double ShearStressOnTarget = norm_2(p_static_cast_potts_mesh->GetTractionOnLattice(targetNodeIndex));

        if(rCellPopulation.IsPottsSimulationPeriodic())
        {
             unsigned Sister_Element = rCellPopulation.GetSister(target_element);
             PottsElement<DIM>* pSisterElement = p_static_cast_potts_mesh->GetElement(Sister_Element);
             SS_0 += norm_2(p_static_cast_potts_mesh->GetTractionOnElement(Sister_Element));
             NumberOfLattices += pSisterElement->GetNumNodes();
        }

        AverageSS_0  = SS_0/NumberOfLattices;
        
        double AverageSS_1 = (SS_0 - ShearStressOnTarget)/(NumberOfLattices -1);
        delta_H -= mShearMinimisationCorrelationParameter * (AverageSS_1 - AverageSS_0);
    }
    return delta_H ;
}


template <unsigned DIM>
double MinimiseShearStressUpdateRule<DIM>::GetShearMinimisationCorrelationParameter()
{
    // TRACE("GetShearMinimisationCorrelationParameter ")
    return mShearMinimisationCorrelationParameter;
}

template <unsigned DIM>
void MinimiseShearStressUpdateRule<DIM>::SetShearMinimisationCorrelationParameter(double ShearMinimisationCorrelationParameter)
{
    // TRACE("SetShearMinimisationCorrelationParameter")
    mShearMinimisationCorrelationParameter = ShearMinimisationCorrelationParameter;
}

template <unsigned DIM>
void MinimiseShearStressUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ShearMinimisationCorrelationParameter>" << mShearMinimisationCorrelationParameter << "</ShearMinimisationCorrelationParameter>\n";

    // Call method on direct parent class
    AbstractWrappedPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class MinimiseShearStressUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MinimiseShearStressUpdateRule<3>)
