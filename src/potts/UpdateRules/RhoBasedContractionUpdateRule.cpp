/*
 * RhoBasedContractionUpdateRule.cpp
 *
 *  Created on: 7 May 2015
 *      Author: mobernabeu
 */

#include "RhoBasedContractionUpdateRule.hpp"
#include "ReplicatableVector.hpp"
#include "GTPasePDESystemParameters.hpp"
#include "PetscVecTools.hpp"
#include "CPMGTPaseEventHandler.hpp"

template<unsigned DIM>
RhoBasedContractionUpdateRule<DIM>::RhoBasedContractionUpdateRule()
: AbstractPottsUpdateRule<DIM>(),
     mDeformationEnergyParameter(1),
     mMatureCellTargetRhoContractionForces(1.25)
{
}

template<unsigned DIM>
RhoBasedContractionUpdateRule<DIM>::~RhoBasedContractionUpdateRule()
{
}

template<unsigned DIM>
double RhoBasedContractionUpdateRule<DIM>::GetRhoValueForNode(unsigned nodeIndex, unsigned pottsElementIndex, PottsBasedCellPopulation<DIM>& rCellPopulation)
{
    CellPtr containing_cell = rCellPopulation.GetCellUsingLocationIndex(pottsElementIndex);

    /// \todo: in the first iteration the modifier hasn't yet populated CellVecData. Return Rho threshold for the moment, needs further investigation, i.e. initialise CellVecData with PDE initial condition.
    if (!containing_cell->HasCellVecData())
    {
        return mMatureCellTargetRhoContractionForces;
    }
    else
    {
        Vec indices = containing_cell->GetCellVecData()->GetItem("GTPasePDEElementIndices");

        unsigned offset = 0;
        bool found = false;

        unsigned num_indices = PetscVecTools::GetSize(indices);
        for (unsigned node_number=0; node_number<num_indices; ++node_number)
        {
            if (PetscVecTools::GetElement(indices, node_number) == nodeIndex)
            {
                found = true;
                break;
            }
            offset++;
        }

        if (!found)
        {
            std::cout << "WARNING: no value of Rho available for node added in the current CPM sweep" << std::endl;
            return BasalRho;
        }

        offset *= GTPASE_PROBLEM_DIM;
        offset += RHO;

        Vec GTPaseSolution = containing_cell->GetCellVecData()->GetItem("GTPasePDESolution");
        return PetscVecTools::GetElement(GTPaseSolution, offset);
    }

}

template<unsigned DIM>
double RhoBasedContractionUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                           unsigned targetNodeIndex,
                                                                           PottsBasedCellPopulation<DIM>& rCellPopulation)
{
    CPMGTPaseEventHandler::BeginEvent(CPMGTPaseEventHandler::RHO_CONTRACTION);

    //    The method is implemented such that we
    //     *select a node A and a Neighbour B we then
    //     *evaluate the hamiltonian as it is (H0)
    //     *and the Hamiltonian with Node A changed to the same spin (cell ID) as node B (H1)
    //    Therefore the target_node is the one we copy the spin to and the current node is the one we copy from.
    //    Hence target_node is A and current node is B (the neighbour)

    //  Look at method VolumeConstraintPottsUpdateRule<DIM>::EvaluateHamiltonianContribution for ways of checking whether a node is contained in a cell, which cell, etc.

    // Equation (16) gives you an expression to evaluate delta_h. delta_h = H1 - H0 according to previous comment.
    // delta_h will have two contributions, one from potential cell retraction and another from a potential cell extension.
    // I say potential because it depends on whether currentNodeIndex/targetNodeIndex belong to a given cell or to the substrate.
    double delta_H = 0.0;

    std::set<unsigned> containing_elements = rCellPopulation.GetNode(currentNodeIndex)->rGetContainingElementIndices();
    std::set<unsigned> new_location_containing_elements = rCellPopulation.GetNode(targetNodeIndex)->rGetContainingElementIndices();

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

    // If currentNodeIndex is contained in a cell, there's a cell extension. Add to delta_h the last term in the first line of equation (16)
    if (current_node_contained)
    {
        unsigned current_element = (*containing_elements.begin());
        double rho = GetRhoValueForNode(currentNodeIndex, current_element, rCellPopulation);
        delta_H += mDeformationEnergyParameter*(rho - mMatureCellTargetRhoContractionForces);
    }

    // If targetNodeIndex is contained in a cell, there's a cell retraction. Add to delta_h the last term in the second line of equation (16)
    if (target_node_contained)
    {
        unsigned target_element = (*new_location_containing_elements.begin());
        double rho = GetRhoValueForNode(targetNodeIndex, target_element, rCellPopulation);
        delta_H -= mDeformationEnergyParameter*(rho - mMatureCellTargetRhoContractionForces);
    }

    CPMGTPaseEventHandler::EndEvent(CPMGTPaseEventHandler::RHO_CONTRACTION);

    return delta_H;
}

template<unsigned DIM>
double RhoBasedContractionUpdateRule<DIM>::GetDeformationEnergyParameter()
{
    return mDeformationEnergyParameter;
}

template<unsigned DIM>
void RhoBasedContractionUpdateRule<DIM>::SetDeformationEnergyParameter(double deformationEnergyParameter)
{
    mDeformationEnergyParameter = deformationEnergyParameter;
}

template<unsigned DIM>
double RhoBasedContractionUpdateRule<DIM>::GetMatureCellTargetRhoContractionForces() const
{
    return mMatureCellTargetRhoContractionForces;
}

template<unsigned DIM>
void RhoBasedContractionUpdateRule<DIM>::SetMatureCellTargetRhoContractionForces(double matureCellTargetRhoContractionForces)
{
    assert(matureCellTargetRhoContractionForces >= 0.0);
    mMatureCellTargetRhoContractionForces = matureCellTargetRhoContractionForces;
}

template<unsigned DIM>
void RhoBasedContractionUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
     *rParamsFile << "\t\t\t<DeformationEnergyParameter>" << mDeformationEnergyParameter << "</DeformationEnergyParameter>\n";
     *rParamsFile << "\t\t\t<MatureCellTargetRhoContractionForces>" << mMatureCellTargetRhoContractionForces << "</MatureCellTargetRhoContractionForcesdx>\n";

     // Call method on direct parent class
     AbstractPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class RhoBasedContractionUpdateRule<1>;
template class RhoBasedContractionUpdateRule<2>;
template class RhoBasedContractionUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RhoBasedContractionUpdateRule)
