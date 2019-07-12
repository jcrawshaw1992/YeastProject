/*
 * RhoBasedContractionUpdateRule.cpp
 *
 *  Created on: 7 May 2015
 *      Author: mobernabeu
 */

#include "BarbedEndsUpdateRule.hpp"
#include "PetscVecTools.hpp"
#include "GTPasePDESystemParameters.hpp"
#include "UblasCustomFunctions.hpp"
#include "CPMGTPaseEventHandler.hpp"

template<unsigned DIM>
BarbedEndsUpdateRule<DIM>::BarbedEndsUpdateRule()
: AbstractPottsUpdateRule<DIM>(),
     mDeformationEnergyParameter(1)
{
    for (unsigned direction = P0; direction <= P4; ++direction)
    {
        unsigned direction_indexed_from_0 = direction - P0;
        mFilamentDirectionVector.push_back(Create_c_vector(cos(ANGLES[direction_indexed_from_0]),
                                                           sin(ANGLES[direction_indexed_from_0])));
    }
}

template<unsigned DIM>
BarbedEndsUpdateRule<DIM>::~BarbedEndsUpdateRule()
{
}

template<unsigned DIM>
double BarbedEndsUpdateRule<DIM>::GetBarbedEndsPushingAgainstInterfaceForNode(unsigned nodeIndex, unsigned neignhourIndex, unsigned nodePottsElementIndex, PottsBasedCellPopulation<DIM>& rCellPopulation)
{
    /// Compute interface normal from nodeIndex and neighbourIndex
    c_vector<double, DIM> interface_outward_normal = rCellPopulation.GetNode(neignhourIndex)->rGetLocation() - rCellPopulation.GetNode(nodeIndex)->rGetLocation();
    interface_outward_normal /= norm_2(interface_outward_normal);

    CellPtr containing_cell = rCellPopulation.GetCellUsingLocationIndex(nodePottsElementIndex);

    /// \todo: in the first iteration the modifier hasn't yet populated CellVecData. Return 0 for the moment, needs further investigation, i.e. initialise CellVecData with PDE initial condition.
    if (!containing_cell->HasCellVecData())
    {
        return 0.0;
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
            // No value of Pushing Barbed Ends (Pi) available means that the node was added in the current CPM sweep. Pi will go up as Bi grows.
            std::cout << "WARNING: no value of Pi available for node added in the current CPM sweep" << std::endl;
            return 0.0;
        }

        offset *= GTPASE_PROBLEM_DIM;

        Vec GTPaseSolution = containing_cell->GetCellVecData()->GetItem("GTPasePDESolution");

        double pushing_barbed_ends=0;

        for (unsigned direction = P0; direction <= P4; ++direction)
        {
            // Compute dot product of interface_outward_normal and corresponding GTPasePDESystemSolver<ELEMENT_DIM, SPACE_DIM>::ANGLES[]
            unsigned direction_indexed_from_0 = direction - P0;
            double interface_barbed_direction_dot = inner_prod(interface_outward_normal,
                                                               mFilamentDirectionVector.at(direction_indexed_from_0));

            // If the dot product is positive, barbed ends in the current direction are pushing against the interface
            if (interface_barbed_direction_dot > 0.0)
            {
                double pushing_barbed_ends_in_current_direction = PetscVecTools::GetElement(GTPaseSolution, offset+direction);
                assert(pushing_barbed_ends_in_current_direction >= -1e16);

                // Weight contribution by dot product to take into account the angle of incidence
                assert(interface_barbed_direction_dot <= 1.0);
                pushing_barbed_ends += interface_barbed_direction_dot * pushing_barbed_ends_in_current_direction;
            }
        }

        return pushing_barbed_ends;
    }
}

template<unsigned DIM>
double BarbedEndsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                  unsigned targetNodeIndex,
                                                                  PottsBasedCellPopulation<DIM>& rCellPopulation)
{
    CPMGTPaseEventHandler::BeginEvent(CPMGTPaseEventHandler::BARBED_ENDS);

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

    // If currentNodeIndex is contained in a cell, there's a cell extension. Add to delta_h the second term in the first line of equation (16)
    if (current_node_contained)
    {
        unsigned current_element = (*containing_elements.begin());
        double barbedEnd = GetBarbedEndsPushingAgainstInterfaceForNode(currentNodeIndex, targetNodeIndex, current_element, rCellPopulation);
        delta_H -= barbedEnd;
    }

    // If targetNodeIndex is contained in a cell, there's a cell retraction. Add to delta_h the second term in the first line of equation (16)
    if (target_node_contained)
    {
        unsigned target_element = (*new_location_containing_elements.begin());
        double barbedEnd = GetBarbedEndsPushingAgainstInterfaceForNode(targetNodeIndex, currentNodeIndex, target_element, rCellPopulation);
        delta_H += barbedEnd;
   }

    CPMGTPaseEventHandler::EndEvent(CPMGTPaseEventHandler::BARBED_ENDS);

   return mDeformationEnergyParameter * delta_H;
}

template<unsigned DIM>
double BarbedEndsUpdateRule<DIM>::GetDeformationEnergyParameter()
{
    return mDeformationEnergyParameter;
}

template<unsigned DIM>
void BarbedEndsUpdateRule<DIM>::SetDeformationEnergyParameter(double deformationEnergyParameter)
{
    mDeformationEnergyParameter = deformationEnergyParameter;
}

template<unsigned DIM>
void BarbedEndsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<DeformationEnergyParameter>" << mDeformationEnergyParameter << "</DeformationEnergyParameter>\n";

    // Call method on direct parent class
    AbstractPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class BarbedEndsUpdateRule<1>;
template class BarbedEndsUpdateRule<2>;
template class BarbedEndsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BarbedEndsUpdateRule)
