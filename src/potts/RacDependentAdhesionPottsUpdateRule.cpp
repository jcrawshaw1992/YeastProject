/*

Copyright (c) 2005-2014, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "PottsArbitrarySurfaceIn3DMesh.hpp"
#include "GTPasePDESystemParameters.hpp"
#include "ReplicatableVector.hpp"
#include "PottsElement.hpp"
#include "RandomNumberGenerator.hpp"
#include "Debug.hpp"
#include "RacDependentAdhesionPottsUpdateRule.hpp"

template<unsigned DIM>
RacDependentAdhesionPottsUpdateRule<DIM>::RacDependentAdhesionPottsUpdateRule()
    : AbstractPottsUpdateRule<DIM>(),
	  mLowLowRacAdhesionEnergyParameter(0.02),
	  mHighLowRacAdhesionEnergyParameter(0.02),
	  mHighHighRacAdhesionEnergyParameter(0.02),
	  mCellBoundaryAdhesionEnergyParameter(0.16)
{
}


template<unsigned DIM>
RacDependentAdhesionPottsUpdateRule<DIM>::~RacDependentAdhesionPottsUpdateRule()
{
}

template<unsigned DIM>
double RacDependentAdhesionPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                unsigned targetNodeIndex,
                                                                PottsBasedCellPopulation<DIM>& rCellPopulation)
{
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

    // Iterate over nodes neighbouring the target node to work out the contact energy contribution
    double delta_H = 0.0;
    std::set<unsigned> target_neighbouring_node_indices = rCellPopulation.rGetMesh().GetVonNeumannNeighbouringNodeIndices(targetNodeIndex);
    for (std::set<unsigned>::iterator iter = target_neighbouring_node_indices.begin();
         iter != target_neighbouring_node_indices.end();
         ++iter)
    {
        std::set<unsigned> neighbouring_node_containing_elements = rCellPopulation.rGetMesh().GetNode(*iter)->rGetContainingElementIndices();

        // Every node must each be in at most one element
        assert(neighbouring_node_containing_elements.size() < 2);

        // Calculate the contact area between the target lattice site and the neighbouring site
        double target_lattice_site_contact_area = static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*>(&rCellPopulation.rGetMesh())->GetContactAreaBetweenLatticeSite(targetNodeIndex,*iter);

        bool neighbouring_node_contained = !neighbouring_node_containing_elements.empty();

        /**
         * Before the move, we have a negative contribution (H_0) to the Hamiltonian if:
         * the target node and neighbouring node are NOT contained in the same Potts element;
         * the neighbouring node is contained in a Potts element, but the target node is not; or
         * the target node is contained in a Potts element, but the neighbouring node is not.
         */
        if (neighbouring_node_contained && target_node_contained)
        {
            unsigned neighbour_element = (*neighbouring_node_containing_elements.begin());
            unsigned target_element = (*new_location_containing_elements.begin());
            if (target_element != neighbour_element)
            {
            	double targetNodeRacDependentAdhesionType = GetNodeAdhesionType(targetNodeIndex, target_element, rCellPopulation);
                double neighbourNodeRacDependentAdhesionType = GetNodeAdhesionType(*iter, neighbour_element, rCellPopulation);

            	// The nodes are currently contained in different elements
            	delta_H -= target_lattice_site_contact_area*GetNodeNodeAdhesionEnergy(targetNodeRacDependentAdhesionType, neighbourNodeRacDependentAdhesionType,rCellPopulation);
            }
        }
        else if (neighbouring_node_contained && !target_node_contained)
        {
            // The neighbouring node is contained in a Potts element, but the target node is not
            unsigned neighbour_element = (*neighbouring_node_containing_elements.begin());
            delta_H -= target_lattice_site_contact_area*GetCellBoundaryAdhesionEnergyParameter();
        }
        else if (!neighbouring_node_contained && target_node_contained)
        {
            // The target node is contained in a Potts element, but the neighbouring node is not
            unsigned target_element = (*new_location_containing_elements.begin());
            delta_H -= target_lattice_site_contact_area*GetCellBoundaryAdhesionEnergyParameter();
       }

        /**
         * After the move, we have a positive contribution (H_1) to the Hamiltonian if:
         * the current node and neighbouring node are contained in different Potts elements;
         * the neighbouring node is contained in a Potts element, but the current node is not; or
         * the current node is contained in a Potts element, but the neighbouring node is not.
         */
        if (neighbouring_node_contained && current_node_contained)
        {
            unsigned neighbour_element = (*neighbouring_node_containing_elements.begin());
            unsigned current_element = (*containing_elements.begin());
            if (current_element != neighbour_element)
            {
                double currentNodeRacDependentAdhesionType = GetNodeAdhesionType(currentNodeIndex, current_element, rCellPopulation);
                double neighbourNodeRacDependentAdhesionType = GetNodeAdhesionType(*iter, neighbour_element, rCellPopulation);

            	// The nodes are currently contained in different elements
            	delta_H += target_lattice_site_contact_area*GetNodeNodeAdhesionEnergy(currentNodeRacDependentAdhesionType,neighbourNodeRacDependentAdhesionType,rCellPopulation);
            }
        }
        else if (neighbouring_node_contained && !current_node_contained)
        {
            // The neighbouring node is contained in a Potts element, but the current node is not
            unsigned neighbour_element = (*neighbouring_node_containing_elements.begin());
            delta_H += target_lattice_site_contact_area*GetCellBoundaryAdhesionEnergyParameter();
        }
        else if (!neighbouring_node_contained && current_node_contained)
        {
            // The current node is contained in a Potts element, but the neighbouring node is not
            unsigned current_element = (*containing_elements.begin());
            delta_H += target_lattice_site_contact_area*GetCellBoundaryAdhesionEnergyParameter();
        }
    }

    return delta_H;
}

template<unsigned DIM>
double RacDependentAdhesionPottsUpdateRule<DIM>::GetNodeNodeAdhesionEnergy(double nodeAdhesionTypeA, unsigned nodeAdhesionTypeB, PottsBasedCellPopulation<DIM>& rCellPopulation)
{
	if (nodeAdhesionTypeA+nodeAdhesionTypeB == 0)
	{
		return GetLowLowRacAdhesionEnergyParameter();
	}
	else if (nodeAdhesionTypeA+nodeAdhesionTypeB == 1)
	{
	    return GetHighLowRacAdhesionEnergyParameter();
	}
	else if (nodeAdhesionTypeA+nodeAdhesionTypeB >= 2)
	{
	    return GetHighHighRacAdhesionEnergyParameter();
	}
	else
	{
	    NEVER_REACHED;
	}
}

template<unsigned DIM>
double RacDependentAdhesionPottsUpdateRule<DIM>::GetNodeAdhesionType(unsigned nodeIndex, unsigned pottsElementIndex, PottsBasedCellPopulation<DIM>& rCellPopulation)
{
	//Return the node adhesion type
	double type;

    // Obtain Node Rac concentration
	double rac = GetRacValueForNode(nodeIndex, pottsElementIndex, rCellPopulation);

	//Rac-level-Dependent classification / Either low or high concentration / Threshold level: BasalRac
	if (rac <= BasalRac)
	{
		type=0; // Low Adhesiveness
	}
	else
	{
		type = 1; // High Adhesiveness
	}
	return type;
}

template<unsigned DIM>
double RacDependentAdhesionPottsUpdateRule<DIM>::GetRacValueForNode(unsigned nodeIndex, unsigned pottsElementIndex, PottsBasedCellPopulation<DIM>& rCellPopulation)
{
    CellPtr containing_cell = rCellPopulation.GetCellUsingLocationIndex(pottsElementIndex);

    /// \todo: in the first iteration the modifier hasn't yet populated CellVecData. Return basal value for the moment, needs further investigation, i.e. initialise CellVecData with PDE initial condition.
    if (!containing_cell->HasCellVecData())
    {
        return BasalRac;
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
            std::cout << "WARNING: no value of RAC available for node added in the current CPM sweep" << std::endl;
            return BasalRac;
        }

        offset *= GTPASE_PROBLEM_DIM;
        offset += RAC;

        Vec GTPaseSolution = containing_cell->GetCellVecData()->GetItem("GTPasePDESolution");
        return PetscVecTools::GetElement(GTPaseSolution, offset);
    }

}

template<unsigned DIM>
double RacDependentAdhesionPottsUpdateRule<DIM>::GetLowLowRacAdhesionEnergyParameter()
{
    return mLowLowRacAdhesionEnergyParameter;
}

template<unsigned DIM>
double RacDependentAdhesionPottsUpdateRule<DIM>::GetHighLowRacAdhesionEnergyParameter()
{
    return mHighLowRacAdhesionEnergyParameter;
}

template<unsigned DIM>
double RacDependentAdhesionPottsUpdateRule<DIM>::GetHighHighRacAdhesionEnergyParameter()
{
    return mHighHighRacAdhesionEnergyParameter;
}

template<unsigned DIM>
double RacDependentAdhesionPottsUpdateRule<DIM>::GetCellBoundaryAdhesionEnergyParameter()
{
    return mCellBoundaryAdhesionEnergyParameter;
}


template<unsigned DIM>
void RacDependentAdhesionPottsUpdateRule<DIM>::SetLowLowRacAdhesionEnergyParameter(double lowLowRacAdhesionEnergyParameter)
{
    mLowLowRacAdhesionEnergyParameter = lowLowRacAdhesionEnergyParameter;
}

template<unsigned DIM>
void RacDependentAdhesionPottsUpdateRule<DIM>::SetHighLowRacAdhesionEnergyParameter(double highLowRacAdhesionEnergyParameter)
{
    mHighLowRacAdhesionEnergyParameter = highLowRacAdhesionEnergyParameter;
}

template<unsigned DIM>
void RacDependentAdhesionPottsUpdateRule<DIM>::SetHighHighRacAdhesionEnergyParameter(double highHighRacAdhesionEnergyParameter)
{
    mHighHighRacAdhesionEnergyParameter = highHighRacAdhesionEnergyParameter;
}

template<unsigned DIM>
void RacDependentAdhesionPottsUpdateRule<DIM>::SetCellBoundaryAdhesionEnergyParameter(double cellBoundaryAdhesionEnergyParameter)
{
    mCellBoundaryAdhesionEnergyParameter = cellBoundaryAdhesionEnergyParameter;
}




template<unsigned DIM>
void RacDependentAdhesionPottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
	 *rParamsFile << "\t\t\t<LowLowRacAdhesionEnergyParameter>" << mLowLowRacAdhesionEnergyParameter << "</LowLowRacAdhesionEnergyParameter>\n";
	 *rParamsFile << "\t\t\t<HighLowRacAdhesionEnergyParameter>" << mHighLowRacAdhesionEnergyParameter << "</HighLowRacAdhesionEnergyParameter>\n";
	 *rParamsFile << "\t\t\t<HighHighRacAdhesionEnergyParameter>" << mHighHighRacAdhesionEnergyParameter << "</HighHighRacAdhesionEnergyParameter>\n";
	 *rParamsFile << "\t\t\t<CellBoundaryAdhesionEnergyParameter>" << mCellBoundaryAdhesionEnergyParameter << "</CellBoundaryAdhesionEnergyParameter>\n";


    // Call method on direct parent class
    AbstractPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class RacDependentAdhesionPottsUpdateRule<1>;
template class RacDependentAdhesionPottsUpdateRule<2>;
template class RacDependentAdhesionPottsUpdateRule<3>;

// Serialization for Boost >= 1.36
//#include "SerializationExportWrapper.hpp"
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RacDependentAdhesionPottsUpdateRule)


