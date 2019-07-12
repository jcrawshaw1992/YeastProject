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

#include "ArbitraryAdhesionPottsUpdateRule.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"
#include "GTPasePDESystemParameters.hpp"
#include "ReplicatableVector.hpp"

template<unsigned DIM>
ArbitraryRacDependentAdhesionPottsUpdateRule<DIM>::ArbitraryRacDependentAdhesionPottsUpdateRule()
    : AbstractPottsUpdateRule<DIM>(),
	  mNodeNode_RacHighHigh_AdhesionEnergyParameter (0.02),
	  mLabelledNodeNode_OR_RacHighLow_AdhesionEnergyParameter (0.02),
	  mLabelledNodeLabelledNode_OR_RacLowLow_AdhesionEnergyParameter (0.02),
	  mOne_LabelledNode_OR_RacLow_AdhesionEnergyParameter (0.02),
	  mLabelledNodeLabelledNode_AND_RacLowLow_AdhesionEnergyParameter (0.02)
{
}

template<unsigned DIM>
ArbitraryRacDependentAdhesionPottsUpdateRule<DIM>::~ArbitraryRacDependentAdhesionPottsUpdateRule()
{
}

template<unsigned DIM>
double ArbitraryRacDependentAdhesionPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
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
    	//Need target and current node adhesion types
    	    	double targetNodeRacDependentAdhesionType;
    	    	double currentNodeRacDependentAdhesionType;

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
            	targetNodeRacDependentAdhesionType=GetNodeAdhesionType(targetNodeIndex,rCellPopulation);
            	// The nodes are currently contained in different elements
            	delta_H -= target_lattice_site_contact_area*GetNodeNodeAdhesionEnergy(targetNodeRacDpendentAdhesionType,*iter,rCellPopulation);
            }
        }
        else if (neighbouring_node_contained && !target_node_contained)
        {
            // The neighbouring node is contained in a Potts element, but the target node is not
            unsigned neighbour_element = (*neighbouring_node_containing_elements.begin());
            delta_H -= target_lattice_site_contact_area*GetCellBoundaryAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(neighbour_element));
        }
        else if (!neighbouring_node_contained && target_node_contained)
        {
            // The target node is contained in a Potts element, but the neighbouring node is not
            unsigned target_element = (*new_location_containing_elements.begin());
            delta_H -= target_lattice_site_contact_area*GetCellBoundaryAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(target_element));
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
            	currentNodeRacDependentAdhesionType=GetNodeAdhesionType(targetNodeIndex,rCellPopulation);
            	// The nodes are currently contained in different elements
            	delta_H += target_lattice_site_contact_area*GetNodeNodeAdhesionEnergy(currentNodeRacDpendentAdhesionType,*iter,rCellPopulation);
            }
        }
        else if (neighbouring_node_contained && !current_node_contained)
        {
            // The neighbouring node is contained in a Potts element, but the current node is not
            unsigned neighbour_element = (*neighbouring_node_containing_elements.begin());
            delta_H += target_lattice_site_contact_area*GetCellBoundaryAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(neighbour_element));
        }
        else if (!neighbouring_node_contained && current_node_contained)
        {
            // The current node is contained in a Potts element, but the neighbouring node is not
            unsigned current_element = (*containing_elements.begin());
            delta_H += target_lattice_site_contact_area*GetCellBoundaryAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(current_element));
        }
    }

    return delta_H;
}


template<unsigned DIM>
double ArbitraryRacDependentAdhesionPottsUpdateRule<DIM>::GetNodeNodeAdhesionEnergy(double nodeAdhesionTypeA, unsigned nodeIndex, PottsBasedCellPopulation<DIM>& rCellPopulation)
{
	double nodeAdhesionTypeB=GetNodeAdhesionType(nodeIndex,rCellPopulation);

	double adhesionEnergy=100;

	if (nodeAdhesionTypeA+nodeAdhesionTypeB == 0)
	{
		adhesionEnergy=Get_NodeNode_RacHighHigh_AdhesionEnergyParameter();
	}
	else if (nodeAdhesionTypeA+nodeAdhesionTypeB == 1)
	{
		adhesionEnergy=Get_LabelledNodeNode_OR_RacHighLow_AdhesionEnergyParameter();
	}
	else if (nodeAdhesionTypeA+nodeAdhesionTypeB >= 2)
	{
		adhesionEnergy=Get_LabelledNodeLabelledNode_OR_RacLowLow_AdhesionEnergyParameter();
	}
	else if (nodeAdhesionTypeA+nodeAdhesionTypeB == 3)
	{
		adhesionEnergy=Get_One_LabelledNode_OR_RacLow_AdhesionEnergyParameter();
	}
	else if (nodeAdhesionTypeA+nodeAdhesionTypeB == 4)
	{
		adhesionEnergy=Get_LabelledNodeLabelledNode_AND_RacLowLow_AdhesionEnergyParameter();
	}
	return adhesionEnergy;
}

template<unsigned DIM>
double ArbitraryRacDependentAdhesionPottsUpdateRule<DIM>::GetNodeAdhesionType(unsigned nodeIndex, PottsBasedCellPopulation<DIM>& rCellPopulation)
{
	//Return the node adhesion type
	double type=100;

	//Get the element that contains the node
	std::set<unsigned> containing_elements = rCellPopulation.GetNode(nodeIndex)->rGetContainingElementIndices();
	unsigned int x= (*containing_elements.begin());
	PottsElement<DIM>* p_element = rCellPopulation.GetElement(*containing_elements.begin());

    // Define 2 possible cases for Rac-Concentration: Low and High
	double rac = GetRacValueForNode(currentNodeIndex, rCellPopulation);


	// NECESITO ENTENDER PRIMERO COMO SE PONDERAN LOS PUNTOS EN EL MAPA CARTESIANO DE STEPHANI. Y COMBINAR EL RAC IF SOLO EN LA ZONA QUE ME INTERESA COMPARAR.
	//Get the aspect ratio of the element
	double aspectRatio=(*p_element).GetAspectRatio();

	if(aspectRatio==1)
	{
		type=0;
	}
	else
	{
		//Critical angle is arctan(1/aspectRatio)
		double criticalAngle= atan(1.0/double(aspectRatio));
		//Get the longest eigenvector of the element
		c_vector<double, DIM> eigenvector=(*p_element).GetLongestEigvec();

		//Get the centroid of the element
		c_vector<double, DIM> centroid=rCellPopulation.rGetMesh().GetCentroidOfElement((*containing_elements.begin()));
		//Get location of node in mesh
		Node<DIM>* p_node = rCellPopulation.GetNode(nodeIndex);
		c_vector<double, DIM> nodeLocation = p_node->rGetLocation();

		//Vector between centroid and node
		c_vector<double, DIM> vector=rCellPopulation.rGetMesh().GetVectorFromAtoB(centroid,nodeLocation);

		//Angle between vector and longest eigenvector
		double numerator=vector[0]*eigenvector[0]+vector[1]*eigenvector[1];
		double denominator= sqrt((pow(vector[0],2)+pow(vector[1],2))*(pow(eigenvector[0],2)+pow(eigenvector[1],2)));
		double angle;
		double cosTheta=double(numerator)/double(denominator);

		//If the nodeLocation is in the same place as the centroid,then node is in the center of the shape,
		//so it doesn't matter what the type is, return type 0.
		if (fabs(nodeLocation[0]-centroid[0])<0.01 &&fabs(nodeLocation[1]-centroid[1])<0.01)
		{
			type=0;
		}
		else
		{
			//If ratio is 1 then acos returns nan, but angle should be 0 acos(1)=acos(-1)=0.
			if(fabs(cosTheta-1)<0.001||fabs(cosTheta+1)<0.001)
			{
				angle=0;
			}
			else
			{
				angle = acos(double(cosTheta));
			}

			//We want the angle (currently in the range [0,pi]) mod pi/2.
			if (angle>1.5707)
			{
				angle=3.14-angle;
			}

			//Exit if the angle=nan.
			assert(!std::isnan(angle));

			//Classify the type by the angle.
			if (fabs(angle) < criticalAngle ) //High Adhesiveness
			{
				//Rac-level-Dependent classification / Either low or high concentration / Threshold level: BasalRac
				if (rac <= 3.0) //BasalRac
				{
					type=0;
				}
				else if (rac > 3.0)
				{
					type = 2; //Regardless of the Opposite Node
				}
			}
			else if (fabs(angle) >= criticalAngle) //Low Adhesiveness
			{
				//Rac-level-Dependent classification / Either low or high concentration / Threshold level: BasalRac
				if (rac <= 3.0) //BasalRac
				{
					type=0;
				}
				else if ()
				{
					type = 1;
				}
			}
		}
	}
		return type;
}


template<unsigned DIM>
double ArbitraryRacDependentAdhesionPottsUpdateRule<DIM>::GetRacValueForNode(unsigned nodeIndex, PottsBasedCellPopulation<DIM>& rCellPopulation)
{
    std::set<unsigned> containing_elements = rCellPopulation.GetNode(nodeIndex)->rGetContainingElementIndices();

    /// \todo: this if must go once proper checks are implemented in the calling method
    if (containing_elements.size() == 0)
    {
        return 0.0;// Replace with mRhoThreshold once is define, that way delta_h will evaluate to 0
    }

    assert(containing_elements.size() == 1);
    unsigned current_element = (*containing_elements.begin());

    CellPtr containing_cell = rCellPopulation.GetCellUsingLocationIndex(current_element);

    /// \todo: in the first iteration the modifier hasn't yet populated CellVecData. Return 0 for the moment, needs further investigation, i.e. initialise CellVecData with PDE initial condition.
    if (!containing_cell->HasCellVecData())
    {
        return 0.0;// Replace with mRacThreshold once is define, that way delta_h will evaluate to 0
    }
    else
    {
        PottsElement<DIM>* potts_element = rCellPopulation.GetElement(current_element);

        unsigned offset = 0;
        bool found = false;
        for (unsigned node_number=0; node_number<potts_element->GetNumNodes(); ++node_number)
        {
            if (potts_element->GetNode(node_number)->GetIndex() == nodeIndex)
            {
                found = true;
                break;
            }
            offset++;
        }
        assert(found);

        offset *= GTPASE_PROBLEM_DIM;
        offset += RAC;

        Vec GTPaseSolution = containing_cell->GetCellVecData()->GetItem("GTPasePDESolution");

        PetscInt vec_current_size;
        VecGetSize(GTPaseSolution, &vec_current_size);
        if(offset < (unsigned) vec_current_size)
        {
            /// \todo Using is very inefficient here, get the value directly from the Vec
            ReplicatableVector GTPaseSolutionReplicated(GTPaseSolution);
            return GTPaseSolutionReplicated[offset];
        }
        else
        {
            std::cout << "WARNING: no value of RAC available for node added in the current CPM sweep" << std::endl;
            return BasalRac;
        }
    }

}


template<unsigned DIM>
double ArbitraryRacDependentAdhesionPottsUpdateRule<DIM>::GetNodeNode_RacHighHigh_AdhesionEnergyParameter()
{
    return mNodeNode_RacHighHigh_AdhesionEnergyParameter;
}

template<unsigned DIM>
double ArbitraryRacDependentAdhesionPottsUpdateRule<DIM>::GetLabelledNodeNode_OR_RacHighLow_AdhesionEnergyParameter()
{
    return mLabelledNodeNode_OR_RacHighLow_AdhesionEnergyParameter;
}

template<unsigned DIM>
double ArbitraryRacDependentAdhesionPottsUpdateRule<DIM>::GetLabelledNodeLabelledNode_OR_RacLowLow_AdhesionEnergyParameter()
{
    return mLabelledNodeLabelledNode_OR_RacLowLow_AdhesionEnergyParameter;
}

template<unsigned DIM>
double ArbitraryRacDependentAdhesionPottsUpdateRule<DIM>::GetOne_LabelledNode_OR_RacLow_AdhesionEnergyParameter()
{
    return mOne_LabelledNode_OR_RacLow_AdhesionEnergyParameter;
}

template<unsigned DIM>
double ArbitraryRacDependentAdhesionPottsUpdateRule<DIM>::GetLabelledNodeLabelledNode_AND_RacLowLow_AdhesionEnergyParameter()
{
    return mLabelledNodeLabelledNode_AND_RacLowLow_AdhesionEnergyParameter;
}



template<unsigned DIM>
void ArbitraryRacDependentAdhesionPottsUpdateRule<DIM>::SetNodeNode_RacHighHigh_AdhesionEnergyParameter(double nodeNode_RacHighHigh_AdhesionEnergyParameter)
{
    mNodeNode_RacHighHigh_AdhesionEnergyParameter = nodeNode_RacHighHigh_AdhesionEnergyParameter;
}

template<unsigned DIM>
void ArbitraryRacDependentAdhesionPottsUpdateRule<DIM>::SetLabelledNodeNode_OR_RacHighLow_AdhesionEnergyParameter(double labelledNodeNode_OR_RacHighLow_AdhesionEnergyParameter)
{
    mLabelledNodeNode_OR_RacHighLow_AdhesionEnergyParameter = labelledNodeNode_OR_RacHighLow_AdhesionEnergyParameter;
}

template<unsigned DIM>
void ArbitraryRacDependentAdhesionPottsUpdateRule<DIM>::SetLabelledNodeLabelledNode_OR_RacLowLow_AdhesionEnergyParameter(double labelledNodeLabelledNode_OR_RacLowLow_AdhesionEnergyParameter)
{
    mLabelledNodeLabelledNode_OR_RacLowLow_AdhesionEnergyParameter = labelledNodeLabelledNode_OR_RacLowLow_AdhesionEnergyParameter;
}

template<unsigned DIM>
void ArbitraryRacDependentAdhesionPottsUpdateRule<DIM>::SetOne_LabelledNode_OR_RacLow_AdhesionEnergyParameter(double one_LabelledNode_OR_RacLow_AdhesionEnergyParameter)
{
    mOne_LabelledNode_OR_RacLow_AdhesionEnergyParameter = one_LabelledNode_OR_RacLow_AdhesionEnergyParameter;
}

template<unsigned DIM>
void ArbitraryRacDependentAdhesionPottsUpdateRule<DIM>::SetLabelledNodeLabelledNode_AND_RacLowLow_AdhesionEnergyParameter(double labelledNodeLabelledNode_AND_RacLowLow_AdhesionEnergyParameter)
{
    mLabelledNodeLabelledNode_AND_RacLowLow_AdhesionEnergyParameter = labelledNodeLabelledNode_AND_RacLowLow_AdhesionEnergyParameter;
}



template<unsigned DIM>
void ArbitraryRacDependentAdhesionPottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<NodeNode_RacHighHigh_AdhesionEnergyParameter>" << mNodeNode_RacHighHigh_AdhesionEnergyParameter << "</NodeNode_RacHighHigh_AdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<LabelledNodeNode_OR_RacHighLow_AdhesionEnergyParameter>" << mLabelledNodeNode_OR_RacHighLow_AdhesionEnergyParameter << "</LabelledNodeNode_OR_RacHighLow_AdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<LabelledNodeLabelledNode_OR_RacLowLow_AdhesionEnergyParameter>" << mLabelledNodeLabelledNode_OR_RacLowLow_AdhesionEnergyParameter << "</LabelledNodeLabelledNode_OR_RacLowLow_AdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<One_LabelledNode_OR_RacLow_AdhesionEnergyParameter>" << mOne_LabelledNode_OR_RacLow_AdhesionEnergyParameter << "</One_LabelledNode_OR_RacLow_AdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<LabelledNodeLabelledNode_AND_RacLowLow_AdhesionEnergyParameter>" << mLabelledNodeLabelledNode_AND_RacLowLow_AdhesionEnergyParameter << "</LabelledNodeLabelledNode_AND_RacLowLow_AdhesionEnergyParameter>\n";


    // Call method on direct parent class
    AbstractPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class ArbitraryRacDependentAdhesionPottsUpdateRule<1>;
template class ArbitraryRacDependentAdhesionPottsUpdateRule<2>;
template class ArbitraryRacDependentAdhesionPottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ArbitraryRacDependentAdhesionPottsUpdateRule)
