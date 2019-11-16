

#include "ArbitraryAdhesionPottsUpdateRule.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"
#include "CPMGTPaseEventHandler.hpp"
#include "Debug.hpp"


template<unsigned DIM>
ArbitraryAdhesionPottsUpdateRule<DIM>::ArbitraryAdhesionPottsUpdateRule()
    : AbstractPottsUpdateRule<DIM>(),
      mCellCellAdhesionEnergyParameter (0.1), // Educated guess
      mCellBoundaryAdhesionEnergyParameter (0.2) // Educated guess
{
}

template<unsigned DIM>
ArbitraryAdhesionPottsUpdateRule<DIM>::~ArbitraryAdhesionPottsUpdateRule()
{
}

template<unsigned DIM>
double ArbitraryAdhesionPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                unsigned targetNodeIndex,
                                                                PottsBasedCellPopulation<DIM>& rCellPopulation)
{
    CPMGTPaseEventHandler::BeginEvent(CPMGTPaseEventHandler::ADHESION);
  
    std::set<unsigned> Element_Of_Source_Lattice = rCellPopulation.GetNode(currentNodeIndex)->rGetContainingElementIndices(); // CPM element that the Lattice site belongs to 
    std::set<unsigned> containing_elements = rCellPopulation.GetNode(currentNodeIndex)->rGetContainingElementIndices(); // CPM element that the Lattice site belongs to 
    std::vector<unsigned> ContainingElementVectors(containing_elements.begin(), containing_elements.end());
    // PRINT_VECTOR(ContainingElementVectors);

    std::set<unsigned> Element_Of_Target_Lattice = rCellPopulation.GetNode(targetNodeIndex)->rGetContainingElementIndices(); // The CPM element that the target lattice site belongs to 
    std::set<unsigned> new_location_containing_elements = rCellPopulation.GetNode(targetNodeIndex)->rGetContainingElementIndices(); // The CPM element that the target lattice site belongs to 

    std::vector<unsigned> Element_Of_Target_LatticeVectors(Element_Of_Target_Lattice.begin(), Element_Of_Target_Lattice.end());
    // PRINT_VECTOR(Element_Of_Target_LatticeVectors);


    bool current_node_contained = !containing_elements.empty(); // Is the current lattice site in a CPM element
    bool target_node_contained = !new_location_containing_elements.empty(); // Is the target lattice site in a CPM element
    // PRINT_2_VARIABLES(current_node_contained,target_node_contained );


    // Every node must each be in at most one element
    assert(new_location_containing_elements.size() < 2);

    if (!current_node_contained && !target_node_contained)
    {
        EXCEPTION("At least one of the current node or target node must be in an element.");
    }

    if (current_node_contained && target_node_contained) // Target node and source node in the same element 
    {
        if (*(new_location_containing_elements.begin()) == *(containing_elements.begin()))
        {
            EXCEPTION("The current node and target node must not be in the same element.");
        }
    }

    // Iterate over nodes neighbouring the target node to work out the contact energy contribution
    double delta_H = 0.0;

    // Get the neighbours??? 
    std::set<unsigned> target_neighbouring_node_indices = rCellPopulation.rGetMesh().GetVonNeumannNeighbouringNodeIndices(targetNodeIndex);
    // std::vector<unsigned> target_neighbouring_node_indicesVectors(target_neighbouring_node_indices.begin(), target_neighbouring_node_indices.end());
    // // PRINT_VECTOR(target_neighbouring_node_indicesVectors);


    // TRACE("Check Get Von Neumann Neighbouring Nodes indices ");

    // iterating over the neighbours of the target lattice site 
    for (std::set<unsigned>::iterator iter = target_neighbouring_node_indices.begin();
         iter != target_neighbouring_node_indices.end();
         ++iter)
    {
        // Get the element of the neighbouring lattice to the target lattice 
        std::set<unsigned> neighbouring_node_containing_elements = rCellPopulation.rGetMesh().GetNode(*iter)->rGetContainingElementIndices();

        // TRACE("Neighbouring lattice's element (nieghbours to the target lattice ");
        // std::vector<unsigned> neighbouring_node_containing_elementsVectors(neighbouring_node_containing_elements.begin(), neighbouring_node_containing_elements.end());
        // PRINT_VECTOR(neighbouring_node_containing_elementsVectors);

        // Every node must each be in at most one  - Check neighbour not belong to two elements
        assert(neighbouring_node_containing_elements.size() < 2);

        // Calculate the contact area between the target lattice site and the neighbouring site

        // This funciton will return the distance between the pair of nodes boarding the neighbouring elements of interest. 
        // If there are no common nodes between the two elements it will return 0 

        /// XXX why would there be no shared nodes between the two elements if the elements are neigjbours                                                // Target Lattice, The neighbour
        double target_lattice_site_contact_area = static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*>(&rCellPopulation.rGetMesh())->GetContactAreaBetweenLatticeSite(targetNodeIndex,*iter);

        // Check neighbour is in an element of not in an element//cell???
        bool neighbouring_node_contained = !neighbouring_node_containing_elements.empty();
        // PRINT_2_VARIABLES(neighbouring_node_contained, target_node_contained);

        /**
         * Before the move, we have a negative contribution (H_0) to the Hamiltonian if:
         * - The target node and neighbouring node are NOT contained in the same Potts element;
         * - The neighbouring node is contained in a Potts element, but the target node is not; or
         * - The target node is contained in a Potts element, but the neighbouring node is not.
         * 
         * delta_H += ContactArea* AdhesionEnergyBetweenCellAandCellB 
         * 
         */

        
        if (neighbouring_node_contained && target_node_contained)
        {
            unsigned neighbour_element = (*neighbouring_node_containing_elements.begin());
            unsigned target_element = (*new_location_containing_elements.begin());
            // PRINT_2_VARIABLES(neighbour_element, target_element);
            if (target_element != neighbour_element)
            {
                // The nodes are currently contained in different elements
                delta_H -= target_lattice_site_contact_area*GetCellCellAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(target_element), rCellPopulation.GetCellUsingLocationIndex(neighbour_element));
                // TRACE("GetCellCellAdhesionEnergy needs updating, I was there to be a difference between cell and matrix")
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
            // PRINT_2_VARIABLES(neighbour_element, current_element);
            if (current_element != neighbour_element)
            {
                // The nodes are currently contained in different elements
                delta_H += target_lattice_site_contact_area*GetCellCellAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(current_element),rCellPopulation.GetCellUsingLocationIndex(neighbour_element));
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

    CPMGTPaseEventHandler::EndEvent(CPMGTPaseEventHandler::ADHESION);

    return delta_H;
}

template<unsigned DIM>
double ArbitraryAdhesionPottsUpdateRule<DIM>::GetCellCellAdhesionEnergy(CellPtr pCellA, CellPtr pCellB)
{
    return GetCellCellAdhesionEnergyParameter();
}

template<unsigned DIM>
double ArbitraryAdhesionPottsUpdateRule<DIM>::GetCellBoundaryAdhesionEnergy(CellPtr pCell)
{
    return GetCellBoundaryAdhesionEnergyParameter();
}

template<unsigned DIM>
double ArbitraryAdhesionPottsUpdateRule<DIM>::GetCellCellAdhesionEnergyParameter()
{
    return mCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
double ArbitraryAdhesionPottsUpdateRule<DIM>::GetCellBoundaryAdhesionEnergyParameter()
{
    return mCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void ArbitraryAdhesionPottsUpdateRule<DIM>::SetCellCellAdhesionEnergyParameter(double cellCellAdhesionEnergyParameter)
{
    mCellCellAdhesionEnergyParameter = cellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void ArbitraryAdhesionPottsUpdateRule<DIM>::SetCellBoundaryAdhesionEnergyParameter(double cellBoundaryAdhesionEnergyParameter)
{
    mCellBoundaryAdhesionEnergyParameter = cellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void ArbitraryAdhesionPottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CellCellAdhesionEnergyParameter>" << mCellCellAdhesionEnergyParameter << "</CellCellAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<CellBoundaryAdhesionEnergyParameter>" << mCellBoundaryAdhesionEnergyParameter << "</CellBoundaryAdhesionEnergyParameter>\n";

    // Call method on direct parent class
    AbstractPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

//template class ArbitraryAdhesionPottsUpdateRule<1>;
template class ArbitraryAdhesionPottsUpdateRule<2>;
template class ArbitraryAdhesionPottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
//EXPORT_TEMPLATE_CLASS_SAME_DIMS(ArbitraryAdhesionPottsUpdateRule)
EXPORT_TEMPLATE_CLASS1(ArbitraryAdhesionPottsUpdateRule, 2)
EXPORT_TEMPLATE_CLASS1(ArbitraryAdhesionPottsUpdateRule, 3)
