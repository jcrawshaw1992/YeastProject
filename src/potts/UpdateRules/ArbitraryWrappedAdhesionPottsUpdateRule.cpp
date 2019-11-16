
// By Jess 

#include "ArbitraryWrappedAdhesionPottsUpdateRule.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"
#include "PottsMesh.hpp"
#include "CPMGTPaseEventHandler.hpp"
#include "Debug.hpp"


template<unsigned DIM>
ArbitraryWrappedAdhesionPottsUpdateRule<DIM>::ArbitraryWrappedAdhesionPottsUpdateRule()
    : AbstractWrappedPottsUpdateRule<DIM>(),
      mCellCellAdhesionEnergyParameter (0.1), // James' Educated guess
      mCellBoundaryAdhesionEnergyParameter (0.2) // James' Educated guess
{
}

template<unsigned DIM>
ArbitraryWrappedAdhesionPottsUpdateRule<DIM>::~ArbitraryWrappedAdhesionPottsUpdateRule()
{
}

template<unsigned DIM>
double ArbitraryWrappedAdhesionPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                unsigned targetNodeIndex,
                                                                WrappedPottsBasedCellPopulation<DIM>& rCellPopulation)
{
     
     
    CPMGTPaseEventHandler::BeginEvent(CPMGTPaseEventHandler::ADHESION);
    PottsArbitrarySurfaceIn3DMesh<DIM>* p_static_cast_potts_mesh = static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*> (&(rCellPopulation.rGetMesh())); // Jess needs to edit to be  PottsArbitrarySurfaceIn3DMesh<DIM>
 
    std::set<unsigned> containing_elements = rCellPopulation.GetNode(currentNodeIndex)->rGetContainingElementIndices(); // CPM element that the Lattice site belongs to 
    std::set<unsigned> new_location_containing_elements = rCellPopulation.GetNode(targetNodeIndex)->rGetContainingElementIndices(); // The CPM element that the target lattice site belongs to 

    bool current_node_contained = !containing_elements.empty(); // Is the current lattice site in a CPM element
    bool target_node_contained = !new_location_containing_elements.empty(); // Is the target lattice site in a CPM element

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
            EXCEPTION("The current node and target node must not be in the same element. Shouldnt reach this line");
        }
    }

    // Iterate over nodes neighbouring the target node to work out the contact energy contribution
    double delta_H = 0.0;

     unsigned target_element = (*new_location_containing_elements.begin());
     unsigned current_element = (*containing_elements.begin());

     unsigned Targets_Sister = NAN; 
     unsigned Current_Sister = NAN; 
 
     if(target_node_contained ) 
     {
        Targets_Sister = rCellPopulation.GetSister(target_element);
     }
     
     if(current_node_contained)
     {
        Current_Sister = rCellPopulation.GetSister(current_element);
     }


     if (current_element == Targets_Sister && rCellPopulation.IsPottsSimulationPeriodic())
     {
         delta_H -= 1000000;
         TRACE("Shouldnt be here ")
     }else{

    
    std::set<unsigned> target_neighbouring_node_indices = rCellPopulation.rGetMesh().GetVonNeumannNeighbouringNodeIndices(targetNodeIndex);
     // TRACE("Check Get Von Neumann Neighbouring Nodes indices ");

    // iterating over the neighbours of the target lattice site 
    for (std::set<unsigned>::iterator iter = target_neighbouring_node_indices.begin();
         iter != target_neighbouring_node_indices.end();
         ++iter)
    {
            std::set<unsigned> neighbouring_node_containing_elements = rCellPopulation.rGetMesh().GetNode(*iter)->rGetContainingElementIndices();
            unsigned neighbouring_elements = (*neighbouring_node_containing_elements.begin());
         
            // Every node must each be in at most one  - Check neighbour not belong to two elements
            assert(neighbouring_node_containing_elements.size() < 2);

            // Calculate the contact perimeter between the target lattice site and the neighbouring site
            double target_lattice_site_contact_area = static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*>(&rCellPopulation.rGetMesh())->GetContactAreaBetweenLatticeSite(targetNodeIndex,*iter);

            // Check neighbour is in an element of not in an element//cell???
            bool neighbouring_node_contained = !neighbouring_node_containing_elements.empty();
            
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
                
                if( target_element != neighbour_element && Targets_Sister != neighbour_element )
                {
                    // The nodes are currently contained in different elements
                    delta_H -= target_lattice_site_contact_area*GetCellCellAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(target_element), rCellPopulation.GetCellUsingLocationIndex(neighbour_element));
                }
                // if(Targets_Sister == neighbour_element )
                // {
                //         delta_H += 100000;
                // }
        
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
                
                if (current_element != neighbour_element  && Current_Sister != neighbour_element  )
                {
                    // The nodes are currently contained in different elements
                    delta_H += target_lattice_site_contact_area*GetCellCellAdhesionEnergy(rCellPopulation.GetCellUsingLocationIndex(current_element),rCellPopulation.GetCellUsingLocationIndex(neighbour_element));
                }
                // if(Current_Sister == neighbour_element )
                // {
                //         delta_H -= 100000;
                // }

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
}

    CPMGTPaseEventHandler::EndEvent(CPMGTPaseEventHandler::ADHESION);

    return delta_H;
}

template<unsigned DIM>
double ArbitraryWrappedAdhesionPottsUpdateRule<DIM>::GetCellCellAdhesionEnergy(CellPtr pCellA, CellPtr pCellB)
{
    return GetCellCellAdhesionEnergyParameter();
}

template<unsigned DIM>
double ArbitraryWrappedAdhesionPottsUpdateRule<DIM>::GetCellBoundaryAdhesionEnergy(CellPtr pCell)
{
    return GetCellBoundaryAdhesionEnergyParameter();
}

template<unsigned DIM>
double ArbitraryWrappedAdhesionPottsUpdateRule<DIM>::GetCellCellAdhesionEnergyParameter()
{
    return mCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
double ArbitraryWrappedAdhesionPottsUpdateRule<DIM>::GetCellBoundaryAdhesionEnergyParameter()
{
    return mCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void ArbitraryWrappedAdhesionPottsUpdateRule<DIM>::SetCellCellAdhesionEnergyParameter(double cellCellAdhesionEnergyParameter)
{
    mCellCellAdhesionEnergyParameter = cellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void ArbitraryWrappedAdhesionPottsUpdateRule<DIM>::SetCellBoundaryAdhesionEnergyParameter(double cellBoundaryAdhesionEnergyParameter)
{
    mCellBoundaryAdhesionEnergyParameter = cellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void ArbitraryWrappedAdhesionPottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CellCellAdhesionEnergyParameter>" << mCellCellAdhesionEnergyParameter << "</CellCellAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<CellBoundaryAdhesionEnergyParameter>" << mCellBoundaryAdhesionEnergyParameter << "</CellBoundaryAdhesionEnergyParameter>\n";

    // Call method on direct parent class
    AbstractWrappedPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

//template class ArbitraryWrappedAdhesionPottsUpdateRule<1>;
template class ArbitraryWrappedAdhesionPottsUpdateRule<2>;
template class ArbitraryWrappedAdhesionPottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
//EXPORT_TEMPLATE_CLASS_SAME_DIMS(ArbitraryWrappedAdhesionPottsUpdateRule)
EXPORT_TEMPLATE_CLASS1(ArbitraryWrappedAdhesionPottsUpdateRule, 2)
EXPORT_TEMPLATE_CLASS1(ArbitraryWrappedAdhesionPottsUpdateRule, 3)
