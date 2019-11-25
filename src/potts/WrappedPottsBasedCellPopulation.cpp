/*

Jess's WrappedPottsBasedCellPopulation

unsigned Sister = rCellPopulation.GetSister(current_element)
*/

#include "WrappedPottsBasedCellPopulation.hpp"
#include "AbstractPottsUpdateRule.hpp"
#include "AbstractWrappedPottsUpdateRule.hpp"
#include "CellIdWriter.hpp"
#include "CellPopulationElementWriter.hpp"
#include "Debug.hpp"
#include "NodesOnlyMesh.hpp"
#include "RandomNumberGenerator.hpp"

// Needed to convert mesh in order to write nodes to VTK (visualize as glyphs)
#include "VtkMeshWriter.hpp"




// Constructor for the periodic mesh 
template <unsigned DIM>
WrappedPottsBasedCellPopulation<DIM>::WrappedPottsBasedCellPopulation(PottsMesh<DIM>& rMesh,
                                                                      std::vector<CellPtr>& rCells,
                                                                      std::map<unsigned, unsigned> ElementPairing,
                                                                      bool deleteMesh,
                                                                      bool validate,
                                                                      const std::vector<unsigned> locationIndices)
        : PottsBasedCellPopulation<DIM>(rMesh, rCells, deleteMesh, validate, locationIndices),
          mElementPairing(ElementPairing)
{
    assert(DIM == 3); 
    int counter = 0;
    double PairingNumber = 1;
    mIsPeriodic = 1;

    for (std::map<unsigned, unsigned>::iterator it = ElementPairing.begin(); it != ElementPairing.end(); it++) // The way the map is constructed is such that it is doubled every second so i connected cell x cell y and i+1 connected cell y cell x
    {
        //  std::cout << it->first << '\t' << it->second << std::endl;

        if (counter % 2 == 0)
        {

            PottsElement<DIM>* p_element_1 = this->GetElement(it->first);

            PottsElement<DIM>* p_element_2 = this->GetElement(it->second);

            c_vector<double, DIM> Location1 = p_element_1->GetNodeLocation(0);
            c_vector<double, DIM> Location2 = p_element_2->GetNodeLocation(0);
            // I want the left one to be the first entry
            if (Location1[0] < Location2[0])
            {
                mElementPairingVector.push_back(Create_c_vector(it->first, it->second));
            }
            else
            {
                mElementPairingVector.push_back(Create_c_vector(it->second, it->first));
            }

            CellPtr pCell1 = this->GetCellUsingLocationIndex(it->first);
            CellPtr pCell2 = this->GetCellUsingLocationIndex(it->second);

            pCell1->GetCellData()->SetItem("CellPair", PairingNumber);
            pCell2->GetCellData()->SetItem("CellPair", PairingNumber);
            PairingNumber += 1;
        }
        counter += 1;
    }

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {

        unsigned elem_index = this->GetLocationIndexUsingCell(*cell_iter);
        cell_iter->GetCellData()->SetItem("Sister", ElementPairing[elem_index]);

        CellPtr p_cell = this->GetCellUsingLocationIndex(ElementPairing[elem_index]);
        PottsElement<DIM>* p_element1 = this->GetElementCorrespondingToCell(*cell_iter);
        PottsElement<DIM>* p_element2 = this->GetElementCorrespondingToCell(p_cell);

        c_vector<double, DIM> Midpoint = (p_element1->GetNodeLocation(0) + p_element2->GetNodeLocation(0)) / 2;

        cell_iter->GetCellData()->SetItem("Center", Midpoint[0]);
        p_cell->GetCellData()->SetItem("Center", Midpoint[0]);

        if (p_element1->GetNodeLocation(0)[0] < Midpoint[0])
        {
            // On the left
            cell_iter->GetCellData()->SetItem("DivisionStatus", 0);
        }
        else
        {
            // On the right
            cell_iter->GetCellData()->SetItem("DivisionStatus", 1);
        }
        // need to get element
    }
}







// Constructor for the non periodic mesh 
template <unsigned DIM>
WrappedPottsBasedCellPopulation<DIM>::WrappedPottsBasedCellPopulation(PottsMesh<DIM>& rMesh,
                                                                      std::vector<CellPtr>& rCells,
                                                                      std::vector<unsigned> BoundaryVector,
                                                                      bool deleteMesh,
                                                                      bool validate,
                                                                      const std::vector<unsigned> locationIndices)
        : PottsBasedCellPopulation<DIM>(rMesh, rCells, deleteMesh, validate, locationIndices),
        mBoundaryVector(BoundaryVector)
{
     mIsPeriodic = 0;
     assert(DIM == 3); 
    //      for (std::vector<unsigned>::iterator i = mBoundaryVector.begin(); i != mBoundaryVector.end(); ++i)
    // {
    //     PRINT_VARIABLE(*i)
    // }
}


template <unsigned DIM>
bool WrappedPottsBasedCellPopulation<DIM>::IsPottsSimulationPeriodic()
{
    return mIsPeriodic;
}

template <unsigned DIM>
WrappedPottsBasedCellPopulation<DIM>::~WrappedPottsBasedCellPopulation()
{
    // delete mpElementTessellation;

    // delete mpMutableMesh;

    if (this->mDeleteMesh)
    {
        delete &this->mrMesh;
    }
}

template <unsigned DIM>
CellPtr WrappedPottsBasedCellPopulation<DIM>::DivideCell(CellPtr pNewCell, CellPtr pParentCell)
{
    // Get the element associated with this cell
    // PottsArbitrarySurfaceIn3DMesh M;

    PottsElement<DIM>* p_element = this->GetElementCorrespondingToCell(pParentCell);

    //-- base this on curvature for the 2d in 3d simulations

    c_vector<double, DIM> Midpoint = zero_vector<double>(DIM);

    for (unsigned i = 0; i < p_element->GetNumNodes(); i++)
    {
        Midpoint += p_element->GetNodeLocation(i);
    }
    Midpoint /= p_element->GetNumNodes();

    // now need to iterate over the Lattice sites. If they are on the right they need to be saved in a vector
    std::vector<unsigned> RightLatticeSites;
    std::vector<Node<DIM>*> nodes_elem;
    for (unsigned i = p_element->GetNumNodes() - 1; i < -1; i--) // for some reason this needs to go backwards
    {
        if (p_element->GetNodeLocation(i)[0] > Midpoint[0])
        {
            // Checks if the lattice site is on the right, if so this lattice will be added to the daugter cell, and deleted from the mother in the comming for loop
            RightLatticeSites.push_back(i);
            nodes_elem.push_back(p_element->GetNode(i));
        }
    }

    unsigned new_element_index = this->GetNumElements();
    this->rGetMesh().AddElement(new PottsElement<DIM>(new_element_index, nodes_elem));

    // Now I need to go through and remove everything in this vector from the old cell and add it to the new one
    for (std::vector<unsigned>::iterator i = RightLatticeSites.begin(); i != RightLatticeSites.end(); ++i)
    {
        p_element->DeleteNode(*i);
    }

    // Associate the new cell with the element
    this->mCells.push_back(pNewCell);

    // Update location cell map
    CellPtr p_created_cell = this->mCells.back();
    this->SetCellUsingLocationIndex(new_element_index, p_created_cell);
    pNewCell->GetCellData()->SetItem("Center", Midpoint[0]);
    pParentCell->GetCellData()->SetItem("Center", Midpoint[0]);

    //  Link the sister cells
    unsigned Mother_element_index = p_element->GetIndex();

    pNewCell->GetCellData()->SetItem("Sister", Mother_element_index);
    pParentCell->GetCellData()->SetItem("Sister", new_element_index);

    // 2 for daughter,
    // 1 for mother
    pNewCell->GetCellData()->SetItem("DivisionStatus", 2);
    pParentCell->GetCellData()->SetItem("DivisionStatus", 1);

    TRACE("A success")
    return p_created_cell;
}

// Jess wrote this
template <unsigned DIM>
unsigned WrappedPottsBasedCellPopulation<DIM>::GetSister(unsigned Element_Index)
{
    // Take an element index, returns element index
    return mElementPairing[Element_Index];
}

// Jess wrote this
template <unsigned DIM>
std::map<unsigned, unsigned> WrappedPottsBasedCellPopulation<DIM>::GetElementPairingMap()
{
    // Returns the mapping between the paired elements
    return mElementPairing;
}

template <unsigned DIM>
std::vector<c_vector<unsigned, 2> > WrappedPottsBasedCellPopulation<DIM>::GetElementPairingVector()
{
    // Returns the mapping between the paired elements
    return mElementPairingVector;
}

// Jess wrote this
// Returns whether is divided or not
template <unsigned DIM>
unsigned WrappedPottsBasedCellPopulation<DIM>::GetDivisionStatus(unsigned Element_Index)
{
    return this->GetCellUsingLocationIndex(Element_Index)->GetCellData()->GetItem("DivisionStatus");
}

template <unsigned DIM>
double WrappedPottsBasedCellPopulation<DIM>::GetNumSweepsPerTimestep()
{
    return mNumSweepsPerTimestep;
}

template <unsigned DIM>
void WrappedPottsBasedCellPopulation<DIM>::SetNumSweepsPerTimestep(double numSweepsPerTimestep)
{
    mNumSweepsPerTimestep = numSweepsPerTimestep;
}



template <unsigned DIM>
void WrappedPottsBasedCellPopulation<DIM>::UpdateCellLocations(double dt)
{
    assert(DIM == 3); 
    /*
     * This method implements a Monte Carlo method to update the cell population.
     * We sample randomly from all nodes in the mesh. Once we have selected a target
     * node we randomly select a neighbour. The Hamiltonian is evaluated in the
     * current configuration (H_0) and with the target node added to the same
     * element as the neighbour (H_1). Based on the vale of deltaH = H_1 - H_0,
     * the switch is either made or not.
     *
     * For each time step (i.e. each time this method is called) we sample
     * mrMesh.GetNumNodes() nodes. This is known as a Monte Carlo Step (MCS).
     * 
     * If I need to find this class again here is a bookmark to find 
     * Potts Core
     */
    // TRACE("A")
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    unsigned num_nodes = this->mrMesh.GetNumNodes();
    unsigned num_elems = this->GetNumElements();
//  TRACE("B")
    // Randomly permute mUpdateRuleCollection if specified
    if (this->mIterateRandomlyOverUpdateRuleCollection)
    {
        // Randomly permute mUpdateRuleCollection
        p_gen->Shuffle(this->mUpdateRuleCollection);
    }

    for (unsigned i = 0; i < num_nodes * this->GetNumSweepsPerTimestep(); i++)
    {
        // unsigned NodesInElement =0;
        unsigned node_index = p_gen->randMod(num_nodes);
        Node<DIM>* p_node = this->mrMesh.GetNode(node_index);
        

        std::vector<unsigned>::iterator Current_iter = mBoundaryVector.begin();
        std::advance(Current_iter , node_index);
      
        if(*Current_iter==0)
        {
        // Each node in the mesh must be in at most one element
        assert(p_node->GetNumContainingElements() <= 1);
        
        // Find a random available neighbouring node to overwrite current site
        std::set<unsigned> neighbouring_node_indices = this->rGetMesh().GetMooreNeighbouringNodeIndices(node_index);
        
        unsigned neighbour_location_index;
        std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();

        //  NumberOfNodesInElement = this->GetElement(*containing_elements.begin())->GetNumNodes(); // If there is only one node in the element, I dont want to try and eat it up into the neighbour

        if ((!containing_elements.empty() && this->GetElement(*containing_elements.begin())->GetNumNodes() > 1) || containing_elements.empty() == 1)
        { // if not empty and has more than one lattice site OR is empty

            unsigned num_neighbours = neighbouring_node_indices.size();
            unsigned chosen_neighbour = p_gen->randMod(num_neighbours);

            std::set<unsigned>::iterator neighbour_iter = neighbouring_node_indices.begin();

            std::advance(neighbour_iter, chosen_neighbour);

            neighbour_location_index = *neighbour_iter;
            
            std::set<unsigned> neighbour_containing_elements = this->GetNode(*neighbour_iter)->rGetContainingElementIndices();

            unsigned Current_Sister = NAN; 
             if(!containing_elements.empty() && mIsPeriodic)
            {
                Current_Sister = this->GetSister(*containing_elements.begin());
            }

            // Only calculate Hamiltonian and update elements if the nodes are from different elements, or one is from the medium
            if ((!containing_elements.empty() && neighbour_containing_elements.empty())
                || (containing_elements.empty() && !neighbour_containing_elements.empty())
                || (!containing_elements.empty() && !neighbour_containing_elements.empty() && *containing_elements.begin() != *neighbour_containing_elements.begin() && Current_Sister != *neighbour_containing_elements.begin()))
            {
                //=====
                unsigned current_Target_element = *(neighbour_containing_elements.begin());
            
              
                
                bool RunPeriodicTrialSwitch =1;
                bool NotOnBoundary =1;
                if ( mIsPeriodic ==1)
                {   
                    unsigned Sister_Element = this->GetSister(current_Target_element);
                    PottsElement<DIM>* pSisterElement = this->GetElement(Sister_Element);
                    PRINT_VARIABLE(Sister_Element);
                    if ( *(containing_elements.begin()) == Sister_Element && pSisterElement->GetNumNodes() < 2 )
                    {
                        RunPeriodicTrialSwitch =0;
                    }
                    
                }if ( mIsPeriodic ==0)
                {
                    std::vector<unsigned>::iterator Current_iter = mBoundaryVector.begin();
                    std::advance(Current_iter , node_index);   
                    NotOnBoundary =*Current_iter;
                }

       

                 //                 This prevents me hitting a boundary in the non periodic mesh -- there seems to be an issure at the boundaries and the cells have a lot of random switches... I think it is due to the area
                if ( (mIsPeriodic ==0 && NotOnBoundary ==0)|| (mIsPeriodic ==1 && RunPeriodicTrialSwitch ==1) )
                {
                    // PRINT_VARIABLE(NotOnBoundary)
                    // This will prevent cell being taken over by the sister if there is only one element in the sister or something
                    //   neighbour_location_index.Sister ~= p_node->GetIndex() when SisterSize <2

                    double delta_H = 0.0; // This is H_1-H_0.

                    // Now add contributions to the Hamiltonian from each AbstractPottsUpdateRule
                    for (typename std::vector<boost::shared_ptr<AbstractUpdateRule<DIM> > >::iterator iter = this->mUpdateRuleCollection.begin();
                         iter != this->mUpdateRuleCollection.end();
                         ++iter)
                    {
                        // This static cast is fine, since we assert the update rule must be a Potts update rule in AddUpdateRule()
                        double dH = (boost::static_pointer_cast<AbstractWrappedPottsUpdateRule<DIM> >(*iter))->EvaluateHamiltonianContribution(neighbour_location_index, p_node->GetIndex(), *this);
                        delta_H += dH;
                    }

                    // Generate a uniform random number to do the random motion
                    double random_number = p_gen->ranf();

                    double p = exp(-delta_H / this->GetTemperature());
                    if (delta_H <= 0 || random_number < p)
                    {
                        // Do swap
                        for (std::set<unsigned>::iterator iter = containing_elements.begin();
                             iter != containing_elements.end();
                             ++iter)
                        {
                          
                            // Remove lattice from the old element
                            this->GetElement(*iter)->DeleteNode(this->GetElement(*iter)->GetNodeLocalIndex(node_index));
                        }
                        if(!neighbour_containing_elements.empty())
                        {
                       
                            this->GetElement(*neighbour_containing_elements.begin())->AddNode(this->mrMesh.GetNode(node_index));
                        }   
                    }      
                    }
                } 
            }

        }
    }
}

template <unsigned DIM>
void WrappedPottsBasedCellPopulation<DIM>::AddUpdateRule(boost::shared_ptr<AbstractUpdateRule<DIM> > pUpdateRule)
{
    assert(bool(dynamic_cast<AbstractWrappedPottsUpdateRule<DIM>*>(pUpdateRule.get())));
    this->mUpdateRuleCollection.push_back(pUpdateRule);
}

template <unsigned DIM>
void WrappedPottsBasedCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    // *rParamsFile << "\t\t<Temperature>" << mTemperature << "</Temperature>\n";
    //  *rParamsFile << "\t\t<NumSweepsPerTimestep>" << mNumSweepsPerTimestep << "</NumSweepsPerTimestep>\n";

    // Call method on direct parent class
    PottsBasedCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

// Explicit instantiation
// template class WrappedPottsBasedCellPopulation<1>;
// template class WrappedPottsBasedCellPopulation<2>;
template class WrappedPottsBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(WrappedPottsBasedCellPopulation)

// template <unsigned DIM>
// void WrappedPottsBasedCellPopulation<DIM>::DividePopulation()
// {

//     // std::list<CellPtr> CellsToDivide = this->rGetCells();

//     // // Save the cells in a list, then do the dividing

//     // for (std::list<CellPtr>::iterator cell_iter = CellsToDivide.begin();
//     //          cell_iter != CellsToDivide.end();
//     //          ++cell_iter)
//     //     {
//     //     PottsElement<DIM>* p_element = this->GetElementCorrespondingToCell(*cell_iter);
//     //     // CellPtr pParentCell = this->GetCellUsingLocationIndex(*cell_iter);

//     //     //-- base this on curvature for the 2d in 3d simulations

//     //     c_vector<double, DIM> Midpoint = zero_vector<double>(DIM);

//     //     for (unsigned i = 0; i < p_element->GetNumNodes(); i++)
//     //     {
//     //         Midpoint += p_element->GetNodeLocation(i);
//     //     }
//     //     Midpoint /= p_element->GetNumNodes();

//     //     // now need to iterate over the Lattice sites. If they are on the right they need to be saved in a vector
//     //     std::vector<unsigned> RightLatticeSites;
//     //     std::vector<Node<DIM>*> nodes_elem;
//     //     for (unsigned i = p_element->GetNumNodes() - 1; i < -1; i--) // for some reason this needs to go backwards
//     //     {
//     //         if (p_element->GetNodeLocation(i)[0] > Midpoint[0])
//     //         {
//     //             RightLatticeSites.push_back(i);
//     //             nodes_elem.push_back(p_element->GetNode(i));
//     //         }
//     //     }

//     //     unsigned new_element_index = this->GetNumElements();
//     //     this->rGetMesh().AddElement(new PottsElement<DIM>(new_element_index, nodes_elem));

//     //     // Now I need to go through and remove everything in this vector from the old cell and add it to the new one
//     //     for (std::vector<unsigned>::iterator i = RightLatticeSites.begin(); i != RightLatticeSites.end(); ++i)
//     //     {
//     //         p_element->DeleteNode(*i);
//     //     }

//     //     // Associate the new cell with the element

//     //     // CellPtr pNewCell(new Cell(p_state, p_model));
//     //     CellPtr pNewCell(new Cell((*cell_iter)->GetMutationState(), (*cell_iter)->GetCellCycleModel()));

//     //     this->mCells.push_back(pNewCell);

//     //     // Update location cell map
//     //     CellPtr p_created_cell = this->mCells.back();
//     //     this->SetCellUsingLocationIndex(new_element_index, p_created_cell);
//     //     pNewCell->GetCellData()->SetItem("CenterLine", Midpoint[0]);
//     //     (*cell_iter)->GetCellData()->SetItem("CenterLine", Midpoint[0]);

//     //     //  Link the sister cells
//     //     unsigned Mother_element_index = p_element->GetIndex();

//     //     pNewCell->GetCellData()->SetItem("Sister", Mother_element_index);
//     //     (*cell_iter)->GetCellData()->SetItem("Sister", new_element_index);

//     //     // 1 for mother -- left
//     //     // 2 for daughter -- right
//     //     pNewCell->GetCellData()->SetItem("DivisionStatus", 2);
//     //     (*cell_iter)->GetCellData()->SetItem("DivisionStatus", 1);

//     //     TRACE("hih");
//     // }

// }

// template <unsigned DIM>
// CellPtr WrappedPottsBasedCellPopulation<DIM>::DivideCell(CellPtr pNewCell, CellPtr pParentCell)
// {
//     // Get the element associated with this cell
//     // PottsArbitrarySurfaceIn3DMesh M;

//     PottsElement<DIM>* p_element = this->GetElementCorrespondingToCell(pParentCell);

//     //-- base this on curvature for the 2d in 3d simulations

//     c_vector<double, DIM> Midpoint = zero_vector<double>(DIM);

//     for (unsigned i = 0; i < p_element->GetNumNodes(); i++)
//     {
//         Midpoint += p_element->GetNodeLocation(i);
//     }
//     Midpoint /= p_element->GetNumNodes();

//     // now need to iterate over the Lattice sites. If they are on the right they need to be saved in a vector
//     std::vector<unsigned> RightLatticeSites;
//     std::vector<Node<DIM>*> nodes_elem;
//     for (unsigned i = p_element->GetNumNodes() - 1; i < -1; i--) // for some reason this needs to go backwards
//     {
//         if (p_element->GetNodeLocation(i)[0] > Midpoint[0])
//         {
//             // Checks if the lattice site is on the right, if so this lattice will be added to the daugter cell, and deleted from the mother in the comming for loop
//             RightLatticeSites.push_back(i);
//             nodes_elem.push_back(p_element->GetNode(i));
//         }
//     }

//     unsigned new_element_index = this->GetNumElements();

//     this->rGetMesh().AddElement(new PottsElement<DIM>(new_element_index, nodes_elem));

//     // Now I need to go through and remove everything in this vector from the old cell and add it to the new one
//     for (std::vector<unsigned>::iterator i = RightLatticeSites.begin(); i != RightLatticeSites.end(); ++i)
//     {
//         p_element->DeleteNode(*i);
//     }

//     // Associate the new cell with the element
//     this->mCells.push_back(pNewCell);

//     // Update location cell map
//     CellPtr p_created_cell = this->mCells.back();
//     this->SetCellUsingLocationIndex(new_element_index, p_created_cell);

//     pNewCell->GetCellData()->SetItem("CenterLine", Midpoint[0]);
//     pParentCell->GetCellData()->SetItem("CenterLine", Midpoint[0]);

//     //  Link the sister cells
//     unsigned Mother_element_index = p_element->GetIndex();

//     pNewCell->GetCellData()->SetItem("Sister", Mother_element_index);
//     pParentCell->GetCellData()->SetItem("Sister", new_element_index);

//     // 2 for daughter,
//     // 1 for mother
//     pNewCell->GetCellData()->SetItem("DivisionStatus", 2);
//     pParentCell->GetCellData()->SetItem("DivisionStatus", 1);

//     TRACE("A success")
//     return p_created_cell;
// }

// Want to keep these sites clear so the lattice sites in the cells have something to switch to as the cell is moving
// double X = 5193;
// double Y = 2020;

// if (node_index == X || node_index == Y)
// {
//     // Need to break here
//     // TRACE("Need to break here 2")
//     ExcuteMCS = 0;
// }
// if (ExcuteMCS == 1)
// {