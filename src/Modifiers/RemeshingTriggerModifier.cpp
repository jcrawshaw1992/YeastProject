
/*

This code controles the membrame properties. 

At the start of a simulation, the intial conditions delievered at 20(?)X less than original starting configuration of the vessel 
The system then develops to an equlibirum (that will be the realisitc size of a network). After this point the memnbrane properties will 
be variable.

For the inital impletmentation of this code, the membrane properties in the center will be changed to decrease the vessel size to 2X the intial 
condition, with some smoothing on the vessel properties just outside the region. 

After this, the membrane constants will be dependent on the number of cells in the region and the denoted mutation state. 

*/

#include "RemeshingTriggerModifier.hpp"
#include <algorithm>
#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"

#include "OutsideFLuidSimulationMutation.hpp"

#include <math.h>
#include "Debug.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RemeshingTriggerModifier<ELEMENT_DIM, SPACE_DIM>::RemeshingTriggerModifier()
        : AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RemeshingTriggerModifier<ELEMENT_DIM, SPACE_DIM>::~RemeshingTriggerModifier()
{
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerModifier<ELEMENT_DIM, SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
     HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation = static_cast<HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);

     mExecute =0;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerModifier<ELEMENT_DIM, SPACE_DIM>::SetRemeshingInterval(int RemeshingInterval)
{
    mRemeshingInterval = RemeshingInterval;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerModifier<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
     HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation = static_cast<HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);
 
    if (mExecute>mRemeshingInterval)
    {
        pCellPopulation->ExecuteHistoryDependentRemeshing();
        mExecute = 0;
    } 
    mExecute +=1;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class RemeshingTriggerModifier<1, 1>;
template class RemeshingTriggerModifier<1, 2>;
template class RemeshingTriggerModifier<2, 2>;
template class RemeshingTriggerModifier<1, 3>;
template class RemeshingTriggerModifier<2, 3>;
template class RemeshingTriggerModifier<3, 3>;
// Serialization for Boost >= 1.36

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(RemeshingTriggerModifier)



// Old ways to get boundaries 

    // Want to see what happens if i identify the boundaries by the size of the local elements 

//     std::map<unsigned, double> AreaOfCells;
//     double AverageCellArea = 0;
//     double NumberOfCells = 0;
//     MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);
//     for ( typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
//             cell_iter != rCellPopulation.End();
//             ++cell_iter)
//     {

//         unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
//         Node<SPACE_DIM>* p_node = rCellPopulation.GetNode(node_index);

//         std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();

//         double Area =0;
//         assert(containing_elements.size() > 0);
//         for (std::set<unsigned>::iterator iter = containing_elements.begin();
//                 iter != containing_elements.end();
//                 ++iter)
//         {
//             Node<SPACE_DIM>* pNode0 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(0));
//             Node<SPACE_DIM>* pNode1= p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(1));
//             Node<SPACE_DIM>* pNode2 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(2));

//             c_vector<double, SPACE_DIM> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
//             c_vector<double, SPACE_DIM> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3

//             c_vector<double, SPACE_DIM> normalVector = VectorProduct(vector_12, vector_13);
//             Area+= 0.5*norm_2(normalVector)/3;
//         }
//         NumberOfCells +=1;
//         AreaOfCells[node_index] = Area;
//         AverageCellArea +=Area;
//     }

//     AverageCellArea/= NumberOfCells;

//    for ( typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
//             cell_iter != rCellPopulation.End();
//             ++cell_iter)
//         {
            
            
//             unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
//             if (AreaOfCells[node_index] < 1e-1*AverageCellArea)
//             {
//                 TRACE("Trip")
//                 cell_iter->GetCellData()->SetItem("SecondBoundary", 1);

//             }
//             else{
//                 cell_iter->GetCellData()->SetItem("SecondBoundary", 0);

//             }

//         }
    