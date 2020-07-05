
/*

    History dependence in the mutable mesh -- For every node in the new mesh we find the 3 local/closest nodes in the old mesh.
    Make a triangle out of these old nodes. 
    Translate and rotate/map this triangle back to the origin
    Make two of the edges the new basis vectors 
    Describe the Node in the new Mesh P using the basis vectors (needs to be translated toward the origina as the triangle was)
    Find the shape funcitons of this triangle 
    Figure out where these nodes where in the initial configuaration and make an inital triangle 
    project/map this triangle back to the orign
    Find the displacement vector
    Using the dispacement equation V4 = P - po = N1*V1 + N2*V2 + N3*V3, find po, the inital position of the new node in the inital configuarion
    Use two edges of the inital triangle as the basis vectors to discribe po, then unmap the inital triangle to its proper/unmapped configuration

*/

#include "HistoryDependentMutableMeshModifier.hpp"
#include <algorithm>
#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"

#include <math.h>
#include "Debug.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HistoryDependentMutableMeshModifier<ELEMENT_DIM, SPACE_DIM>::HistoryDependentMutableMeshModifier()
        : AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HistoryDependentMutableMeshModifier<ELEMENT_DIM, SPACE_DIM>::~HistoryDependentMutableMeshModifier()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDependentMutableMeshModifier<ELEMENT_DIM, SPACE_DIM>::RemeshingWithHistoryDepenance(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, MutableMesh<2,3>& rOldMesh,MutableMesh<2,3>& rNewMesh, MutableMesh<2,3>& rInitalOldMesh,)
{
    // Need - The new node from the inital mesh P
    //      - Triangle made up of old nodes 
    //      - Triangle made up of old nodes initial configuration     
}




template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDependentMutableMeshModifier<ELEMENT_DIM, SPACE_DIM>::NewNodeInInitalConfiguration(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
    // Need - The new node from the inital mesh P
    //      - Triangle made up of old nodes 
    //      - Triangle made up of old nodes initial configuration     
}

    

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDependentMutableMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{


    TRACE("SET UpSolve")
    // assert(ELEMENT_DIM == 2);
    // assert(SPACE_DIM == 3);
    // // std::map<unsigned, c_vector<unsigned, 5> > mNearestNodesMap;

    // for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
    //      cell_iter != rCellPopulation.End();
    //      ++cell_iter)
    // {
    //     unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
    //     Node<SPACE_DIM>* pNode = rCellPopulation.rGetMesh().GetNode(node_index);
    //     c_vector<long double, SPACE_DIM> CellLocation = pNode->rGetLocation();

    //     c_vector<unsigned, 5> NearestNodes;
    //     double Distance5 = 5.44; // Upper bound
    //     double Distance4 = 5.44;
    //     double Distance3 = 5.44;
    //     double Distance2 = 5.44;
    //     double Distance1 = 5.44;

    //     //find the distance to the nearest neighbour
    //     // std::set<unsigned> Neighbours  = rCellPopulation.GetNeighbouringNodeIndices(node_index);
      
    //     // Node<SPACE_DIM>* pNeighbour_Node = rCellPopulation.rGetMesh().GetNode(*Neighbours.begin());
    //     // c_vector<long double, SPACE_DIM> NeighbourLocation = pNeighbour_Node->rGetLocation();
    //     // double LatticeSpacing = norm_2(NeighbourLocation - CellLocation);

    //     // double MinimimNeighbourDistance=20 * LatticeSpacing;

    //     // PRINT_VECTOR(CellLocation)

    //     if (cell_iter->GetCellData()->GetItem("Boundary") == 1)
    //     {
    //         // 
    //         for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter_Search = rCellPopulation.Begin();
    //              cell_iter_Search != rCellPopulation.End();
    //              ++cell_iter_Search)
    //         {
    //             unsigned node_index_Search = rCellPopulation.GetLocationIndexUsingCell(*cell_iter_Search);

    //             if (cell_iter_Search->GetCellData()->GetItem("Boundary") == 0 && node_index_Search != node_index)
    //             {
    //                 // Need location of each

    //                 Node<SPACE_DIM>* pNode1 = rCellPopulation.rGetMesh().GetNode(node_index_Search);

    //                 c_vector<long double, SPACE_DIM> LocationOfNode = pNode1->rGetLocation();
    //                 //  PRINT_VECTOR(LocationOfNode)

    //                 double Distance = norm_2(LocationOfNode - CellLocation);
    //                 // if (Distance > MinimimNeighbourDistance) 
    //                 // {
    //                     //   PRINT_VARIABLE(Distance)
    //                     if (Distance <= Distance1)
    //                     {
    //                         NearestNodes[0] = node_index_Search;
    //                         Distance1 = Distance;

    //                         // TRACE("D1")
    //                     }
    //                     else if (Distance <= Distance2)
    //                     {
    //                         NearestNodes[1] = node_index_Search;
    //                         Distance2 = Distance;
    //                         // TRACE("D2")
    //                     }
    //                     else if (Distance <= Distance3)
    //                     {
    //                         NearestNodes[2] = node_index_Search;
    //                         Distance3 = Distance;
    //                         // TRACE("D3")
    //                     }
    //                     else if (Distance <= Distance4)
    //                     {
    //                         NearestNodes[3] = node_index_Search;
    //                         Distance4 = Distance;
    //                         // TRACE("D4")
    //                     }
    //                     else if (Distance <= Distance5)
    //                     {
    //                         NearestNodes[4] = node_index_Search;
    //                         Distance5 = Distance;
    //                         // TRACE("D5")
    //                     }
    //                 // }
    //             }
    //         }
    //         mNearestNodesMap[node_index] = NearestNodes;
    //     }
    //     else{
    //         NearestNodes[0] = 50;
    //         NearestNodes[1] = 50;
    //         NearestNodes[2] = 50;
    //         NearestNodes[3] = 50;
    //         NearestNodes[4] = 50;
    //         mNearestNodesMap[node_index] = NearestNodes;
    //     }
    // }

    // TRACE("done SetUpSolve");
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDependentMutableMeshModifier<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    // TRACE("Not Doing this yet");
    // UpdateCellData( rCellPopulation);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDependentMutableMeshModifier<ELEMENT_DIM, SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HistoryDependentMutableMeshModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class HistoryDependentMutableMeshModifier<1, 1>;
template class HistoryDependentMutableMeshModifier<1, 2>;
template class HistoryDependentMutableMeshModifier<2, 2>;
template class HistoryDependentMutableMeshModifier<1, 3>;
template class HistoryDependentMutableMeshModifier<2, 3>;
template class HistoryDependentMutableMeshModifier<3, 3>;
// Serialization for Boost >= 1.36

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(HistoryDependentMutableMeshModifier)
