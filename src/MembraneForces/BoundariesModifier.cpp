
/*

This code controles the membrame properties. 

At the start of a simulation, the intial conditions delievered at 20(?)X less than original starting configuration of the vessel 
The system then develops to an equlibirum (that will be the realisitc size of a network). After this point the memnbrane properties will 
be variable.

For the inital impletmentation of this code, the membrane properties in the center will be changed to decrease the vessel size to 2X the intial 
condition, with some smoothing on the vessel properties just outside the region. 

After this, the membrane constants will be dependent on the number of cells in the region and the denoted mutation state. 

*/

#include "BoundariesModifier.hpp"
#include <algorithm>
#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"

#include <math.h>
#include "Debug.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BoundariesModifier<ELEMENT_DIM, SPACE_DIM>::BoundariesModifier()
        : AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BoundariesModifier<ELEMENT_DIM, SPACE_DIM>::~BoundariesModifier()
{
}

// Need boundaires and cells
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BoundariesModifier<ELEMENT_DIM, SPACE_DIM>::CreateBoundaryNodes(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::vector<c_vector<double, 3> > boundary_plane_normals, std::vector<c_vector<double, 3> > boundary_plane_points)
{

    // std::vector<c_vector<double,3> > boundary_plane_normals, std::vector<c_vector<double,3> > boundary_plane_points
    double Counter = 0;
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        Node<SPACE_DIM>* pNode = rCellPopulation.rGetMesh().GetNode(node_index);
        c_vector<double, SPACE_DIM> node_location = pNode->rGetLocation();
        std::vector<c_vector<double, 3> >::iterator Normal_iter = boundary_plane_normals.begin();
        std::vector<c_vector<double, 3> >::iterator Point_iter = boundary_plane_points.begin();
        for (int i = 0; i < boundary_plane_normals.size(); i++)
        {
            double signed_distance = inner_prod(node_location - *Point_iter, *Normal_iter);
            if (signed_distance > 0.0)
            {
                cell_iter->GetCellData()->SetItem("Boundary", 1);
                break;
            }
            else
            {
                cell_iter->GetCellData()->SetItem("Boundary", 0);
            }
            advance(Normal_iter, 1);
            advance(Point_iter, 1);
        }
    }

    TRACE("Finished assigning boundaries")
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BoundariesModifier<ELEMENT_DIM, SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{

    TRACE("SET UpSolve")
    assert(ELEMENT_DIM == 2);
    assert(SPACE_DIM == 3);
    // std::map<unsigned, c_vector<unsigned, 5> > mNearestNodesMap;

    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        Node<SPACE_DIM>* pNode = rCellPopulation.rGetMesh().GetNode(node_index);
        c_vector<long double, SPACE_DIM> CellLocation = pNode->rGetLocation();

        c_vector<unsigned, 5> NearestNodes;
        double Distance5 = 5.44; // Upper bound
        double Distance4 = 5.44;
        double Distance3 = 5.44;
        double Distance2 = 5.44;
        double Distance1 = 5.44;

        //find the distance to the nearest neighbour
        // std::set<unsigned> Neighbours  = rCellPopulation.GetNeighbouringNodeIndices(node_index);
      
        // Node<SPACE_DIM>* pNeighbour_Node = rCellPopulation.rGetMesh().GetNode(*Neighbours.begin());
        // c_vector<long double, SPACE_DIM> NeighbourLocation = pNeighbour_Node->rGetLocation();
        // double LatticeSpacing = norm_2(NeighbourLocation - CellLocation);

        // double MinimimNeighbourDistance=20 * LatticeSpacing;

        // PRINT_VECTOR(CellLocation)

        if (cell_iter->GetCellData()->GetItem("Boundary") == 1)
        {
            // 
            for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter_Search = rCellPopulation.Begin();
                 cell_iter_Search != rCellPopulation.End();
                 ++cell_iter_Search)
            {
                unsigned node_index_Search = rCellPopulation.GetLocationIndexUsingCell(*cell_iter_Search);

                if (cell_iter_Search->GetCellData()->GetItem("Boundary") == 0 && node_index_Search != node_index)
                {
                    // Need location of each

                    Node<SPACE_DIM>* pNode1 = rCellPopulation.rGetMesh().GetNode(node_index_Search);

                    c_vector<long double, SPACE_DIM> LocationOfNode = pNode1->rGetLocation();
                    //  PRINT_VECTOR(LocationOfNode)

                    double Distance = norm_2(LocationOfNode - CellLocation);
                    // if (Distance > MinimimNeighbourDistance) 
                    // {
                        //   PRINT_VARIABLE(Distance)
                        if (Distance <= Distance1)
                        {
                            NearestNodes[0] = node_index_Search;
                            Distance1 = Distance;

                            // TRACE("D1")
                        }
                        else if (Distance <= Distance2)
                        {
                            NearestNodes[1] = node_index_Search;
                            Distance2 = Distance;
                            // TRACE("D2")
                        }
                        else if (Distance <= Distance3)
                        {
                            NearestNodes[2] = node_index_Search;
                            Distance3 = Distance;
                            // TRACE("D3")
                        }
                        else if (Distance <= Distance4)
                        {
                            NearestNodes[3] = node_index_Search;
                            Distance4 = Distance;
                            // TRACE("D4")
                        }
                        else if (Distance <= Distance5)
                        {
                            NearestNodes[4] = node_index_Search;
                            Distance5 = Distance;
                            // TRACE("D5")
                        }
                    // }
                }
            }
            mNearestNodesMap[node_index] = NearestNodes;
        }
        else{
            NearestNodes[0] = 50;
            NearestNodes[1] = 50;
            NearestNodes[2] = 50;
            NearestNodes[3] = 50;
            NearestNodes[4] = 50;
            mNearestNodesMap[node_index] = NearestNodes;
        }
    }

    TRACE("done SetUpSolve");
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<unsigned, 5> BoundariesModifier<ELEMENT_DIM, SPACE_DIM>::GetNeighbouringNodes(unsigned node_index)
{
    // assert(NodeIsInSet)
    return mNearestNodesMap[node_index];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::map<unsigned, c_vector<unsigned, 5> > BoundariesModifier<ELEMENT_DIM, SPACE_DIM>::GetNeighbouringNodesMap()
{
    PRINT_4_VARIABLES(mNearestNodesMap[171][0], mNearestNodesMap[171][1], mNearestNodesMap[171][2], mNearestNodesMap[171][3])
    return mNearestNodesMap;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BoundariesModifier<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    // TRACE("Not Doing this yet");
    // UpdateCellData( rCellPopulation);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BoundariesModifier<ELEMENT_DIM, SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{

    // assert(ELEMENT_DIM ==2);
    // TRACE("Switharoo")
    // assert(SPACE_DIM == 3);
    // // std::map<unsigned, c_vector<unsigned, 5> > mNearestNodesMap;

    //      for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
    //             cell_iter != rCellPopulation.End();
    //             ++cell_iter)
    //     {
    //         unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
    //         Node<SPACE_DIM>* pNode = rCellPopulation.rGetMesh().GetNode(node_index);
    //         c_vector<unsigned, 5>  NearestNodes = mNearestNodesMap[node_index];
    //         // p_new_node->ClearAppliedForce();
    //         // pNode->AddAppliedForceContribution(MembraneForceMap[node_index] ); // Add the new force
    //        cell_iter->GetCellData()->SetItem("MembraneForce", -350 );

    //         //             Node<3>* pReferenceNode = p_cell_population->rGetMesh().GetNode(ReferenceNode);
    //       c_vector<long double, 3> ForceOnNode = pNode->rGetAppliedForce();
    //       pNode->ClearAppliedForce(); // remove any applied force, this stops there begin two expanding forces at this node
    //       pNode->AddAppliedForceContribution(Create_c_vector(-1,-1,0));
    //     }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BoundariesModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class BoundariesModifier<1, 1>;
template class BoundariesModifier<1, 2>;
template class BoundariesModifier<2, 2>;
template class BoundariesModifier<1, 3>;
template class BoundariesModifier<2, 3>;
template class BoundariesModifier<3, 3>;
// Serialization for Boost >= 1.36

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(BoundariesModifier)
