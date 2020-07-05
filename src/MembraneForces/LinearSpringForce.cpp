

#include "LinearSpringForce.hpp"
#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
LinearSpringForce<ELEMENT_DIM,SPACE_DIM>::LinearSpringForce()
   : AbstractForce<ELEMENT_DIM,SPACE_DIM>()
{
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForce<ELEMENT_DIM,SPACE_DIM>::AddForceContribution(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    // Throw an exception message if not using a subclass of AbstractCentreBasedCellPopulation
    if (dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("Subclasses of LinearSpringForce are to be used with subclasses of AbstractCentreBasedCellPopulation only");
    }

        MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_static_cast_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation);

        std::map<unsigned, c_vector<double, 3> > SpringForceOnNodeMap;

         // Iterate over all springs and add force contributions

        for (typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator node_iter = p_static_cast_cell_population->rGetMesh().GetNodeIteratorBegin();
            node_iter != p_static_cast_cell_population->rGetMesh().GetNodeIteratorEnd();
            ++node_iter)
        {   
            unsigned node_index = node_iter->GetIndex();
            SpringForceOnNodeMap[node_index] = Create_c_vector(0,0,0);
        }



        // Iterate over all springs and add force contributions
        for (typename MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringIterator spring_iterator = p_static_cast_cell_population->SpringsBegin();
             spring_iterator != p_static_cast_cell_population->SpringsEnd();
             ++spring_iterator)
        {
            unsigned NodeA_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned NodeB_index = spring_iterator.GetNodeB()->GetIndex();

        // Find the pointer for each of the nodes
            CellPtr pCellA =  rCellPopulation.GetCellUsingLocationIndex(NodeA_index); // This gives us the pointer based on the index.
            CellPtr pCellB =  rCellPopulation.GetCellUsingLocationIndex(NodeB_index);


           // Calculate the force between nodes
            double volumeA =rCellPopulation.GetVolumeOfCell(pCellA);
            double volumeB = rCellPopulation.GetVolumeOfCell(pCellB);
            // PRINT_2_VARIABLES(volumeA, volumeB)

            // Calculate the force between nodes
            c_vector<double, SPACE_DIM> force = CalculateForceBetweenNodes(NodeA_index, NodeB_index, rCellPopulation);

            // Add the force contribution to each node
            // c_vector<double, SPACE_DIM> negative_force = -1.0*force;
            // spring_iterator.GetNodeB()->AddAppliedForceContribution(negative_force/volumeB);
            // spring_iterator.GetNodeA()->AddAppliedForceContribution(force/volumeA);
            // PRINT_VECTOR(SpringForceOnNodeMap[NodeA_index])
            SpringForceOnNodeMap[NodeA_index] += force/volumeA;
            SpringForceOnNodeMap[NodeB_index] -= force/volumeB;

           
            // 
            // PRINT_VECTOR(SpringForceOnNodeMap[NodeA_index])
            // PRINT_VARIABLE(volumeB)
     
        }
        // Add the forces and do what I need to messing with stuff
         for (typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator node_iter = p_static_cast_cell_population->rGetMesh().GetNodeIteratorBegin();
            node_iter != p_static_cast_cell_population->rGetMesh().GetNodeIteratorEnd();
            ++node_iter)
        {   
            unsigned node_index = node_iter->GetIndex();
            // Get the cell corresponding to this node
            CellPtr p_cell = p_static_cast_cell_population->GetCellUsingLocationIndex(node_index);

            if (p_cell->GetCellData()->GetItem("Boundary") ==1)
            {
                c_vector<double, 3> ForceOnEdgeNode = Create_c_vector(0,0,0);
                c_vector<unsigned, 2>  LocalNodes = mNearestNodesMap[node_index];
                // PRINT_3_VARIABLES(LocalNodses[0], LocalNodes[1], LocalNodes[2])
                for(int i =0 ; i<2; i++)
                {   
                    ForceOnEdgeNode += SpringForceOnNodeMap[LocalNodes[i]];
                }
                ForceOnEdgeNode/=2;


                // Now I want to get the normal 

                std::set<unsigned>& containing_elements = node_iter->rGetContainingElementIndices();
                c_vector<double, SPACE_DIM> Normal = zero_vector<double>(SPACE_DIM);
                MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation);
                assert(containing_elements.size() > 0);
                for (std::set<unsigned>::iterator iter = containing_elements.begin();
                    iter != containing_elements.end();
                    ++iter)
                {
                    Node<SPACE_DIM>* pNode0 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(0));
                    Node<SPACE_DIM>* pNode1 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(1));
                    Node<SPACE_DIM>* pNode2 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(2));

                    c_vector<double, SPACE_DIM> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
                    c_vector<double, SPACE_DIM> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3

                    c_vector<double, SPACE_DIM> normalVector = VectorProduct(vector_12, vector_13);
                    
                    Normal += normalVector;///norm_2(normalVector);

                }
                // Normal[2] =0;

                assert(norm_2(Normal) <10e10);
                assert(norm_2(Normal) !=0);

                Normal /=norm_2(Normal);
                // PRINT_VARIABLE(norm_2(Normal))
                if (ForceOnEdgeNode[0] >10e10 || ForceOnEdgeNode[1] >10e10 || ForceOnEdgeNode[2] >10e10)
                {
                    PRINT_VECTOR(ForceOnEdgeNode)
                     PRINT_VECTOR(Normal)

                      for(int i =0 ; i<2; i++)
                        {   
                            PRINT_VECTOR(SpringForceOnNodeMap[LocalNodes[i]]);
                        }

                }


                
                ForceOnEdgeNode = -norm_2(ForceOnEdgeNode)* Normal;
                node_iter->AddAppliedForceContribution(ForceOnEdgeNode);
                

                 assert(ForceOnEdgeNode[0] <10e10);
                 assert(ForceOnEdgeNode[1] <10e10);
                 assert(ForceOnEdgeNode[2] <10e10);
                 
                
                
                p_cell->GetCellData()->SetItem("SpringForce", norm_2(ForceOnEdgeNode) );
                p_cell->GetCellData()->SetItem("SpringForceX", ForceOnEdgeNode[0] );
                p_cell->GetCellData()->SetItem("SpringForceY", ForceOnEdgeNode[1] );
                p_cell->GetCellData()->SetItem("SpringForceZ", ForceOnEdgeNode[2] );
            }
            else
            {
                // Add applied force from the map 
                // SpringForceOnNodeMap[node_index][2]=0;

        

             if (SpringForceOnNodeMap[node_index][0] >10e10 || SpringForceOnNodeMap[node_index][1] >10e10 || SpringForceOnNodeMap[node_index][2] >10e10)
                {
                    PRINT_VECTOR(SpringForceOnNodeMap[node_index])
                }



                 assert(SpringForceOnNodeMap[node_index][0] <10e10);
                 assert(SpringForceOnNodeMap[node_index][1] <10e10);
                 assert(SpringForceOnNodeMap[node_index][2] <10e10);

                node_iter->AddAppliedForceContribution(SpringForceOnNodeMap[node_index]);
                p_cell->GetCellData()->SetItem("SpringForce", norm_2(SpringForceOnNodeMap[node_index]) );
                p_cell->GetCellData()->SetItem("SpringForceX", SpringForceOnNodeMap[node_index][0] );
                p_cell->GetCellData()->SetItem("SpringForceY", SpringForceOnNodeMap[node_index][1] );
                p_cell->GetCellData()->SetItem("SpringForceZ", SpringForceOnNodeMap[node_index][2] );
                // PRINT_VECTOR(SpringForceOnNodeMap[node_index])

            }
        }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> LinearSpringForce<ELEMENT_DIM,SPACE_DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                                                    unsigned nodeBGlobalIndex,
                                                                                    AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    Node<SPACE_DIM>* p_node_a = rCellPopulation.GetNode(nodeAGlobalIndex);
    Node<SPACE_DIM>* p_node_b = rCellPopulation.GetNode(nodeBGlobalIndex);

    // Get the node locations
    const c_vector<double, SPACE_DIM>& r_node_a_location = p_node_a->rGetLocation();
    const c_vector<double, SPACE_DIM>& r_node_b_location = p_node_b->rGetLocation();

    /*
     * Get the unit vector parallel to the line joining the two nodes
     * We use the mesh method GetVectorFromAtoB() to compute the direction of the
     * unit vector along the line joining the two nodes, rather than simply subtract
     * their positions, because this method can be overloaded (e.g. to enforce a
     * periodic boundary in Cylindrical2dMesh).
     */
    c_vector<double, SPACE_DIM> UnitSpring = rCellPopulation.rGetMesh().GetVectorFromAtoB(r_node_a_location, r_node_b_location); // Not yet normalised 

    // Calculate the distance between the two nodes
    double distance_between_nodes = norm_2(UnitSpring);
    UnitSpring /= distance_between_nodes;
    assert(distance_between_nodes > 0);
    assert(!std::isnan(distance_between_nodes));

    /*
     * Calculate the rest length of the spring connecting the two nodes with a default
     */

    double rest_length = static_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation)->GetRestLength(nodeAGlobalIndex, nodeBGlobalIndex);
   
    return mMeinekeSpringStiffness * UnitSpring * (distance_between_nodes - rest_length);

}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForce<ELEMENT_DIM,SPACE_DIM>::SetNearestNeighboursMap(std::map<unsigned, c_vector<unsigned, 2> > NearestNodesMap)
{
    mNearestNodesMap = NearestNodesMap;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForce<ELEMENT_DIM,SPACE_DIM>::SetMeinekeSpringStiffness(double springStiffness)
{
    assert(springStiffness > 0.0);
    mMeinekeSpringStiffness = springStiffness;
}




template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{

    // Call method on direct parent class
    AbstractForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(rParamsFile);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void LinearSpringForce<ELEMENT_DIM,SPACE_DIM>::WriteDataToVisualizerSetupFile(out_stream& pVizSetupFile)
{
    // *pVizSetupFile << "Cutoff\t" << mMechanicsCutOffLength << "\n";
}





// Explicit instantiation
template class LinearSpringForce<1,1>;
template class LinearSpringForce<1,2>;
template class LinearSpringForce<2,2>;
template class LinearSpringForce<1,3>;
template class LinearSpringForce<2,3>;
template class LinearSpringForce<3,3>;
