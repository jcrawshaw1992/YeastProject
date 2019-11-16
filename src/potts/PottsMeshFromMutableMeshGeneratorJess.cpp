

#include "PottsMeshFromMutableMeshGeneratorJess.hpp"
#include "Debug.hpp"
#include "UblasCustomFunctions.hpp"

template<unsigned SPACE_DIM>
// Send in the Mutable mesh 
PottsMeshFromMutableMeshGeneratorJess<SPACE_DIM>::PottsMeshFromMutableMeshGeneratorJess(MutableMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
{

	assert( SPACE_DIM==3);

    std::vector<Node<SPACE_DIM>*> nodes; // Object definition for the vector to hold the node objects
    std::vector<std::set<unsigned> > node_neighbours; // Object definition for the vector to hold the neighbours to each node

    for (typename MutableMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator node_iter = rMesh.GetNodeIteratorBegin();
           node_iter != rMesh.GetNodeIteratorEnd();
           ++node_iter)
           {    
         
                std::vector<unsigned> NeighboursToNode;
                unsigned node_index =  node_iter->GetIndex(); //rCellPopulation.GetLocationIndexUsingCell(*cell_iter); //p_cell_population->GetLocationIndexUsingCell(*cell_iter);
                // Node<3>* p_node = node_iter;
            
                // This piece records all the elements associated with a mutant cell
                std::set<unsigned>& containing_elements = node_iter->rGetContainingElementIndices();
                assert(containing_elements.size() > 0);
                // LabledElements.push_back(containing_elements);
                // Loops over the Elements connected to this node
                for (std::set<unsigned>::iterator iter = containing_elements.begin();
                    iter != containing_elements.end();
                    ++iter)
                {
                  // Loops over the nodes in the element and saves them as neighbours
                    for (int i=0; i<3; i++)
                    {
                          Node<3>* pNode = rMesh.GetNode(rMesh.GetElement(*iter)->GetNodeGlobalIndex(i));
                          unsigned node_index_i  = pNode->GetIndex();
                          NeighboursToNode.push_back(node_index_i);
                    }
                    
                }

                // Finished iterating over each of the elements for this node, and have everything saved in a long vector, 
                // Now I need to loop over the neighbours vector and remove anything repeting and in node in question 
 
                sort(NeighboursToNode.begin(), NeighboursToNode.end());
                NeighboursToNode.erase(unique(NeighboursToNode.begin(), NeighboursToNode.end()), NeighboursToNode.end());
               int  NeighbourNodeLenght = NeighboursToNode.size();

                // remove the cell of interest from its own neighbour list -- this will find the element in the vector that has the same index as the node in question, and then erase it 
                // previously had this code running using a int to loop over vector, but caused seg 11 fault down the line because apparently i< NeighboursToNode.size() doesnt return an int
                for(std::vector<unsigned>::iterator it = NeighboursToNode.begin(); it != NeighboursToNode.end(); ++it)
                   {
                       if (*it == node_index)
                        {   
                            NeighboursToNode.erase(it);
                            break;
                        }
                    }
            // mNeighbours[node_index] =  NeighboursToNode;
            std::set<unsigned> neighbouring_node_indices(NeighboursToNode.begin(), NeighboursToNode.end());

            // Save the neighbouring nodes as a set, which can be pushed back into the node_neighbour vector
            node_neighbours.push_back(neighbouring_node_indices);

            // Check if node is on the boundary
            bool is_boundary_node = node_iter->IsBoundaryNode();

            // Get the node location 
            c_vector<double, SPACE_DIM> NodeLocation = node_iter->rGetLocation() ;

            // Push the node back into the set 
            nodes.push_back(new Node<SPACE_DIM>(node_index, NodeLocation, is_boundary_node)); 


        }


        std::vector<PottsElement<SPACE_DIM>*> elements;
        // TRACE("Construct the new Potts Mesh");
        mpMesh = new PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>(nodes, elements, node_neighbours, node_neighbours, &rMesh);          
 }

template<unsigned SPACE_DIM>
PottsMeshFromMutableMeshGeneratorJess<SPACE_DIM>::~PottsMeshFromMutableMeshGeneratorJess()
{
    delete mpMesh;
}

template<unsigned SPACE_DIM>
PottsArbitrarySurfaceIn3DMesh<SPACE_DIM>* PottsMeshFromMutableMeshGeneratorJess<SPACE_DIM>::GetMesh()
{
    return mpMesh;
}

// template class PottsMeshFromMutableMeshGeneratorJess<2>;
template class PottsMeshFromMutableMeshGeneratorJess<3>;

// #include "SerializationExportWrapperForCpp.hpp"
// EXPORT_TEMPLATE_CLASS1(PottsMeshFromMutableMeshGeneratorJess, 3);
// EXPORT_TEMPLATE_CLASS1(PottsMeshFromMutableMeshGeneratorJess, 2);

    // for (typename MutableMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator node_iter = rMesh.GetNodeIteratorBegin();
    //             node_iter != rMesh.GetNodeIteratorEnd();
    //             ++node_iter)
    //         {
    //             // node_indices.insert(node_iter->GetIndex());
    //             unsigned node_index2 = node_iter->GetIndex();

    //             std::set<unsigned> SOmeneighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index2);
    //             node_neighbours.push_back(SOmeneighbouring_node_indices );

    //             // Check is boundary node 
    //             bool is_boundary_node2 = GetNode(node_index2)->IsBoundaryNode();

    //             c_vector<double, SPACE_DIM> NodeLocation = node_iter->rGetLocation() ;


    //             // nodes2.push_back(new Node<SPACE_DIM>(iter->GetIndex(), iter->CalculateCentroid(), is_element_boundary)); //(iter->GetIndex(), iter->CalculateCentroid(), is_element_boundary)); 
       
            
    //         }




        // I THINK I WANT SOMETHING TO REPRESENT THE ELEMENTS, ITS A VORINOI THING 
        // Iterate over all of the elements and save the center of the element in a map 
   
            // for (typename MutableMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator elem_iter = rMesh.GetElementIteratorBegin();
            // elem_iter != rMesh.GetElementIteratorEnd();
            // ++elem_iter)
            // {
        
            //     // Calculate the center point for each node 
            //     unsigned Element_index = elem_iter->GetIndex();
            //     c_vector<double , 3> CenterPoint = Create_c_vector(0,0,0);
            //     c_vector<double , 3> NodeIncides;

            //     for (int i=0; i<3; i++)
            //     {
            //         Node<SPACE_DIM>* pNode = rMesh.GetNode(rMesh.GetElement(Element_index)->GetNodeGlobalIndex(i));
            //         CenterPoint += pNode->rGetLocation();
            //         NodeIncides[i] = pNode->GetIndex();  
            //     }
            //     MeshELementCenterPoints[Element_index] = CenterPoint /3;
            // }

        