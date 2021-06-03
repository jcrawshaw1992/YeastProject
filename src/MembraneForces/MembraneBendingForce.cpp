
#include "MembraneBendingForce.hpp"
#include "UblasCustomFunctions.hpp"
#include "Debug.hpp"
#include "HistoryDepMeshBasedCellPopulation.hpp"

MembraneBendingForce::MembraneBendingForce()
   : AbstractForce<2,3>(),
     mMembraneStiffness(1)
{
}

void MembraneBendingForce::SetMembraneStiffness(double membraneStiffnes)
{
    mMembraneStiffness = membraneStiffnes;
}

// double MembraneBendingForce::GetMembraneStiffness() const
// {
// 	return mMembraneStiffness;
// }

// double MembraneBendingForce::GetOriginalAngle(std::pair<Node<3>*, Node<3>*> edge)
// {
//     return mOriginalAngles.at(std::pair<unsigned,unsigned>(edge.first->GetIndex(),edge.second->GetIndex()));
// }

// double MembraneBendingForce::SetOriginalAngle(std::pair<Node<3>*, Node<3>*> edge, double angle)
// {
//     mOriginalAngles[std::pair<unsigned,unsigned>(edge.first->GetIndex(),edge.second->GetIndex())] = angle;
// }

bool MembraneBendingForce::CalculateElementNormals(MutableMesh<2, 3>& rMesh, std::pair<Node<3>*, Node<3>*> edge,
                                                     std::pair<c_vector<double, 3>, c_vector<double, 3> >& nonUnitNormals,
                                                     std::pair<Node<3>*,  Node<3>*>& otherNodes)
{
    Node<3>* pNode1 = edge.first;
    Node<3>* pNode3 = edge.second;

    /*
     *  Find common triangles
     */
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_node1 = pNode1->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_node3 = pNode3->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_node1.begin(),
                          elements_containing_node1.end(),
                          elements_containing_node3.begin(),
                          elements_containing_node3.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    switch (shared_elements.size())
    {
    case 1:
        // We found a boundary edge and we finish here
        return true;
    case 2:
        break;
    default:
        NEVER_REACHED;
    }

    std::set<unsigned>::iterator set_iter = shared_elements.begin();
    Element<2,3>* pElement1 = rMesh.GetElement(*set_iter);
    ++set_iter;
    Element<2,3>* pElement2 = rMesh.GetElement(*set_iter);

    // Find additional nodes
    Node<3>* pNode2 = NULL;
    Node<3>* pNode4 = NULL;
    for (unsigned local_index = 0; local_index < 3; ++local_index)
    {
        unsigned index_for_node2 = pElement1->GetNodeGlobalIndex(local_index);
        unsigned index_for_node4 = pElement2->GetNodeGlobalIndex(local_index);

        if ((index_for_node2 != pNode1->GetIndex()) && (index_for_node2 != pNode3->GetIndex()))
        {
            pNode2 = pElement1->GetNode(local_index);
        }

        if ((index_for_node4 != pNode1->GetIndex()) && (index_for_node4 != pNode3->GetIndex()))
        {
            pNode4 = pElement2->GetNode(local_index);
        }
    }
    assert(pNode2 != NULL);
    assert(pNode4 != NULL);

    // Calculate the force acting on each node
    c_vector<double, 3> vector_A = pNode1->rGetLocation() - pNode3->rGetLocation();
    c_vector<double, 3> vector_B = pNode2->rGetLocation() - pNode3->rGetLocation();
    c_vector<double, 3> normal_1 = VectorProduct(vector_A,vector_B);

    vector_A = pNode4->rGetLocation() - pNode3->rGetLocation();
    vector_B = pNode1->rGetLocation() - pNode3->rGetLocation();
    c_vector<double, 3> normal_2 = VectorProduct(vector_A,vector_B);

    nonUnitNormals = std::pair<c_vector<double, 3>, c_vector<double, 3> >(normal_1, normal_2);
    otherNodes = std::pair<Node<3>*,  Node<3>*>(pNode2, pNode4);

    return false;
}



void MembraneBendingForce::AddForceContribution(AbstractCellPopulation<2,3>& rCellPopulation)
{
    HistoryDepMeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<HistoryDepMeshBasedCellPopulation<2, 3>*>(&rCellPopulation);

    // Iterate over all springs and add force contributions
    for (MeshBasedCellPopulation<2,3>::SpringIterator spring_iterator = p_cell_population->SpringsBegin();
         spring_iterator != p_cell_population->SpringsEnd();
         ++spring_iterator)
    {
        Node<3>* pNode1 = spring_iterator.GetNodeA();
        Node<3>* pNode3 = spring_iterator.GetNodeB();
       
        std::pair<Node<3>*, Node<3>*> edge = std::pair<Node<3>*, Node<3>*>(pNode1, pNode3);

        std::pair<c_vector<double, 3>, c_vector<double, 3> > nonUnitNormals;
        std::pair<Node<3>*, Node<3>*> otherNodes;

        bool boundary_edge_found = CalculateElementNormals(p_cell_population->rGetMesh(),
                                                           edge, nonUnitNormals, otherNodes);

        if (boundary_edge_found)
        {
            continue;
        }

        Node<3>* pNode2 = otherNodes.first;
        Node<3>* pNode4 = otherNodes.second;


        unsigned node_index1 = pNode1->GetIndex();
        unsigned node_index2 = pNode2->GetIndex();
        unsigned node_index3 = pNode3->GetIndex();
        unsigned node_index4 = pNode4->GetIndex();
        

        // unsigned node_index;
        CellPtr p_cell1 = p_cell_population->GetCellUsingLocationIndex(node_index1);  //MembraneStiffness += p_cell1->GetCellData()->GetItem("BendingConstant");
        CellPtr p_cell2 = p_cell_population->GetCellUsingLocationIndex(node_index2);   //MembraneStiffness += p_cell2->GetCellData()->GetItem("BendingConstant");
        CellPtr p_cell3 = p_cell_population->GetCellUsingLocationIndex(node_index3); // MembraneStiffness += p_cell3->GetCellData()->GetItem("BendingConstant");
        CellPtr p_cell4 = p_cell_population->GetCellUsingLocationIndex(node_index4);  //MembraneStiffness += p_cell4->GetCellData()->GetItem("BendingConstant");
 
        double  MembraneStiffness = p_cell1->GetCellData()->GetItem("BendingConstant"); 
        MembraneStiffness += p_cell2->GetCellData()->GetItem("BendingConstant");
        MembraneStiffness += p_cell3->GetCellData()->GetItem("BendingConstant");
        MembraneStiffness += p_cell4->GetCellData()->GetItem("BendingConstant"); 
        MembraneStiffness/=4;

        c_vector<double, 3> normal_1 = nonUnitNormals.first;
        double area_1 = 0.5 * norm_2(normal_1);
        normal_1 /= norm_2(normal_1);
   
        c_vector<double, 3> normal_2 = nonUnitNormals.second;
        double area_2 = 0.5 * norm_2(normal_2);
        normal_2 /= norm_2(normal_2);

 
        c_vector<double, 3> projection_21 = normal_2 - inner_prod(normal_1, normal_2) * normal_1;
        c_vector<double, 3> projection_12 = normal_1 - inner_prod(normal_1, normal_2) * normal_2;
        if (norm_2(projection_21) >1e-12)// If normals are parallel then orthogonal projections are zero.
        {
            projection_21 /= norm_2(projection_21);// Not sure why this is here???? Its not in my calculations in the confirmation
        }
        
        
        if (norm_2(projection_12) >1e-12)
        {
            projection_12 /= norm_2(projection_12);
        }

        c_vector<double, 3> vector_23 = pNode2->rGetLocation() - pNode3->rGetLocation();
        c_vector<double, 3> vector_34 = pNode3->rGetLocation() - pNode4->rGetLocation();
        c_vector<double, 3> node1_contribution = 1 / (2*area_1) * VectorProduct(vector_23, projection_21) +
                                                 1 / (2*area_2) * VectorProduct(vector_34, projection_12);

        c_vector<double, 3> vector_31 = pNode3->rGetLocation() - pNode1->rGetLocation();
        c_vector<double, 3> node2_contribution = 1 / (2*area_1) * VectorProduct(vector_31, projection_21);

        c_vector<double, 3> vector_12 = pNode1->rGetLocation() - pNode2->rGetLocation();
        c_vector<double, 3> vector_41 = pNode4->rGetLocation() - pNode1->rGetLocation();
        c_vector<double, 3> node3_contribution = 1 / (2*area_1) * VectorProduct(vector_12, projection_21) +
                                                 1 / (2*area_2) * VectorProduct(vector_41, projection_12);

        c_vector<double, 3> vector_13 = pNode1->rGetLocation() - pNode3->rGetLocation();
        c_vector<double, 3> node4_contribution = 1 / (2*area_2) * VectorProduct(vector_13, projection_12);

        double OriginalAngle = p_cell_population->GetOriginalAngle(edge);
        assert(OriginalAngle != DOUBLE_UNSET);
        double CurrentAngle = acos(inner_prod(normal_1, normal_2));
        
        
        double force_coefficient = MembraneStiffness * (acos(inner_prod(normal_1, normal_2)) - OriginalAngle);
        // force_coefficient  /=sqrt(1-inner_prod(normal_1, normal_2)*inner_prod(normal_1, normal_2));


        node1_contribution *= force_coefficient;///rCellPopulation.GetVolumeOfCell(p_cell2);
        node2_contribution *= force_coefficient;////rCellPopulation.GetVolumeOfCell(p_cell2);
        node3_contribution *= force_coefficient;////rCellPopulation.GetVolumeOfCell(p_cell3);
        node4_contribution *= force_coefficient;////rCellPopulation.GetVolumeOfCell(p_cell4);


        double Boundary1 = p_cell1->GetCellData()->GetItem("Boundary"); double Boundary2 = p_cell2->GetCellData()->GetItem("Boundary");
        double Boundary3 = p_cell3->GetCellData()->GetItem("Boundary"); double Boundary4 = p_cell4->GetCellData()->GetItem("Boundary");
        // Add the force contribution to each node
        if (Boundary1 ==0 && Boundary2 ==0 && Boundary3 ==0 && Boundary4 ==0)
        {
            // assert(norm_2(node1_contribution)<10000);
            // assert(norm_2(node2_contribution)<10000);
            // assert(norm_2(node3_contribution)<10000);
            // assert(norm_2(node4_contribution)<10000);

            p_cell1->GetCellData()->SetItem("BendingForce", norm_2(node1_contribution));
            p_cell2->GetCellData()->SetItem("BendingForce", norm_2(node2_contribution));
            p_cell3->GetCellData()->SetItem("BendingForce", norm_2(node3_contribution));
            p_cell4->GetCellData()->SetItem("BendingForce", norm_2(node4_contribution));


            if (norm_2(node1_contribution)>1000 || norm_2(node2_contribution)>1000 ||norm_2(node3_contribution)>1000 || norm_2(node4_contribution)>1000  )
            {
                PRINT_VECTOR(node1_contribution);
                PRINT_2_VARIABLES(force_coefficient,(acos(inner_prod(normal_1, normal_2)) - OriginalAngle));
                

                // PRINT_4_VARIABLES(rCellPopulation.GetVolumeOfCell(p_cell1),rCellPopulation.GetVolumeOfCell(p_cell2),rCellPopulation.GetVolumeOfCell(p_cell3),rCellPopulation.GetVolumeOfCell(p_cell4) )
                    node1_contribution =Create_c_vector(0,0,0);
                    node2_contribution= Create_c_vector(0,0,0);
                    node3_contribution= Create_c_vector(0,0,0);
                    node4_contribution= Create_c_vector(0,0,0);
                    p_cell1->GetCellData()->SetItem("BendingForce", -100);
                    p_cell2->GetCellData()->SetItem("BendingForce", -100);
                    p_cell3->GetCellData()->SetItem("BendingForce", -100);
                    p_cell4->GetCellData()->SetItem("BendingForce", -100);
            }

            pNode1->AddAppliedForceContribution(node1_contribution);
            pNode2->AddAppliedForceContribution(node2_contribution);
            pNode3->AddAppliedForceContribution(node3_contribution);
            pNode4->AddAppliedForceContribution(node4_contribution);
            
        }
        else{
        p_cell1->GetCellData()->SetItem("BendingForce", 0);
        p_cell2->GetCellData()->SetItem("BendingForce", 0);
        p_cell3->GetCellData()->SetItem("BendingForce", 0);
        p_cell4->GetCellData()->SetItem("BendingForce", 0);

        }
        

        
        // if (Boundary2 ==0)
        // {
        //     pNode2->AddAppliedForceContribution(node2_contribution);
        // }
        
        // if (Boundary3 ==0)
        // {
        //     pNode3->AddAppliedForceContribution(node3_contribution);
        // }
        
        // if (Boundary4 ==0)
        // {
        //     pNode4->AddAppliedForceContribution(node4_contribution);
        // }
        
        // pNode2->AddAppliedForceContribution(node2_contribution);
        // pNode3->AddAppliedForceContribution(node3_contribution);
        // pNode4->AddAppliedForceContribution(node4_contribution);
    }
}


void MembraneBendingForce::OutputForceParameters(out_stream& rParamsFile)
{
//    *rParamsFile << "\t\t\t<UseCutOffLength>" << mUseCutOffLength << "</UseCutOffLength>\n";
//    *rParamsFile << "\t\t\t<CutOffLength>" << mMechanicsCutOffLength << "</CutOffLength>\n";

    // Call method on direct parent class
    AbstractForce<2,3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MembraneBendingForce)






// -------------------


// void MembraneBendingForce::SetupInitialMembrane(MutableMesh<2,3>& rMesh, AbstractCellPopulation<2, 3>& rCellPopulation)
// {
//      MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);
//     for (MutableMesh<2,3>::EdgeIterator edge_iterator = rMesh.EdgesBegin();
//         edge_iterator != rMesh.EdgesEnd();
//         ++edge_iterator)
//     {
//         std::pair<c_vector<double, 3>, c_vector<double, 3> > nonUnitNormals;
//         std::pair<Node<3>*,  Node<3>*> otherNodes;

//         std::pair<Node<3>*, Node<3>*> edge = std::pair<Node<3>*, Node<3>*>(edge_iterator.GetNodeA(), edge_iterator.GetNodeB());

//         bool boundary_edge_found = CalculateElementNormalsInital(rMesh, edge, nonUnitNormals, otherNodes, p_cell_population);

//         if (boundary_edge_found)
//         {
//             SetOriginalAngle(edge, DOUBLE_UNSET);
//         }
//         else
//         {

//             c_vector<double, 3> normal_1 =  nonUnitNormals.first;
//             normal_1 /= norm_2(normal_1);
//             c_vector<double, 3> normal_2 = nonUnitNormals.second;
//             normal_2 /= norm_2(normal_2);

//             SetOriginalAngle(edge, acos(inner_prod(normal_1, normal_2)));
//         }
//     }
// }

// bool MembraneBendingForce::CalculateElementNormalsInital(MutableMesh<2, 3>& rMesh, std::pair<Node<3>*, Node<3>*> edge,
//                                                      std::pair<c_vector<double, 3>, c_vector<double, 3> >& nonUnitNormals,
//                                                      std::pair<Node<3>*,  Node<3>*>& otherNodes,
//                                                      MeshBasedCellPopulation<2, 3>* p_cell_population)
// {
//     Node<3>* pNode1 = edge.first;
//     Node<3>* pNode3 = edge.second;

//     /*
//      *  Find common triangles
//      */
//     // Find the indices of the elements owned by each node
//     std::set<unsigned> elements_containing_node1 = pNode1->rGetContainingElementIndices();
//     std::set<unsigned> elements_containing_node3 = pNode3->rGetContainingElementIndices();

//     // Find common elements
//     std::set<unsigned> shared_elements;
//     std::set_intersection(elements_containing_node1.begin(),
//                           elements_containing_node1.end(),
//                           elements_containing_node3.begin(),
//                           elements_containing_node3.end(),
//                           std::inserter(shared_elements, shared_elements.begin()));

//     // Check that the nodes have a common edge
//     switch (shared_elements.size())
//     {
//     case 1:
//         // We found a boundary edge and we finish here
//         return true;
//     case 2:
//         break;
//     default:
//         NEVER_REACHED;
//     }

//     std::set<unsigned>::iterator set_iter = shared_elements.begin();
//     Element<2,3>* pElement1 = rMesh.GetElement(*set_iter);
//     ++set_iter;
//     Element<2,3>* pElement2 = rMesh.GetElement(*set_iter);

//     // Find additional nodes
//     Node<3>* pNode2 = NULL;
//     Node<3>* pNode4 = NULL;
//     for (unsigned local_index = 0; local_index < 3; ++local_index)
//     {
//         unsigned index_for_node2 = pElement1->GetNodeGlobalIndex(local_index);
//         unsigned index_for_node4 = pElement2->GetNodeGlobalIndex(local_index);

//         if ((index_for_node2 != pNode1->GetIndex()) && (index_for_node2 != pNode3->GetIndex()))
//         {
//             pNode2 = pElement1->GetNode(local_index);
//         }

//         if ((index_for_node4 != pNode1->GetIndex()) && (index_for_node4 != pNode3->GetIndex()))
//         {
//             pNode4 = pElement2->GetNode(local_index);
//         }
//     }
//     assert(pNode2 != NULL);
//     assert(pNode4 != NULL);

//     // Calculate the force acting on each node


//             CellPtr p_cell1 = p_cell_population->GetCellUsingLocationIndex( pNode1->GetIndex());
//             CellPtr p_cell2 = p_cell_population->GetCellUsingLocationIndex( pNode2->GetIndex());
//             CellPtr p_cell3 = p_cell_population->GetCellUsingLocationIndex( pNode3->GetIndex());
//             CellPtr p_cell4 = p_cell_population->GetCellUsingLocationIndex( pNode4->GetIndex());

//             c_vector<double, 3> Location_1 = pNode1->rGetLocation();
//             c_vector<double, 3> Location_2 = pNode2->rGetLocation();
//             c_vector<double, 3> Location_3 = pNode3->rGetLocation();
//             c_vector<double, 3> Location_4 = pNode4->rGetLocation();


//             // Location_1[0] = p_cell1->GetCellData()->GetItem("Initial_Location_X");  Location_1[1] = p_cell1->GetCellData()->GetItem("Initial_Location_Y");  Location_1[2] = p_cell1->GetCellData()->GetItem("Initial_Location_Z");

//             // c_vector<double, 3> Location_2;
//             // Location_2[0] = p_cell2->GetCellData()->GetItem("Initial_Location_X");  Location_2[1] = p_cell2->GetCellData()->GetItem("Initial_Location_Y");  Location_2[2] = p_cell2->GetCellData()->GetItem("Initial_Location_Z");

//             // c_vector<double, 3> Location_3;
//             // Location_3[0] = p_cell3->GetCellData()->GetItem("Initial_Location_X");  Location_3[1] = p_cell3->GetCellData()->GetItem("Initial_Location_Y");  Location_3[2] = p_cell3->GetCellData()->GetItem("Initial_Location_Z");

//             // c_vector<double, 3> Location_4;
//             // Location_4[0] = p_cell4->GetCellData()->GetItem("Initial_Location_X");  Location_4[1] = p_cell4->GetCellData()->GetItem("Initial_Location_Y");  Location_4[2] = p_cell4->GetCellData()->GetItem("Initial_Location_Z");

//             c_vector<double, 3> vector_A = Location_1 - Location_3;
//             c_vector<double, 3> vector_B = Location_2 - Location_3;
//             c_vector<double, 3> normal_1 = VectorProduct(vector_A,vector_B);

//             vector_A = Location_4 - Location_3;
//             vector_B = Location_1 - Location_3;
//             c_vector<double, 3> normal_2 = VectorProduct(vector_A,vector_B);

//             nonUnitNormals = std::pair<c_vector<double, 3>, c_vector<double, 3> >(normal_1, normal_2);
//             otherNodes = std::pair<Node<3>*,  Node<3>*>(pNode2, pNode4);

//             return false;
// }
