
#include "MembraneBendingForce0TargetAngle.hpp"
#include "UblasCustomFunctions.hpp"
#include "Debug.hpp"
#include "HistoryDepMeshBasedCellPopulation.hpp"

MembraneBendingForce0TargetAngle::MembraneBendingForce0TargetAngle()
   : AbstractForce<2,3>(),
     mMembraneStiffness(1)
{
}

void MembraneBendingForce0TargetAngle::SetMembraneStiffness(double membraneStiffnes)
{
    mMembraneStiffness = membraneStiffnes;
}

bool MembraneBendingForce0TargetAngle::CalculateElementNormals(MutableMesh<2, 3>& rMesh, std::pair<Node<3>*, Node<3>*> edge,
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



void MembraneBendingForce0TargetAngle::AddForceContribution(AbstractCellPopulation<2,3>& rCellPopulation)
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


        double  MembraneStiffness = mMembraneStiffness;// THis needs updating 
        // unsigned node_index;
        CellPtr p_cell1 = p_cell_population->GetCellUsingLocationIndex(node_index1);  //MembraneStiffness += p_cell1->GetCellData()->GetItem("BendingConstant");
        CellPtr p_cell2 = p_cell_population->GetCellUsingLocationIndex(node_index2);   //MembraneStiffness += p_cell2->GetCellData()->GetItem("BendingConstant");
        CellPtr p_cell3 = p_cell_population->GetCellUsingLocationIndex(node_index3); // MembraneStiffness += p_cell3->GetCellData()->GetItem("BendingConstant");
        CellPtr p_cell4 = p_cell_population->GetCellUsingLocationIndex(node_index4);  //MembraneStiffness += p_cell4->GetCellData()->GetItem("BendingConstant");

        c_vector<double, 3> normal_1 = nonUnitNormals.first;
        assert( std::isnan(normal_1[0])==0 & std::isnan(normal_1[1]) ==0 &std::isnan(normal_1[2])==0  );
        double area_1 = 0.5 * norm_2(normal_1);
        assert( std::isnan(area_1) ==0);
        normal_1 /= norm_2(normal_1);
        assert(std::isnan(norm_2(normal_1))==0);

        c_vector<double, 3> normal_2 = nonUnitNormals.second;
        double area_2 = 0.5 * norm_2(normal_2);
        normal_2 /= norm_2(normal_2);
        assert( std::isnan(normal_2[0])==0 & std::isnan(normal_2[1]) ==0 &std::isnan(normal_2[2])==0  );
        assert( std::isnan(area_2) ==0);
        assert(std::isnan(norm_2(normal_2))==0);

        c_vector<double, 3> projection_21 = normal_2 - inner_prod(normal_1, normal_2) * normal_1;
        c_vector<double, 3> projection_12 = normal_1 - inner_prod(normal_1, normal_2) * normal_2;
        assert(std::isnan(norm_2(projection_21))==0);
        assert(std::isnan(norm_2(projection_12))==0);

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
        assert(std::isnan(norm_2(node1_contribution))==0);
        c_vector<double, 3> vector_31 = pNode3->rGetLocation() - pNode1->rGetLocation();
        c_vector<double, 3> node2_contribution = 1 / (2*area_1) * VectorProduct(vector_31, projection_21);
         assert(std::isnan(norm_2(node2_contribution))==0);
        c_vector<double, 3> vector_12 = pNode1->rGetLocation() - pNode2->rGetLocation();
        c_vector<double, 3> vector_41 = pNode4->rGetLocation() - pNode1->rGetLocation();
        c_vector<double, 3> node3_contribution = 1 / (2*area_1) * VectorProduct(vector_12, projection_21) +
                                                 1 / (2*area_2) * VectorProduct(vector_41, projection_12);

        c_vector<double, 3> vector_13 = pNode1->rGetLocation() - pNode3->rGetLocation();
        c_vector<double, 3> node4_contribution = 1 / (2*area_2) * VectorProduct(vector_13, projection_12);


        assert(std::isnan(norm_2(node3_contribution))==0);
        assert(std::isnan(norm_2(node4_contribution))==0);
        double OriginalAngle = 0;
        double CurrentAngle = acos(inner_prod(normal_1, normal_2));
        if ( std::isnan(CurrentAngle)  == 1)
        {
            CurrentAngle = 0; // they are parellel
        }
        assert(std::isnan(CurrentAngle)  == 0);

        
        double force_coefficient = MembraneStiffness * (CurrentAngle - OriginalAngle);
  
        node1_contribution *= force_coefficient;
        node2_contribution *= force_coefficient;
        node3_contribution *= force_coefficient;
        node4_contribution *= force_coefficient;



        pNode1->AddAppliedForceContribution(node1_contribution);
        pNode2->AddAppliedForceContribution(node2_contribution);
        pNode3->AddAppliedForceContribution(node3_contribution);
        pNode4->AddAppliedForceContribution(node4_contribution);
         
          }
}


void MembraneBendingForce0TargetAngle::OutputForceParameters(out_stream& rParamsFile)
{
//    *rParamsFile << "\t\t\t<UseCutOffLength>" << mUseCutOffLength << "</UseCutOffLength>\n";
//    *rParamsFile << "\t\t\t<CutOffLength>" << mMechanicsCutOffLength << "</CutOffLength>\n";

    // Call method on direct parent class
    AbstractForce<2,3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MembraneBendingForce0TargetAngle)




