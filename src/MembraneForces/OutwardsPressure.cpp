#include "OutwardsPressure.hpp"

#include "UblasCustomFunctions.hpp"
#include "MathsFunctions.hpp"
#include "HistoryDepMeshBasedCellPopulation.hpp"


OutwardsPressure::OutwardsPressure()
        : AbstractForce<2, 3>()
{
}

void OutwardsPressure::SetPressure(double Pressure)
{
      mStrength =  Pressure;
}

void OutwardsPressure::SetNearestNeighboursMap(std::map<unsigned, c_vector<unsigned, 2> > NearestNodesMap)
{
    mNearestNodesMap = NearestNodesMap;
}


void OutwardsPressure::AddForceContribution(AbstractCellPopulation<2, 3>& rCellPopulation)
{
    HistoryDepMeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<HistoryDepMeshBasedCellPopulation<2, 3>*>(&rCellPopulation);        
    std::map<unsigned, c_vector<double, 3> > ForceMap;
    for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
    {
        cell_iter->GetCellData()->SetItem("FixedBoundary", 0);
        if (cell_iter->GetCellData()->GetItem("Boundary") ==0)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<3>* p_node = rCellPopulation.GetNode(node_index);

            c_vector<double, 3> Normal = zero_vector<double>(3);

            std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();

            assert(containing_elements.size() > 0);
            for (std::set<unsigned>::iterator iter = containing_elements.begin();
                    iter != containing_elements.end();
                    ++iter)
            {
                Node<3>* pNode0 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(0));
                Node<3>* pNode1 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(1));
                Node<3>* pNode2 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(2));

                c_vector<double, 3> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
                c_vector<double, 3> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3

                c_vector<double, 3> normalVector = VectorProduct(vector_12, vector_13);
                Normal += normalVector;
            }
            Normal /=norm_2(Normal);
        
            c_vector<double, 3> force = mStrength *Normal; // / norm_2(cell_location);
            
            // rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(force);
            // cell_iter->GetCellData()->SetItem("OutwardForce", norm_2(force) );
            ForceMap[node_index] = force;
        }
    }
    for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
        cell_iter != rCellPopulation.End();
        ++cell_iter)
    {
        // TRACE("PRES")
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
         if (cell_iter->GetCellData()->GetItem("Boundary") ==1)
        {
            // PRINT_VARIABLE(cell_iter->GetCellData()->GetItem("Boundary"))
            Node<3>* p_node = rCellPopulation.GetNode(node_index);

            c_vector<double, 3> ForceOnEdgeNode = Create_c_vector(0,0,0);
            c_vector<unsigned, 3>  LocalNodes = p_cell_population->GetNearestInternalNodes(node_index);
            for(int i =0 ; i<3; i++)
            {   
                ForceOnEdgeNode += ForceMap[LocalNodes[i]];
            }
            ForceOnEdgeNode/=3;
            
            ForceOnEdgeNode/=norm_2(ForceOnEdgeNode);
            ForceOnEdgeNode*=mStrength;    
        
            assert(norm_2(ForceOnEdgeNode) <10e10); 
            p_node->AddAppliedForceContribution(ForceOnEdgeNode);
            cell_iter->GetCellData()->SetItem("OutwardForce", norm_2(ForceOnEdgeNode) );

        }else 
        {
            // PRINT_VARIABLE(cell_iter->GetCellData()->GetItem("Boundary"))
            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(ForceMap[node_index]);
            cell_iter->GetCellData()->SetItem("OutwardForce", norm_2(ForceMap[node_index]) );
        }
    }
}
 

void OutwardsPressure::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<2, 3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(OutwardsPressure)


// Code grave yard


    // }



    //  if (bool(dynamic_cast<MeshBasedCellPopulation<2,3>*>(&rCellPopulation)))
    // {
    //     std::map<unsigned, c_vector<double, 3> > ForceMap;
    //     MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);        
    //     for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
    //          cell_iter != rCellPopulation.End();
    //          ++cell_iter)
    //     {
    //         ForceMap[rCellPopulation.GetLocationIndexUsingCell(*cell_iter)] = Create_c_vector(0,0,0);
    //     }
     
    //     for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
    //          cell_iter != rCellPopulation.End();
    //          ++cell_iter)
    //     {

    //         unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
    //         Node<3>* p_node = rCellPopulation.GetNode(node_index);

    //         c_vector<double, 3> Normal = zero_vector<double>(3);
 
    //         std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();
    //         // PRINT_2_VARIABLES(node_index,containing_elements.size())
    //         // cell_iter->GetCellData()->SetItem("Node", node_index);

    //     // TRACE("do we get here")

    //         assert(containing_elements.size() > 0);
    //         for (std::set<unsigned>::iterator iter = containing_elements.begin();
    //              iter != containing_elements.end();
    //              ++iter)
    //         {
    //             Node<3>* pNode0 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(0));
    //             Node<3>* pNode1= p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(1));
    //             Node<3>* pNode2 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(2));

    //             c_vector<double, 3> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
    //             c_vector<double, 3> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3

    //             c_vector<double, 3> normalVector = VectorProduct(vector_12, vector_13);
                
    //             Normal += normalVector;///norm_2(normalVector);

    //         }
    //         Normal /=norm_2(Normal);
        
    //         c_vector<double, 3> force = mStrength *Normal; // / norm_2(cell_location);
    //         // rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(force);
    //         ForceMap[node_index] = force;
    
    //     }


    //     for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
    //         cell_iter != rCellPopulation.End();
    //         ++cell_iter)
    //     {   
            
    //         unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
    //         Node<3>* p_node = rCellPopulation.GetNode(node_index);
    //         if (cell_iter->GetCellData()->GetItem("Boundary") ==1)
    //         {
    //             c_vector<double, 3> ForceOnEdgeNode = Create_c_vector(0,0,0);
    //             c_vector<unsigned, 2>  LocalNodes = mNearestNodesMap[node_index];
    //             // PRINT_3_VARIABLES(LocalNodses[0], LocalNodes[1], LocalNodes[2])
    //             for(int i =0 ; i<2; i++)
    //             {   
    //                 ForceOnEdgeNode += ForceMap[LocalNodes[i]];
    //             }
    //             ForceOnEdgeNode/=2;
    //             ForceOnEdgeNode = norm_2(ForceOnEdgeNode) * ForceMap[node_index] / norm_2(ForceMap[node_index]);


    //             if (norm_2(ForceOnEdgeNode) >10e10) 
    //              {
    //                  PRINT_VECTOR(ForceOnEdgeNode);
    //                  PRINT_VECTOR(ForceMap[node_index]);
    //                   for(int i =0 ; i<2; i++)
    //                     {   
    //                         PRINT_VECTOR(ForceMap[LocalNodes[i]]);
    //                     }

    //              }

    //             assert(ForceOnEdgeNode[0] <10e10);
    //             assert(ForceOnEdgeNode[1] <10e10);
    //              assert(ForceOnEdgeNode[2] <10e10);
                 


    //             // PRINT_VARIABLE(norm_2(ForceMap[node_index]))
    //             p_node->AddAppliedForceContribution(ForceOnEdgeNode);
    //             cell_iter->GetCellData()->SetItem("OutwardForce", norm_2(ForceOnEdgeNode) );
    //             // PRINT_VECTOR(ForceOnEdgeNode)
    //             // Need to mess around with shit 
    //         }
    //         else
    //         {
    //             assert(ForceMap[node_index][0] <10e10);
    //             assert(ForceMap[node_index][1] <10e10);
    //             assert(ForceMap[node_index][2] <10e10);
    //             p_node->AddAppliedForceContribution(ForceMap[node_index]);
    //             cell_iter->GetCellData()->SetItem("OutwardForce", norm_2(ForceMap[node_index]) );
    //         }
    // //     }
    // // }
    // //   else if (bool(dynamic_cast<HistoryDepMeshBasedCellPopulation<2,3>*>(&rCellPopulation)))
    // // {
    //     std::map<unsigned, c_vector<double, 3> > ForceMap;
    //     HistoryDepMeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<HistoryDepMeshBasedCellPopulation<2, 3>*>(&rCellPopulation);        
    //     for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
    //          cell_iter != rCellPopulation.End();
    //          ++cell_iter)
    //     {
    //         ForceMap[rCellPopulation.GetLocationIndexUsingCell(*cell_iter)] = Create_c_vector(0,0,0);
    //     }