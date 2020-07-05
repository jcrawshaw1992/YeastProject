#include "MembraneBendingForce.hpp"

// THis is the bending force!!!!
// The bending force will not work if the normals are in the same direction, and
MembraneBendingForce::MembraneBendingForce()
        : AbstractForce<2, 3>()
{
}

void MembraneBendingForce::SetMembraneStiffness(double membraneStiffness, double Nc, double Nz)
{
    mMembraneStiffness = membraneStiffness;
    mNc = Nc;
    mNz = Nz;
}

void MembraneBendingForce::AddForceContribution(AbstractCellPopulation<2, 3>& rCellPopulation)
{
     std::map<unsigned, c_vector<double, 3> > BendingForceMap;
     HistoryDepMeshBasedCellPopulation<2, 3>* p_static_cast_cell_population = static_cast<HistoryDepMeshBasedCellPopulation<2, 3>*>(&rCellPopulation);

    for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {

        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        BendingForceMap[node_index] = Create_c_vector(0, 0, 0);
        cell_iter->GetCellData()->SetItem("BendingForce", 0);
    }

    // Throw an exception message if not using a subclass of MeshBasedCellPopulation
    if (dynamic_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("MembraneBendingForce is to be used with subclasses of MeshBasedCellPopulation<2,3> only");
    }
    // This loop saves the area for each cell so then we can divide the force by the area, therby normalising the drag to the nodes in the geometry
  
    // Iterate over all springs and add force contributions
    for (typename MeshBasedCellPopulation<2, 3>::SpringIterator spring_iterator = p_static_cast_cell_population->SpringsBegin();
         spring_iterator != p_static_cast_cell_population->SpringsEnd();
         ++spring_iterator)
    {

        
        Node<3>* pNode1 = spring_iterator.GetNodeA();
        Node<3>* pNode3 = spring_iterator.GetNodeB();

        std::pair<Node<3>*, Node<3>*> edge = std::pair<Node<3>*, Node<3>*>(pNode1, pNode3);
        double OriginalAngle = p_static_cast_cell_population->GetOriginalAngle(edge);

        std::pair<c_vector<double, 3>, c_vector<double, 3> > nonUnitNormals;
        
        std::pair<Node<3>*, Node<3>*> otherNodes;

        bool boundary_edge_found = p_static_cast_cell_population->CalculateElementNormals( edge, nonUnitNormals, otherNodes);
         
        if (boundary_edge_found)
        {
            continue;
        }
        

        Node<3>* pNode2 = otherNodes.first;
        Node<3>* pNode4 = otherNodes.second;

        c_vector<double, 3> normal_1 = nonUnitNormals.first;
        double area_1 = 0.5 * norm_2(normal_1);
        normal_1 /= norm_2(normal_1);

        c_vector<double, 3> normal_2 = nonUnitNormals.second;
        double area_2 = 0.5 * norm_2(normal_2);
        normal_2 /= norm_2(normal_2);
        double NormalsDot = inner_prod(normal_1, normal_2);
        if (std::isnan(NormalsDot))
        {
            PRINT_VECTOR(normal_1)
            PRINT_VECTOR(normal_2)
            PRINT_2_VARIABLES(area_1, area_2)
        }
        
        double angle = acos(NormalsDot);
        if (NormalsDot == 1)
        {
            angle = 0;
        }
        if (std::isnan(angle))
        {
            angle = 0;
        }
        double On = 0;
        if (abs(OriginalAngle- angle) < 1e-10)
        {
            On = 1;
            continue;
        }
        if (On == 1)
        {
            TRACE("Didint Jump")
        }

        unsigned node_index1 = pNode1->GetIndex();
        unsigned node_index2 = pNode2->GetIndex();
        unsigned node_index3 = pNode3->GetIndex();
        unsigned node_index4 = pNode4->GetIndex();

        // // unsigned node_index;
        CellPtr p_cell1 = p_static_cast_cell_population->GetCellUsingLocationIndex(node_index1);
        CellPtr p_cell2 = p_static_cast_cell_population->GetCellUsingLocationIndex(node_index2);
        CellPtr p_cell3 = p_static_cast_cell_population->GetCellUsingLocationIndex(node_index3);
        CellPtr p_cell4 = p_static_cast_cell_population->GetCellUsingLocationIndex(node_index4);

        double MembraneStiffness = mMembraneStiffness;//p_cell1->GetCellData()->GetItem("BendingConstant");

        c_vector<double, 3> projection_21 = normal_2 - inner_prod(normal_1, normal_2) * normal_1;

        if (norm_2(projection_21) > 1e-12) // If normals are parallel then orthogonal projections are zero.
        {
            projection_21 /= norm_2(projection_21);
        }

        c_vector<double, 3> projection_12 = normal_1 - inner_prod(normal_1, normal_2) * normal_2;
        if (norm_2(projection_12) > 1e-12)
        {
            projection_12 /= norm_2(projection_12);
        }

        c_vector<double, 3> vector_23 = pNode2->rGetLocation() - pNode3->rGetLocation();
        c_vector<double, 3> vector_34 = pNode3->rGetLocation() - pNode4->rGetLocation();

        c_vector<double, 3> node1_contribution = 1 / (2 * area_1) * VectorProduct(vector_23, projection_21) + 1 / (2 * area_2) * VectorProduct(vector_34, projection_12);

        c_vector<double, 3> vector_31 = pNode3->rGetLocation() - pNode1->rGetLocation();
        c_vector<double, 3> node2_contribution = 1 / (2 * area_1) * VectorProduct(vector_31, projection_21);

        c_vector<double, 3> vector_12 = pNode1->rGetLocation() - pNode2->rGetLocation();
        c_vector<double, 3> vector_41 = pNode4->rGetLocation() - pNode1->rGetLocation();

        c_vector<double, 3> node3_contribution = 1 / (2 * area_1) * VectorProduct(vector_12, projection_21) + 1 / (2 * area_2) * VectorProduct(vector_41, projection_12);
        // PRINT_2_VARIABLES(area_1, area_2);
        c_vector<double, 3> vector_13 = pNode1->rGetLocation() - pNode3->rGetLocation();
        c_vector<double, 3> node4_contribution = 1 / (2 * area_2) * VectorProduct(vector_13, projection_12);

        assert(OriginalAngle != DOUBLE_UNSET);
        double force_coefficient = MembraneStiffness * (angle - OriginalAngle) / sqrt(1 - inner_prod(normal_1, normal_2) * inner_prod(normal_1, normal_2));
        // PRINT_3_VARIABLES(angle,OriginalAngle, (angle - OriginalAngle))
        // correct for the area of each cell, can do this for each element contribution to the cell here rather than later as

        node1_contribution *= (force_coefficient/ rCellPopulation.GetVolumeOfCell(p_cell1));
        node2_contribution *= (force_coefficient/ rCellPopulation.GetVolumeOfCell(p_cell2));
        node3_contribution *= (force_coefficient/ rCellPopulation.GetVolumeOfCell(p_cell3));
        node4_contribution *= (force_coefficient/ rCellPopulation.GetVolumeOfCell(p_cell4));


        if (abs(norm_2(node1_contribution)) > 1000 || std::isnan(norm_2(node1_contribution)))
        {
            node1_contribution = Create_c_vector(0, 0, 0);
        }
        if (abs(norm_2(node2_contribution)) > 1000 || std::isnan(norm_2(node2_contribution)))
        {
            node2_contribution = Create_c_vector(0, 0, 0);
        }
        if (abs(norm_2(node3_contribution)) > 1000 || std::isnan(norm_2(node3_contribution)))
        {
            node3_contribution = Create_c_vector(0, 0, 0);
        }
        if (abs(norm_2(node4_contribution)) > 1000 || std::isnan(norm_2(node4_contribution)))
        {
            node4_contribution = Create_c_vector(0, 0, 0);
        }

        if (p_cell1->GetCellData()->GetItem("Boundary") == 0)
        {
            pNode1->AddAppliedForceContribution(node1_contribution ); // Add the new force
            double BendingForce = p_cell1->GetCellData()->GetItem("BendingForce");
            p_cell1->GetCellData()->SetItem("BendingForce", norm_2(node1_contribution) + BendingForce);
        }
          
        if (p_cell2->GetCellData()->GetItem("Boundary") == 0)
        {
             pNode2->AddAppliedForceContribution(node2_contribution ); // Add the new force
             double BendingForce = p_cell2->GetCellData()->GetItem("BendingForce");
             p_cell2->GetCellData()->SetItem("BendingForce", norm_2(node2_contribution) + BendingForce);

        }if (p_cell3->GetCellData()->GetItem("Boundary") == 0)
        {
             pNode3->AddAppliedForceContribution(node3_contribution ); // Add the new force
             double BendingForce = p_cell3->GetCellData()->GetItem("BendingForce");
             p_cell3->GetCellData()->SetItem("BendingForce", norm_2(node3_contribution) + BendingForce);

        }
        if (p_cell4->GetCellData()->GetItem("Boundary") == 0)
        {
             pNode4->AddAppliedForceContribution(node4_contribution ); // Add the new force
             double BendingForce = p_cell4->GetCellData()->GetItem("BendingForce");
             p_cell4->GetCellData()->SetItem("BendingForce", norm_2(node4_contribution) + BendingForce);
        }
         

        // BendingForceMap[node_index1] += node1_contribution;
        // BendingForceMap[node_index2] += node2_contribution;
        // BendingForceMap[node_index3] += node3_contribution;
        // // BendingForceMap[node_index4] += node4_contribution;
        // pNode1->AddAppliedForceContribution(node1_contribution ); // Add the new force
        // pNode2->AddAppliedForceContribution(node2_contribution ); // Add the new force
        // pNode3->AddAppliedForceContribution(node3_contribution ); // Add the new force
        // pNode4->AddAppliedForceContribution(node4_contribution ); // Add the new force
    //         

        



    }


    
}

void MembraneBendingForce::SetNearestNodesForBoundaryNodesBending(std::map<unsigned, c_vector<unsigned, 5> > NearestNodesMap)
{

    mNearestNodesMap = NearestNodesMap;
     
}

void MembraneBendingForce::OutputForceParameters(out_stream & rParamsFile)
{
    //    *rParamsFile << "\t\t\t<UseCutOffLength>" << mUseCutOffLength << "</UseCutOffLength>\n";
    //    *rParamsFile << "\t\t\t<CutOffLength>" << mMechanicsCutOffLength << "</CutOffLength>\n";

    // Call method on direct parent class
    AbstractForce<2, 3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MembraneBendingForce)
