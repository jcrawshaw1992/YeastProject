#include "MembraneStiffnessForce.hpp"
#include "Debug.hpp"
#include "UblasCustomFunctions.hpp"

#include "BetaCateninOneHitCellMutationState.hpp"
// THis is the bending force!!!!
// The bending force will not work if the normals are in the same direction, and
MembraneStiffnessForce::MembraneStiffnessForce()
        : AbstractForce<2, 3>()
{
}

void MembraneStiffnessForce::SetMembraneStiffness(double membraneStiffnes)
{
    mMembraneStiffness = membraneStiffnes;
}

// double MembraneStiffnessForce::GetMembraneStiffness() const
// {
//     return mMembraneStiffness;
// }

double MembraneStiffnessForce::GetOriginalAngle(std::pair<Node<3>*, Node<3>*> edge)
{
    return mOriginalAngles.at(std::pair<unsigned, unsigned>(edge.first->GetIndex(), edge.second->GetIndex()));
}

double MembraneStiffnessForce::SetOriginalAngle(std::pair<Node<3>*, Node<3>*> edge, double angle)
{
    mOriginalAngles[std::pair<unsigned, unsigned>(edge.first->GetIndex(), edge.second->GetIndex())] = angle;
}

void MembraneStiffnessForce::SetScallingBending(double Scalling)
{
    mScalling = Scalling;
}

void MembraneStiffnessForce::SetupInitialMembrane(MutableMesh<2, 3>& rMesh, AbstractCellPopulation<2, 3>& rCellPopulation)
{
    // TRACE("Add Stiffness Force");
    MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);
    // Loop over all nodes and come up with a new inital position
    for (typename AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter); //p_cell_population->GetLocationIndexUsingCell(*cell_iter);
        Node<3>* p_node = rCellPopulation.GetNode(node_index);

        c_vector<double, 3> NodeLocation = p_node->rGetLocation();

        c_vector<long double, 3> normal = zero_vector<long double>(3);

        std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();
        assert(containing_elements.size() > 0);
        for (std::set<unsigned>::iterator iter = containing_elements.begin();
                iter != containing_elements.end();
                ++iter)
        {
            // Negative as normals point inwards for these surface meshes
            normal += - p_cell_population->rGetMesh().GetElement(*iter)->CalculateNormal();
        }
        normal /= norm_2(normal);

      // Need to walk backward into the mesh by the scalling factor 
        c_vector<double, 3> PositionVector = NodeLocation - mScalling * normal;

        (cell_iter)->GetCellData()->SetItem("Initial_Location_X", PositionVector[0]);
        (cell_iter)->GetCellData()->SetItem("Initial_Location_Y", PositionVector[1]);
        (cell_iter)->GetCellData()->SetItem("Initial_Location_Z", PositionVector[2]);
    }

    for (MutableMesh<2, 3>::EdgeIterator edge_iterator = rMesh.EdgesBegin();
         edge_iterator != rMesh.EdgesEnd();
         ++edge_iterator)
    {
        std::pair<c_vector<double, 3>, c_vector<double, 3> > nonUnitNormals;
        std::pair<Node<3>*, Node<3>*> otherNodes;

        std::pair<Node<3>*, Node<3>*> edge = std::pair<Node<3>*, Node<3>*>(edge_iterator.GetNodeA(), edge_iterator.GetNodeB());

        bool boundary_edge_found = CalculateElementNormals(rMesh, edge, nonUnitNormals, otherNodes);

        if (boundary_edge_found)
        {
            SetOriginalAngle(edge, DOUBLE_UNSET);
        }
        else
        { // This is where the angle is set
            Node<3>* pNode1 = edge_iterator.GetNodeA();
            Node<3>* pNode3 = edge_iterator.GetNodeB();
            Node<3>* pNode2 = otherNodes.first;
            Node<3>* pNode4 = otherNodes.second;

            std::vector<double> Node_Indices;
            Node_Indices.push_back(pNode1->GetIndex());
            Node_Indices.push_back(pNode2->GetIndex());
            Node_Indices.push_back(pNode3->GetIndex());
            Node_Indices.push_back(pNode4->GetIndex());

            unsigned Node_1_index = pNode1->GetIndex();
            unsigned Node_2_index = pNode3->GetIndex();

            unsigned Shared_Node_1_index = pNode2->GetIndex();
            unsigned Shared_Node_2_index = pNode4->GetIndex();

            c_vector<c_vector<double, 3>, 4> PositionVectors;
            //  Location of each of the nodes
            for (double i = 0; i < Node_Indices.size(); i++)
            {
                CellPtr p_cell = p_cell_population->GetCellUsingLocationIndex(Node_Indices[i]);
                PositionVectors[i][0] = p_cell->GetCellData()->GetItem("Initial_Location_X");
                PositionVectors[i][1] = p_cell->GetCellData()->GetItem("Initial_Location_Y");
                PositionVectors[i][2] = p_cell->GetCellData()->GetItem("Initial_Location_Z");
            }

            // Element 1
            c_vector<double, 3> Element_1_vector_12 = PositionVectors[0] - PositionVectors[3]; // Vector 1 to 2
            c_vector<double, 3> Element_1_vector_13 = PositionVectors[2] - PositionVectors[3]; // Vector 1 to 3

            c_vector<double, 3> normal_1 = VectorProduct(Element_1_vector_12, Element_1_vector_13);

            // Element 2
            c_vector<double, 3> Element_2_vector_12 = PositionVectors[1] - PositionVectors[0]; // Vector 1 to 2
            c_vector<double, 3> Element_2_vector_13 = PositionVectors[2] - PositionVectors[0]; // Vector 1 to 3

            c_vector<double, 3> normal_2 = VectorProduct(Element_2_vector_12, Element_2_vector_13);

            double DirectionOf_2 = norm_2(normal_2);
            double DirectionOf_1 = norm_2(normal_1);

            c_vector<double, 2> NodeLocation = Create_c_vector(PositionVectors[2][0], PositionVectors[2][2]);
            c_vector<double, 2> Projection = Create_c_vector(normal_1[0], normal_1[1]);

            normal_2 = normal_2 / DirectionOf_2;
            normal_1 = normal_1 / DirectionOf_1;

            // if (DirectionOf_2 < 0)
            // {
            //     normal_2 = -normal_2;
            // }
            // if (DirectionOf_1 < 0)
            // {
            //     normal_1 = -normal_1;
            // }
            double Angle = acos(inner_prod(normal_1, normal_2));
            // double OriginalAngle = GetOriginalAngle(edge);
            // double DifferenceInAngel = OriginalAngle - angle;

            // PRINT_3_VARIABLES(OriginalAngle, angle, DifferenceInAngel );
            SetOriginalAngle(edge, Angle);

            // if (abs(PositionVectors[0][2]) < 2e-3 || abs(PositionVectors[1][2]) < 2e-3 || abs(PositionVectors[2][2]) < 2e-3 || abs(PositionVectors[3][2]) < 2e-3)
            // {
            //     mMembraneStiffnessMap[edge] = mMembraneStiffness*2000;
            // }
            // else
            // {
               mMembraneStiffnessMap[edge] = mMembraneStiffness ;
            // }

            // Check that if mScalling =1 there in no difference in the angle
            //  normal_1 = nonUnitNormals.first;
            // normal_1 /= norm_2(normal_1);
            // normal_2 = nonUnitNormals.second;
            // normal_2 /= norm_2(normal_2);

            // double angle = acos(inner_prod(normal_1, normal_2));
            // SetOriginalAngle(edge, angle);
            // double Difference = abs(Angle - angle);
            // if (Difference>1e-10)
            // {
            //     TRACE("PROBLEM");
            //     PRINT_2_VARIABLES(Angle, angle);
            // }
        }
    }
}

bool MembraneStiffnessForce::CalculateElementNormals(MutableMesh<2, 3>& rMesh, std::pair<Node<3>*, Node<3>*> edge,
                                                     std::pair<c_vector<double, 3>, c_vector<double, 3> >& nonUnitNormals,
                                                     std::pair<Node<3>*, Node<3>*>& otherNodes)
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
    Element<2, 3>* pElement1 = rMesh.GetElement(*set_iter);
    ++set_iter;
    Element<2, 3>* pElement2 = rMesh.GetElement(*set_iter);

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
    c_vector<double, 3> normal_1 = VectorProduct(vector_A, vector_B);

    vector_A = pNode4->rGetLocation() - pNode3->rGetLocation();
    vector_B = pNode1->rGetLocation() - pNode3->rGetLocation();
    c_vector<double, 3> normal_2 = VectorProduct(vector_A, vector_B);

    nonUnitNormals = std::pair<c_vector<double, 3>, c_vector<double, 3> >(normal_1, normal_2);
    otherNodes = std::pair<Node<3>*, Node<3>*>(pNode2, pNode4);

    return false;
}

// ---------------------------
void MembraneStiffnessForce::AddForceContribution(AbstractCellPopulation<2, 3>& rCellPopulation)
{

    // Throw an exception message if not using a subclass of MeshBasedCellPopulation
    if (dynamic_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("MembraneStiffnessForce is to be used with subclasses of MeshBasedCellPopulation<2,3> only");
    }
    // This loop saves the area for each cell so then we can divide the force by the area, therby normalising the drag to the nodes in the geometry
    MeshBasedCellPopulation<2, 3>* p_static_cast_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);

    // Iterate over all springs and add force contributions
    for (typename MeshBasedCellPopulation<2, 3>::SpringIterator spring_iterator = p_static_cast_cell_population->SpringsBegin();
         spring_iterator != p_static_cast_cell_population->SpringsEnd();
         ++spring_iterator)
    {

        unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
        unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
        Node<3>* pNode1 = spring_iterator.GetNodeA();
        Node<3>* pNode3 = spring_iterator.GetNodeB();

        std::pair<Node<3>*, Node<3>*> edge = std::pair<Node<3>*, Node<3>*>(pNode1, pNode3);

        std::pair<c_vector<double, 3>, c_vector<double, 3> > nonUnitNormals;
        std::pair<Node<3>*, Node<3>*> otherNodes;

        bool boundary_edge_found = CalculateElementNormals(p_static_cast_cell_population->rGetMesh(),
                                                           edge, nonUnitNormals, otherNodes);

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
        double angle = acos(NormalsDot);
        if (NormalsDot > 1)
        {
            angle = 0;
        }

        // Calculate the membrane constant here
        // Take the average of the Membrane constants at each of the nodes of the two elements
        // Me = 1/4(M1+M2+M3+M4)

        double MembraneStiffness = 0;
        unsigned node_index1 = pNode1->GetIndex();
        unsigned node_index2 = pNode2->GetIndex();
        unsigned node_index3 = pNode3->GetIndex();
        unsigned node_index4 = pNode4->GetIndex();

        // unsigned node_index;
        CellPtr p_cell1 = p_static_cast_cell_population->GetCellUsingLocationIndex(node_index1);
        // MembraneStiffness += p_cell1->GetCellData()->GetItem("BendingModulus");

        CellPtr p_cell2 = p_static_cast_cell_population->GetCellUsingLocationIndex(node_index2);
        // MembraneStiffness += p_cell2->GetCellData()->GetItem("BendingModulus");

        CellPtr p_cell3 = p_static_cast_cell_population->GetCellUsingLocationIndex(node_index3);
        // MembraneStiffness += p_cell3->GetCellData()->GetItem("BendingModulus");

        CellPtr p_cell4 = p_static_cast_cell_population->GetCellUsingLocationIndex(node_index4);
        // MembraneStiffness += p_cell4->GetCellData()->GetItem("BendingModulus");

        // MembraneStiffness /= 4;

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

        // c_vector<double, 3> VectorA = VectorProduct(vector_23, projection_21) ;

        c_vector<double, 3> node1_contribution = 1 / (2 * area_1) * VectorProduct(vector_23, projection_21) + 1 / (2 * area_2) * VectorProduct(vector_34, projection_12);
        // PRINT_2_VARIABLES(area_1, area_2)

        // PRINT_VECTOR(node1_contribution);
        c_vector<double, 3> vector_31 = pNode3->rGetLocation() - pNode1->rGetLocation();

        c_vector<double, 3> node2_contribution = 1 / (2 * area_1) * VectorProduct(vector_31, projection_21);

        c_vector<double, 3> vector_12 = pNode1->rGetLocation() - pNode2->rGetLocation();
        c_vector<double, 3> vector_41 = pNode4->rGetLocation() - pNode1->rGetLocation();

        c_vector<double, 3> node3_contribution = 1 / (2 * area_1) * VectorProduct(vector_12, projection_21) + 1 / (2 * area_2) * VectorProduct(vector_41, projection_12);
        // PRINT_2_VARIABLES(area_1, area_2);
        c_vector<double, 3> vector_13 = pNode1->rGetLocation() - pNode3->rGetLocation();

        c_vector<double, 3> node4_contribution = 1 / (2 * area_2) * VectorProduct(vector_13, projection_12);

        assert(GetOriginalAngle(edge) != DOUBLE_UNSET);
        double force_coefficient = mMembraneStiffnessMap[edge] * (angle - GetOriginalAngle(edge));

        // correct for the area of each cell, can do this for each element contribution to the cell here rather than later as
        // (F1+F2+F3)/area = F1/area + F2/area + F3/area
        node1_contribution *= force_coefficient / rCellPopulation.GetVolumeOfCell(p_cell1);
        node2_contribution *= force_coefficient / rCellPopulation.GetVolumeOfCell(p_cell2);
        node3_contribution *= force_coefficient / rCellPopulation.GetVolumeOfCell(p_cell3);
        node4_contribution *= force_coefficient / rCellPopulation.GetVolumeOfCell(p_cell4);

        if (angle == 0)
        {
            node1_contribution = Create_c_vector(0, 0, 0);
            node2_contribution = node1_contribution;
            node3_contribution = node1_contribution;
            node4_contribution = node1_contribution;
        }

        // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        pNode1->AddAppliedForceContribution(node1_contribution);
        pNode2->AddAppliedForceContribution(node2_contribution);
        pNode3->AddAppliedForceContribution(node3_contribution);
        pNode4->AddAppliedForceContribution(node4_contribution);
    }
}

void MembraneStiffnessForce::UpdateMembraneStiffnessProperties(AbstractCellPopulation<2, 3>& rCellPopulation)
{

    CellPtr p_cell;
    Node<3>* pNode;

    MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);

    for (typename MeshBasedCellPopulation<2, 3>::SpringIterator spring_iterator = p_cell_population->SpringsBegin();
         spring_iterator != p_cell_population->SpringsEnd();
         ++spring_iterator)
    {
        // double ElementPair = 0;
        // unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
        // unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();
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

        std::vector<double> Node_Indices;
        Node_Indices.push_back(pNode1->GetIndex());
        Node_Indices.push_back(pNode2->GetIndex());
        Node_Indices.push_back(pNode3->GetIndex());
        Node_Indices.push_back(pNode4->GetIndex());

        unsigned Node_1_index = pNode1->GetIndex();
        unsigned Node_2_index = pNode3->GetIndex();

        unsigned Shared_Node_1_index = pNode2->GetIndex();
        unsigned Shared_Node_2_index = pNode4->GetIndex();
        std::map<double, bool> IsMutantMap;
        double ElementPair = 0;
        // PRINT_VARIABLE(Node_Indices.size());
        for (double i = 0; i < Node_Indices.size(); i++)
        {
            p_cell = p_cell_population->GetCellUsingLocationIndex(Node_Indices[i]);
            if (p_cell->GetMutationState()->IsType<LostEndothelialCell>()) // If on edge
            {
                ElementPair += 1;
                IsMutantMap[i + 1] = 1;
            }
            else
            {
                IsMutantMap[i + 1] = 0;
            }
        }

        // PRINT_VARIABLE(ElementPair);
        c_vector<c_vector<double, 3>, 4> PositionVectors;
        if (ElementPair > 0)
        {
            // std::pair<c_vector<double, 3>, c_vector<double, 3> > nonUnitNormals;
            // std::pair<Node<3>*, Node<3>*> edge = std::pair<Node<3>*, Node<3>*>(edge_iterator.GetNodeA(), edge_iterator.GetNodeB());

            //  Location of each of the nodes
            for (double i = 0; i < Node_Indices.size(); i++)
            {
                p_cell = p_cell_population->GetCellUsingLocationIndex(Node_Indices[i]);
                PositionVectors[i][0] = p_cell->GetCellData()->GetItem("Initial_Location_X");
                PositionVectors[i][1] = p_cell->GetCellData()->GetItem("Initial_Location_Y");
                PositionVectors[i][2] = p_cell->GetCellData()->GetItem("Initial_Location_Z");
                // PRINT_VECTOR(PositionVectors[i]);
            }

            // Element 1

            c_vector<double, 3> Element_1_vector_12 = PositionVectors[0] - PositionVectors[3]; // Vector 1 to 2
            c_vector<double, 3> Element_1_vector_13 = PositionVectors[2] - PositionVectors[3]; // Vector 1 to 3

            c_vector<double, 3> normal_1 = VectorProduct(Element_1_vector_12, Element_1_vector_13);

            // Element 2

            c_vector<double, 3> Element_2_vector_12 = PositionVectors[1] - PositionVectors[0]; // Vector 1 to 2
            c_vector<double, 3> Element_2_vector_13 = PositionVectors[2] - PositionVectors[0]; // Vector 1 to 3

            c_vector<double, 3> normal_2 = VectorProduct(Element_2_vector_12, Element_2_vector_13);

            double DirectionOf_2 = norm_2(normal_2);
            double DirectionOf_1 = norm_2(normal_1);

            c_vector<double, 2> NodeLocation = Create_c_vector(PositionVectors[2][0], PositionVectors[2][2]);
            c_vector<double, 2> Projection = Create_c_vector(normal_1[0], normal_1[1]);

            normal_2 = normal_2 / DirectionOf_2;
            normal_1 = normal_1 / DirectionOf_1;

            // if (DirectionOf_2 < 0)
            // {
            //     normal_2 = -normal_2;
            // }
            // if (DirectionOf_1 < 0)
            // {
            //     normal_1 = -normal_1;
            // }
            double angle = acos(inner_prod(normal_1, normal_2));
            double OriginalAngle = GetOriginalAngle(edge);
            double DifferenceInAngel = OriginalAngle - angle;

            // PRINT_3_VARIABLES(OriginalAngle, angle, DifferenceInAngel );
            SetOriginalAngle(edge, angle);
            // mMembraneStiffnessMap[edge] = mMembraneStiffness /10000;
        }

        //           4
        //        -     -
        //      -         -
        //    1 ----------- 3
        //      -         -
        //        -     -
        //           2

        if (ElementPair == 4)
        {
            //           0
            //        -     -
            //      -         -
            //    0 ----------- 0
            //      -         -
            //        -     -
            //           0
            mMembraneStiffnessMap[edge] = mMembraneStiffness / 1000;
        }
        else if (ElementPair == 1)
        {
            if (IsMutantMap[4] | IsMutantMap[2])
            {
                //           *
                //        -     -
                //      -         -
                //    * ----------- *
                //      -         -
                //        -     -
                //           0
                mMembraneStiffnessMap[edge] = mMembraneStiffness / 1; // Trial and error found this was best
            }
            // Dont need to consider the other options, because it should effect the bending stiffness between the two elements
        }
        else if (ElementPair == 2)
        {
            if (IsMutantMap[1] & IsMutantMap[3])
            { // Dont want this happening, it will mean cells are breaking apart
                //           *
                //        -     -
                //      -         -
                //    0 ----------- 0
                //      -         -
                //        -     -
                //           *
                mMembraneStiffnessMap[edge] = mMembraneStiffness / 10;
            }
            else
            {
                //           *
                //        -     -
                //      -         -
                //    0 ----------- *
                //      -         -
                //        -     -
                //           0
                mMembraneStiffnessMap[edge] = mMembraneStiffness / 1; // Trial and error found this was best
            }
        }
        else if (ElementPair == 3)
        {
            if (IsMutantMap[1] == 0 | IsMutantMap[3] == 0) // Meaning there is a cell at node 1 or 3
            {
                //          0                      0
                //        -     -                -     -
                //      -         -            -         -
                //    0 ----------- *  or    * ----------- 0
                //      -         -            -         -
                //        -     -                -     -
                //           0                      0
                mMembraneStiffnessMap[edge] = mMembraneStiffness / 100;
            }
            if (IsMutantMap[2] == 0 | IsMutantMap[4] == 0) // Meaning there is a cell at node 1 or 3
            {
                //          *                      0
                //        -     -                -     -
                //      -         -            -         -
                //    0 ----------- 0  or    0 ----------- 0
                //      -         -            -         -
                //        -     -                -     -
                //           0                      *
                mMembraneStiffnessMap[edge] = mMembraneStiffness / 100;
            }
        }

        // // PRINT_2_VARIABLES(IsMutantMap[1],IsMutantMap[3]);

        // else if (IsMutantMap[1] & IsMutantMap[3])
        // {
        //      //          M                      M
        //     //        -     -                -     -
        //     //      -         -            -         -
        //     //    * ----------- *  or    M ----------- M
        //     //      -         -            -         -
        //     //        -     -                -     -
        //     //           M                      M
        //     mMembraneStiffnessMap[edge] = mMembraneStiffness /100;
        // }
        // else if (IsMutantMap[1] & IsMutantMap[3])
        // {
        //      //          M                      M
        //     //        -     -                -     -
        //     //      -         -            -         -
        //     //    * ----------- *  or    M ----------- M
        //     //      -         -            -         -
        //     //        -     -                -     -
        //     //           M                      M
        //     mMembraneStiffnessMap[edge] = mMembraneStiffness /100;
        // }
        // else if (IsMutantMap[1] | IsMutantMap[3])
        // {
        //     if (IsMutantMap[2] & IsMutantMap[4])
        //     {
        //     //           *
        //     //        -     -
        //     //      -         -
        //     //    M ----------- M
        //     //      -         -
        //     //        -     -
        //     //           M
        //             mMembraneStiffnessMap[edge] = mMembraneStiffness/1;
        //     } else if (IsMutantMap[2] | IsMutantMap[4])
        //     {
        //     //           *
        //     //        -     -
        //     //      -         -
        //     //    M ----------- *
        //     //      -         -
        //     //        -     -
        //     //           M
        //          mMembraneStiffnessMap[edge] = mMembraneStiffness/1;
        //     } else
        //     {
        //     //           *
        //     //        -     -
        //     //      -         -
        //     //    * ----------- *
        //     //      -         -
        //     //        -     -
        //     //           M
        //         //   mMembraneStiffnessMap[edge] = mMembraneStiffness /2;
        //     }
        // }

        // else if (IsMutantMap[1] & IsMutantMap[3])
        // {
        //     mMembraneStiffnessMap[edge] = mMembraneStiffness /100;
        // }

        // if (ElementPair == 1 )
        // {
        //     mMembraneStiffnessMap[edge] = mMembraneStiffness/6;
        // } else if (ElementPair == 2 )
        // {
        //     mMembraneStiffnessMap[edge] = mMembraneStiffness/10;
        // } else if (ElementPair > 2 )
        // {
        //     mMembraneStiffnessMap[edge] = mMembraneStiffness/100;
        // }
        // PRINT_VARIABLE(mMembraneStiffnessMap[edge]);
    }
}

// void MembraneStiffnessForce::UpdateMembraneProperties(AbstractCellPopulation<2, 3>& rCellPopulation)
// {
//     MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);

//     for (typename AbstractTetrahedralMesh<2, 3>::ElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
//          elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
//          ++elem_iter)
//     {
//         unsigned elem_index = elem_iter->GetIndex();

//         // Check the mutant state of each of the cells in the element, then set this number to decide the
//         // inital area and the membrane constants of the element
//         unsigned node_index;
//         CellPtr p_cell ;
//         double NumberOfMutantCells = 0;
//         // double EleElasticShearModulus =  mElasticShearModulus;
//         // double EleAreaDilationModulus = mAreaDilationModulus ;
//          // Iterate over the cells in this element to check if they are mutant
//          // want to know how many mutant cells we have
//          for (int i = 0; i < 3; i++)
//         {
//             node_index = elem_iter->GetNodeGlobalIndex(i);
//             p_cell = p_cell_population->GetCellUsingLocationIndex(node_index);
//             if (p_cell->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>()) // If on edge
//             {
//                 NumberOfMutantCells+=1;
//             }
//         }

//         // Membrane constants ....
//         // if 0 mutants then ElasticShearModulus = mElasticShearModulus
//         if (NumberOfMutantCells == 0)
//         {
//             mMembraneStiffnessMap[elem_index]  = mMembraneStiffness ;
//         }
//         else
//         {
//             // mMembraneStiffnessMap[elem_index]  = 100 * mMembraneStiffness ;
//         }

//     }

// }

void MembraneStiffnessForce::OutputForceParameters(out_stream& rParamsFile)
{
    //    *rParamsFile << "\t\t\t<UseCutOffLength>" << mUseCutOffLength << "</UseCutOffLength>\n";
    //    *rParamsFile << "\t\t\t<CutOffLength>" << mMechanicsCutOffLength << "</CutOffLength>\n";

    // Call method on direct parent class
    AbstractForce<2, 3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MembraneStiffnessForce)
