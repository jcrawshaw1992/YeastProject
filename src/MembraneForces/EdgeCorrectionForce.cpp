#include "EdgeCorrectionForce.hpp"

#include "LostEndothelialCell.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"


// #include "MathsFunctions.hpp"

EdgeCorrectionForce::EdgeCorrectionForce()
        : AbstractForce<2, 3>()
{
}



void EdgeCorrectionForce::AddForceContribution(AbstractCellPopulation<2, 3>& rCellPopulation)
{
    // TRACE("Add Edge correction Force");
    double counter = 0;
    MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);

   for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
     cell_iter != rCellPopulation.End();
     ++cell_iter)
    {
    
        unsigned ReferenceNode = 0;

        if ((cell_iter)->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>())// | (is_boundary_node &&  node_index > (mNc *mNz)- mNc))
        {      
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        Node<3>* pNode = p_cell_population->rGetMesh().GetNode(node_index);

            // PRINT_2_VARIABLES(counter, node_index);
            if (node_index < mNc + 1) // if on lower edge
            {
                ReferenceNode = node_index + (2 * mNc); // select node from two rows up
            }
            else if (node_index > mNc) // if on upper edge
            {
                ReferenceNode = node_index - (2 * mNc); // select node from two rows down
            }
            Node<3>* pReferenceNode = p_cell_population->rGetMesh().GetNode(ReferenceNode);

            // TRACE("clear the force");
            pNode->ClearAppliedForce(); // remove the already present force at this node
            pNode->AddAppliedForceContribution(pReferenceNode->rGetAppliedForce()); // Add the new force

    }

    }

}
      
      
void EdgeCorrectionForce::SetMeshType(bool RegularCylinder, unsigned Nc, unsigned Nz)
{
    mNc = Nc;
    mNz = Nz;
    mRegularCylinder = RegularCylinder;

}

void EdgeCorrectionForce::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<2, 3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(EdgeCorrectionForce)

// 
//     }
//   }

//   if ( mCorrectDrag == true) // Need to write a code the hangles not regular meshes
//   {
// for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
//      cell_iter != rCellPopulation.End();
//      ++cell_iter)
// {
//     unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
//     Node<3>* pNode = p_cell_population->rGetMesh().GetNode(node_index);
//     double voronoi_cell_area = rCellPopulation.GetVolumeOfCell(*cell_iter);
//     ForceMap[node_index] = ForceMap[node_index];/// voronoi_cell_area;
//     // PRINT_VARIABLE(voronoi_cell_area);

//     pNode->AddAppliedForceContribution(ForceMap[node_index]);

//     if (node_index ==100)
//     {
//         PRINT_VECTOR(ForceMap[node_index]);
//     }

// }

// void EdgeCorrectionForce::RegularBoundCylinder(bool Cylinder, double Nc , double Nz)
// {

//     // this will be used in an if statment in the AddForceContribution to prevent boundary effects on a basic cylinder
//     // (boundary effects arrising due to the different number of associated elements on nodes at the edge. This particual option is for a strctured mesh)
//     mCylinder = Cylinder;
//     mNc = Nc;
//     mNz = Nz;
// }

// void EdgeCorrectionForce::DragCorrection(bool CorrectDrag)
// {
//      // this will be used in an if statment in the AddForceContribution to prevent boundary effects on a basic cylinder
//     // (boundary effects arrising due to the different number of associated elements on nodes at the edge. This particual option is for a unstructured mesh)

//     mCorrectDrag = CorrectDrag;
// }

//    for (int i = 0; i < 3; i++)
//         {
//             EleElasticShearModulus   EleAreaDilationModulus
//             dedvX = mElasticShearModulus * (I1 + 1) * (Dxx * a_i[i] + Dxy * b_i[i]) + (-mElasticShearModulus + mAreaDilationModulus * I2) * ((Dxx * Dyy * Dyy - Dxy * Dyx * Dyy) * a_i[i] + ( Dxy*Dyx*Dyx - Dxx * Dyx * Dyy) * b_i[i]);
//             dedvY = mElasticShearModulus * (I1 + 1) * (Dyx * a_i[i] + Dyy * b_i[i]) + (-mElasticShearModulus + mAreaDilationModulus * I2) * ((Dxx * Dxx * Dyy - Dxx * Dxy * Dyx) * b_i[i] + (Dyx*Dxy*Dxy - Dxx * Dxy * Dyy ) * a_i[i]);
//             RotatedForceOnNode[i] = -mArea0[elem_index] / 3 * Create_c_vector(dedvX, dedvY, 0);
//             RotatedMag[i] = norm_2(RotatedForceOnNode[i]);
//             // This is the force for each node  for the triangle situated at the origin
//             // Matrix with each row containing the force on the corresponding element
//         }

// CellPtr p_cell_0 = p_cell_population->GetCellUsingLocationIndex(elem_iter->GetNodeGlobalIndex(0));
// CellPtr p_cell_1 = p_cell_population->GetCellUsingLocationIndex(elem_iter->GetNodeGlobalIndex(1));
// CellPtr p_cell_2 = p_cell_population->GetCellUsingLocationIndex(elem_iter->GetNodeGlobalIndex(2));

// // This could definitly be written in a loop to make more efficent

// c_vector<long double, 3> Vector0;
// Vector0[0] = p_cell_0->GetCellData()->GetItem("Initial_Location_X");
// Vector0[1] = p_cell_0->GetCellData()->GetItem("Initial_Location_Y");
// Vector0[2] = p_cell_0->GetCellData()->GetItem("Initial_Location_Z");

// double Radius_0 = norm_2(Create_c_vector(Vector0[0],Vector0[1]));
// if (Vector0[0]< 1e-19)
// {
//     theta_0 = M_PI/2;
// } else
// {
//     theta_0 = atan(Vector0[1]/Vector0[0]);
// }

// c_vector<long double, 3> Vector1;
// Vector1[0] = p_cell_1->GetCellData()->GetItem("Initial_Location_X");
// Vector1[1] = p_cell_1->GetCellData()->GetItem("Initial_Location_Y");
// Vector1[2] = p_cell_1->GetCellData()->GetItem("Initial_Location_Z");

// double Radius_1 = norm_2(Create_c_vector(Vector1[0],Vector1[1]));
// if (Vector1[0]< 1e-19)
// {
//    theta_1 = M_PI/2;
// } else
// {
//    theta_1 = atan(Vector1[1]/Vector1[0]);
// }

//  c_vector<long double, 3> Vector2;
// Vector2[0] = p_cell_2->GetCellData()->GetItem("Initial_Location_X");
// Vector2[1] = p_cell_2->GetCellData()->GetItem("Initial_Location_Y");
// Vector2[2] = p_cell_2->GetCellData()->GetItem("Initial_Location_Z");

// double Radius_2 = norm_2(Create_c_vector(Vector2[0],Vector2[1]));
// if (Vector2[0]< 1e-19)
// {
//    theta_2 = M_PI/2;
// } else
// {
//   theta_2 = atan(Vector2[1]/Vector2[0]);
// }
// // PRINT_3_VARIABLES(theta_0,theta_1, theta_2);

// // The radius of the new nodes
// double Adapted_Radius_0 = Radius_0/ScallingFactor;
// double Adapted_Radius_1 = Radius_1/ScallingFactor;
// double Adapted_Radius_2 = Radius_2/ScallingFactor;

// c_vector<long double, 3> Adapted_Node_0 = Create_c_vector(Adapted_Radius_0* cos(theta_0),Adapted_Radius_0* sin(theta_0),Vector0[2] );
// c_vector<long double, 3> Adapted_Node_1 = Create_c_vector(Adapted_Radius_1* cos(theta_1),Adapted_Radius_1* sin(theta_1),Vector1[2] );
// c_vector<long double, 3> Adapted_Node_2 = Create_c_vector(Adapted_Radius_2* cos(theta_2),Adapted_Radius_2* sin(theta_2),Vector2[2] );

// PRINT_VARIABLE(NodeLocation[1]);

// if (abs(NodeLocation[0]) < 1e-20) // along the y axis of the unit circle
// {
//     PRINT_VARIABLE(theta_0 );
//     if (NodeLocation[1] >0)
//     {
//         theta_0 = M_PI/2;
//     }else if  (NodeLocation[1] < 0)
//     {
//         theta_0 = 3 *M_PI/2;
//     }
//     // TRACE("NodeLocation[0]) ==0");
//     // theta_0 = 10;//; M_PI/2;//std::copysignf(1.0, NodeLocation[1])* M_PI/2; //std::copysignf(1.0, NodeLocation[1]) returns the sign of NodeLocation
// } else if (abs(NodeLocation[1]) <1e-19) //  along the x axis of the unit circle
// {
//     TRACE("NodeLocation[1] ==0");
//     theta_0 = signbit(NodeLocation[0])* M_PI; // signbit returns 0 for positive numbers and 1 for negative numbers
// }

//    PRINT_VARIABLE(theta_0);



 // if ((cell_iter)->GetMutationState()->IsType<LostEndothelialCell>()) // If on edge
        // {
        //     unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter); //p_cell_population->GetLocationIndexUsingCell(*cell_iter);
        //     Node<3>* p_node = rCellPopulation.GetNode(node_index);
        //     // TRACE("Have a mutant");
        //     c_vector<double, 3> NodeLocation;
        //     NodeLocation[0] = cell_iter->GetCellData()->GetItem("Initial_Location_X");
        //     NodeLocation[1] = cell_iter->GetCellData()->GetItem("Initial_Location_Y");
        //     NodeLocation[2] = cell_iter->GetCellData()->GetItem("Initial_Location_Z");

        //     double Original_Radius = norm_2(Create_c_vector(NodeLocation[0], NodeLocation[1]));
        //     double New_Radius = Original_Radius / ScallingFactor2;
        //     theta_0 = atan(NodeLocation[1] / NodeLocation[0]);

        //     double x = NodeLocation[0];
        //     double y = NodeLocation[1];

        //     if (x < 0 && y > 0)
        //     {
        //         theta_0 = theta_0 + M_PI;
        //     }
        //     else if (x < 0 && y < 0)
        //     {
        //         theta_0 = theta_0 + M_PI;
        //     }

        //     c_vector<double, 3> PositionVector = Create_c_vector(New_Radius * cos(theta_0), New_Radius * sin(theta_0), NodeLocation[2]);
        //     mAdjustedInitalNodeLocation[node_index] = PositionVector;
        //     (cell_iter)->GetCellData()->SetItem("Initial_Location_X", PositionVector[0]);
        //     (cell_iter)->GetCellData()->SetItem("Initial_Location_Y", PositionVector[1]);
        //     (cell_iter)->GetCellData()->SetItem("Initial_Location_Z", PositionVector[2]);
        //     mHaveAdjustedPositionMap[node_index] = 1;

        //     // This piece records all the elements associated with a mutant cell
        //     std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();
        //     assert(containing_elements.size() > 0);
        //     // LabledElements.push_back(containing_elements);
        //     for (std::set<unsigned>::iterator iter = containing_elements.begin();
        //          iter != containing_elements.end();
        //          ++iter)
        //     {
        //         LabledElements.push_back(*iter);
        //     }
        //  }  else 



    // for (typename AbstractTetrahedralMesh<2, 3>::ElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
    //      elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
    //      ++elem_iter)
    // {
    //     // TRACE("In modifering loop");

    //     unsigned elem_index = elem_iter->GetIndex();
    //     unsigned node_index;
    //     CellPtr p_cell;
    //     double NumberOfMutantCells = 0;

    //     // Iterate over the 3 cells in this element to check if they are mutant
    //     // want to know how many mutant cells we have
    //     c_vector<c_vector<double, 3>, 3> PositionVectors;

    //     for (int i = 0; i < 3; i++)
    //     {
    //         // TRACE("In the for loop");
    //         node_index = elem_iter->GetNodeGlobalIndex(i);

    //         p_cell = p_cell_population->GetCellUsingLocationIndex(node_index);

    //         if (p_cell->GetMutationState()->IsType<LostEndothelialCell>()) // If on edge
    //         {
    //             NumberOfMutantCells += 1;

    //             node_index = elem_iter->GetNodeGlobalIndex(i);
    //             // PositionVectors[i] =  mAdjustedInitalNodeLocation[node_index];
    //             // PRINT_VECTOR(mAdjustedInitalNodeLocation[node_index+1] );

    //         }
    //     }

    //     // Membrane constants ....
    //     // if 0 mutants then ElasticShearModulus = mElasticShearModulus
    //     if (NumberOfMutantCells ==3 )
    //     {

    //         for (int i = 0; i < 3; i++)
    //         {
    //             // TRACE("In the for loop");
    //             node_index = elem_iter->GetNodeGlobalIndex(i);
    //             p_cell = p_cell_population->GetCellUsingLocationIndex(node_index);

    //             PositionVectors[i][0] = p_cell->GetCellData()->GetItem("Initial_Location_X");
    //             PositionVectors[i][1] = p_cell->GetCellData()->GetItem("Initial_Location_Y");
    //             PositionVectors[i][2] = p_cell->GetCellData()->GetItem("Initial_Location_Z");
    //         }

    //         c_vector<double, 3> Adapted_vector_12 = PositionVectors[1] - PositionVectors[0]; // Vector 1 to 2
    //         c_vector<double, 3> Adapted_vector_13 = PositionVectors[2] - PositionVectors[0] ; // Vector 1 to 3

    //         // Find the intial area, A0 = 0.5*norm(normal)
    //         c_vector<long double, 3> normalVector = VectorProduct(Adapted_vector_12, Adapted_vector_13);
    //         // PRINT_VECTOR(normalVector);
    //         long double Area = 0.5 * norm_2(normalVector);

    //         long double a = norm_2(Adapted_vector_12); // Lenght a -> edge connecting P1 and P2
    //         long double b = norm_2(Adapted_vector_13); // Lenght b -> edge connecting P1 and P3
    //         long double alpha = acos(inner_prod(Adapted_vector_12, Adapted_vector_13) / (a * b));
    //         double Y = b * sin(alpha);

    //         mArea0[elem_index] = Area;

    //         c_vector<long double, 2> x1 = Create_c_vector(0, 0);
    //         c_vector<long double, 2> x2 = Create_c_vector(a, 0);
    //         c_vector<long double, 2> x3 = Create_c_vector(b * cos(alpha), b * sin(alpha));

    //         //Save the intial position vectors
    //         mInitalVectors[elem_index][0] = x1;
    //         mInitalVectors[elem_index][1] = x2;
    //         mInitalVectors[elem_index][2] = x3;

    //         // FInd the 6 shape function constants
    //         aVector[0] = (x2[1] - x3[1]) / (2 * Area);
    //         aVector[1] = (x3[1] - x1[1]) / (2 * Area);
    //         aVector[2] = (x1[1] - x2[1]) / (2 * Area);

    //         bVector[0] = (x3[0] - x2[0]) / (2 * Area);
    //         bVector[1] = (x1[0] - x3[0]) / (2 * Area);
    //         bVector[2] = (x2[0] - x1[0]) / (2 * Area);

    //         // PRINT_VECTOR(mACoefficients[elem_index]);
    //         // PRINT_VECTOR(mBCoefficients[elem_index]);

    //         mACoefficients[elem_index] = aVector;
    //         mBCoefficients[elem_index] = bVector;
    //     }

    //  }