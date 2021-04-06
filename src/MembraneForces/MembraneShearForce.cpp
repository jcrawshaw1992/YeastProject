#include "MembraneShearForce.hpp"


// #include "MathsFunctions.hpp"

MembraneShearForce::MembraneShearForce()
        : AbstractForce<2, 3>()
{
}

void MembraneShearForce::SetElasticShearModulus(double ElasticShearModulus)
{
    mElasticShearModulus = ElasticShearModulus;
}

void MembraneShearForce::SetAreaDilationModulus(double AreaDilationModulus)
{
    mAreaDilationModulus = AreaDilationModulus;
}

void MembraneShearForce::AddForceContribution(AbstractCellPopulation<2, 3>& rCellPopulation)
{
    // TRACE("Add Shear Force");
    MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);

    for (typename AbstractTetrahedralMesh<2, 3>::ElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        
        unsigned elem_index = elem_iter->GetIndex();

        Node<3>* pNode0 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(0));
        Node<3>* pNode1 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(1));
        Node<3>* pNode2 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(2));

        c_vector<long double, 3> Vector1 = pNode1->rGetLocation();
        c_vector<long double, 3> Vector2 = pNode2->rGetLocation();
        c_vector<long double, 3> Vector0 = pNode0->rGetLocation();
        // PRINT_VECTOR(Vector0 );

        c_vector<long double, 3> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
        c_vector<long double, 3> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 33);

        unsigned node_index;
        CellPtr p_cell;

        // Get side lengths and angle to create an equal triangle at the origin
        long double a = norm_2(vector_12); // Lenght a -> edge connecting P0 and P1
        long double b = norm_2(vector_13); // Lenght b -> edge connecting P0 and P2
        long double alpha = acos(inner_prod(vector_12, vector_13) / (a * b));

        // This will create an equal triangle at the origin
        c_vector<long double, 2> x1 = Create_c_vector(0, 0);
        c_vector<long double, 2> x2 = Create_c_vector(a, 0);
        c_vector<long double, 2> x3 = Create_c_vector(b * cos(alpha), b * sin(alpha));

        // Get the original triangle
        c_vector<c_vector<long double, 2>, 3> InitialVectors = mInitalVectors[elem_index];

        //Displacement vectors.
        c_vector<long double, 2> V1 = x1 - InitialVectors[0];
        c_vector<long double, 2> V2 = x2 - InitialVectors[1];
        c_vector<long double, 2> V3 = x3 - InitialVectors[2];

        // Get the shape function coefficents for this element
        c_vector<long double, 3> a_i = mACoefficients[elem_index];
        c_vector<long double, 3> b_i = mBCoefficients[elem_index];

        // Deformation tensor
        long double Dxx = 1 + a_i[0] * V1[0] + a_i[1] * V2[0] + a_i[2] * V3[0];
        long double Dxy = b_i[0] * V1[0] + b_i[1] * V2[0] + b_i[2] * V3[0];
        long double Dyx = a_i[0] * V1[1] + a_i[1] * V2[1] + a_i[2] * V3[1];
        long double Dyy = 1 + b_i[0] * V1[1] + b_i[1] * V2[1] + b_i[2] * V3[1];

        c_vector<c_vector<long double, 2>, 2> G;

        // G =DTD  -- Caughy green 
        G[0][0] = Dxx * Dxx + Dyx * Dyx;
        G[0][1] = Dxx * Dxy + Dyx * Dyy;
        G[1][0] = Dxx * Dxy + Dyy * Dyx;
        G[1][1] = Dxy * Dxy + Dyy * Dyy;

        // Strain invarients 
        long double I1 = tr(G) - 2;
        long double I2 = det(G) - 1;


        c_vector<c_vector<long double, 3>, 3> RotatedForceOnNode;
        c_vector<long double, 3> RotatedMag;
        c_vector<long double, 3> OrientatedMag;
        long double dedvX;
        long double dedvY;

// double S = log10(mElasticShearModulusMap[elem_index] );
// double Alpha = log10(mAreaDilationModulusMap[elem_index] );
// PRINT_2_VARIABLES(S, Alpha);

        for (int i = 0; i < 3; i++)
        {
            dedvX = mElasticShearModulusMap[elem_index] * (I1 + 1) * (Dxx * a_i[i] + Dxy * b_i[i]) + (-mElasticShearModulusMap[elem_index] + mAreaDilationModulusMap[elem_index] * I2) * ((Dxx * Dyy * Dyy - Dxy * Dyx * Dyy) * a_i[i] + (Dxy * Dyx * Dyx - Dxx * Dyx * Dyy) * b_i[i]);
            dedvY = mElasticShearModulusMap[elem_index] * (I1 + 1) * (Dyx * a_i[i] + Dyy * b_i[i]) + (-mElasticShearModulusMap[elem_index] + mAreaDilationModulusMap[elem_index] * I2) * ((Dxx * Dxx * Dyy - Dxx * Dxy * Dyx) * b_i[i] + (Dyx * Dxy * Dxy - Dxx * Dxy * Dyy) * a_i[i]);
            RotatedForceOnNode[i] = -mArea0[elem_index] / 3 * Create_c_vector(dedvX, dedvY, 0);
            RotatedMag[i] = norm_2(RotatedForceOnNode[i]);
            // This is the force for each node  for the triangle situated at the origin
            // Matrix with each row containing the force on the corresponding element
        }

        // Need the transpose of the RotatedForceOnNode for the correct order for matrix multipication
        // c_vector<c_vector<long double, 3>, 3> TRotatedForceOnNode = MatrixTranspose(RotatedForceOnNode);

        // Roate the RotatedForceOnNode matrix back to the orientation of the orignial triangle. The MappingMatrix function is defined below at line 573
        // and produces the inverse of the roation matrix needed to put the original triangle at the origin for the given element
        // c_vector<c_vector<long double, 3>, 3> ForceOnNode = MatrixTranspose(MatrixMultiplication(MappingMatrix(p_cell_population, elem_iter, a, b, alpha), TRotatedForceOnNode)); // Transpose the force vector to get the right order for the matrix mult

        // c_vector<c_vector<long double, 3>, 3> ForceOnNode;
        c_vector<c_vector<long double, 3>, 3> X;
        c_vector<c_vector<long double, 3>, 3> F_rp;
        c_vector<c_vector<long double, 3>, 3> ForceOnNode;

        X[0] = Create_c_vector(a, 0, 0);
        X[1] = Create_c_vector(b * cos(alpha), b * sin(alpha), 0);
        X[2] = Create_c_vector(0, 0, 1);
        X = MatrixTranspose(X);

        // Rotate the force to the original configuretion

        c_vector<long double, 3> F0_rp = MatrixMultiplication(Inverse(X), RotatedForceOnNode[0]);
        c_vector<long double, 3> F1_rp = MatrixMultiplication(Inverse(X), RotatedForceOnNode[1]);
        c_vector<long double, 3> F2_rp = MatrixMultiplication(Inverse(X), RotatedForceOnNode[2]);

        c_vector<c_vector<long double, 3>, 3> Ident = MatrixMultiplication(Inverse(X), X);
        ForceOnNode[0] = F0_rp[0] * vector_12 + F0_rp[1] * vector_13 + F0_rp[2] * X[2];
        ForceOnNode[1] = F1_rp[0] * vector_12 + F1_rp[1] * vector_13 + F1_rp[2] * X[2];
        ForceOnNode[2] = F2_rp[0] * vector_12 + F2_rp[1] * vector_13 + F2_rp[2] * X[2];

        long double CellArea;

        // Correct the force for the cell area so drag is area dependent 
        for (int i = 0; i < 3; i++)
        {
            node_index = elem_iter->GetNodeGlobalIndex(i);
            CellArea = rCellPopulation.GetVolumeOfCell(rCellPopulation.GetCellUsingLocationIndex(node_index));
            ForceOnNode[i] /= CellArea;
            // PRINT_VARIABLE(CellArea);
        }


        //Add the forces to the node
        pNode0->AddAppliedForceContribution(ForceOnNode[0]);
        pNode1->AddAppliedForceContribution(ForceOnNode[1]);
        pNode2->AddAppliedForceContribution(ForceOnNode[2]);
    }

// TRACE("DONE");
    // for (typename AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
    //      cell_iter != rCellPopulation.End();
    //      ++cell_iter)
    // {
    //     unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
    //     Node<3>* pNode = p_cell_population->rGetMesh().GetNode(node_index);
    //     c_vector<double, 3> Force = pNode->rGetAppliedForce();

    //     double RadialForce = norm_2(Create_c_vector(Force[0], Force[1]));
    //     // double Beta = atan(Force[1] / Force[0]);
    //     // double Z = Force[2];

    //     cell_iter->GetCellData()->SetItem("RadialStrainForce", RadialForce);
    //     cell_iter->GetCellData()->SetItem("StrainForceX", Force[0]);
    //     cell_iter->GetCellData()->SetItem("StrainForceY", Force[1]);
    //     cell_iter->GetCellData()->SetItem("StrainForceZ", Force[2]);
    // }
}

/*
 *  Update the membrane properties based on mutation states
*/

// This piece of code needs to be updated such that it only changes the mechanical properties of the vessel, but
// the original areas are already set to be 1/10 of what they originally where. 


void MembraneShearForce::UpdateMembraneProperties(AbstractCellPopulation<2, 3>& rCellPopulation)
{
    //     TRACE("Modifiing membrane properties");

    //     TRACE("Currently assume only if all three nodes of the element are mutant, then the element will be altered, however I might want to think about this down the track. It might be more physically relevant to set conditions for if 1 or 2 or 3 nodes are mutant");

    //     c_vector<c_vector<double, 2>, 3> InitalVectors;
    //     c_vector<long double, 3> aVector;
    //     c_vector<long double, 3> bVector;
    //     double theta_0;
    //     double ScallingFactor = 15;
    //     std::vector<unsigned> LabledElements;
    //     double NumberOfNeighbouringECs;

    //     PRINT_VARIABLE(LabledElements.size());


    // //  mNeighbours[node_index] 
        
    //     // LabledElements.push_back(1);
    //     MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);
    //     // Loop over all the mutant elements and
    //     for (typename AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
    //          cell_iter != rCellPopulation.End();
    //          ++cell_iter)
    //     {
    //         // if mutant
    //          if ((cell_iter)->GetMutationState()->IsType<LostEndothelialCell>()) // If on edge
    //         {
    //             unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter); //p_cell_population->GetLocationIndexUsingCell(*cell_iter);
    //             Node<3>* p_node = rCellPopulation.GetNode(node_index);
    //             PRINT_VARIABLE( node_index );
    //             // Need to check how many of the neighbours have ECs -> the more EC surrounding the empty node, the slower it is going to pull
    //             // Into the lumen. Conversly an emtpy node with no surrounding ECs is going to fall quicker. 

    //             std::vector<unsigned> NeighboursVector = mNeighbours[node_index] ;
    //             PRINT_VECTOR(NeighboursVector)
    //             NumberOfNeighbouringECs = 0;

    //             // PRINT_VECTOR(NeighboursVector);
    //             // itterate over the nighbours and count the number with cells
    //             for (int i=0; i < NeighboursVector.size(); i++)
    //             {
    //                 // PRINT_VARIABLE(i);
    //                 // get the cell iterator from the node index
    //                 CellPtr p_cell = p_cell_population->GetCellUsingLocationIndex(NeighboursVector[i]);
    //                 if (p_cell->GetMutationState()->IsType<HasEndothelialCell>())
    //                     {
    //                         NumberOfNeighbouringECs +=1;
    //                         PRINT_VARIABLE( node_index );
    //                     } 
    //             }
    //             PRINT_VARIABLE(NumberOfNeighbouringECs);
            
    //             // // This will set if has two neighbours with nodes scalling factor is halved
    //             double Scalling = ScallingFactor ;
    //             PRINT_VARIABLE(Scalling );
    //             // double Scalling = 1+ (ScallingFactor -1)/(exp(NumberOfNeighbouringECs));
    //             // PRINT_VARIABLE(Scalling);

    //             (cell_iter)->GetCellData()->SetItem("Scalling", Scalling);
        

    //             c_vector<double, 3> NodeLocation;
    //             NodeLocation[0] = cell_iter->GetCellData()->GetItem("Initial_Location_X");
    //             NodeLocation[1] = cell_iter->GetCellData()->GetItem("Initial_Location_Y");
    //             NodeLocation[2] = cell_iter->GetCellData()->GetItem("Initial_Location_Z");
    //             PRINT_VECTOR(NodeLocation);

    //             double Original_Radius = norm_2(Create_c_vector(NodeLocation[0], NodeLocation[1]));
    //             double New_Radius = (Original_Radius) / Scalling;
    //             theta_0 = atan(NodeLocation[1] / NodeLocation[0]);

    //             double x = NodeLocation[0];
    //             double y = NodeLocation[1];

    //             if (x < 0 && y > 0)
    //             {
    //                 theta_0 = theta_0 + M_PI;
    //             }
    //             else if (x < 0 && y < 0)
    //             {
    //                 theta_0 = theta_0 + M_PI;
    //             }
    //             PRINT_3_VARIABLES(New_Radius, Original_Radius, Scalling);
                

    //             c_vector<double, 3> PositionVector = Create_c_vector(New_Radius * cos(theta_0), New_Radius * sin(theta_0), NodeLocation[2]);
    //             // mAdjustedInitalNodeLocation[node_index] = PositionVector;
    //             (cell_iter)->GetCellData()->SetItem("Initial_Location_X", PositionVector[0]);
    //             (cell_iter)->GetCellData()->SetItem("Initial_Location_Y", PositionVector[1]);
    //             (cell_iter)->GetCellData()->SetItem("Initial_Location_Z", PositionVector[2]);
    //             mHaveAdjustedPositionMap[node_index] = 1;

    //             // This piece records all the elements associated with a mutant cell
    //             std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();
    //             assert(containing_elements.size() > 0);
    //             // LabledElements.push_back(containing_elements);
    //             for (std::set<unsigned>::iterator iter = containing_elements.begin();
    //                  iter != containing_elements.end();
    //                  ++iter)
    //             {
    //                 LabledElements.push_back(*iter);
    //             }
    //         }
    //     }
    //     sort(LabledElements.begin(), LabledElements.end());
    //     LabledElements.erase(unique(LabledElements.begin(), LabledElements.end()), LabledElements.end());
    //     PRINT_VARIABLE(LabledElements.size());
    //     if (LabledElements.size()>0)
    //     {
    //         for (int i = 0; i < LabledElements.size(); i++) //
    //         {
    //         unsigned elem_index = LabledElements[i];
    //             c_vector<c_vector<double, 3>, 3> PositionVectors;
    //             Node<3>* pNode;
    //             unsigned node_index;
    //             CellPtr p_cell;
    //             double NumberOfMutants=0;

    //             for (int j = 0; j < 3; j++) //
    //             {
    //                 pNode = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(LabledElements[i])->GetNodeGlobalIndex(j));
    //                 node_index = pNode->GetIndex();
    //                 p_cell = p_cell_population->GetCellUsingLocationIndex(node_index);
    //                 PositionVectors[j][0] = p_cell->GetCellData()->GetItem("Initial_Location_X");
    //                 PositionVectors[j][1] = p_cell->GetCellData()->GetItem("Initial_Location_Y");
    //                 PositionVectors[j][2] = p_cell->GetCellData()->GetItem("Initial_Location_Z");

    //             if (( p_cell)->GetMutationState()->IsType<LostEndothelialCell>()) // If on edge
    //             {
    //                 NumberOfMutants+=1;
    //             }
                
    //             }

    //             c_vector<double, 3> Adapted_vector_12 = PositionVectors[1] - PositionVectors[0]; // Vector 1 to 2
    //             c_vector<double, 3> Adapted_vector_13 = PositionVectors[2] - PositionVectors[0]; // Vector 1 to 3

    //             // Find the intial area, A0 = 0.5*norm(normal)
    //             c_vector<long double, 3> normalVector = VectorProduct(Adapted_vector_12, Adapted_vector_13);
    //             // PRINT_VECTOR(normalVector);
    //             long double Area = 0.5 * norm_2(normalVector);

    //             long double a = norm_2(Adapted_vector_12); // Lenght a -> edge connecting P1 and P2
    //             long double b = norm_2(Adapted_vector_13); // Lenght b -> edge connecting P1 and P3
    //             long double alpha = acos(inner_prod(Adapted_vector_12, Adapted_vector_13) / (a * b));
    //             double Y = b * sin(alpha);

    //             mArea0[elem_index] = Area;

    //             c_vector<long double, 2> x1 = Create_c_vector(0, 0);
    //             c_vector<long double, 2> x2 = Create_c_vector(a, 0);
    //             c_vector<long double, 2> x3 = Create_c_vector(b * cos(alpha), b * sin(alpha));

    //             //Save the intial position vectors
    //             mInitalVectors[elem_index][0] = x1;
    //             mInitalVectors[elem_index][1] = x2;
    //             mInitalVectors[elem_index][2] = x3;

    //             // FInd the 6 shape function constants
    //             aVector[0] = (x2[1] - x3[1]) / (2 * Area);
    //             aVector[1] = (x3[1] - x1[1]) / (2 * Area);
    //             aVector[2] = (x1[1] - x2[1]) / (2 * Area);

    //             bVector[0] = (x3[0] - x2[0]) / (2 * Area);
    //             bVector[1] = (x1[0] - x3[0]) / (2 * Area);
    //             bVector[2] = (x2[0] - x1[0]) / (2 * Area);


    //             mACoefficients[elem_index] = aVector;
    //             mBCoefficients[elem_index] = bVector;

    //             /// ADAPT THE MEMBRNE CONSTANTS BASED ON THE NUMBER OF MUTANTS ON THE ELEMENT 

                    
            
    //             if (NumberOfMutants==1 ) // 6 - NumberOfNeighbourinEcs = Basement matrix sites 
    //             {

    //                 mElasticShearModulusMap[elem_index] /=25;
    //                 mAreaDilationModulusMap[elem_index] /=50;

    //             } else if(NumberOfMutants== 2)
    //             {
    //                 mElasticShearModulusMap[elem_index] /=100;
    //                 mAreaDilationModulusMap[elem_index] /=250;

    //             } else if(NumberOfMutants== 3)
    //             {
    //                 mElasticShearModulusMap[elem_index] /=200;
    //                 mAreaDilationModulusMap[elem_index] /=300;
    //             }

    //         }
    //     }

}

// double MembraneShearForce::GetoriginalElementArea(unsigned elem_index)
    // {
    //     TRACE("Here I am ");
    //     double OriginalArea = mArea0[elem_index];
    //     return OriginalArea;
    // }

    // double MembraneShearForce::SetoriginalElementArea(unsigned elem_index, double ElementArea)
    // {
    //     double OriginalArea = mArea0[elem_index];
    //     return OriginalArea;

    // }

    /*
    * Set up the inital configuration -> find the inital position vectors, the inital area, and the shape function 
    * 
    * The inital poisitions of all the nodes have been updated to be closer to the radius by a scalling factor (S) 
    * This reflects state where a vessel is deflated with no fluid flow to expand it, leaving only the compressive forces 
    * from the tissue, naturally pushing it inwards. 
    */


void MembraneShearForce::SetScallingShear(double Scalling)
{
       mScalling=Scalling;
}


void MembraneShearForce::SetupMembraneConfiguration(AbstractCellPopulation<2, 3>& rCellPopulation)
{
    MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);
    std::vector<unsigned> LabledElements;
    double NumberOfNeighbouringECs;
    // PRINT_VARIABLE(LabledElements.size());

   // Loop over all nodes and come up with a new inital position 

   
    for (typename AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter); //p_cell_population->GetLocationIndexUsingCell(*cell_iter);
        Node<3>* p_node = rCellPopulation.GetNode(node_index);
          
        c_vector<double, 3> NodeLocation =  p_node->rGetLocation();

        c_vector<long double, 3> normal = zero_vector<long double>(3);

        std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();
        assert(containing_elements.size() > 0);
        for (std::set<unsigned>::iterator iter = containing_elements.begin();
                iter != containing_elements.end();
                ++iter)
        {
            // Negative as normals point inwards for these surface meshes

            c_vector<long double, 3> LocationNormal = - p_cell_population->rGetMesh().GetElement(*iter)->CalculateNormal();
            // LocationNormal /=norm_2(LocationNormal );
            normal +=  LocationNormal;
        }
        normal /= norm_2(normal);
        // PRINT_VECTOR(normal);

      // Need to walk backward into the mesh by the scalling factor 
        // c_vector<double, 3> PositionVector = NodeLocation ;// + mScalling * normal;

        // (cell_iter)->GetCellData()->SetItem("Initial_Location_X", PositionVector[0]);
        // (cell_iter)->GetCellData()->SetItem("Initial_Location_Y", PositionVector[1]);
        // (cell_iter)->GetCellData()->SetItem("Initial_Location_Z", PositionVector[2]);

        (cell_iter)->GetCellData()->SetItem("Initial_Location_X",  NodeLocation[0]);
        (cell_iter)->GetCellData()->SetItem("Initial_Location_Y",  NodeLocation[1]);
        (cell_iter)->GetCellData()->SetItem("Initial_Location_Z",  NodeLocation[2]);
    }

    // Determine the inital element shape (shape function, internal angle, and area) for all elements
    for (typename AbstractTetrahedralMesh<2, 3>::ElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        // TRACE("Inital set up ");
        // define the necessary objects to be used in the loop 
        c_vector<c_vector<double, 3>, 3> PositionVector;
        unsigned elem_index = elem_iter->GetIndex();
        c_vector<double, 3> aVector;
        c_vector<double, 3> bVector;

        unsigned node_index;
        CellPtr p_cell;

        for (int i = 0; i < 3; i++) // Loop over the three cells and get their intial vessel locations
        {
            node_index = elem_iter->GetNodeGlobalIndex(i);
            p_cell = p_cell_population->GetCellUsingLocationIndex(node_index);
        
            // Get the inital locations for node i
            PositionVector[i][0] = p_cell->GetCellData()->GetItem("Initial_Location_X");
            PositionVector[i][1] = p_cell->GetCellData()->GetItem("Initial_Location_Y");
            PositionVector[i][2] = p_cell->GetCellData()->GetItem("Initial_Location_Z");
        }


        // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter); //p_cell_population->GetLocationIndexUsingCell(*cell_iter);
        // Node<3>* p_node = rCellPopulation.GetNode(node_index);
        // c_vector<double, 3> PositionVector =  p_node->rGetLocation();


        // Vectors connecting the nodes 
        c_vector<double, 3> vector_12 = PositionVector[1] - PositionVector[0]; // Vector 1 to 2
        c_vector<double, 3> vector_13 = PositionVector[2] - PositionVector[0]; // Vector 1 to 3

        // Find the intial area, A0 = 0.5*norm(normal)
        c_vector<double, 3> normalVector = VectorProduct(vector_12, vector_13);
        // PRINT_VECTOR(normalVector);
        double Area = 0.5 * norm_2(normalVector);

        // Save intial area in a map with the element as the key
        mArea0[elem_index] = Area;
        double a = norm_2(vector_12); // Lenght a -> edge connecting P1 and P2
        double b = norm_2(vector_13); // Lenght b -> edge connecting P1 and P3

        double alpha = acos(inner_prod(vector_12, vector_13) / (a * b));

        c_vector<double, 2> x1 = Create_c_vector(0, 0);
        c_vector<double, 2> x2 = Create_c_vector(a, 0);
        c_vector<double, 2> x3 = Create_c_vector(b * cos(alpha), b * sin(alpha));

        //Save the intial position vectors
        mInitalVectors[elem_index][0] = x1;
        mInitalVectors[elem_index][1] = x2;
        mInitalVectors[elem_index][2] = x3;

        // FInd the 6 shape function constants
        aVector[0] = (x2[1] - x3[1]) / (2 * Area);
        aVector[1] = (x3[1] - x1[1]) / (2 * Area);
        aVector[2] = (x1[1] - x2[1]) / (2 * Area);

        bVector[0] = (x3[0] - x2[0]) / (2 * Area);
        bVector[1] = (x1[0] - x3[0]) / (2 * Area);
        bVector[2] = (x2[0] - x1[0]) / (2 * Area);

        mACoefficients[elem_index] = aVector;
        mBCoefficients[elem_index] = bVector;

        // if (std::abs(PositionVector[0][2]) < 1.6e-3 || std::abs(PositionVector[1][2]) < 1.5e-3 || std::abs(PositionVector[2][2]) < 1.5e-3 )
        // {
        //     mElasticShearModulusMap[elem_index] = mElasticShearModulus*40;
        //     mAreaDilationModulusMap[elem_index] = mAreaDilationModulus*20;
        // }else 
        // {
        mElasticShearModulusMap[elem_index] = mElasticShearModulus;
        mAreaDilationModulusMap[elem_index] = mAreaDilationModulus;
        // }
        
        mMutantMap[elem_index] = 0;
    }


       // Want to save the neighbouring nodes
       // Loop over each cell and save all the nodes in a vector, then remove the double ups in the vector, and remove the cell of interest from its own neighbour list
     for (typename AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
            {      
                std::vector<unsigned> NeighboursToNode;
                unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter); //p_cell_population->GetLocationIndexUsingCell(*cell_iter);
                Node<3>* p_node = rCellPopulation.GetNode(node_index);
                // This piece records all the elements associated with a mutant cell
                std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();
                assert(containing_elements.size() > 0);
                // LabledElements.push_back(containing_elements);
                for (std::set<unsigned>::iterator iter = containing_elements.begin();
                    iter != containing_elements.end();
                    ++iter)
                {
                    // Loops over the neighbours to this cell and save them 
                    for (int i=0; i<3; i++)
                    {
                        Node<3>* pNode = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(i));
                        unsigned node_index_i  = pNode->GetIndex();
                        NeighboursToNode.push_back(node_index_i);

                    }
                    
                }
            
            //  
            sort(NeighboursToNode.begin(), NeighboursToNode.end());
            NeighboursToNode.erase(unique(NeighboursToNode.begin(), NeighboursToNode.end()), NeighboursToNode.end());

            // remove the cell of interest from its own neighbour list 
            for (int i= 0; i <  NeighboursToNode.size(); i++) //myvector.begin()+5
            {
                if (NeighboursToNode[i] == node_index)
                {
                    NeighboursToNode.erase(NeighboursToNode.begin()+i);
                    break;
                }
            }
            mNeighbours[node_index] =  NeighboursToNode;
        }
}




       // // Get the Nodes in this element
        // Node<3>* pNode0 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(0));
        // Node<3>* pNode1 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(1));
        // Node<3>* pNode2 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(2));










c_vector<c_vector<long double, 3>, 3> MembraneShearForce::MatrixMultiplication(c_vector<c_vector<long double, 3>, 3> Matrix1, c_vector<c_vector<long double, 3>, 3> Matrix2)
{

    c_vector<c_vector<long double, 3>, 3> MatrixTranspose;
    c_vector<c_vector<long double, 3>, 3> Answer;

    // This will give us the determinat
    for (int i = 0; i < 3; i++)
    {
        MatrixTranspose[i] = Create_c_vector(Matrix2[0][i], Matrix2[1][i], Matrix2[2][i]);
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Answer[i][j] = inner_prod(Matrix1[i], MatrixTranspose[j]);
        }
    }

    return Answer;
}

c_vector<long double, 3> MembraneShearForce::MatrixMultiplication(c_vector<c_vector<long double, 3>, 3> Matrix, c_vector<long double, 3> Vector)
{

    c_vector<long double, 3> Answer;

    for (int i = 0; i < 3; i++)
    {
        Answer[i] = inner_prod(Matrix[i], Vector);
    }

    return Answer;
}

c_vector<c_vector<long double, 3>, 3> MembraneShearForce::MatrixTranspose(c_vector<c_vector<long double, 3>, 3> Matrix)
{
    c_vector<c_vector<long double, 3>, 3> MatrixTranspose;
    for (int i = 0; i < 3; i++)
    { // Transporse of the matrix
        MatrixTranspose[i] = Create_c_vector(Matrix[0][i], Matrix[1][i], Matrix[2][i]);
    }
    return MatrixTranspose;
}

c_vector<c_vector<long double, 3>, 3> MembraneShearForce::Inverse(c_vector<c_vector<long double, 3>, 3> Matrix)
{

    c_vector<c_vector<long double, 3>, 3> MTranspose = MatrixTranspose(Matrix);

    long double det = Matrix[0][0] * (Matrix[1][1] * Matrix[2][2] - Matrix[1][2] * Matrix[2][1]) - Matrix[0][1] * (Matrix[1][0] * Matrix[2][2] - Matrix[1][2] * Matrix[2][0]) + Matrix[0][2] * (Matrix[1][0] * Matrix[2][1] - Matrix[1][1] * Matrix[2][0]);

    double l = 0, m = 1;
    c_vector<c_vector<long double, 3>, 3> InverseMatrix;
    for (int k = 0; k < 3; k++)
    {
        l += 1, m += 1;
        if (l == 3)
        {
            l = 0;
        }
        if (m == 3)
        {
            m = 0;
        }
        InverseMatrix[k] = (1 / det) * VectorProduct(MTranspose[l], MTranspose[m]);
    }
    return InverseMatrix;
}

c_vector<c_vector<long double, 3>, 3> MembraneShearForce::Inverse(c_vector<c_vector<long double, 3>, 3> Matrix, double elem)
{

    long double det = Matrix[0][0] * (Matrix[1][1] * Matrix[2][2] - Matrix[1][2] * Matrix[2][1]) - Matrix[0][1] * (Matrix[1][0] * Matrix[2][2] - Matrix[1][2] * Matrix[2][0]) + Matrix[0][2] * (Matrix[1][0] * Matrix[2][1] - Matrix[1][1] * Matrix[2][0]);

    c_vector<c_vector<long double, 3>, 3> MTranspose = MatrixTranspose(Matrix);

    double l = 0, m = 1;
    c_vector<c_vector<long double, 3>, 3> InverseMatrix;
    for (int k = 0; k < 3; k++)
    {
        l += 1, m += 1;
        if (l == 3)
        {
            l = 0;
        }
        if (m == 3)
        {
            m = 0;
        }
        InverseMatrix[k] = (1 / det) * VectorProduct(MTranspose[l], MTranspose[m]);
    }
    return InverseMatrix;
}

c_vector<c_vector<long double, 3>, 3> MembraneShearForce::RowReduction(c_vector<c_vector<long double, 3>, 3> Matrix)
{
    // row reduction to reduced row eshioln form
    c_vector<c_vector<long double, 3>, 3> Identity;
    Identity[0] = Create_c_vector(1, 0, 0);
    Identity[1] = Create_c_vector(0, 1, 0);
    Identity[2] = Create_c_vector(0, 0, 1);
    for (int i = 0; i < 3; i++)
    {

        Identity[i] /= Matrix[i][i];
        Matrix[i] /= Matrix[i][i]; // sets all the diagonal elements to 1

        for (int j = 0; j < 3; j++)
        {
            if (j != i) //prevents taking out the t row that we just make 1
            {
                Identity[j] -= (Identity[i] * Matrix[j][i]);
            }
            if (j != i) //prevents taking out the t row that we just make 1
            {
                Matrix[j] -= (Matrix[i] * Matrix[j][i]);
            }
        }
    }

    // by this point the oridinal identity matrix is now the inverse
    return Identity;
}

c_vector<c_vector<long double, 3>, 3> MembraneShearForce::MappingMatrix(MeshBasedCellPopulation<2, 3>* p_cell_population, typename AbstractTetrahedralMesh<2, 3>::ElementIterator elem_iter, double a, double b, double theta)
{
    // Inverse Mapping Matrix
    c_vector<c_vector<long double, 3>, 3> NewPositionVector;
    c_vector<c_vector<long double, 3>, 3> PositionVector;

    Node<3>* pNode0 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(0));
    Node<3>* pNode1 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(1));
    Node<3>* pNode2 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(2));

    int SpConst = 450; // prevents matrix begin singular
    PositionVector[0] = Create_c_vector(0, 0, SpConst);
    PositionVector[1] = pNode1->rGetLocation() - pNode0->rGetLocation() + Create_c_vector(0, 0, SpConst);
    PositionVector[2] = pNode2->rGetLocation() - pNode0->rGetLocation() + Create_c_vector(0, 0, SpConst);
    PositionVector = MatrixTranspose(PositionVector);

    NewPositionVector[0] = Create_c_vector(0, 0, SpConst);
    NewPositionVector[1] = Create_c_vector(a, 0, SpConst);
    NewPositionVector[2] = Create_c_vector(b * cos(theta), b * sin(theta), SpConst);
    NewPositionVector = MatrixTranspose(NewPositionVector);

    c_vector<c_vector<long double, 3>, 3> InverseE = MatrixMultiplication(PositionVector, Inverse(NewPositionVector));

    return InverseE;
}

long double MembraneShearForce::det(c_vector<c_vector<long double, 2>, 2> Matrix)
{

    long double Determinate = Matrix[0][0] * Matrix[1][1] - Matrix[0][1] * Matrix[1][0];
    return Determinate;
}

long double MembraneShearForce::tr(c_vector<c_vector<long double, 2>, 2> Matrix)
{

    long double Trace = Matrix[0][0] + Matrix[1][1];
    return Trace;
}

void MembraneShearForce::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<2, 3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MembraneShearForce)

// c_vector<long double, 3> NewForce = ForceMap[node_index] ;
// double NormNewForce = norm_2(NewForce);

// long double Area = 0;
// std::set<unsigned>& containing_elements = pNode->rGetContainingElementIndices();
// assert(containing_elements.size() > 0);
// for (std::set<unsigned>::iterator iter = containing_elements.begin();
//     iter != containing_elements.end();
//     ++iter)
// {
//     Node<3>* pNode0 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(0));
//     Node<3>* pNode1= p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(1));
//     Node<3>* pNode2 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(2));

//     c_vector<long double, 3> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
//     c_vector<long double, 3> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3

//     c_vector<long double, 3> normalVector = VectorProduct(vector_12, vector_13);
//     Area+= norm_2(normalVector)/6;
// }

// }

// UNCOMMENT FOR CYLINDER
// COMMENT FOR NON REGULAR MESH

//  Edge Conditions ::
//   if ( mCylinder == true) // Need to write a code the handles not regular meshes
//   {
//

//         // The nodes around the edge are marked with a mutation, if the node is a mutation need to select
//         // the forces from a node two rows up, this node will have the same x,y location

//             // unsigned ReferenceNode = 0;
//             // if (cell_iter->GetMutationState()->IsType<LostEndothelialCell>()) // If on edge
//             // {
//             //     if (node_index < mNc + 1) // if on lower edge
//             //     {
//             //         ReferenceNode = node_index + (2 * mNc); // select node from two rows up
//             //     }
//             //     else if (node_index > mNc) // if on upper edge
//             //     {
//             //         ReferenceNode = node_index - (2 * mNc); // select node from two rows down
//             //     }
//             //     Node<3>* pReferenceNode = p_cell_population->rGetMesh().GetNode(ReferenceNode);
//             //     // swapped the force over

//             //     if ( mCorrectDrag == false)
//             //     {
//             //     pNode->ClearAppliedForce(); // remove the already present force at this node
//             //     pNode->AddAppliedForceContribution(pReferenceNode->rGetAppliedForce()); // Add the new force
//             //     }

//             //    //PRINT_VECTOR(pReferenceNode->rGetAppliedForce());
//             //     ForceMap[node_index] = ForceMap[ReferenceNode]; // Do the same with the redundant force map
//             // }
//         // pNode->AddAppliedForceContribution(ForceMap[node_index]);
//         // Force contribution is first saved in a map and then added to the simulation.
//         // The inial way to deal with edges was to clear the force on the edge node and replace with the equivialant node two rows
//         // down, but when combinded with the other forces this messes with the additional forces

//         // //  Save the force contribution in a visualiable manner
//         // cell_iter->GetCellData()->SetItem("mesh_element", ForceMap[node_index][0]);
//         // cell_iter->GetCellData()->SetItem("Strain_force_y", ForceMap[node_index][1]);
//         // cell_iter->GetCellData()->SetItem("Strain_force_z", ForceMap[node_index][2]);

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

// void MembraneShearForce::RegularBoundCylinder(bool Cylinder, double Nc , double Nz)
// {

//     // this will be used in an if statment in the AddForceContribution to prevent boundary effects on a basic cylinder
//     // (boundary effects arrising due to the different number of associated elements on nodes at the edge. This particual option is for a strctured mesh)
//     mCylinder = Cylinder;
//     mNc = Nc;
//     mNz = Nz;
// }

// void MembraneShearForce::DragCorrection(bool CorrectDrag)
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

// if (std::abs(NodeLocation[0]) < 1e-20) // along the y axis of the unit circle
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
// } else if (std::abs(NodeLocation[1]) <1e-19) //  along the x axis of the unit circle
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