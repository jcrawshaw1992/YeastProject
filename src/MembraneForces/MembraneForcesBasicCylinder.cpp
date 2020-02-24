#include "MembraneForcesBasicCylinder.hpp"

#include "UblasCustomFunctions.hpp"
#include "MathsFunctions.hpp"


// #include "MathsFunctions.hpp"

MembraneForcesBasicCylinder::MembraneForcesBasicCylinder()
        : AbstractForce<2, 3>()
{
}

// void MembraneForcesBasicCylinder::SetElasticShearModulus(double ElasticShearModulus)
// {
//     mElasticShearModulus = ElasticShearModulus;
// }

// void MembraneForcesBasicCylinder::SetAreaDilationModulus(double AreaDilationModulus)
// {
//     mAreaDilationModulus = AreaDilationModulus;
// }

// void SetMembraneStiffness(double AreaConstant);

// void MembraneForcesBasicCylinder::SetMembraneStiffness(double MembraneSurface)
// {
//     mMembraneSurface = MembraneSurface;
// }



void MembraneForcesBasicCylinder::AddForceContribution(AbstractCellPopulation<2, 3>& rCellPopulation)
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

        double Kalpha;
        double KA;
        double KS;

        // FInd the min node and apply that one 
        c_vector<double, 3> Node0 = pNode0->rGetLocation();
            double Z0 = Node0[2];
            double Z1 = pNode1 ->rGetLocation()[2];
            double Z2 = pNode2 ->rGetLocation()[2];
    

     double MinZ = 1e-3;
     CellPtr p_cell1 = p_cell_population->GetCellUsingLocationIndex(elem_iter->GetNodeGlobalIndex(1));
     CellPtr p_cell2 = p_cell_population->GetCellUsingLocationIndex(elem_iter->GetNodeGlobalIndex(1));
     CellPtr p_cell3 = p_cell_population->GetCellUsingLocationIndex(elem_iter->GetNodeGlobalIndex(1));

     
        if (p_cell1->GetMutationState()->template IsType<HasEndothelialCell>() && p_cell2->GetMutationState()->template IsType<HasEndothelialCell>() && p_cell3->GetMutationState()->template IsType<HasEndothelialCell>() )
        {
                    
            // Kalpha =  pow(10, -8.2459);    
            // KA = pow(10, -6.9) ;
            // KS = pow(10, -9) ;

            p_cell = p_cell_population->GetCellUsingLocationIndex(0);
            Kalpha =p_cell->GetCellData()->GetItem("AreaDilationModulus");
            KA =p_cell->GetCellData()->GetItem("AreaConstant");
            KS =p_cell->GetCellData()->GetItem("ShearModulus");


        }
        else 
        {   
            // Just take the smallest modulus
            // Kalpha = std::min(std::min(p_cell1->GetCellData()->GetItem("AreaDilationModulus"), p_cell2->GetCellData()->GetItem("AreaDilationModulus")),p_cell3->GetCellData()->GetItem("AreaDilationModulus"));
            // KA = std::min(std::min(p_cell1->GetCellData()->GetItem("AreaConstant"), p_cell2->GetCellData()->GetItem("AreaConstant")),p_cell3->GetCellData()->GetItem("AreaConstant"));
            // KS = std::min(std::min(p_cell1->GetCellData()->GetItem("ShearModulus"), p_cell2->GetCellData()->GetItem("ShearModulus")),p_cell3->GetCellData()->GetItem("ShearModulus"));

            if (Node0[2] < 0 )
            { // is on left, we want the max z value 
                if ( Z0 == Z1)
                {
                    if ( Z0 > Z2)
                    {
                        node_index = elem_iter->GetNodeGlobalIndex(0);
                    }
                    else
                    {
                        node_index = elem_iter->GetNodeGlobalIndex(2);
                    }
                } else if ( Z0 == Z2)
                {
                    if ( Z0 > Z1)
                    {
                        node_index = elem_iter->GetNodeGlobalIndex(0);
                    }
                    else
                    {
                        node_index = elem_iter->GetNodeGlobalIndex(1);
                    }
                }else if ( Z1 == Z2)
                {
                    if ( Z0 > Z1)
                    {
                        node_index = elem_iter->GetNodeGlobalIndex(0);
                    }
                    else
                    {
                        node_index = elem_iter->GetNodeGlobalIndex(1);
                    }
                }
            }else if (Node0[2] >= 0 )
            {// is on righ, we want the min z value 
                if ( Z0 == Z1)
                {
                    if ( Z0 < Z2)
                    {
                        node_index = elem_iter->GetNodeGlobalIndex(0);
                    }
                    else
                    {
                        node_index = elem_iter->GetNodeGlobalIndex(2);
                    }
                } else if ( Z0 == Z2)
                {
                    if ( Z0 < Z1)
                    {
                        node_index = elem_iter->GetNodeGlobalIndex(0);
                    }
                    else
                    {
                        node_index = elem_iter->GetNodeGlobalIndex(1);
                    }
                }else if ( Z1 == Z2)
                {
                    if ( Z0 < Z1)
                    {
                        node_index = elem_iter->GetNodeGlobalIndex(0);
                    }
                    else
                    {
                        node_index = elem_iter->GetNodeGlobalIndex(1);
                    }
                }
            }
            p_cell = p_cell_population->GetCellUsingLocationIndex(node_index);
            Kalpha =p_cell->GetCellData()->GetItem("AreaDilationModulus");
            KA =p_cell->GetCellData()->GetItem("AreaConstant");
            KS =p_cell->GetCellData()->GetItem("ShearModulus");

        } 
     


        for (int i = 0; i < 3; i++)
        {
            dedvX = KS * (I1 + 1) * (Dxx * a_i[i] + Dxy * b_i[i]) + (-KS + Kalpha * I2) * ((Dxx * Dyy * Dyy - Dxy * Dyx * Dyy) * a_i[i] + (Dxy * Dyx * Dyx - Dxx * Dyx * Dyy) * b_i[i]);
            dedvY = KS * (I1 + 1) * (Dyx * a_i[i] + Dyy * b_i[i]) + (-KS + Kalpha * I2) * ((Dxx * Dxx * Dyy - Dxx * Dxy * Dyx) * b_i[i] + (Dyx * Dxy * Dxy - Dxx * Dxy * Dyy) * a_i[i]);
            RotatedForceOnNode[i] = -mArea0[elem_index] / 3 * Create_c_vector(dedvX, dedvY, 0);
            RotatedMag[i] = norm_2(RotatedForceOnNode[i]);
            // This is the force for each node  for the triangle situated at the origin
            // Matrix with each row containing the force on the corresponding element
        }




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



        // Area force

        // Calculate the edge vectors of this element
        c_vector<double, 3> vector_A = pNode0->rGetLocation() - pNode2->rGetLocation();
        c_vector<double, 3> vector_B = pNode1->rGetLocation() - pNode2->rGetLocation();

        // Calculate the normal, unit normal and |normal|.
        c_vector<double, 3> normal = VectorProduct(vector_A, vector_B);
        double NormNormal = norm_2(normal);
        c_vector<double, 3> UnitNormal = normal / NormNormal;
      
        // Calculate the area of the element
        double Area = 0.5 * NormNormal;
        double AreaDiff = (Area - mArea0[elem_index])/mArea0[elem_index];
        
          // Force on Node 0
        c_vector<double, 3> vector_12t = pNode2->rGetLocation() - pNode1->rGetLocation();
        ForceOnNode[0] -= 0.5 * KA * AreaDiff * VectorProduct(UnitNormal, vector_12t);
        
        // Force on Node 1
        c_vector<double, 3> vector_20 = pNode0->rGetLocation() - pNode2->rGetLocation();
        ForceOnNode[1] -= 0.5 * KA * AreaDiff * VectorProduct(UnitNormal, vector_20);
//  double A = log10(KA );
//  PRINT_VARIABLE(A);
        // Force on Node 2
        c_vector<double, 3> vector_01 = pNode1->rGetLocation() - pNode0->rGetLocation();
        ForceOnNode[2] -= 0.5 * KA * AreaDiff * VectorProduct(UnitNormal, vector_01);

        double CellArea;
        for (int i = 0; i < 3; i++)
        {
            node_index = elem_iter->GetNodeGlobalIndex(i);
            CellArea= rCellPopulation.GetVolumeOfCell(rCellPopulation.GetCellUsingLocationIndex(node_index));
            ForceOnNode[i] /= CellArea;
            // ForceMap[node_index] += ForceOnNode[i];      
        }
       

        pNode0->AddAppliedForceContribution(ForceOnNode[0]);
        pNode1->AddAppliedForceContribution(ForceOnNode[1]);
        pNode2->AddAppliedForceContribution(ForceOnNode[2]);



    }
}
    /*
    * Set up the inital configuration -> find the inital position vectors, the inital area, and the shape function 
    * 
    * The inital poisitions of all the nodes have been updated to be closer to the radius by a scalling factor (S) 
    * This reflects state where a vessel is deflated with no fluid flow to expand it, leaving only the compressive forces 
    * from the tissue, naturally pushing it inwards. 
    */


void MembraneForcesBasicCylinder::SetupMembraneConfiguration(AbstractCellPopulation<2, 3>& rCellPopulation)
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
    }
}



void MembraneForcesBasicCylinder::FindNeighbours(AbstractCellPopulation<2, 3>& rCellPopulation)
{
    MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);
   
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


c_vector<c_vector<long double, 3>, 3> MembraneForcesBasicCylinder::MatrixMultiplication(c_vector<c_vector<long double, 3>, 3> Matrix1, c_vector<c_vector<long double, 3>, 3> Matrix2)
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

c_vector<long double, 3> MembraneForcesBasicCylinder::MatrixMultiplication(c_vector<c_vector<long double, 3>, 3> Matrix, c_vector<long double, 3> Vector)
{

    c_vector<long double, 3> Answer;

    for (int i = 0; i < 3; i++)
    {
        Answer[i] = inner_prod(Matrix[i], Vector);
    }

    return Answer;
}

c_vector<c_vector<long double, 3>, 3> MembraneForcesBasicCylinder::MatrixTranspose(c_vector<c_vector<long double, 3>, 3> Matrix)
{
    c_vector<c_vector<long double, 3>, 3> MatrixTranspose;
    for (int i = 0; i < 3; i++)
    { // Transporse of the matrix
        MatrixTranspose[i] = Create_c_vector(Matrix[0][i], Matrix[1][i], Matrix[2][i]);
    }
    return MatrixTranspose;
}

c_vector<c_vector<long double, 3>, 3> MembraneForcesBasicCylinder::Inverse(c_vector<c_vector<long double, 3>, 3> Matrix)
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

c_vector<c_vector<long double, 3>, 3> MembraneForcesBasicCylinder::Inverse(c_vector<c_vector<long double, 3>, 3> Matrix, double elem)
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

c_vector<c_vector<long double, 3>, 3> MembraneForcesBasicCylinder::RowReduction(c_vector<c_vector<long double, 3>, 3> Matrix)
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

c_vector<c_vector<long double, 3>, 3> MembraneForcesBasicCylinder::MappingMatrix(MeshBasedCellPopulation<2, 3>* p_cell_population, typename AbstractTetrahedralMesh<2, 3>::ElementIterator elem_iter, double a, double b, double theta)
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

long double MembraneForcesBasicCylinder::det(c_vector<c_vector<long double, 2>, 2> Matrix)
{

    long double Determinate = Matrix[0][0] * Matrix[1][1] - Matrix[0][1] * Matrix[1][0];
    return Determinate;
}

long double MembraneForcesBasicCylinder::tr(c_vector<c_vector<long double, 2>, 2> Matrix)
{

    long double Trace = Matrix[0][0] + Matrix[1][1];
    return Trace;
}

void MembraneForcesBasicCylinder::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<2, 3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MembraneForcesBasicCylinder)
