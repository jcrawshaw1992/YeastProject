#include "MembraneDeformationForce.hpp"

#include "UblasCustomFunctions.hpp"
#include "MathsFunctions.hpp"
#include "BoundariesModifier.hpp"
#include "HistoryDepMeshBasedCellPopulation.hpp"


// #include "MathsFunctions.hpp"

MembraneDeformationForce::MembraneDeformationForce()
        : AbstractForce<2, 3>()
{
}

void MembraneDeformationForce::AddForceContribution(AbstractCellPopulation<2, 3>& rCellPopulation)
{
   HistoryDepMeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<HistoryDepMeshBasedCellPopulation<2, 3>*>(&rCellPopulation);

    std::map<unsigned, c_vector<double, 3> > MembraneForceMap;
    for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
        cell_iter != rCellPopulation.End();
        ++cell_iter)
    {
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        // MembraneForceMap[node_index] = Create_c_vector(0,0,0);
        cell_iter->GetCellData()->SetItem("MembraneForce",0);
    }

    for (typename AbstractTetrahedralMesh<2, 3>::ElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        // THis will be needed later -- going to need to figure out how to stream line this later 
        double Kalpha = 0;
        double KA = 0;
        double KS = 0;
        for (int i = 0; i < 3; i++)
        {
            unsigned node_index = elem_iter->GetNodeGlobalIndex(i);
            CellPtr p_cell = p_cell_population->GetCellUsingLocationIndex(node_index);

            Kalpha +=p_cell->GetCellData()->GetItem("AreaDilationModulus");
            KS +=p_cell->GetCellData()->GetItem("ShearModulus");
            KA += p_cell->GetCellData()->GetItem("AreaConstant"); 
        }
        Kalpha/=3;
        KA/=3;
        KS/=3;
        // PRINT_3_VARIABLES(Kalpha,KA,KS ) 

        unsigned elem_index = elem_iter->GetIndex();

        Node<3>* pNode0 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(0));
        Node<3>* pNode1 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(1));
        Node<3>* pNode2 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(2));

        // c_vector<double, 3> Vector1 = pNode1->rGetLocation();
        // c_vector<double, 3> Vector2 = pNode2->rGetLocation();
        // c_vector<double, 3> Vector0 = pNode0->rGetLocation();
        // PRINT_VECTOR(Vector0 );

        c_vector<double, 3> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
        c_vector<double, 3> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 33);
        


        unsigned node_index;
        // CellPtr p_cell;

        // Get side lengths and angle to create an equal triangle at the origin
        double a = norm_2(vector_12); // Lenght a -> edge connecting P0 and P1
        double b = norm_2(vector_13); // Lenght b -> edge connecting P0 and P2
        double alpha = acos(inner_prod(vector_12, vector_13) / (a * b));

        // This will create an equal triangle at the origin
        c_vector<double, 2> x1 = Create_c_vector(0, 0);
        c_vector<double, 2> x2 = Create_c_vector(a, 0);
        c_vector<double, 2> x3 = Create_c_vector(b * cos(alpha), b * sin(alpha));

        // Get the original triangle
        c_vector<c_vector<double, 2>, 3> InitialVectors = p_cell_population->GetInitalVectors(elem_index);

        //Displacement vectors.
        c_vector<double, 2> V1 = x1 - InitialVectors[0];
        c_vector<double, 2> V2 = x2 - InitialVectors[1];
        c_vector<double, 2> V3 = x3 - InitialVectors[2];

        // Get the shape function coefficents for this element
        c_vector<c_vector<double, 3>, 2>  ShapeFunction = p_cell_population->GetInitalShapeFunction(elem_index);
        c_vector<double, 3> a_i = ShapeFunction[0];
        c_vector<double, 3> b_i = ShapeFunction[1];

        double Area0 = p_cell_population->GetOriginalArea(elem_index);

        // Deformation tensor
        double Dxx = 1 + a_i[0] * V1[0] + a_i[1] * V2[0] + a_i[2] * V3[0];
        double Dxy = b_i[0] * V1[0] + b_i[1] * V2[0] + b_i[2] * V3[0];
        double Dyx = a_i[0] * V1[1] + a_i[1] * V2[1] + a_i[2] * V3[1];
        double Dyy = 1 + b_i[0] * V1[1] + b_i[1] * V2[1] + b_i[2] * V3[1];

        c_vector<c_vector<double, 2>, 2> G;

        // G =DTD  -- Caughy green 
        G[0][0] = Dxx * Dxx + Dyx * Dyx;
        G[0][1] = Dxx * Dxy + Dyx * Dyy;
        G[1][0] = Dxx * Dxy + Dyy * Dyx;
        G[1][1] = Dxy * Dxy + Dyy * Dyy;

        // Strain invarients 
        double I1 = tr(G) - 2;
        double I2 = det(G) - 1;


        c_vector<c_vector<double, 3>, 3> RotatedForceOnNode;
        // c_vector<double, 3> RotatedMag;
    
        double dedvX;
        double dedvY;

        for (int i = 0; i < 3; i++)
        {
            dedvX = KS * (I1 + 1) * (Dxx * a_i[i] + Dxy * b_i[i]) + (-KS + Kalpha * I2) * ((Dxx * Dyy * Dyy - Dxy * Dyx * Dyy) * a_i[i] + (Dxy * Dyx * Dyx - Dxx * Dyx * Dyy) * b_i[i]);
            dedvY = KS * (I1 + 1) * (Dyx * a_i[i] + Dyy * b_i[i]) + (-KS + Kalpha * I2) * ((Dxx * Dxx * Dyy - Dxx * Dxy * Dyx) * b_i[i] + (Dyx * Dxy * Dxy - Dxx * Dxy * Dyy) * a_i[i]);
            RotatedForceOnNode[i] = -Area0 / 3 * Create_c_vector(dedvX, dedvY, 0);
            // RotatedMag[i] = norm_2(RotatedForceOnNode[i]);
            // This is the force for each node  for the triangle situated at the origin
            // Matrix with each row containing the force on the corresponding element
        }

        c_vector<c_vector<double, 3>, 3> X;
        c_vector<c_vector<double, 3>, 3> ForceOnNode;

        X[0] = Create_c_vector(a, 0, 0);
        X[1] = Create_c_vector(b * cos(alpha), b * sin(alpha), 0);
        X[2] = Create_c_vector(0, 0, 1);
        X = MatrixTranspose(X);

        // Rotate the force to the original configuretion

        c_vector<double, 3> F0_rp = MatrixMultiplication(Inverse(X), RotatedForceOnNode[0]);
        c_vector<double, 3> F1_rp = MatrixMultiplication(Inverse(X), RotatedForceOnNode[1]);
        c_vector<double, 3> F2_rp = MatrixMultiplication(Inverse(X), RotatedForceOnNode[2]);

        c_vector<c_vector<double, 3>, 3> Ident = MatrixMultiplication(Inverse(X), X);
        ForceOnNode[0] = F0_rp[0] * vector_12 + F0_rp[1] * vector_13 + F0_rp[2] * X[2];
        ForceOnNode[1] = F1_rp[0] * vector_12 + F1_rp[1] * vector_13 + F1_rp[2] * X[2];
        ForceOnNode[2] = F2_rp[0] * vector_12 + F2_rp[1] * vector_13 + F2_rp[2] * X[2];
        

        // Area force
        
        c_vector<double, 3> vector_A = pNode0->rGetLocation() - pNode2->rGetLocation();
        c_vector<double, 3> vector_B = pNode1->rGetLocation() - pNode2->rGetLocation();

        // Calculate the normal, unit normal and |normal|.
        c_vector<double, 3> normal = VectorProduct(vector_A, vector_B);
        double NormNormal = norm_2(normal);
        c_vector<double, 3> UnitNormal = normal / NormNormal;
      
        // Calculate the area of the element
        double Area = 0.5 * NormNormal;
        double AreaDiff = (Area - Area0)/Area0;
        
          // Force on Node 0
        c_vector<double, 3> vector_12t = pNode2->rGetLocation() - pNode1->rGetLocation();
        ForceOnNode[0] -= 0.5 * KA * AreaDiff * VectorProduct(UnitNormal, vector_12t);
        // Force on Node 1
        c_vector<double, 3> vector_20 = pNode0->rGetLocation() - pNode2->rGetLocation();
        ForceOnNode[1] -= 0.5 * KA * AreaDiff * VectorProduct(UnitNormal, vector_20);
        
        // Force on Node 2
        c_vector<double, 3> vector_01 = pNode1->rGetLocation() - pNode0->rGetLocation();
        ForceOnNode[2] -= 0.5 * KA * AreaDiff * VectorProduct(UnitNormal, vector_01);

        double CellArea;
        double currentForce;
        for (int i = 0; i < 3; i++)
        {
            node_index = elem_iter->GetNodeGlobalIndex(i);
            CellArea= rCellPopulation.GetVolumeOfCell(rCellPopulation.GetCellUsingLocationIndex(node_index));
            ForceOnNode[i] /= CellArea;
            MembraneForceMap[node_index] += ForceOnNode[i];  
            
            CellPtr p_cell = p_cell_population->GetCellUsingLocationIndex(node_index);
            currentForce = p_cell->GetCellData()->GetItem("MembraneForce");
            p_cell->GetCellData()->SetItem("MembraneForce",currentForce + norm_2(ForceOnNode[i]));


        }
    
        pNode0->AddAppliedForceContribution(ForceOnNode[0]);
        pNode1->AddAppliedForceContribution(ForceOnNode[1]);
        pNode2->AddAppliedForceContribution(ForceOnNode[2]);
        
        
       

        
       
    }

    // // // THis bit takes care of the edges ...
    // for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
    //         cell_iter != rCellPopulation.End();
    //         ++cell_iter)
    // {
        
    //     unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
    //     Node<3>* pNode = p_cell_population->rGetMesh().GetNode(node_index);  
    //     if (cell_iter->GetCellData()->GetItem("Boundary") == 1)
    //     {
    //         c_vector<double, 3> AverageForce = Create_c_vector(0,0,0);
    //         c_vector<unsigned, 2> NearestNodes =  p_cell_population->GetNearestInternalNodes(node_index);
    //         PRINT_VECTOR(NearestNodes);

    //         for ( int i = 0; i <2; i++)
    //         {  
    //             AverageForce += MembraneForceMap[NearestNodes[i]];
    //         }
    //         AverageForce/=2;
    //         PRINT_VARIABLE(AverageForce)
    //        c_vector<double, 3>  Normal = MembraneForceMap[node_index] /norm_2(MembraneForceMap[node_index]) ;
    //         // AverageForce/=norm_2(AverageForce);
    //         AverageForce = norm_2(AverageForce) * Normal;
    //         cell_iter->GetCellData()->SetItem("MembraneForce",norm_2(AverageForce) );
    //         pNode->AddAppliedForceContribution(AverageForce); // Add the new force
    //     }
    //     else
    //     {
    //         cell_iter->GetCellData()->SetItem("MembraneForce", norm_2(MembraneForceMap[node_index] ));
    //         pNode->AddAppliedForceContribution(MembraneForceMap[node_index] ); // Add the new force
    //     }
  

    // }
        
}


/*
* Save a map of the closest nodes for the boudary nodes
*/


void MembraneDeformationForce::SetNearestNodesForBoundaryNodes(std::map<unsigned, c_vector<unsigned, 2> > NearestNodesMap)
{

    mNearestNodesMap = NearestNodesMap;
     
}



/*
* Set up the inital configuration -> find the inital position vectors, the inital area, and the shape function 
* 
* The inital poisitions of all the nodes have been updated to be closer to the radius by a scalling factor (S) 
* This reflects state where a vessel is deflated with no fluid flow to expand it, leaving only the compressive forces 
* from the tissue, naturally pushing it inwards. 
*/

c_vector<c_vector<double, 3>, 3> MembraneDeformationForce::MatrixMultiplication(c_vector<c_vector<double, 3>, 3> Matrix1, c_vector<c_vector<double, 3>, 3> Matrix2)
{

    c_vector<c_vector<double, 3>, 3> MatrixTranspose;
    c_vector<c_vector<double, 3>, 3> Answer;

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

c_vector<double, 3> MembraneDeformationForce::MatrixMultiplication(c_vector<c_vector<double, 3>, 3> Matrix, c_vector<double, 3> Vector)
{

    c_vector<double, 3> Answer;

    for (int i = 0; i < 3; i++)
    {
        Answer[i] = inner_prod(Matrix[i], Vector);
    }

    return Answer;
}

c_vector<c_vector<double, 3>, 3> MembraneDeformationForce::MatrixTranspose(c_vector<c_vector<double, 3>, 3> Matrix)
{
    c_vector<c_vector<double, 3>, 3> MatrixTranspose;
    for (int i = 0; i < 3; i++)
    { // Transporse of the matrix
        MatrixTranspose[i] = Create_c_vector(Matrix[0][i], Matrix[1][i], Matrix[2][i]);
    }
    return MatrixTranspose;
}

c_vector<c_vector<double, 3>, 3> MembraneDeformationForce::Inverse(c_vector<c_vector<double, 3>, 3> Matrix)
{

    c_vector<c_vector<double, 3>, 3> MTranspose = MatrixTranspose(Matrix);

    double det = Matrix[0][0] * (Matrix[1][1] * Matrix[2][2] - Matrix[1][2] * Matrix[2][1]) - Matrix[0][1] * (Matrix[1][0] * Matrix[2][2] - Matrix[1][2] * Matrix[2][0]) + Matrix[0][2] * (Matrix[1][0] * Matrix[2][1] - Matrix[1][1] * Matrix[2][0]);

    double l = 0, m = 1;
    c_vector<c_vector<double, 3>, 3> InverseMatrix;
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

c_vector<c_vector<double, 3>, 3> MembraneDeformationForce::Inverse(c_vector<c_vector<double, 3>, 3> Matrix, double elem)
{

    double det = Matrix[0][0] * (Matrix[1][1] * Matrix[2][2] - Matrix[1][2] * Matrix[2][1]) - Matrix[0][1] * (Matrix[1][0] * Matrix[2][2] - Matrix[1][2] * Matrix[2][0]) + Matrix[0][2] * (Matrix[1][0] * Matrix[2][1] - Matrix[1][1] * Matrix[2][0]);

    c_vector<c_vector<double, 3>, 3> MTranspose = MatrixTranspose(Matrix);

    double l = 0, m = 1;
    c_vector<c_vector<double, 3>, 3> InverseMatrix;
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

double MembraneDeformationForce::det(c_vector<c_vector<double, 2>, 2> Matrix)
{

    double Determinate = Matrix[0][0] * Matrix[1][1] - Matrix[0][1] * Matrix[1][0];
    return Determinate;
}

double MembraneDeformationForce::tr(c_vector<c_vector<double, 2>, 2> Matrix)
{

    double Trace = Matrix[0][0] + Matrix[1][1];
    return Trace;
}

void MembraneDeformationForce::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<2, 3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MembraneDeformationForce)


  // double AppliedPressure = norm_2(AverageForce);        
            // // Loop over neighbouring elements to get normal 

            // c_vector<long double, 3> Normal = zero_vector<long double>(3);
            // std::set<unsigned>& containing_elements = pNode->rGetContainingElementIndices();
            // for (std::set<unsigned>::iterator iter = containing_elements.begin();
            //     iter != containing_elements.end();
            //     ++iter)
            // {
            //     Node<3>* pNode0 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(0));
            //     Node<3>* pNode1= p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(1));
            //     Node<3>* pNode2 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(2));

            //     c_vector<long double, 3> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
            //     c_vector<long double, 3> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3

            //     Normal += VectorProduct(vector_12, vector_13);
            //     // Dont know if always pointing the right way
            // }
            
            // if (inner_prod(AverageForce, Normal)<0)
            // {
            //     Normal/=-norm_2(Normal);
            // }else{
            // Normal/=norm_2(Normal);
            // }