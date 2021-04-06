#include "MembraneSurfaceForce.hpp"
#include "Debug.hpp"
#include "UblasCustomFunctions.hpp"
// Area constrant on a membrane

MembraneSurfaceForce::MembraneSurfaceForce()
        : AbstractForce<2, 3>()
{
}

void MembraneSurfaceForce::SetMembraneStiffness(double membrane_constant) // Might need later so just scilence for now 
{
    mMembraneSurface = membrane_constant;
}

void MembraneSurfaceForce::AddForceContribution(AbstractCellPopulation<2, 3>& rCellPopulation)
{
    // TRACE("Add SurfaceArea Force");
  MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);
	
/*
 * Cast the cell population to the MeshBasedCellPopulation so that I can loop over all of the elements, calulcate the area
 * of each element, then caclulate the area difference and add the area force contribution of this element to each of the 
 * of this element.
*/

    for (typename AbstractTetrahedralMesh<2, 3>::ElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();

        // Get the nodes of this element
        Node<3>* pNode0 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(0));
        Node<3>* pNode1 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(1));
        Node<3>* pNode2 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(2));

       
        
        // Calculate the edge vectors of this element
        c_vector<double, 3> vector_A = pNode0->rGetLocation() - pNode2->rGetLocation();
        c_vector<double, 3> vector_B = pNode1->rGetLocation() - pNode2->rGetLocation();

        // Calculate the normal, unit normal and |normal|.
        c_vector<double, 3> normal = VectorProduct(vector_A, vector_B);
        double NormNormal = norm_2(normal);
        c_vector<double, 3> UnitNormal = normal / NormNormal;
      
        // Calculate the area of the element
        double Area = 0.5 * NormNormal;
        double AreaDiff = (Area - mOriginalAreas[elem_index]) / mOriginalAreas[elem_index];
        
        c_vector<c_vector<long double, 3>, 3> ForceOnNode;
        
        // Force on Node 0
        c_vector<double, 3> vector_12 = pNode2->rGetLocation() - pNode1->rGetLocation();
        ForceOnNode[0] = -0.5 * mMembraneSurfaceMaps[elem_index] * AreaDiff * VectorProduct(UnitNormal, vector_12);
        
        // Force on Node 1
        c_vector<double, 3> vector_20 = pNode0->rGetLocation() - pNode2->rGetLocation();
        ForceOnNode[1] = -0.5 * mMembraneSurfaceMaps[elem_index] * AreaDiff * VectorProduct(UnitNormal, vector_20);
//  double A = log10(mMembraneSurfaceMaps[elem_index] );
//  PRINT_VARIABLE(A);
        // Force on Node 2
        c_vector<double, 3> vector_01 = pNode1->rGetLocation() - pNode0->rGetLocation();
        ForceOnNode[2] = -0.5 * mMembraneSurfaceMaps[elem_index] * AreaDiff * VectorProduct(UnitNormal, vector_01);


        unsigned node_index;
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



void MembraneSurfaceForce::SetScallingArea(double Scalling)
{
       mScalling=Scalling;
}


/*
 * This method is called at the begining, specifically from the primary simulation to calculate the areas of each element. This 
 * is then saved as a protected member variable (a map) with the element index as the key, and the area as the vaule. This member 
 * can be easily accessed from the other methods in this class. 
*/
void MembraneSurfaceForce::SetupInitialAreas(AbstractCellPopulation<2, 3>& rCellPopulation)
{

    //std::map<unsigned, double> InitialAreaMap;
    MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);
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
            normal += - p_cell_population->rGetMesh().GetElement(*iter)->CalculateNormal();
        }
        normal /= norm_2(normal);

        // Need to walk backward into the mesh by the scalling factor 
        c_vector<double, 3> PositionVector = NodeLocation ;//- mScalling * normal;


        (cell_iter)->GetCellData()->SetItem("Initial_Location_X", NodeLocation[0]);
        (cell_iter)->GetCellData()->SetItem("Initial_Location_Y", NodeLocation[1]);
        (cell_iter)->GetCellData()->SetItem("Initial_Location_Z", NodeLocation[2]);

    }



    for (typename AbstractTetrahedralMesh<2, 3>::ElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        // Get location of each node in the element
        unsigned elem_index = elem_iter->GetIndex();
        c_vector<c_vector<double, 3>, 3> PositionVector;
        
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
        c_vector<double, 3> vector_A = PositionVector[0] - PositionVector[2] ;
        c_vector<double, 3> vector_B = PositionVector[1]  - PositionVector[2] ;

        c_vector<double, 3> normal = VectorProduct(vector_A, vector_B);

        // Find the inital area
        double Area = 0.5 * norm_2(normal);

        // Save the inital area in the Original area map
         mOriginalAreas[elem_index] = Area;
         mMembraneSurfaceMaps[elem_index]  = mMembraneSurface;


        //  if (std::abs(PositionVector[0][2]) < 1.6e-3 || std::abs(PositionVector[1][2]) < 1.5e-3 || std::abs(PositionVector[2][2]) < 1.5e-3 )
        // {
        //     mMembraneSurfaceMaps[elem_index]  = mMembraneSurface*10;
        // }else 
        // {
            mMembraneSurfaceMaps[elem_index]  = mMembraneSurface;
        // }
        



    }
}




// modify this such that it is only altering the membrane forces 
void MembraneSurfaceForce::UpdateMembraneSurfaceForceProperties(AbstractCellPopulation<2, 3>& rCellPopulation)
{
    // TRACE("Modify area properties");
    // MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);

    //  c_vector<c_vector<double, 3>, 3> PositionVectors;
    //  std::vector<unsigned> LabledElements; 

    //      for (typename AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
    //      cell_iter != rCellPopulation.End();
    //      ++cell_iter)
    // {
    //     // if mutant
    //       if ((cell_iter)->GetMutationState()->IsType<LostEndothelialCell>()) // If on edge
    //         {
    //             unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter); //p_cell_population->GetLocationIndexUsingCell(*cell_iter);
    //             Node<3>* p_node = rCellPopulation.GetNode(node_index);
    //             // This piece records all the elements associated with a mutant cell
    //             std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();
    //             assert(containing_elements.size() > 0);
    //             // LabledElements.push_back(containing_elements);
    //             for (std::set<unsigned>::iterator iter = containing_elements.begin();
    //                 iter != containing_elements.end();
    //                 ++iter)
    //             {
    //                 LabledElements.push_back(*iter);
    //             }
    //         }
    // }
    // sort(LabledElements.begin(), LabledElements.end());
    // LabledElements.erase(unique(LabledElements.begin(), LabledElements.end()), LabledElements.end());
    
    // for (int i = 0; i < LabledElements.size(); i++) //
    // {
    //     unsigned elem_index = LabledElements[i];

    //         unsigned node_index;
    //         CellPtr p_cell ; 
    //         Node<3>* pNode;
    //         double NumberOfMutants=0;
    
    //         for (int j = 0; j < 3; j++) //
    //             {
    //                 pNode = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(LabledElements[i])->GetNodeGlobalIndex(j));
    //                 node_index = pNode->GetIndex();
    //                 p_cell = p_cell_population->GetCellUsingLocationIndex(node_index);
    //                 PositionVectors[j][0] = p_cell->GetCellData()->GetItem("Initial_Location_X");
    //                 PositionVectors[j][1] = p_cell->GetCellData()->GetItem("Initial_Location_Y");
    //                 PositionVectors[j][2] = p_cell->GetCellData()->GetItem("Initial_Location_Z");

    //                 if (( p_cell)->GetMutationState()->IsType<LostEndothelialCell>()) // If on edge
    //                 {
    //                     NumberOfMutants+=1;
    //                 }
    //             }

    //         c_vector<double, 3> Adapted_vector_12 = PositionVectors[1] - PositionVectors[0]; // Vector 1 to 2
    //         c_vector<double, 3> Adapted_vector_13 = PositionVectors[2] - PositionVectors[0] ; // Vector 1 to 3

    //         // Find the intial area, A0 = 0.5*norm(normal)
    //         c_vector<long double, 3> normalVector = VectorProduct(Adapted_vector_12, Adapted_vector_13);
    //         // PRINT_VECTOR(normalVector);
    //         mOriginalAreas[elem_index] =  0.5 * norm_2(normalVector);    

    //         if (NumberOfMutants ==1 ) // 6 - NumberOfNeighbourinEcs = Basement matrix sites 
    //         {
    //          mMembraneSurfaceMaps[elem_index] /=10;
    //         } else if(NumberOfMutants == 2)
    //         {
    //          mMembraneSurfaceMaps[elem_index] /=25;
    //         } else if(NumberOfMutants == 3)
    //         {
    //          mMembraneSurfaceMaps[elem_index] /=50;
    //         }
    // } 
}






void MembraneSurfaceForce::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<2, 3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MembraneSurfaceForce)







    // for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
    //         cell_iter != rCellPopulation.End();
    //         ++cell_iter)
    //     {
    //         unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
    //          Node<3>* pNode = p_cell_population->rGetMesh().GetNode(node_index);
    //         double CellArea = rCellPopulation.GetVolumeOfCell(*cell_iter);
    //         ForceMap[node_index]/=CellArea ;
    //         pNode->AddAppliedForceContribution(ForceMap[node_index]);
    //     }

//     for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
//          cell_iter != rCellPopulation.End();
//          ++cell_iter)
//     {
//         unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
//         Node<3>* pNode = p_cell_population->rGetMesh().GetNode(node_index);

//         long double Area = 0;
//         std::set<unsigned>& containing_elements = pNode->rGetContainingElementIndices();
//         assert(containing_elements.size() > 0);
//         for (std::set<unsigned>::iterator iter = containing_elements.begin();
//             iter != containing_elements.end();
//             ++iter)
//         {
//             Node<3>* pNode0 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(0));
//             Node<3>* pNode1 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(1));
//             Node<3>* pNode2 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(2));

//             c_vector<long double, 3> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
//             c_vector<long double, 3> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3

//             c_vector<long double, 3> normalVector = VectorProduct(vector_12, vector_13);
//             Area+= norm_2(normalVector)/6;
//         }
//         // Now get the force map
//         // c_vector<long double, 3> OldForce = ForceMap[node_index] ;
//         // double NormOldForce = norm_2(OldForce);
//         ForceMap[node_index] = ForceMap[node_index] / Area;
//         c_vector<long double, 3> NewForce = ForceMap[node_index] ;
//         // double NormNewForce = norm_2(NewForce);
//         pNode->AddAppliedForceContribution(ForceMap[node_index]);

//         // TRACE("Drag corrected");

//   }

// UNCOMMENT FOR CYLINDER
// COMMENT FOR NON REGULAR MESH

//Edge Conditions ::

//  if ( mCylinder == true) // Need to write a code the handles not regular meshes
//     {
//         for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
//             cell_iter != rCellPopulation.End();
//             ++cell_iter)
//         {
//             unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
//             Node<3>* pNode = p_cell_population->rGetMesh().GetNode(node_index);

//             // The nodes around the edge are marked with a mutation, if the node is a mutation need to select
//             // the forces from a node two rows up, this node will have the same x,y location
//             unsigned ReferenceNode = 0;
//             if (cell_iter->GetMutationState()->IsType<BetaCateninOneHitCellMutationState>()) // If on edge
//             {
//                 if (node_index < Nc + 1) // if on lower edge
//                 {
//                     ReferenceNode = node_index + (2 * Nc); // select node from two rows up
//                 }
//                 else if (node_index > Nc) // if on upper edge
//                 {
//                     ReferenceNode = node_index - (2 * Nc); // select node from two rows down
//                 }
//                 Node<3>* pReferenceNode = p_cell_population->rGetMesh().GetNode(ReferenceNode);
//                 // swapped the force over

//                 pNode->ClearAppliedForce(); // remove the already present force at this node
//                 pNode->AddAppliedForceContribution(pReferenceNode->rGetAppliedForce()); // Add the new force
//             // PRINT_VECTOR(pReferenceNode->rGetAppliedForce());
//                 //ForceMap[node_index] = ForceMap[ReferenceNode]; // Do the same with the redundant force map
//             }

//         }
//     }

//    unsigned node_index = elem_iter->GetNodeGlobalIndex(0);
//     ForceMap[node_index] += ForceNode0;

//    node_index = elem_iter->GetNodeGlobalIndex(1);
//    ForceMap[node_index] += ForceNode1;

//    node_index = elem_iter->GetNodeGlobalIndex(2);
//    ForceMap[node_index] += ForceNode2;



        // Check you have an outwards facing normal.
        //  c_vector<double, 3> PositionVector  = pNode0->rGetLocation();
        //  c_vector<double, 3> DirectionVector = UnitNormal + PositionVector;
        //  double RadialDiretion = sqrt(pow(DirectionVector[0],2) + pow(DirectionVector[1],2) );
        //  double RadialPosition = sqrt(pow(PositionVector [0],2) + pow(PositionVector[1],2) );

        //  if (RadialDiretion < RadialPosition)
        //  {
        //     UnitNormal = -UnitNormal;
        //     TRACE("Maintained outwards pointing normal");
        //  };



         // for (typename AbstractTetrahedralMesh<2, 3>::ElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
    //      elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
    //      ++elem_iter)
    // {
        
    //     unsigned elem_index = elem_iter->GetIndex();

    //     // Check the mutant state of each of the cells in the element, then set this number to decide the 
    //     // inital area and the membrane constants of the element 
    //     unsigned node_index;
    //     CellPtr p_cell ;
    //     double NumberOfMutantCells = 0;  
    //     //  // Iterate over the cells in this element to check if they are mutant
    //     //  want to know how many mutant cells we have
    //      for (int i = 0; i < 3; i++)
    //     {    
    //         node_index = elem_iter->GetNodeGlobalIndex(i);
    //         p_cell = p_cell_population->GetCellUsingLocationIndex(node_index);
    //         if (p_cell->GetMutationState()->IsType<LostEndothelialCell>()) // If on edge
    //         {
    //             NumberOfMutantCells+=1;
    //         }


    //     }

    //     if (NumberOfMutantCells ==3 )
    //     {   

    //       for (int i = 0; i < 3; i++)
    //         {    
    //             node_index = elem_iter->GetNodeGlobalIndex(i);
    //             p_cell = p_cell_population->GetCellUsingLocationIndex(node_index);
    //             // The adapted node locations are set in Shear strain force. 
    //             PositionVectors[i][0] =  p_cell->GetCellData()->GetItem("Initial_Location_X");
    //             PositionVectors[i][1] =  p_cell->GetCellData()->GetItem("Initial_Location_Y");
    //             PositionVectors[i][2] =  p_cell->GetCellData()->GetItem("Initial_Location_Z");
    //         }

    //     c_vector<double, 3> Adapted_vector_12 = PositionVectors[1] - PositionVectors[0]; // Vector 1 to 2
    //     c_vector<double, 3> Adapted_vector_13 = PositionVectors[2] - PositionVectors[0] ; // Vector 1 to 3


    //     // Find the intial area, A0 = 0.5*norm(normal)
    //     c_vector<long double, 3> normalVector = VectorProduct(Adapted_vector_12, Adapted_vector_13);
    //     // PRINT_VECTOR(normalVector);
    //     mOriginalAreas[elem_index] =  0.5 * norm_2(normalVector);
    //     }     
    // }