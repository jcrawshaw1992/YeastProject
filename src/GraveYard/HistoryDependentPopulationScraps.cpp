
// VtkMeshReader<ELEMENT_DIM, SPACE_DIM> mesh_reader(mChasteOutputDirectory + "RemeshedGeometry.vtu");
// mNew_mesh.ConstructFromMeshReader(mesh_reader);

//  TRACE("History dependent remeshing")

// static_cast<HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>&>((this->mrMesh)).AddANewNodeBehindBoundary();
//     VtkMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer("CheckingremeshsedMesh", "config", false);
// mesh_writer.WriteFilesUsingMesh(static_cast<AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>&>((this->mrMesh)));

//     // Now look at adding a new cell for the new node
//     // Option 1, follow the line of constructors to fine where I need to clear/delete cells and add new cells on

// old, when using vmtk and stls
/*
	 * Here I will remesh the current chaste geometery
     * 1) Take the latest .vtu file convert it into an .stl 
     * 2) Remesh this stl
     * 3) Convert this stl into a .vtu so I can read it into chaste and use it  
	 */

// // Convert the .vtu into an .stl

// std::string stlfile = mChasteOutputDirectory + "CurrentMesh.stl";
// TRACE("Create stl")
// PRINT_VARIABLE(mChasteOutputDirectory)
// std::string vtu2stlCommand = "python projects/VascularRemodelling/apps/ConvertCurrentGeometryToSTL.py -ChasteOutput " + mChasteOutputDirectory + "  -stlOutput " + stlfile;
// std::system(vtu2stlCommand.c_str()); // system only takes char *

// // Remesh the .stl
// std::string Remeshedstl = mChasteOutputDirectory + "RemeshedGeometry.stl";
// TRACE("Remeshing geometry")

// // More info for the remeshing options found at https://sourceforge.net/p/vmtk/mailman/message/29391820/ and  http://www.vmtk.org/vmtkscripts/vmtksurfaceremeshing.html
// std::string RemeshCommand = "vmtksurfaceremeshing -ifile " + stlfile + " -iterations " + std::to_string(mIterations) + " -area " + std::to_string(mTargetRemeshingElementArea) + " -maxarea "+ std::to_string(1.2*mTargetRemeshingElementArea) +" -ofile " + Remeshedstl;
// std::system(RemeshCommand.c_str());

// //Finally conver the Remeshed stl into a vtu -- meshio is your friend
// std::string Remeshedvtu = mChasteOutputDirectory + "RemeshedGeometry.vtu";
// std::string stl2vtuCommand = " meshio-convert " + Remeshedstl + " " + Remeshedvtu;
// std::system(stl2vtuCommand.c_str());

//         // Read in the new mesh
// VtkMeshReader<ELEMENT_DIM, SPACE_DIM> mesh_reader(mChasteOutputDirectory + "RemeshedGeometry.vtu");
// mNew_mesh.ConstructFromMeshReader(mesh_reader);
// TRACE("Have the new mesh ;) ");

// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::CheckCurvature()
// {
//     // Loop over all springs and check the angles between each
//     for (typename MeshBasedCellPopulation<2, 3>::SpringIterator spring_iterator = p_static_cast_cell_population->SpringsBegin();
//          spring_iterator != p_static_cast_cell_population->SpringsEnd();
//          ++spring_iterator)
//     {

//         Node<3>* pNode1 = spring_iterator.GetNodeA();
//         Node<3>* pNode3 = spring_iterator.GetNodeB();

//         std::pair<Node<3>*, Node<3>*> edge = std::pair<Node<3>*, Node<3>*>(pNode1, pNode3);
//         double OriginalAngle = p_static_cast_cell_population->GetOriginalAngle(edge);

//         std::pair<c_vector<double, 3>, c_vector<double, 3> > nonUnitNormals;

//         std::pair<Node<3>*, Node<3>*> otherNodes;

//         bool boundary_edge_found = p_static_cast_cell_population->CalculateElementNormals( edge, nonUnitNormals, otherNodes);

// }

// ----------------------------------------------------------------------

/*

        //  c_vector<double,3> Bin(xStep, yStep, zStep); //This is the key in the map, then the value is some vector that this is added to
       

        // // Now I need to get the nearest centroid
        // for (typename std::map<unsigned, c_vector<double, SPACE_DIM> >::iterator Centroid_iter = mCentroidMap.begin();
        //      Centroid_iter != mCentroidMap.end();
        //      ++Centroid_iter)
        // {
        //     // PRINT_2_VARIABLES(Centroid_iter->first, Centroid_iter->second)
        //     // if (abs(norm_2(NewNodeLocation - Centroid_iter->second)) < distance && PointInTriangle3D(NewNodeLocation,  Centroid_iter->first))
        //     // {
        //      if (PointInTriangle3D(NewNodeLocation,  Centroid_iter->first) ==1)
        //     {
        //         // now test if it is in the element %
        //         ClosestElement = Centroid_iter->first;
        //         break;

        //     }
        // }

        //  if (ClosestElement == -10)  //
        // {
        //     TRACE("If this fails the new mesh might not be a good representation of the old mesh, will need to loop over everything again")
        //     TRACE("first try 2D edge checking .....")
        //     for (typename std::map<unsigned, c_vector<double, SPACE_DIM> >::iterator Centroid_iter = mCentroidMap.begin();
        //      Centroid_iter != mCentroidMap.end();
        //      ++Centroid_iter)
        //     {
        //         // PRINT_2_VARIABLES(Centroid_iter->first, Centroid_iter->second)
        //           if (PointInTriangle3D2(NewNodeLocation,  Centroid_iter->first) ==1)
        //             {
        //                 // now test if it is in the element %
        //                 ClosestElement = Centroid_iter->first;
        //                 break;

        //             }
        //     }
        // }

        // if (ClosestElement == -10)  //
        // {
        // TRACE(" 2D edge checking didnt work, now need to try cloeset center point -- there is something wrong")
        for (typename std::map<unsigned, c_vector<double, SPACE_DIM> >::iterator Centroid_iter = mCentroidMap.begin();
             Centroid_iter != mCentroidMap.end();
             ++Centroid_iter)
        {
            // PRINT_2_VARIABLES(Centroid_iter->first, Centroid_iter->second)
            if (abs(norm_2(NewNodeLocation - Centroid_iter->second)) < distance)
            {
                // now test if it is in the element %
                ClosestElement = Centroid_iter->first;
                distance = abs(norm_2(NewNodeLocation - Centroid_iter->second));
            }
        }
       

        // PRINT_VARIABLE(InTriangle)
        distance *= 10;
        unsigned WrongClosestElement = ClosestElement;
        if (InTriangle == 0)
        {

            // this is the wrong triangle, try the next one
            for (typename std::map<unsigned, c_vector<double, SPACE_DIM> >::iterator Centroid_iter = mCentroidMap.begin();
                 Centroid_iter != mCentroidMap.end();
                 ++Centroid_iter)
            {
                // PRINT_2_VARIABLES(Centroid_iter->first, Centroid_iter->second)
                if (abs(norm_2(NewNodeLocation - Centroid_iter->second)) < distance && Centroid_iter->first != WrongClosestElement)
                {
                    // now test if it is in the element %
                    ClosestElement = Centroid_iter->first;
                    distance = abs(norm_2(NewNodeLocation - Centroid_iter->second));
                }
            }
        }
        InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
        distance *= 10;
        unsigned WrongClosestElement2 = ClosestElement;
        if (InTriangle == 0)
        {
            TRACE("This is the third try")
            // this is the wrong triangle, try the next one
            for (typename std::map<unsigned, c_vector<double, SPACE_DIM> >::iterator Centroid_iter = mCentroidMap.begin();
                 Centroid_iter != mCentroidMap.end();
                 ++Centroid_iter)
            {
                // PRINT_2_VARIABLES(Centroid_iter->first, Centroid_iter->second)
                if (abs(norm_2(NewNodeLocation - Centroid_iter->second)) < distance && Centroid_iter->first != WrongClosestElement && Centroid_iter->first != WrongClosestElement2)
                {
                    // now test if it is in the element %
                    ClosestElement = Centroid_iter->first;
                    distance = abs(norm_2(NewNodeLocation - Centroid_iter->second));
                }
            }
        }
        InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
        distance *= 10;
        unsigned WrongClosestElement3 = ClosestElement;

        if (InTriangle == 0)
        {
            TRACE("This is the fourth try")
            // this is the wrong triangle, try the next one
            for (typename std::map<unsigned, c_vector<double, SPACE_DIM> >::iterator Centroid_iter = mCentroidMap.begin();
                 Centroid_iter != mCentroidMap.end();
                 ++Centroid_iter)
            {
                // PRINT_2_VARIABLES(Centroid_iter->first, Centroid_iter->second)
                if (abs(norm_2(NewNodeLocation - Centroid_iter->second)) < distance && Centroid_iter->first != WrongClosestElement && Centroid_iter->first != WrongClosestElement2 && Centroid_iter->first != WrongClosestElement3)
                {
                    // now test if it is in the element %
                    ClosestElement = Centroid_iter->first;
                    distance = abs(norm_2(NewNodeLocation - Centroid_iter->second));
                }
            }
        }

        InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
        distance *= 10;
        unsigned WrongClosestElement4 = ClosestElement;

        if (InTriangle == 0)
        {
            TRACE("This is the fith try")
            // this is the wrong triangle, try the next one
            for (typename std::map<unsigned, c_vector<double, SPACE_DIM> >::iterator Centroid_iter = mCentroidMap.begin();
                 Centroid_iter != mCentroidMap.end();
                 ++Centroid_iter)
            {
                // PRINT_2_VARIABLES(Centroid_iter->first, Centroid_iter->second)
                if (abs(norm_2(NewNodeLocation - Centroid_iter->second)) < distance && Centroid_iter->first != WrongClosestElement4 && Centroid_iter->first != WrongClosestElement && Centroid_iter->first != WrongClosestElement2 && Centroid_iter->first != WrongClosestElement3)
                {
                    // now test if it is in the element %
                    ClosestElement = Centroid_iter->first;
                    distance = abs(norm_2(NewNodeLocation - Centroid_iter->second));
                }
            }
        }
        InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
        distance *= 10;
        unsigned WrongClosestElement5 = ClosestElement;

        if (InTriangle == 0)
        {
            TRACE("This is the sixth try")
            // this is the wrong triangle, try the next one
            for (typename std::map<unsigned, c_vector<double, SPACE_DIM> >::iterator Centroid_iter = mCentroidMap.begin();
                 Centroid_iter != mCentroidMap.end();
                 ++Centroid_iter)
            {
                // PRINT_2_VARIABLES(Centroid_iter->first, Centroid_iter->second)
                if (abs(norm_2(NewNodeLocation - Centroid_iter->second)) < distance && Centroid_iter->first != WrongClosestElement5 && Centroid_iter->first != WrongClosestElement4 && Centroid_iter->first != WrongClosestElement && Centroid_iter->first != WrongClosestElement2 && Centroid_iter->first != WrongClosestElement3)
                {
                    // now test if it is in the element %
                    ClosestElement = Centroid_iter->first;
                    distance = abs(norm_2(NewNodeLocation - Centroid_iter->second));
                }
            }
        }
        InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);

        if (InTriangle == 0)
        {
            ClosestElement = WrongClosestElement;
        }

        // PRINT_VARIABLE(InTriangle)
        // I want to know what this says about being in the same ele
        // }
        */

// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::UpdateCellData()
// {
//      // When outputting any CellData, we assume that the first cell is representative of all cells
//     unsigned num_cell_data_items = this->Begin()->GetCellData()->GetNumItems();
//     std::vector<std::string> cell_data_names = this->Begin()->GetCellData()->GetKeys();

//     std::vector<std::vector<double> > cell_data;
//     for (unsigned var=0; var<num_cell_data_items; var++)
//     {
//         std::vector<double> cell_data_var(num_cells_from_mesh);
//         cell_data.push_back(cell_data_var);
//     }

//     if (mOutputMeshInVtk)
//     {
//         // Create mesh writer for VTK output
//         VtkMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(rDirectory, "mesh_"+time.str(), false);
//         mesh_writer.WriteFilesUsingMesh(rGetMesh());
//     }

// }

// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetBinningRegions()
// {
//     SetMaxDomainDimensions();

//     SetCentroidMap();
//     // Clear the map so I dont have the old elements and the new elements messing the binning up

//     mBinMap.clear();
//     // Need to iterate over the elements and determine which bin each element centroid is in. I look for the cloesest old centeroid for each new node, so it fits that I will have the centroids sorted in the bins
//     for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator elem_iter = this->rGetMesh().GetElementIteratorBegin();
//          elem_iter != this->rGetMesh().GetElementIteratorEnd();
//          ++elem_iter)
//     {
//         unsigned elem_index = elem_iter->GetIndex();
//         c_vector<double, 3> Centroid = mCentroidMap[elem_index];
//         bool FoundXRegion = 0;
//         bool FoundYRegion = 0;
//         bool FoundZRegion = 0;
//         int xStep = 1;
//         int yStep = 1;
//         int zStep = 1;
//         double FudgeX = mFudgeScalling * (mMaxX - mMinX);
//         double FudgeY = mFudgeScalling * (mMaxY - mMinY);
//         double FudgeZ = mFudgeScalling * (mMaxZ - mMinZ);

//         bool FudgeInX = 0;
//         bool FudgeInY = 0;
//         bool FudgeInZ = 0;
//         bool lowerFudgeX = 0;
//         bool lowerFudgeY = 0;
//         bool lowerFudgeZ = 0;
//         while (FoundXRegion == 0 && xStep <= mNx)
//         {
//             if (Centroid[0] - mMinX <= (mMaxX - mMinX) * xStep / mNx) // Look for x bins
//             {
//                 FoundXRegion = 1; // SO we know it should be in this spacial step
//                 while (FoundYRegion == 0 && yStep <= mNy)
//                 {
//                     if (Centroid[1] - mMinY <= (mMaxY - mMinY) * yStep / mNy) // Look for y bins
//                     {
//                         FoundYRegion = 1; // We now also have the y bin
//                         while (FoundZRegion == 0 && zStep <= mNz) // Look for z bins
//                         {
//                             if (Centroid[2] - mMinZ <= (mMaxZ - mMinZ) * zStep / mNz)
//                             {
//                                 FoundZRegion = 1; // We now also have the y bin

//                                 if (Centroid[2] - mMinZ >= (mMaxZ - mMinZ) * zStep / mNz - FudgeZ && zStep + 1 <= mNz) // Look for x bins
//                                 {
//                                     FudgeInZ = 1;
//                                     // SHould also be in the next X bin across
//                                     c_vector<double,3> Bin= Create_c_vector(xStep, yStep, zStep + 1); //This is the key in the map, then the value is some vector that this is added to
//                                     mBinMap[Bin].push_back(elem_index);
//                                 }
//                                 else if (Centroid[2] - mMinZ <= (mMaxZ - mMinZ) * (zStep - 1) / mNz + FudgeZ && zStep > 1) // // What about the region behind?????
//                                 {
//                                     lowerFudgeZ = 1;
//                                     c_vector<double,3> Bin= Create_c_vector(xStep, yStep, zStep - 1); // SHould also be in the next X bin across
//                                     mBinMap[Bin].push_back(elem_index);
//                                 }
//                             }
//                             else
//                             {
//                                 zStep += 1;
//                             }
//                         }

//                         if (Centroid[1] - mMinY >= (mMaxY - mMinY) * yStep / mNy - FudgeY && yStep + 1 <= mNy) // Look for x bins
//                         {
//                             FudgeInY = 1;
//                             // SHould also be in the next X bin across
//                             c_vector<double,3> Bin= Create_c_vector(xStep, yStep + 1, zStep); //This is the key in the map, then the value is some vector that this is added to
//                             mBinMap[Bin].push_back(elem_index);

//                             if (FudgeInZ)
//                             {
//                                 c_vector<double,3> Bin= Create_c_vector(xStep, yStep + 1, zStep + 1); //This is the key in the map, then the value is some vector that this is added to
//                                 mBinMap[Bin].push_back(elem_index);
//                             }
//                         }
//                         else if (Centroid[1] - mMinY <= (mMaxY - mMinY) * (yStep - 1) / mNy + FudgeY && yStep > 1) // // What about the region behind?????
//                         {
//                             lowerFudgeY = 1;
//                             c_vector<double,3> Bin= Create_c_vector(xStep, yStep - 1, zStep); // SHould also be in the next X bin across
//                             mBinMap[Bin].push_back(elem_index);

//                             if (lowerFudgeZ)
//                             {
//                                 c_vector<double,3> Bin= Create_c_vector(xStep, yStep - 1, zStep - 1); //This is the key in the map, then the value is some vector that this is added to
//                                 mBinMap[Bin].push_back(elem_index);
//                             }
//                         }
//                     }
//                     else
//                     {
//                         yStep += 1;
//                     }
//                 }
//                 // Is this element also in the fudge region of the next bin across?
//                 if (Centroid[0] - mMinX >= (mMaxX - mMinX) * xStep / mNx - FudgeX && xStep + 1 <= mNx) // Look for x bins
//                 {
//                     FudgeInX = 1;
//                     c_vector<double,3> Bin= Create_c_vector(xStep + 1, yStep, zStep); // SHould also be in the next X bin across
//                     mBinMap[Bin].push_back(elem_index);
//                     // Also want to account for the fudge in the other directions :)
//                     if (FudgeInY)
//                     {
//                         c_vector<double,3> Bin= Create_c_vector(xStep + 1, yStep + 1, zStep); //This is the key in the map, then the value is some vector that this is added to
//                         mBinMap[Bin].push_back(elem_index);
//                         if (FudgeInZ)
//                         {
//                             c_vector<double,3> Bin= Create_c_vector(xStep + 1, yStep + 1, zStep + 1); //This is the key in the map, then the value is some vector that this is added to
//                             mBinMap[Bin].push_back(elem_index);
//                         }
//                     }
//                     if (FudgeInZ)
//                     {
//                         c_vector<double,3> Bin= Create_c_vector(xStep + 1, yStep, zStep + 1); //This is the key in the map, then the value is some vector that this is added to
//                         mBinMap[Bin].push_back(elem_index);
//                     }
//                 }
//                 else if (Centroid[0] - mMinX <= (mMaxX - mMinX) * (xStep - 1) / mNx + FudgeX && xStep > 1) // What about the region behind?????
//                 {
//                     lowerFudgeX = 1;
//                     c_vector<double,3> Bin= Create_c_vector(xStep - 1, yStep, zStep); // SHould also be in the next X bin across
//                     mBinMap[Bin].push_back(elem_index);

//                     // Also want to account for the fudge in the other directions :)
//                     if (lowerFudgeY)
//                     {
//                         c_vector<double,3> Bin= Create_c_vector(xStep - 1, yStep - 1, zStep); //This is the key in the map, then the value is some vector that this is added to
//                         mBinMap[Bin].push_back(elem_index);
//                         if (lowerFudgeZ)
//                         {
//                             c_vector<double,3> Bin= Create_c_vector(xStep - 1, yStep - 1, zStep - 1); //This is the key in the map, then the value is some vector that this is added to
//                             mBinMap[Bin].push_back(elem_index);
//                         }
//                     }
//                     if (lowerFudgeZ)
//                     {
//                         c_vector<double,3> Bin =Create_c_vector(xStep - 1, yStep, zStep - 1); //This is the key in the map, then the value is some vector that this is added to
//                         mBinMap[Bin].push_back(elem_index);
//                     }
//                 }
//             }
//             else
//             {
//                 xStep += 1;
//             }
//         }

//         c_vector<double,3> Bin= Create_c_vector(xStep, yStep, zStep); //This is the key in the map, then the value is some vector that this is added to

//         // I think the fugde region needs to be implemented else where -- there are still some problems with my nearest element search too --- this could be causing me pain
//         mBinMap[Bin].push_back(elem_index);
//         // TRACE(" In bin mapping")
//         // PRINT_VARIABLE(elem_index)
//         // assert(elem_index < this->mrMesh->mElements.size());

//         // I want a visulisation so I can easily check if my binning is looking okay
//         for (int i = 0; i < 3; i++)
//         {
//             CellPtr p_cell = this->GetCellUsingLocationIndex(elem_iter->GetNodeGlobalIndex(i));
//             p_cell->GetCellData()->SetItem("BinX", xStep);
//             p_cell->GetCellData()->SetItem("BinY", yStep);
//             p_cell->GetCellData()->SetItem("BinZ", zStep);

//             if (lowerFudgeZ == 1 || lowerFudgeX == 1 || lowerFudgeY == 1)
//             {
//                 p_cell->GetCellData()->SetItem("FudgeL", 1);
//             }
//             else
//             {
//                 p_cell->GetCellData()->SetItem("FudgeL", 0);
//             }
//             if (FudgeInZ == 1 || FudgeInX == 1 || FudgeInY == 1)
//             {
//                 p_cell->GetCellData()->SetItem("Fudge", 1);
//             }
//             else
//             {
//                 p_cell->GetCellData()->SetItem("Fudge", 0);
//             }
//         }
//     }
// }

// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// unsigned HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetClosestElementInOldMeshMethod3(unsigned node_index, c_vector<double, SPACE_DIM> NewNodeLocation)
// {

//     mMapOfProbNodes[node_index] = 5; // -5 meaning I cant find the needed element
//     unsigned ClosestElement = -10;
//     double distance =  std::max( (mMaxX - mMinX), (mMaxY - mMinY))/(2*mNx);

//     // Determine which bin I am in, then look in that bin for the element containing this node
//     c_vector<double,3> Bin = DetermineBin(NewNodeLocation);
// std::vector<double> CurrentBin = mBinMap[Bin];
//     bool HaveElement = 0;
//      for (std::vector<double>::iterator ele_iter = CurrentBin.begin(); ele_iter != CurrentBin.end(); ++ele_iter) //  bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
//     {
//         // Now need to find the cloeset element in this bin
//         c_vector<double, SPACE_DIM> Centroid = mCentroidMap[*ele_iter];
//         if (abs(norm_2(NewNodeLocation - Centroid)) < distance)
//         {
//             distance = abs(norm_2(NewNodeLocation - Centroid));
//             ClosestElement = *ele_iter;
//             if ( PointInTriangle3D(NewNodeLocation, ClosestElement))
//             {
//                 HaveElement =1;
//                 mMapOfProbNodes[node_index] = 1;
//                 break;
//             }
//         }
//     }

//     // unsigned num_elements = this->rGetMesh().GetNumElements();

//     // if (HaveElement == 0)
//     // {
//     //     distance = mMaxEdgelength;
//     //     std::pair<unsigned,unsigned> ClosestEdge;
//     //     c_vector<double, SPACE_DIM> ClosestEdgeNorm ;
//     //     c_vector<double, SPACE_DIM> ClosestMidPoint;
//     //      for (typename MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SpringIterator spring_iterator = this->SpringsBegin();
//     //      spring_iterator != this->SpringsEnd();
//     //      ++spring_iterator)
//     //      {
//     //         // mInitalPositionOfRemeshedNodes
//     //         unsigned NodeIndexA = spring_iterator.GetNodeA()->GetIndex();
//     //         unsigned NodeIndexB = spring_iterator.GetNodeB()->GetIndex();

//     //         c_vector<double, SPACE_DIM> P1 = this->GetNode(NodeIndexA)->rGetLocation();
//     //         c_vector<double, SPACE_DIM> P2 = this->GetNode(NodeIndexB)->rGetLocation();
//     //         c_vector<double, SPACE_DIM> EdgeVector = P2- P1;
//     //         c_vector<double, SPACE_DIM> MidPoint = (P2+ P1)/2;
//     //         EdgeVector/=norm_2(EdgeVector);
//     //         c_vector<double, SPACE_DIM> V = NewNodeLocation - P1;

//     //         // Normal to the plane containing the edge and the new node
//     //         c_vector<double, SPACE_DIM> N = VectorProduct(EdgeVector, V)/norm_2(V);
//     //         N/=norm_2(N);
//     //         // Now can find the normal to the edge that is contaned in the plane defined by the three points

//     //         c_vector<double, SPACE_DIM> normal = VectorProduct(EdgeVector, N);
//     //         normal/=norm_2(normal);
//     //         if (norm_2(MidPoint-normal -NewNodeLocation) < norm_2(MidPoint+normal -NewNodeLocation) )
//     //         {
//     //             normal*=-1;
//     //         }
//     //         // Need to check the direction of the normal :)

//     //         // normal to edge will be normal to the edge and the element normal -- se get the cross product of
//     //         double theta = acos(inner_prod(normal, V/norm_2(V)) );
//     //         double d =  norm_2(V) *cos(theta);
//     //         if (d< distance )
//     //         {
//     //             distance =d;
//     //             // need to think about what to do here,
//     //             ClosestEdge = std::pair<unsigned,unsigned>(NodeIndexA, NodeIndexB);
//     //             ClosestEdgeNorm = normal;
//     //             ClosestMidPoint = MidPoint;
//     //         }

//     //     }

//     //     c_vector<double, SPACE_DIM> P1 = this->GetNode(ClosestEdge.first)->rGetLocation();
//     //     c_vector<double, SPACE_DIM> P2 = this->GetNode(ClosestEdge.second)->rGetLocation();

//     //     //Check that this is close enough // how ??
//     //     if (distance < norm_2(P1 - NewNodeLocation )  || distance < norm_2(P2 - NewNodeLocation )  )
//     //     {
//     //         std::set<unsigned> elements_containing_node1 = this->GetNode(ClosestEdge.first)->rGetContainingElementIndices();
//     //         std::set<unsigned> elements_containing_node3 = this->GetNode(ClosestEdge.second)->rGetContainingElementIndices();

//     //         std::set<unsigned> shared_elements;
//     //         std::set_intersection(elements_containing_node1.begin(),
//     //                                 elements_containing_node1.end(),
//     //                                 elements_containing_node3.begin(),
//     //                                 elements_containing_node3.end(),
//     //                                 std::inserter(shared_elements, shared_elements.begin()));
//     //         typename std::set<unsigned>::iterator CommonElements = shared_elements.begin();
//     //         unsigned Element1 = *CommonElements;
//     //         if (shared_elements.size() ==1)
//     //         {
//     //             // This is a boundary edge
//     //              ClosestElement =Element1;
//     //              HaveElement=1;
//     //         }else
//     //         {
//     //             std::advance(CommonElements, 1);
//     //             unsigned Element2 = *CommonElements;

//     //             if (PointInTriangle3D(NewNodeLocation, Element1))
//     //             {
//     //                 ClosestElement =Element1;
//     //                 HaveElement=1;
//     //                 // std::cout<<" Inside element 1\n";
//     //                 // TRACE(" Inside element 1\n")
//     //             } else if( PointInTriangle3D(NewNodeLocation, Element2) ==1)
//     //             {
//     //                 ClosestElement =Element2;
//     //                 HaveElement=1;

//     //                 // TRACE(" Inside element 2\n")
//     //             }else
//     //             {
//     //                 // TRACE("  Not inside the elements connected to the edge. Dam it\n")
//     //                  // unsigned Node1 = ClosestEdge.first;
//     //             // unsigned Node2 = ClosestEdge.second;

//     //                 c_vector<double, SPACE_DIM> Centroid1 = mCentroidMap[Element1];
//     //                 c_vector<double, SPACE_DIM> Centroid2 = mCentroidMap[Element2];

//     //                 double Distance1 = norm_2(Centroid1 - NewNodeLocation);
//     //                 double Distance2 = norm_2(Centroid2 - NewNodeLocation);
//     //                 if (Distance1 <Distance2)
//     //                 {
//     //                     ClosestElement =Element1;
//     //                 }
//     //                 else{
//     //                     ClosestElement =Element2;

//     //                 }

//     //                 // // FIgure out where it is the same way I used as above in method 2 -- the current way I have doesnt work!!!!!!
//     //                 // c_vector<double, 3> Edge = P2 - P1;
//     //                 // c_vector<double, 3> EdgeMidpoint = (P2 + P1) / 2;

//     //                 // c_vector<double, SPACE_DIM> Normal1 = GetElementNormal(Element1);
//     //                 // c_vector<double, SPACE_DIM> Normal2 = GetElementNormal(Element2);

//     //                 // // Average normal of the two elements
//     //                 // c_vector<double, SPACE_DIM> Normal = (Normal1 + Normal2) / 2;

//     //                 // // The normal should be pointing outwards, so the direction of the plane normal will depend on the direction of the edge vector.
//     //                 // // The edge vector is totally depenedent on the order of the nodes pair, which should be sorted numerically
//     //                 // c_vector<double, SPACE_DIM> PlaneNormal = VectorProduct(Normal, Edge);
//     //                 // PlaneNormal /= norm_2(PlaneNormal);

//     //                 // // Need to determine the direction of the normal

//     //                 // c_vector<double, SPACE_DIM> NormalArm = PlaneNormal + EdgeMidpoint;
//     //                 // double DistanceToA = norm_2(NormalArm - mCentroidMap[Element1]); // A being the center point of element 1
//     //                 // double DistanceToB = norm_2(NormalArm - mCentroidMap[Element2]); // A being the center point of element 2
//     //                 // assert(DistanceToB != DistanceToA);
//     //                 // if (DistanceToA < DistanceToB)
//     //                 // {
//     //                 //     // TRACE("The normal is pointing towards element 1 ")
//     //                 //     // I now have a plane, Dependeing on what side of the plane the point is, will give me the closest element, but I dont know which way the plane is pointing
//     //                 //     double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
//     //                 //     if (side > 0)
//     //                 //     {
//     //                 //         // TRACE("The new node is on the element ONE side of the midline plane ")
//     //                 //         ClosestElement = Element1;
//     //                 //         HaveElement=1;
//     //                 //     }
//     //                 //     else
//     //                 //     {
//     //                 //         // TRACE("The new node is on the element TWO side of the midline plane")
//     //                 //         ClosestElement = Element2;
//     //                 //         HaveElement=1;
//     //                 //     }

//     //                 // }  else if (DistanceToB < DistanceToA)
//     //                 // {
//     //                 // // TRACE("The normal is pointing towards element 2 ")
//     //                 // double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
//     //                 // if (side > 0)
//     //                 // {
//     //                 //     // TRACE("The new node is on the element TWO side of the midline plane")
//     //                 //     ClosestElement = Element2;
//     //                 //     HaveElement=1;
//     //                 // }
//     //                 // else
//     //                 // {
//     //                 //     //   TRACE("The new node is on the element ONE side of the midline plane")
//     //                 //     ClosestElement = Element1;
//     //                 //     HaveElement=1;
//     //                 // }
//     //             }

//     //             // }

//     //         // //////////////////////////////v//////
//     //         //     // Edge has two elements
//     //         //    // Now determine which way the normal is faceing? -- Away from new point or towards?
//     //         //     std::advance(CommonElements, 1);
//     //         //     unsigned Element2 = *CommonElements;
//     //         //     c_vector<double, SPACE_DIM> Centroid1 = mCentroidMap[Element1];
//     //         //     c_vector<double, SPACE_DIM> Centroid2 = mCentroidMap[Element2];

//     //         //      if (norm_2((ClosestMidPoint+ ClosestEdgeNorm) - Centroid1 )< norm_2((ClosestMidPoint+ ClosestEdgeNorm) - Centroid2)  && num_elements<=Element1)
//     //         //     {
//     //         //         // Element 1 is closer//
//     //         //         ClosestElement =Element1;
//     //         //         HaveElement=1;
//     //         //     }else
//     //         //     {
//     //         //         // Element 2 is closer//
//     //         //         ClosestElement =Element2;
//     //         //         HaveElement=1;
//     //         //      }
//     //         }
//     //       mMapOfProbNodes[node_index] =3;
//     //     }
//     //     //So I Have the element this node is closest to.
//     // }
//     // if (HaveElement == 0)
//     // {
//     //      mMapOfProbNodes[node_index] =4;
//     // }
//     assert(ClosestElement !=-10);
//     return ClosestElement;
// }

// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// unsigned HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetClosestElementInOldMeshMethod5(unsigned node_index, c_vector<double, SPACE_DIM> NewNodeLocation)
// {

//     mMapOfProbNodes[node_index] = 5; // -5 meaning I cant find the needed element
//     unsigned ClosestElement = -10;
//     double distance =  std::max( (mMaxX - mMinX), (mMaxY - mMinY))/(2*mNx);

//     // Determine which bin I am in, then look in that bin for the element containing this node
//     c_vector<double,3> Bin = DetermineBin(NewNodeLocation);

// std::vector<double> CurrentBin = mBinMap[Bin];

//     bool HaveElement = 0;
//      for (std::vector<double>::iterator ele_iter = CurrentBin.begin(); ele_iter != CurrentBin.end(); ++ele_iter) //  bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
//     {
//         // Now need to find the cloeset element in this bin
//         c_vector<double, SPACE_DIM> Centroid = mCentroidMap[*ele_iter];
//         if (abs(norm_2(NewNodeLocation - Centroid)) < distance)
//         {
//             distance = abs(norm_2(NewNodeLocation - Centroid));
//             ClosestElement = *ele_iter;
//             if ( PointInTriangle3D(NewNodeLocation, ClosestElement))
//             {
//                 HaveElement =1;
//                 mMapOfProbNodes[node_index] = 1;
//                 break;
//             }
//         }
//     }

//     unsigned num_elements = this->rGetMesh().GetNumElements();

//     if (HaveElement == 0)
//     {
//         distance = mMaxEdgelength;
//         std::pair<unsigned,unsigned> ClosestEdge;
//         c_vector<double, SPACE_DIM> ClosestEdgeNorm ;
//         c_vector<double, SPACE_DIM> ClosestMidPoint;
//          for (typename MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::SpringIterator spring_iterator = this->SpringsBegin();
//          spring_iterator != this->SpringsEnd();
//          ++spring_iterator)
//          {
//             // mInitalPositionOfRemeshedNodes
//             unsigned NodeIndexA = spring_iterator.GetNodeA()->GetIndex();
//             unsigned NodeIndexB = spring_iterator.GetNodeB()->GetIndex();

//             c_vector<double, SPACE_DIM> P1 = this->GetNode(NodeIndexA)->rGetLocation();
//             c_vector<double, SPACE_DIM> P2 = this->GetNode(NodeIndexB)->rGetLocation();
//             c_vector<double, SPACE_DIM> EdgeVector = P2- P1;
//             c_vector<double, SPACE_DIM> MidPoint = (P2+ P1)/2;
//             EdgeVector/=norm_2(EdgeVector);
//             c_vector<double, SPACE_DIM> V = NewNodeLocation - P1;

//             // Normal to the plane containing the edge and the new node
//             c_vector<double, SPACE_DIM> N = VectorProduct(EdgeVector, V)/norm_2(V);
//             N/=norm_2(N);
//             // Now can find the normal to the edge that is contaned in the plane defined by the three points

//             c_vector<double, SPACE_DIM> normal = VectorProduct(EdgeVector, N);
//             normal/=norm_2(normal);
//             if (norm_2(MidPoint-normal -NewNodeLocation) < norm_2(MidPoint+normal -NewNodeLocation) )
//             {
//                 normal*=-1;
//             }
//             // Need to check the direction of the normal :)

//             // normal to edge will be normal to the edge and the element normal -- se get the cross product of
//             double theta = acos(inner_prod(normal, V/norm_2(V)) );
//             double d =  norm_2(V) *cos(theta);
//             if (d< distance )
//             {
//                 distance =d;
//                 // need to think about what to do here,
//                 ClosestEdge = std::pair<unsigned,unsigned>(NodeIndexA, NodeIndexB);
//                 ClosestEdgeNorm = normal;
//                 ClosestMidPoint = MidPoint;
//             }

//         }

//         c_vector<double, SPACE_DIM> P1 = this->GetNode(ClosestEdge.first)->rGetLocation();
//         c_vector<double, SPACE_DIM> P2 = this->GetNode(ClosestEdge.second)->rGetLocation();

//         //Check that this is close enough // how ??  -- distance to edge is smaller than distance to element center
//         if (distance < norm_2(P1 - NewNodeLocation )  || distance < norm_2(P2 - NewNodeLocation )  )
//         {
//             std::set<unsigned> elements_containing_node1 = this->GetNode(ClosestEdge.first)->rGetContainingElementIndices();
//             std::set<unsigned> elements_containing_node3 = this->GetNode(ClosestEdge.second)->rGetContainingElementIndices();

//             std::set<unsigned> shared_elements;
//             std::set_intersection(elements_containing_node1.begin(),
//                                     elements_containing_node1.end(),
//                                     elements_containing_node3.begin(),
//                                     elements_containing_node3.end(),
//                                     std::inserter(shared_elements, shared_elements.begin()));
//             typename std::set<unsigned>::iterator CommonElements = shared_elements.begin();
//             unsigned Element1 = *CommonElements;
//             if (shared_elements.size() ==1)
//             {
//                 // This is a boundary edge
//                  ClosestElement =Element1;
//                  HaveElement=1;
//             }else
//             {
//                 std::advance(CommonElements, 1);
//                 unsigned Element2 = *CommonElements;

//                 if (PointInTriangle3D(NewNodeLocation, Element1)==1)
//                 {
//                     ClosestElement =Element1;
//                     HaveElement=1;
//                     // std::cout<<" Inside element 1\n";
//                     // TRACE(" Inside element 1\n")
//                 } else if( PointInTriangle3D(NewNodeLocation, Element2) ==1)
//                 {
//                     ClosestElement =Element2;
//                     HaveElement=1;
//                     // TRACE(" Inside element 2\n")
//                 }
//                 // The new node is now not in either of the elements, need to find which one it is closest to
//                 else if ( ClosestPointInTriangle(NewNodeLocation, Element1) <= ClosestPointInTriangle(NewNodeLocation, Element2))
//                 {
//                     ClosestElement =Element1;
//                     HaveElement=1;

//                 }else
//                 {
//                     ClosestElement =Element2;
//                     HaveElement=1;

//                 }

//             }
//         }
//     }

//                 // else
//                 // {
//                 //     // TRACE("  Not inside the elements connected to the edge. Dam it\n")
//                 //      // unsigned Node1 = ClosestEdge.first;
//                 // // unsigned Node2 = ClosestEdge.second;

//                 //     c_vector<double, SPACE_DIM> Centroid1 = mCentroidMap[Element1];
//                 //     c_vector<double, SPACE_DIM> Centroid2 = mCentroidMap[Element2];

//                 //     double Distance1 = norm_2(Centroid1 - NewNodeLocation);
//                 //     double Distance2 = norm_2(Centroid2 - NewNodeLocation);
//                 //     if (Distance1 <Distance2)
//                 //     {
//                 //         ClosestElement =Element1;
//                 //     }
//                 //     else{
//                 //         ClosestElement =Element2;

//                 //     }

//                     // // FIgure out where it is the same way I used as above in method 2 -- the current way I have doesnt work!!!!!!
//                     // c_vector<double, 3> Edge = P2 - P1;
//                     // c_vector<double, 3> EdgeMidpoint = (P2 + P1) / 2;

//                     // c_vector<double, SPACE_DIM> Normal1 = GetElementNormal(Element1);
//                     // c_vector<double, SPACE_DIM> Normal2 = GetElementNormal(Element2);

//                     // // Average normal of the two elements
//                     // c_vector<double, SPACE_DIM> Normal = (Normal1 + Normal2) / 2;

//                     // // The normal should be pointing outwards, so the direction of the plane normal will depend on the direction of the edge vector.
//                     // // The edge vector is totally depenedent on the order of the nodes pair, which should be sorted numerically
//                     // c_vector<double, SPACE_DIM> PlaneNormal = VectorProduct(Normal, Edge);
//                     // PlaneNormal /= norm_2(PlaneNormal);

//                     // // Need to determine the direction of the normal

//                     // c_vector<double, SPACE_DIM> NormalArm = PlaneNormal + EdgeMidpoint;
//                     // double DistanceToA = norm_2(NormalArm - mCentroidMap[Element1]); // A being the center point of element 1
//                     // double DistanceToB = norm_2(NormalArm - mCentroidMap[Element2]); // A being the center point of element 2
//                     // assert(DistanceToB != DistanceToA);
//                     // if (DistanceToA < DistanceToB)
//                     // {
//                     //     // TRACE("The normal is pointing towards element 1 ")
//                     //     // I now have a plane, Dependeing on what side of the plane the point is, will give me the closest element, but I dont know which way the plane is pointing
//                     //     double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
//                     //     if (side > 0)
//                     //     {
//                     //         // TRACE("The new node is on the element ONE side of the midline plane ")
//                     //         ClosestElement = Element1;
//                     //         HaveElement=1;
//                     //     }
//                     //     else
//                     //     {
//                     //         // TRACE("The new node is on the element TWO side of the midline plane")
//                     //         ClosestElement = Element2;
//                     //         HaveElement=1;
//                     //     }

//                     // }  else if (DistanceToB < DistanceToA)
//                     // {
//                     // // TRACE("The normal is pointing towards element 2 ")
//                     // double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
//                     // if (side > 0)
//                     // {
//                     //     // TRACE("The new node is on the element TWO side of the midline plane")
//                     //     ClosestElement = Element2;
//                     //     HaveElement=1;
//                     // }
//                     // else
//                     // {
//                     //     //   TRACE("The new node is on the element ONE side of the midline plane")
//                     //     ClosestElement = Element1;
//                     //     HaveElement=1;
//                     // }
//             // }

//                 // }

//             // //////////////////////////////v//////
//             //     // Edge has two elements
//             //    // Now determine which way the normal is faceing? -- Away from new point or towards?
//             //     std::advance(CommonElements, 1);
//             //     unsigned Element2 = *CommonElements;

//             //     c_vector<double, SPACE_DIM> Centroid1 = mCentroidMap[Element1];
//             //     c_vector<double, SPACE_DIM> Centroid2 = mCentroidMap[Element2];

//             //      if (norm_2((ClosestMidPoint+ ClosestEdgeNorm) - Centroid1 )< norm_2((ClosestMidPoint+ ClosestEdgeNorm) - Centroid2)  && num_elements<=Element1)
//             //     {
//             //         // Element 1 is closer//
//             //         ClosestElement =Element1;
//             //         HaveElement=1;

//             //     }else
//             //     {
//             //         // Element 2 is closer//
//             //         ClosestElement =Element2;
//             //         HaveElement=1;
//             //      }
//         //     // }
//         //   mMapOfProbNodes[node_index] =3;
//         // }
//         //So I Have the element this node is closest to.
//     // }

//     if (HaveElement == 0)
//     {
//          mMapOfProbNodes[node_index] =4;
//     }
//     assert(ClosestElement !=-10);
//     return ClosestElement;
// }

// std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > UnitNormals;
// std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> otherNodes;

// std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> edgey = std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*>(this->GetNode(ClosestEdge.first), this->GetNode(ClosestEdge.second));

// bool boundary_edge_found = CalculateElementNormals(edgey, UnitNormals, otherNodes);
// if (boundary_edge_found ==1)
// {
//     // PRINT_VARIABLE(boundary_edge_found)
//     // PRINT_4_VARIABLES(num_elements, Element1, Element2, shared_elements.size())
//     ClosestElement =Element1;
//     HaveElement=1;

// }

// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// unsigned HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetClosestElementInOldMesh(unsigned node_index, c_vector<double, SPACE_DIM> NewNodeLocation)
// {
//     mMapOfProbNodes[node_index] = 5; // -5 meaning I cant find the needed element
//     unsigned ClosestElement = -10;

//     // Determine which bin I am in, then look in that bin for the element containing this node
//     c_vector<double,3> Bin = DetermineBin(NewNodeLocation);
//     //  std::cout << Bin[0] << " " << Bin[1] << " " << Bin[2] << std::endl;
//     unsigned num_elements = this->rGetMesh().GetNumElements();
//     assert(mNx >= Bin[0] && mNy >= Bin[1] && mNz >= Bin[2]);

//     bool HaveElement = 0;

//      std::vector<double> CurrentBin = mBinMap[Bin];
//     for (std::vector<double>::iterator ele_iter = CurrentBin.begin(); ele_iter != CurrentBin.end(); ++ele_iter) //  bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
//     {
//         // Now need to find the element containing this node in this bin
//         bool InTriangle = PointInTriangle3D(NewNodeLocation, *ele_iter);
//         if (InTriangle)
//         {
//             // TRACE("found nearest element by looping over bin and checking if in each element element")
//             ClosestElement = *ele_iter;
//             HaveElement = 1; // If I have the element I can end this search here :) --- I might wrap this search up to its own funciton

//             mMapOfProbNodes[node_index] = 0;
//             break;
//         }
//     }
//     // If I havent found the element by searching the closest bin, then look across whole mesh

//     if (HaveElement == 0)
//     {
//         for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator elem_iter = this->rGetMesh().GetElementIteratorBegin();
//              elem_iter != this->rGetMesh().GetElementIteratorEnd();
//              ++elem_iter)
//         {
//             bool InTriangle = PointInTriangle3D(NewNodeLocation, elem_iter->GetIndex());
//             if (InTriangle)
//             {
//                 ClosestElement = elem_iter->GetIndex(); //TRACE(" Have found the nearest element by looping over every element and checking if it is in the element")

//                 HaveElement = 1;
//                 mMapOfProbNodes[node_index] = 1;
//                 break;
//             }
//         }
//     }

//     if (HaveElement == 0)
//     {

//         double distance = 1e3; //
//         double distance2 = 1e3; //

// std::vector<double> CurrentBin = mBinMap[Bin];
//         unsigned ClosestElement2 = -10;
//         // Now need to loop over this bin to find the nearest element centroid  -- iterate over vector -- getting element centorid indices out -- can go straight to the centroid here
//         for (std::vector<double>::iterator ele_iter = CurrentBin.begin(); ele_iter != CurrentBin.end(); ++ele_iter) //  bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
//         {
//             // Now need to find the cloeset element in this bin
//             c_vector<double, SPACE_DIM> Centroid = mCentroidMap[*ele_iter];
//             if (abs(norm_2(NewNodeLocation - Centroid)) < distance)
//             {
//                 distance = abs(norm_2(NewNodeLocation - Centroid));
//                 ClosestElement = *ele_iter;
//             }
//         }
//         bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
//         if (InTriangle)
//         {
//             mMapOfProbNodes[node_index] = 2;
//         }
//         else if (InTriangle == 0)
//         {
//             // TRACE("Not in triganlge")

//             // Now I have the 'closest element', get the neighbours to this element, and check the node isnt actually in one of them
//             std::set<unsigned> NeighbouringElements = GetNeighbouringElements(ClosestElement);
//             for (typename std::set<unsigned>::iterator neigh_iter = NeighbouringElements.begin();
//                  neigh_iter != NeighbouringElements.end();
//                  ++neigh_iter)
//             {
//                 if (*neigh_iter != ClosestElement)
//                 {
//                     c_vector<double, SPACE_DIM> Centroid = mCentroidMap[*neigh_iter];

//                     if (abs(norm_2(NewNodeLocation - Centroid)) < distance2)
//                     {
//                         distance2 = abs(norm_2(NewNodeLocation - Centroid));
//                         ClosestElement2 = *neigh_iter;
//                     }
//                 }
//             }
//             // TRACE("Not in the second cloeset element either, need to figure out which side of the plane containing the edge")
//             std::vector<unsigned> CommonNodes = GetCommonNodes(ClosestElement, ClosestElement2);
//             assert(CommonNodes.size() < 3);
//             if (CommonNodes.size() == 2)
//             {

//                 typename std::vector<unsigned>::iterator CommonNodeIter = CommonNodes.begin();

//                 Node<SPACE_DIM>* pNode1 = this->rGetMesh().GetNode(*CommonNodeIter);
//                 c_vector<double, SPACE_DIM> P1 = pNode1->rGetLocation();
//                 std::advance(CommonNodeIter, 1);

//                 Node<SPACE_DIM>* pNode2 = this->rGetMesh().GetNode(*CommonNodeIter);
//                 c_vector<double, SPACE_DIM> P2 = pNode2->rGetLocation();

//                 c_vector<double, 3> Edge = P2 - P1;
//                 c_vector<double, 3> EdgeMidpoint = (P2 + P1) / 2;

//                 c_vector<double, SPACE_DIM> Normal1 = GetElementNormal(ClosestElement);
//                 c_vector<double, SPACE_DIM> Normal2 = GetElementNormal(ClosestElement2);

//                 // Average normal of the two elements
//                 c_vector<double, SPACE_DIM> Normal = (Normal1 + Normal2) / 2;

//                 // The normal should be pointing outwards, so the direction of the plane normal will depend on the direction of the edge vector.
//                 // The edge vector is totally depenedent on the order of the nodes pair, which should be sorted numerically
//                 c_vector<double, SPACE_DIM> PlaneNormal = VectorProduct(Normal, Edge);
//                 PlaneNormal /= norm_2(PlaneNormal);

//                 // Need to determine the direction of the normal

//                 c_vector<double, SPACE_DIM> NormalArm = PlaneNormal + EdgeMidpoint;
//                 double DistanceToA = norm_2(NormalArm - mCentroidMap[ClosestElement]); // A being the center point of element 1
//                 double DistanceToB = norm_2(NormalArm - mCentroidMap[ClosestElement2]); // A being the center point of element 2
//                 assert(DistanceToB != DistanceToA);

//                 if (DistanceToA < DistanceToB)
//                 {

//                     // TRACE("The normal is pointing towards element 1 ")
//                     // I now have a plane, Dependeing on what side of the plane the point is, will give me the closest element, but I dont know which way the plane is pointing
//                     double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
//                     if (side > 0)
//                     {
//                         // TRACE("The new node is on the element ONE side of the midline plane ")
//                     }
//                     else
//                     {
//                         // TRACE("The new node is on the element TWO side of the midline plane")
//                         ClosestElement = ClosestElement2;
//                     }

//                 }
//                 else if (DistanceToB < DistanceToA)
//                 {
//                     // TRACE("The normal is pointing towards element 2 ")
//                     double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
//                     if (side > 0)
//                     {
//                         // TRACE("The new node is on the element TWO side of the midline plane")
//                         ClosestElement = ClosestElement2;
//                     }
//                     else
//                     {
//                         //   TRACE("The new node is on the element ONE side of the midline plane")
//                     }

//                 }
//                 mMapOfProbNodes[node_index] = 3;
//                 HaveElement = 1;
//             }
//         }
//     }
//     assert(ClosestElement != -10);
//     if (HaveElement == 0)
//     {
//         // I think this is the cases I cant find
//         mMapOfProbNodes[node_index] = 5; // -5 meaning I cant find the needed element
//         // If I havent been able to find the containing or the closest element, then the defult is the closest
//     }
//     return ClosestElement;
// }

// // Now need to loop over this bin to find the nearest element centroid  -- iterate over vector -- getting element centorid indices out -- can go straight to the centroid here
// for (std::vector<double>::iterator ele_iter = CurrentBin.begin(); ele_iter != mBinMap[Bin].end(); ++ele_iter) //  bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
// {

//     // Now need to find the cloeset element in this bin
//     c_vector<double, SPACE_DIM> Centroid = mCentroidMap[*ele_iter];
//     if (abs(norm_2(NewNodeLocation - Centroid)) < distance)
//     {
//         distance = abs(norm_2(NewNodeLocation - Centroid));
//         ClosestElement = *ele_iter;
//     }
// }
// assert(num_elements >= ClosestElement && ClosestElement != -10);
// bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
// if (InTriangle ==1 )
// {
//     A+=1;
// }
// if (InTriangle ==0 )
// {
//     // TRACE("Not in triganlge")

//     // Now I have the 'closest element', get the neighbours to this element, and check the node isnt actually in one of them
//     std::set<unsigned> NeighbouringElements = GetNeighbouringElements(ClosestElement);
//     for (typename std::set<unsigned>::iterator neigh_iter = NeighbouringElements.begin();
//         neigh_iter != NeighbouringElements.end();
//         ++neigh_iter)
//     {
//         if (*neigh_iter != ClosestElement)
//         {
//             c_vector<double, SPACE_DIM> Centroid = mCentroidMap[*neigh_iter];

//             if (abs(norm_2(NewNodeLocation - Centroid)) < distance2)
//             {
//                 distance2 = abs(norm_2(NewNodeLocation - Centroid));
//                 ClosestElement2 = *neigh_iter;
//             }
//         }
//     }

//     bool InTriangle2 = PointInTriangle3D(NewNodeLocation, ClosestElement);
//     if  (InTriangle2==1)
//     {
//         // TRACE("In the second triangle")// Looks like this never happens?
//         B+=1;
//     }
//     if (InTriangle2==0)
//     {
//         // TRACE("Not in the second cloeset element either, need to figure out which side of the plane containing the edge")
//         std::vector<unsigned> CommonNodes = GetCommonNodes(ClosestElement, ClosestElement2);
//         assert(CommonNodes.size() < 3);
//         if (CommonNodes.size() == 2)
//         {
//             // TRACE("Two common nodes")
//             typename std::vector<unsigned>::iterator CommonNodeIter = CommonNodes.begin();

//             Node<SPACE_DIM>* pNode1 = this->rGetMesh().GetNode(*CommonNodeIter);
//             c_vector<double, SPACE_DIM> P1 = pNode1->rGetLocation();
//             std::advance(CommonNodeIter, 1);

//             Node<SPACE_DIM>* pNode2 = this->rGetMesh().GetNode(*CommonNodeIter);
//             c_vector<double, SPACE_DIM> P2 = pNode2->rGetLocation();

//             c_vector<double, 3> Edge = P2 - P1;
//             c_vector<double, 3> EdgeMidpoint = (P2 + P1) / 2;

//             c_vector<double, SPACE_DIM> Normal1 = GetElementNormal(ClosestElement);
//             c_vector<double, SPACE_DIM> Normal2 = GetElementNormal(ClosestElement2);

//             // Average normal of the two elements
//             c_vector<double, SPACE_DIM> Normal = (Normal1 + Normal2) / 2;

//             // The normal should be pointing outwards, so the direction of the plane normal will depend on the direction of the edge vector.
//             // The edge vector is totally depenedent on the order of the nodes pair, which should be sorted numerically
//             c_vector<double, SPACE_DIM> PlaneNormal = VectorProduct(Normal, Edge);
//             PlaneNormal /= norm_2(PlaneNormal);

//             // Need to determine the direction of the normal

//             c_vector<double, SPACE_DIM> NormalArm = PlaneNormal + EdgeMidpoint;
//             double DistanceToA = norm_2(NormalArm -mCentroidMap[ClosestElement]); // A being the center point of element 1
//             double DistanceToB = norm_2(NormalArm -mCentroidMap[ClosestElement2]); // A being the center point of element 2
//             assert(DistanceToB != DistanceToA);

//             if (DistanceToA < DistanceToB)
//             {

//                 // TRACE("The normal is pointing towards element 1 ")
//                 // I now have a plane, Dependeing on what side of the plane the point is, will give me the closest element, but I dont know which way the plane is pointing
//                 double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
//                 if (side>0)
//                 {
//                     // TRACE("The new node is on the element ONE side of the midline plane ")
//                 }
//                 else
//                 {
//                         // TRACE("The new node is on the element TWO side of the midline plane")
//                         ClosestElement = ClosestElement2;
//                 }
//             C+=1;

//             }else if  (DistanceToB < DistanceToA)
//             {
//                 // TRACE("The normal is pointing towards element 2 ")
//                 double side = inner_prod((NewNodeLocation - EdgeMidpoint), PlaneNormal);
//                 if (side>0)
//                 {
//                     // TRACE("The new node is on the element TWO side of the midline plane")
//                     ClosestElement = ClosestElement2;
//                 }
//                 else
//                 {
//                     //   TRACE("The new node is on the element ONE side of the midline plane")
//                 }
//             C+=1;

//             }

//         } else if (CommonNodes.size() <2)
//         {



    // for (typename std::vector<Element<ELEMENT_DIM, SPACE_DIM>*>::iterator iter =  ElementsToRefine.begin(); iter!= ElementsToRefine.end(); ++iter)
    // {
    //     // PRINT_VARIABLE(iter)
    // //        for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator elem_iter = this->rGetMesh().GetElementIteratorBegin();
    // //      elem_iter != this->rGetMesh().GetElementIteratorEnd();
    // //      ++elem_iter)
    // // {
    //     //    unsigned elem_index = iter->sGetIndex();
    //     // unsigned elem_index = this->rGetMesh().GetIndex(iter);
    //     unsigned temp = iter->GetNodeGlobalIndex(0);
    //         // Node<SPACE_DIM>* pNode0 = this->rGetMesh().GetNode(iter->GetNodeGlobalIndex(0));
    //     Node<SPACE_DIM>* pNode1 = this->rGetMesh().GetNode(iter->GetNodeGlobalIndex(1));
    //     Node<SPACE_DIM>* pNode2 = this->rGetMesh().GetNode(iter->GetNodeGlobalIndex(2));
    // c_vector<double, 3> Centroid = (pNode0->rGetLocation() + pNode1->rGetLocation() + pNode2->rGetLocation()) / 3;

    // Add this as a new node
    //     Element<ELEMENT_DIM, SPACE_DIM>* p_element = this->rGetMesh().GetElement(elem_index);

    // }

    // }


// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// unsigned HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetClosestElementInOldMeshMethod4(unsigned node_index, c_vector<double, SPACE_DIM> NewNodeLocation)
// {

//     unsigned ClosestElement = -10;
//     double distance =  std::max( (mMaxX - mMinX), (mMaxY - mMinY))/(2*mNx);

//     // Determine which bin I am in, then look in that bin for the element containing this node
//     c_vector<double,3> Bin = DetermineBin(NewNodeLocation);

//      for (std::vector<double>::iterator ele_iter = mBinMap[Bin].begin(); ele_iter != mBinMap[Bin].end(); ++ele_iter) //  bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
//     {
//         // Now need to find the cloeset element in this bin
//         c_vector<double, SPACE_DIM> Centroid = mCentroidMap[*ele_iter];
//         if (abs(norm_2(NewNodeLocation - Centroid)) < distance)
//         {
//             distance = abs(norm_2(NewNodeLocation - Centroid));
//             ClosestElement = *ele_iter;
//             if ( PointInTriangle3D(NewNodeLocation, ClosestElement))
//             {

//                 mMapOfProbNodes[node_index] = 1;
//                 break;
//             }
//         }
//     }

//     assert(ClosestElement !=-10);
//     mMapOfProbNodes[node_index] = 3;
//     return ClosestElement;mMinX
// }



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetClosestElementInOldMeshMethod6(unsigned node_index, c_vector<double, SPACE_DIM> NewNodeLocation)
{

    unsigned num_elements = this->rGetMesh().GetNumElements();
    mMapOfProbNodes[node_index] = 5; // -5 meaning I cant find the needed element
    double ClosestElement = -10;
    double distance = 1000 * std::max((mMaxX - mMinX), (mMaxY - mMinY)) / (2 * mNx);

    // Determine which bin I am in, then look in that bin for the element containing this node
    std::vector<int> Bin2 = DetermineBin(NewNodeLocation);
    PRINT_VECTOR(Bin2)
    std::vector<int> Bin = { 1, 1, 1 };

    //

    bool HaveElement = 0;
    for (std::vector<unsigned>::iterator ele_iter = mBinMap[Bin].begin(); ele_iter != mBinMap[Bin].end(); ++ele_iter) //  bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
    {
        // Now need to find the cloeset element in this bin
        c_vector<double, SPACE_DIM> Centroid = mCentroidMap[*ele_iter];
        if (abs(norm_2(NewNodeLocation - Centroid)) < distance)
        {
            distance = abs(norm_2(NewNodeLocation - Centroid));
            ClosestElement = *ele_iter;
        }
    }
    // if ( PointInTriangle3D(NewNodeLocation, ClosestElement))
    // {
    //     HaveElement =1;
    //     mMapOfProbNodes[node_index] = 1;
    // }
    if (HaveElement == 0)
    {
        typename std::vector<int>::iterator iter = Bin.begin();
        int XBin = *iter;
        std::advance(iter, 1);
        int YBin = *iter;
        std::advance(iter, 1);
        int ZBin = *iter;

        std::vector<int> NeighbouringBins;
        for (int i = -1; i <= 1; ++i)
        {
            if (XBin + i > mNx || XBin + i < 1)
            {
                continue;
            }
            if (HaveElement == 1)
            {
                break;
            }
            for (int j = -1; j <= 1; ++j)
            {
                if (YBin + j > mNy || YBin + j < 1 || HaveElement == 1)
                {
                    continue;
                }
                for (int k = -1; k <= 1; ++k)
                {
                    if (ZBin + k > mNz || ZBin + k < 1 || (k == 0 && j == 0 && i == 0) || HaveElement == 1)
                    {
                        continue;
                    }
                    //
                    std::vector<int> NeighBinAddress = { XBin + i, YBin + j, ZBin + k };
                    std::vector<unsigned> NeighBin = mBinMap[NeighBinAddress];
                    // NeighbouringBins.insert(NeighbouringBins.end(), NeighBin.begin(), NeighBin.end());
                    // PRINT_3_VARIABLES(ZBin+k , YBin+j, XBin+i)

                    if (NeighBin.size() < 1)
                    {
                        continue;
                    }
                    for (std::vector<unsigned>::iterator ele_iter = NeighBin.begin(); ele_iter != NeighBin.end(); ++ele_iter) //  bool InTriangle = PointInTriangle3D(NewNodeLocation, ClosestElement);
                    {
                        // Now need to find the cloeset element in this bin
                        c_vector<double, SPACE_DIM> Centroid = mCentroidMap[*ele_iter];
                        if (abs(norm_2(NewNodeLocation - Centroid)) < distance)
                        {
                            distance = abs(norm_2(NewNodeLocation - Centroid));
                            ClosestElement = *ele_iter;
                            if (PointInTriangle3D(NewNodeLocation, ClosestElement))
                            {
                                HaveElement = 1;
                                mMapOfProbNodes[node_index] = 1;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    assert(ClosestElement != -10 && ClosestElement <= num_elements);
    return ClosestElement;
    //
}


 unsigned GetClosestElementInOldMesh(unsigned node_index, c_vector<double, SPACE_DIM> NewNodeLocation);
    unsigned GetClosestElementInOldMeshMethod3(unsigned node_index, c_vector<double, SPACE_DIM> NewNodeLocation);
    unsigned GetClosestElementInOldMeshMethod4(unsigned node_index, c_vector<double, SPACE_DIM> NewNodeLocation);

    unsigned GetClosestElementInOldMeshMethod5(unsigned node_index, c_vector<double, SPACE_DIM> NewNodeLocation);
    unsigned GetClosestElementInOldMeshMethod6(unsigned node_index, c_vector<double, SPACE_DIM> NewNodeLocation);