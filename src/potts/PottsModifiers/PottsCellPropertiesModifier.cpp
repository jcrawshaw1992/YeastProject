

#include "PottsCellPropertiesModifier.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"
#include "WrappedPottsBasedCellPopulation.hpp"

#include <cxxtest/TestSuite.h>

template <unsigned DIM>
PottsCellPropertiesModifier<DIM>::PottsCellPropertiesModifier()
        : AbstractCellBasedSimulationModifier<DIM>()
{
    // friend class AbstractCellPopulation;
}

template <unsigned DIM>
PottsCellPropertiesModifier<DIM>::~PottsCellPropertiesModifier()
{
}

template <unsigned DIM>
void PottsCellPropertiesModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template <unsigned DIM>
void PottsCellPropertiesModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */

    WrappedPottsBasedCellPopulation<DIM>* PottsPopulation = static_cast<WrappedPottsBasedCellPopulation<DIM>*>(&(rCellPopulation)); // Jess needs to edit to be  PottsArbitrarySurfaceIn3DMesh<DIM>
    PottsArbitrarySurfaceIn3DMesh<DIM>* p_static_cast_potts_mesh = static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*>(&(rCellPopulation.rGetMesh())); // Jess needs to edit to be  PottsArbitrarySurfaceIn3DMesh<DIM>

    // UpdateCellData(rCellPopulation);

    /*
          Need to record the store the data as cell data so the writer can take it and write it out
        */

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned Cell_Index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        double Perimeter = p_static_cast_potts_mesh->GetPerimeterOfElement(Cell_Index);
        double Area = p_static_cast_potts_mesh->GetVolumeOfElement(Cell_Index);

        cell_iter->GetCellData()->SetItem("Area", Area);
        cell_iter->GetCellData()->SetItem("Perimeter", Perimeter);

        PottsElement<DIM>* p_element = PottsPopulation->GetElement(Cell_Index);

        // Get center point -- this will let me track the center point
        c_vector<double, DIM> CenterPoint;
        for (unsigned node_index = 0; node_index < p_element->GetNumNodes(); node_index++) // for some reason this needs to go backwards
        {
            CenterPoint += p_element->GetNodeLocation(node_index);
        }
        CenterPoint /= p_element->GetNumNodes();

        cell_iter->GetCellData()->SetItem("CenterX", CenterPoint[0]);
        cell_iter->GetCellData()->SetItem("CenterY", CenterPoint[1]);
        cell_iter->GetCellData()->SetItem("CenterZ", CenterPoint[2]);


        double AspectRatio = p_static_cast_potts_mesh->GetAspectRatio(Cell_Index);
        c_vector<double, 2> MajorAxis = p_static_cast_potts_mesh->GetMajorAxisVector(Cell_Index);
        double MajorAxisAngle = atan(MajorAxis[1]/MajorAxis[0]);
        MajorAxis/=norm_2( MajorAxis);

        c_vector<double, 2> FluidDirection = Create_c_vector(0,1);

        double Orientation = inner_prod(MajorAxis,FluidDirection);


        cell_iter->GetCellData()->SetItem("AspectRatio", AspectRatio);
        cell_iter->GetCellData()->SetItem("Orientation", Orientation);
        cell_iter->GetCellData()->SetItem("MajorAxisAngle", MajorAxisAngle);



    }
}

template <unsigned DIM>
void PottsCellPropertiesModifier<DIM>::SetMeshDimensions(double N_D, double N_Z, double Width, double Length)
{
    mN_D = N_D;
    mN_Z = N_Z;
    mWidth = Width, mLength = Length;
}

template <unsigned DIM>
void PottsCellPropertiesModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{

    /* 
    Periodic Potts population 
    If the Potts population is on a periodic domain, I need to reshuffle the center line 

    Not Perioidic Potts
    Only need to get the stats for the population to save for a writer :) 
    */

    assert(DIM == 3);
    PottsArbitrarySurfaceIn3DMesh<DIM>* p_static_cast_potts_mesh = static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*>(&(rCellPopulation.rGetMesh())); // Jess needs to edit to be  PottsArbitrarySurfaceIn3DMesh<DIM>
    WrappedPottsBasedCellPopulation<DIM>* PottsPopulation = static_cast<WrappedPottsBasedCellPopulation<DIM>*>(&(rCellPopulation)); // Jess needs to edit to be  PottsArbitrarySurfaceIn3DMesh<DIM>
    bool Periodic = PottsPopulation->IsPottsSimulationPeriodic();
    if (PottsPopulation->IsPottsSimulationPeriodic())
    {
        // Periodic potts population -- need to move center lines

        std::vector<c_vector<unsigned, 2> > ElementPairing = PottsPopulation->GetElementPairingVector();
        //
        double G =15;
        double H =7;
        double UnitStep = mWidth / (2*(mN_D - 1)); // *1/2 because the top of a triangle is 1/2 between the two bottom nodes
        for (std::vector<c_vector<unsigned, 2> >::iterator it = ElementPairing.begin(); it != ElementPairing.end(); it++)
        {
            //  Get both elements stored in this vector of the vector

            unsigned Element_1 = (*it)[0];
            unsigned Element_2 = (*it)[1];

            PottsElement<DIM>* p_element_1 = PottsPopulation->GetElement(Element_1);
            PottsElement<DIM>* p_element_2 = PottsPopulation->GetElement(Element_2);

            CellPtr p_cell_1 = PottsPopulation->GetCellUsingLocationIndex(Element_1);
            CellPtr p_cell_2 = PottsPopulation->GetCellUsingLocationIndex(Element_2);
            double XCenter = p_cell_2->GetCellData()->GetItem("Center");
            assert(p_cell_2->GetCellData()->GetItem("Center") == p_cell_1->GetCellData()->GetItem("Center"));

            double TotalNodesInCell = p_element_1->GetNumNodes() + p_element_2->GetNumNodes();
            double NodesInLeftCell = p_element_1->GetNumNodes();
            double NodesInRightCell = p_element_2->GetNumNodes();
            if(TotalNodesInCell < 3)
            {
                continue;
            }
        
            if (NodesInLeftCell > 0.51* TotalNodesInCell || NodesInRightCell > 0.51* TotalNodesInCell)
            {
                double MoveDirection = 0;
                // PRINT_2_VARIABLES(NodesInLeftCell, NodesInRightCell);
                if (NodesInLeftCell > 0.51* TotalNodesInCell)
                {
                    // There are too mamy nodes in the left cell, need to move the midline left
                    // TRACE("Moving to the left")
                    MoveDirection = -1;
                }
                else if (NodesInRightCell  > 0.51* TotalNodesInCell)
                {
                    // There are too mamy nodes in the right cell, need to move the midline right
                    // TRACE("Moving to the right")
                    MoveDirection = 1;
                }
                double LatticeSitesOnLeft = NodesInLeftCell;  double LatticeSitesOnRight = NodesInRightCell;
                // Need to move the center line one unit
                XCenter += MoveDirection * UnitStep;

                // Check if I have fallen off either edge
                if (XCenter < 0)
                { // TRACE("FALLEN OF THE LEFT EDGE");//  std::cout<<" "<<endl;
                    XCenter = mWidth-1.1*UnitStep;
                }
                else if (XCenter > mWidth)
                {
                    XCenter =1.1*UnitStep;
                    // TRACE("FALLEN OF THE RIGHT EDGE") ;// std::cout<<" "<<endl;
                }   
                 for (unsigned node_index = 0; node_index < p_element_1->GetNumNodes(); node_index++)
                {
                    c_vector<double, DIM> LatticeLocation = p_element_1->GetNodeLocation(node_index);
                    // If lattice is within certain distance of center and on the right
                    if (XCenter <= LatticeLocation[0] && LatticeLocation[0] < XCenter +G* UnitStep)
                    {
                        LatticeSitesOnRight +=1;
                        LatticeSitesOnLeft -=1;  
                    }
                    else if( XCenter > 0.8* mWidth && LatticeLocation[0] < H* UnitStep)
                        {// Was 0.8 
                            LatticeSitesOnRight +=1;
                            LatticeSitesOnLeft -=1;
                        }
                }
                // Right cell 
                for (unsigned node_index = 0; node_index < p_element_2->GetNumNodes(); node_index++) 
                {
                    c_vector<double, DIM> LatticeLocation = p_element_2->GetNodeLocation(node_index);
                    if (XCenter > LatticeLocation[0]  && LatticeLocation[0] > XCenter - G * UnitStep)
                    {
                        // TRACE("Removing a cell from the right, and adding it to the left")
                        LatticeSitesOnLeft +=1; 
                        LatticeSitesOnRight -=1; 
                    } 
                    else if (XCenter > 0.2* mWidth && LatticeLocation[0] > XCenter  + H * UnitStep)
                    {
                        // TRACE(" Lattice is in right cell, but is on the left side of the centerlines (periodically)")
                        LatticeSitesOnLeft +=1; 
                        LatticeSitesOnRight -=1;  
                    }
                }
                assert(NodesInLeftCell + NodesInRightCell == LatticeSitesOnLeft + LatticeSitesOnRight);
                if (LatticeSitesOnRight ==0 || LatticeSitesOnLeft==0)
                {continue;
                }
                assert(LatticeSitesOnRight > 0); assert(LatticeSitesOnLeft > 0);

                while( ( (LatticeSitesOnRight > 0.51* TotalNodesInCell ) && MoveDirection==1 ) || ( (LatticeSitesOnLeft > 0.51* TotalNodesInCell) && MoveDirection==-1  ) ) 
                {
                    // Now check where things are on each side
                    LatticeSitesOnLeft = NodesInLeftCell; LatticeSitesOnRight = NodesInRightCell;

                    //TRACE("Midline is still to out of place")
                    XCenter += MoveDirection * UnitStep;
                    // Check if I have fallen off either edge
                    if (XCenter < 0)
                    { XCenter = mWidth - 1.1*UnitStep;
                    }
                    else if (XCenter > mWidth)
                    {XCenter =1.1*UnitStep;
                    }

                    // Left cell 
                    LatticeSitesOnLeft = NodesInLeftCell;  LatticeSitesOnRight = NodesInRightCell;
                    for (unsigned node_index = 0; node_index < p_element_1->GetNumNodes(); node_index++)
                    {
                        c_vector<double, DIM> LatticeLocation = p_element_1->GetNodeLocation(node_index);
                        // If lattice is within certain distance of center and on the right
                        if (XCenter <= LatticeLocation[0] && LatticeLocation[0] < XCenter + G* UnitStep)
                        {
                            LatticeSitesOnRight +=1;
                            LatticeSitesOnLeft -=1;  
                        }
                        else if( XCenter > 0.8* mWidth && LatticeLocation[0] < H* UnitStep)
                        {// Was 0.8 
                            LatticeSitesOnRight +=1;
                            LatticeSitesOnLeft -=1;
                        }
                    }
                    // right cell
                    for (unsigned node_index = 0; node_index < p_element_2->GetNumNodes(); node_index++) 
                    {
                        c_vector<double, DIM> LatticeLocation = p_element_2->GetNodeLocation(node_index);
                        if (XCenter > LatticeLocation[0]  && LatticeLocation[0] > XCenter - G* UnitStep)
                        {
                            // TRACE("Removing a cell from the right, and adding it to the left")
                            LatticeSitesOnLeft +=1; 
                            LatticeSitesOnRight -=1; 
                        }    //( XCenter < 0.2* mWidth && LatticeLocation[0] > mWidth - 3* UnitStep)// && LatticeLocation[0] > XCenter  + 6 * UnitStep) 
                        else if (XCenter < 0.2* mWidth && LatticeLocation[0] > mWidth - H* UnitStep)
                        {
                            // TRACE(" Lattice is in right cell, but is on the left side of the centerlines (periodically)")
                            LatticeSitesOnLeft +=1; 
                            LatticeSitesOnRight -=1;  
                        } 
                    }
                    assert(NodesInLeftCell + NodesInRightCell == LatticeSitesOnLeft + LatticeSitesOnRight);
                    if (LatticeSitesOnRight <=0 || LatticeSitesOnLeft<=0)
                    {
                        // PRINT_2_VARIABLES(LatticeSitesOnRight,LatticeSitesOnLeft)
                        XCenter = p_cell_2->GetCellData()->GetItem("Center") + MoveDirection * UnitStep;
                        break;
                    }
                    assert(LatticeSitesOnRight > 0 && LatticeSitesOnLeft > 0);                
                }
                // TRACE("Finished shuffeling midline")
                // Now I  need to loop over each of the nodes in each and if they are on the wrong side of the line they need to be switched --- I will need to remember which is left and right buhhh
                // Assumed i is left and i+1 is right
                // I want to fist add nodes to elements, then take the relevant ones away 

                // Cell om the left getting added to right


              //-------------------------
                double NodesToRemoveOnLeft = 0 ;
                double NodesToRemoveOnRight = 0 ;
                    for (unsigned node_index = 0; node_index < p_element_1->GetNumNodes(); node_index++)
                    {
                        c_vector<double, DIM> LatticeLocation = p_element_1->GetNodeLocation(node_index);
                        // If lattice is within certain distance of center and on the right
                        if (XCenter <= LatticeLocation[0] && LatticeLocation[0] < XCenter + G* UnitStep)
                        {
                            NodesToRemoveOnLeft+=1;
                        }
                        else if( XCenter > 0.8* mWidth && LatticeLocation[0] < H* UnitStep)
                        { 
                            NodesToRemoveOnLeft+=1;
                        }
                    }
                    for (unsigned node_index = 0; node_index < p_element_2->GetNumNodes(); node_index++) 
                    {
                        c_vector<double, DIM> LatticeLocation = p_element_2->GetNodeLocation(node_index);
                        if (XCenter > LatticeLocation[0]  && LatticeLocation[0] > XCenter - G * UnitStep)
                        { 
                           NodesToRemoveOnRight+=1;
                        }
                        else if (XCenter < 0.2* mWidth && LatticeLocation[0]> mWidth - H* UnitStep)
                        {
                            NodesToRemoveOnRight+=1;
                        }  
                    }

                // PRINT_2_VARIABLES(NodesInLeftCell, NodesInRightCell);
                // PRINT_2_VARIABLES(NodesToRemoveOnLeft, NodesToRemoveOnRight)


                //--------------------
                unsigned Offset;

                if (NodesToRemoveOnLeft!=NodesInLeftCell && NodesToRemoveOnLeft !=0 )
                {
                     for (unsigned node_index = 0; node_index < p_element_1->GetNumNodes(); node_index++)
                    {
                        c_vector<double, DIM> LatticeLocation = p_element_1->GetNodeLocation(node_index);
                        // If lattice is within certain distance of center and on the right
                        if (XCenter <= LatticeLocation[0] && LatticeLocation[0] < XCenter + G* UnitStep)
                        {
                            p_element_2->AddNode(p_element_1->GetNode(node_index)); 
                        } 
                        else if( XCenter > 0.8* mWidth && LatticeLocation[0] < H* UnitStep)
                        {// Was 0.8 
                            // TRACE(" Lattice is in left cell, but is on the right side of the centerlines (periodically)")
                            p_element_2->AddNode(p_element_1->GetNode(node_index));// record this node/lattice site by local identifier to delete
                        }
                    }
                }
               
               if (NodesToRemoveOnRight!= NodesInRightCell && NodesToRemoveOnRight !=0 )
                {
                    // Cell on the right getting added to left
                    for (unsigned node_index = 0; node_index < p_element_2->GetNumNodes(); node_index++) 
                    {
                        c_vector<double, DIM> LatticeLocation = p_element_2->GetNodeLocation(node_index);
                        if (XCenter > LatticeLocation[0]  && LatticeLocation[0] > XCenter - G * UnitStep)
                        {  // TRACE("Removing a cell from the right, and adding it to the left")
                            p_element_1->AddNode(p_element_2->GetNode(node_index));
                        }
                        else if (XCenter < 0.2* mWidth && LatticeLocation[0] > mWidth - H* UnitStep)
                        {
                            // TRACE(" Lattice is in right cell, but is on the left side of the centerlines (periodically)")
                            p_element_1->AddNode(p_element_2->GetNode(node_index)); 
                        }  
                    }
                }

                // Cell on the left getting added to right
                
                if ( NodesToRemoveOnLeft!=NodesInLeftCell && NodesToRemoveOnLeft !=0 )
                {
                    std::vector<unsigned> NodesToRemoveFromElement1;
                    for (unsigned node_index = 0; node_index < p_element_1->GetNumNodes(); node_index++)
                    {
                        c_vector<double, DIM> LatticeLocation = p_element_1->GetNodeLocation(node_index);
                        // If lattice is within certain distance of center and on the right
                        if (XCenter <= LatticeLocation[0] && LatticeLocation[0] < XCenter + G* UnitStep)
                        {
                            NodesToRemoveFromElement1.push_back(node_index); // record this node/lattice site by local identifier to delete
                        }
                        else if( XCenter > 0.8* mWidth && LatticeLocation[0] < H* UnitStep)
                        { 
                            NodesToRemoveFromElement1.push_back(node_index); // record this node/lattice site by local identifier to delete
                        }
                    }

                    Offset=0;
                    assert(p_element_1->GetNumNodes() >0 && p_element_2->GetNumNodes() >0);
                    if (p_element_1->GetNumNodes() > NodesToRemoveFromElement1.size())
                    {
                        for (std::vector<unsigned>::iterator iter = NodesToRemoveFromElement1.begin(); iter != NodesToRemoveFromElement1.end(); iter++)
                        {    
                            // TRACE("Deleting On Left")
                            // PRINT_2_VARIABLES(p_element_1->GetNumNodes(), *iter - Offset)
                            p_element_1->DeleteNode(*iter - Offset);
                            Offset += 1;
                        }
                    }

                }

                // Cell on the right getting added to left

                if (NodesToRemoveOnRight!= NodesInRightCell && NodesToRemoveOnRight !=0 )
                {
                    std::vector<unsigned> NodesToRemoveFromElement2;
                    for (unsigned node_index = 0; node_index < p_element_2->GetNumNodes(); node_index++) 
                    {
                        c_vector<double, DIM> LatticeLocation = p_element_2->GetNodeLocation(node_index);
                        if (XCenter > LatticeLocation[0]  && LatticeLocation[0] > XCenter - G * UnitStep)
                        {  // TRACE("Removing a cell from the right, and adding it to the left")
                            NodesToRemoveFromElement2.push_back(node_index);
                        }
                        else if (XCenter < 0.2* mWidth && LatticeLocation[0]> mWidth - H* UnitStep)
                        {
                            // TRACE(" Lattice is in right cell, but is on the left side of the centerlines (periodically)")
                            NodesToRemoveFromElement2.push_back(node_index);
                        }  
                    }

                    Offset = 0;
                    if (p_element_2->GetNumNodes() > NodesToRemoveFromElement2.size())
                    {

                        for (std::vector<unsigned>::iterator iter = NodesToRemoveFromElement2.begin(); iter != NodesToRemoveFromElement2.end(); iter++)
                        { 
                            // TRACE("Deleting On Right")
                            // PRINT_3_VARIABLES(p_element_2->GetNumNodes(), *iter,*iter - Offset)
                            p_element_2->DeleteNode(*iter - Offset);
                            Offset += 1;
                        }
                    }

                }
                p_cell_1->GetCellData()->SetItem("Center", XCenter);
                p_cell_2->GetCellData()->SetItem("Center", XCenter);
            }
        }
        
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned Cell_Index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            unsigned Current_Sister = PottsPopulation->GetSister(Cell_Index);
        
            PottsElement<DIM>* p_element = PottsPopulation->GetElement(Cell_Index);
            PottsElement<DIM>* p_Sister_element = PottsPopulation->GetElement(Current_Sister);

            // Get center point -- this will let me track the center point
            c_vector<double, DIM> CenterPoint;
            for (unsigned node_index = p_element->GetNumNodes() - 1; node_index < -1; node_index--) // for some reason this needs to go backwards
            {
                CenterPoint += p_element->GetNodeLocation(node_index);
            }
            for (unsigned node_index = p_Sister_element->GetNumNodes() - 1; node_index < -1; node_index--) // for some reason this needs to go backwards
            {
                CenterPoint += p_Sister_element->GetNodeLocation(node_index);
            }
            CenterPoint /= (p_element->GetNumNodes() + p_Sister_element->GetNumNodes() );

            cell_iter->GetCellData()->SetItem("CenterX", CenterPoint[0]);
            cell_iter->GetCellData()->SetItem("CenterY", CenterPoint[1]);
            cell_iter->GetCellData()->SetItem("CenterZ", CenterPoint[2]);
            cell_iter->GetCellData()->SetItem("LatticeCount", p_element->GetNumNodes());
            double Xcenter = cell_iter->GetCellData()->GetItem("Center");

            double Perimeter = p_static_cast_potts_mesh->GetPerimeterOfCoupledElements(Cell_Index, Current_Sister, Xcenter);
            double Area = p_static_cast_potts_mesh->GetVolumeOfElement(Cell_Index) + p_static_cast_potts_mesh->GetVolumeOfElement(Current_Sister);
            double SS = norm_2(p_static_cast_potts_mesh->GetTractionOnElement(Cell_Index)+p_static_cast_potts_mesh->GetTractionOnElement(Current_Sister) )/2 ;

            cell_iter->GetCellData()->SetItem("Area", Area);
            cell_iter->GetCellData()->SetItem("Perimeter", Perimeter);
            cell_iter->GetCellData()->SetItem("ShearStress", SS);

            double AspectRatio = p_static_cast_potts_mesh->GetAspectRatio(Cell_Index, Current_Sister);
             c_vector<double, 2> MajorAxis = p_static_cast_potts_mesh->GetMajorAxisVector(Cell_Index);
             double MajorAxisAngle = atan(MajorAxis[1]/MajorAxis[0]);
              MajorAxis/=norm_2( MajorAxis);
              
             c_vector<double, 2> FluidDirection = Create_c_vector(0,1);

             double Orientation = inner_prod(MajorAxis,FluidDirection);
            
            cell_iter->GetCellData()->SetItem("AspectRatio", AspectRatio);
            cell_iter->GetCellData()->SetItem("Orientation", Orientation);
            cell_iter->GetCellData()->SetItem("MajorAxisAngle", MajorAxisAngle);
        }
    }
    else
    { // Not periodic

        // I need the center point of each cell, and the area and perimeter

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned Cell_Index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            double Perimeter = p_static_cast_potts_mesh->GetPerimeterOfElement(Cell_Index);
            double Area = p_static_cast_potts_mesh->GetVolumeOfElement(Cell_Index);
            double SS = norm_2(p_static_cast_potts_mesh->GetTractionOnElement(Cell_Index));

            cell_iter->GetCellData()->SetItem("Area", Area);
            cell_iter->GetCellData()->SetItem("Perimeter", Perimeter);
            cell_iter->GetCellData()->SetItem("ShearStress", SS);

            PottsElement<DIM>* p_element = PottsPopulation->GetElement(Cell_Index);

            // Get center point -- this will let me track the center point
            c_vector<double, DIM> CenterPoint;
            for (unsigned node_index = p_element->GetNumNodes() - 1; node_index < -1; node_index--) // for some reason this needs to go backwards
            {
                CenterPoint += p_element->GetNodeLocation(node_index);
            }
            CenterPoint /= p_element->GetNumNodes();

            cell_iter->GetCellData()->SetItem("CenterX", CenterPoint[0]);
            cell_iter->GetCellData()->SetItem("CenterY", CenterPoint[1]);
            cell_iter->GetCellData()->SetItem("CenterZ", CenterPoint[2]);


             double AspectRatio = p_static_cast_potts_mesh->GetAspectRatio(Cell_Index);
             c_vector<double, 2> MajorAxis = p_static_cast_potts_mesh->GetMajorAxisVector(Cell_Index);
              MajorAxis/=norm_2( MajorAxis);
             c_vector<double, 2> FluidDirection = Create_c_vector(0,1);

             double Orientation = inner_prod(MajorAxis,FluidDirection);

            cell_iter->GetCellData()->SetItem("AspectRatio", AspectRatio);
            cell_iter->GetCellData()->SetItem("Orientation", Orientation);

        }
    }

    /*
            Now need to record the store the stuff I want to output in cell data so the writer can take it and use it 
        */

    // Now need to loop over the map

    // I need to loop over all the cell pairings -- for all the ones that are not where they should  be, i need to adjhust
    //     // p_static_cast_potts_mesh->CalculateCurvature();
    //     // p_static_cast_potts_mesh->CalculateTraction();
    //   for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
    //          cell_iter != rCellPopulation.End();
    //          ++cell_iter)
    //     {
    //         unsigned Cell_Index= rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

    // 		double ApsectRatio = p_static_cast_potts_mesh->GetAspectRatio(Cell_Index);
    // 		double MajorAxisAngle = p_static_cast_potts_mesh->GetMajorAxisAngle(Cell_Index);
    //         double WallShearStress = norm_2(p_static_cast_potts_mesh->GetTractionOnElement(Cell_Index));

    //         cell_iter->GetCellData()->SetItem("AspectRatio", ApsectRatio);
    // 		cell_iter->GetCellData()->SetItem("MajorAxisAngle", MajorAxisAngle);

    // 		cell_iter->GetCellData()->SetItem("WallShearStress", WallShearStress);

    //     }
    // TRACE("This code works with the non-periodic version, but not the periodic version")
}

template <unsigned DIM>
void PottsCellPropertiesModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class PottsCellPropertiesModifier<1>;
template class PottsCellPropertiesModifier<2>;
template class PottsCellPropertiesModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PottsCellPropertiesModifier)












//
//         // if (MoveDirection != 0)
//         // {
//             double success = 0;  double NumberOfTrys = 1;
//             // while (success == 0)
//             // {
//                 // if (NumberOfTrys > 6)
//                 // {
//                 //     if (NumberOfTrys > 7)
//                 //     {
//                 //         success = 1;
//                 //         TRACE("Midline not moving right, Number of trys")
//                 //         PRINT_VARIABLE(NumberOfTrys)
//                 //     }

//                 // }

//                 // Need to move the center line one unit
//                 XCenter += MoveDirection * UnitStep;

//                 // if (XCenter < 0)
//                 // {
//                 //     TRACE("WE HAVE FALLEND OF THE LEFT EDGE")
//                 //     XCenter = mWidth/2;
//                 // }
//                 // else if (XCenter > mWidth)
//                 // {
//                 //     XCenter = mWidth/2;
//                 //     TRACE("WE HAVE FALLEND OF THE RIGHT EDGE")
//                 // }

//                 // // Now check where things are on each side
//                 // double LatticeSitesOnLeft = 0;
//                 // double LatticeSitesOnRight = 0;
//                 // for (unsigned node_index = 0; node_index < p_element_1->GetNumNodes(); node_index++) // for some reason this needs to go backwards
//                 // {
//                 //     c_vector<double, DIM> LatticeLocation = p_element_1->GetNodeLocation(node_index);
//                 //     if (LatticeLocation[0] <= XCenter)
//                 //     {
//                 //         LatticeSitesOnLeft += 1;
//                 //     }
//                 //     else
//                 //     {
//                 //         LatticeSitesOnRight += 1;
//                 //     }
//                 // }

//                 // // Now check where things are on each side
//                 // for (unsigned node_index = 0; node_index < p_element_2->GetNumNodes(); node_index++)  // for some reason this needs to go backwards
//                 // {
//                 //     c_vector<double, DIM> LatticeLocation = p_element_2->GetNodeLocation(node_index);
//                 //     if (LatticeLocation[0] <= XCenter)
//                 //     {
//                 //         LatticeSitesOnLeft += 1;
//                 //     }
//                 //     else
//                 //     {
//                 //         LatticeSitesOnRight += 1;
//                 //     }
//                 // }

//                 // double TotalNumberOfNodesAfterSwitch = LatticeSitesOnRight + LatticeSitesOnLeft;
//                 // assert(TotalNumberOfNodesAfterSwitch == TotalNodesInCell);
//                 // //moving left because too much in left

//                 // if (MoveDirection == -1 && LatticeSitesOnRight > 0.53 * TotalNumberOfNodesAfterSwitch)
//                 // { //Trying to move right because too much on right
//                 //     success = 1;
//                 // }
//                 // else if (MoveDirection == 1 && LatticeSitesOnLeft > 0.53 * TotalNumberOfNodesAfterSwitch)
//                 // {
//                 //     // Trying to move left because too much on left
//                 //     success = 1;
//                 // }
//                 // NumberOfTrys += 1;
//             // }
//         // }

//         // Now I  need to loop over each of the nodes in each and if they are on the wrong side of the line they need to be switched --- I will need to remember which is left and right buhhh
//         // Assumed i is left and i+1 is right

//         std::vector<unsigned> NodesToRemoveFromElement1;

//         for (unsigned node_index =0;  node_index < p_element_1->GetNumNodes(); node_index++) // - 1; node_index < -1; node_index--) // for some reason this needs to go backwards
//         {
//             c_vector<double, DIM> LatticeLocation = p_element_1->GetNodeLocation(node_index);
//             // If lattice is within certain distance of center and on the right
//             if (XCenter < LatticeLocation[0])// && LatticeLocation[0] < XCenter + 5 * UnitStep)
//             {

//                     TRACE("Removing a cell from the left, and adding it to the right")
//                 // Need to think about changing but also should be checking first
//                 // TRACE("Lattice needs to be removed from the left cell and added to the right")
//                 p_element_2->AddNode(p_element_1->GetNode(node_index)); // Need to record this node/lattice site by its global identifier to add it to the other cell
//                 NodesToRemoveFromElement1.push_back(node_index); // record this node/lattice site by local identifier to delete
//             break;
//             }
//         }
//         for (std::vector<unsigned>::iterator iter = NodesToRemoveFromElement1.begin(); iter != NodesToRemoveFromElement1.end(); iter++)
//         {
//              TRACE("Deleting")
//             p_element_1->DeleteNode(*iter);
//         }

//         std::vector<unsigned> NodesToRemoveFromElement2;
//         for (unsigned node_index = 0; node_index < p_element_2->GetNumNodes();  node_index++)//- 1; node_index < -1; node_index--) // for some reason this needs to go backwards
//         {
//             c_vector<double, DIM> LatticeLocation = p_element_2->GetNodeLocation(node_index);

//             // THis is the right cell, so everthing should be greater than the center point.

//             if (XCenter > LatticeLocation[0])//0 && LatticeLocation[0] > XCenter - 5 * UnitStep)
//             {
//                 TRACE("Removing a cell from the right, and adding it to the left")
//                 // TRACE("Lattice needs to be removed from the right and added to the left ")
//                 p_element_1->AddNode(p_element_2->GetNode(node_index));

//                 NodesToRemoveFromElement2.push_back(node_index);
//                 break;
//             }
//         }

//         for (std::vector<unsigned>::iterator iter = NodesToRemoveFromElement2.begin(); iter != NodesToRemoveFromElement2.end(); iter++)
//         {
//             TRACE("Deleting")

//             p_element_2->DeleteNode(*iter);
//         }

//         TRACE("Have Shuffled")
//         p_cell_1->GetCellData()->SetItem("Center", XCenter);
//         p_cell_2->GetCellData()->SetItem("Center", XCenter);
//         NodesInLeftCell = p_element_1->GetNumNodes();

//         NodesInRightCell = p_element_2->GetNumNodes();
//         PRINT_2_VARIABLES(NodesInLeftCell, NodesInRightCell)
//          TRACE("Finished shuffling the mid-line")
//     }
// }

           // LatticeSitesOnLeft = 0; LatticeSitesOnRight = 0;
                    // for (unsigned node_index = 0; node_index < p_element_1->GetNumNodes(); node_index++) // for some reason this needs to go backwards
                    // {  c_vector<double, DIM> LatticeLocation = p_element_1->GetNodeLocation(node_index);

                    //     if (LatticeLocation[0] <= XCenter && LatticeLocation[0] > XCenter - 5 * UnitStep)
                    //     {    LatticeSitesOnLeft += 1;}
                    //     else
                    //     {  LatticeSitesOnRight += 1;  }
                    // }

                    // // Now check where things are on each side
                    // for (unsigned node_index = 0; node_index < p_element_2->GetNumNodes(); node_index++) // for some reason this needs to go backwards
                    // {  c_vector<double, DIM> LatticeLocation = p_element_2->GetNodeLocation(node_index);
                    //     if (LatticeLocation[0] <= XCenter && LatticeLocation[0] > XCenter - 5 * UnitStep)
                    //     { LatticeSitesOnLeft += 1; }
                    //     else {LatticeSitesOnRight += 1; }
                    // }
























































