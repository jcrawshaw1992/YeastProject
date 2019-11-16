

#include "PottsCellPropertiesModifier.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "PottsArbitrarySurfaceIn3DMesh.hpp"
#include "WrappedPottsBasedCellPopulation.hpp"

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
             unsigned Cell_Index= rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

             double Perimeter = p_static_cast_potts_mesh->GetPerimeterOfElement(Cell_Index);
             double Area = p_static_cast_potts_mesh->GetVolumeOfElement(Cell_Index) ; 


            cell_iter->GetCellData()->SetItem("Area", Area);
    		cell_iter->GetCellData()->SetItem("Perimeter", Perimeter);


            PottsElement<DIM>* p_element = PottsPopulation->GetElement(Cell_Index);
            
            // Get center point -- this will let me track the center point 
            c_vector<double, DIM> CenterPoint; 
            for (unsigned node_index = 0; node_index < p_element->GetNumNodes(); node_index++) // for some reason this needs to go backwards
            {
                CenterPoint += p_element->GetNodeLocation(node_index);
            }
             CenterPoint /=p_element->GetNumNodes();

             cell_iter->GetCellData()->SetItem("CenterX", CenterPoint[0]);
             cell_iter->GetCellData()->SetItem("CenterY", CenterPoint[1]);
             cell_iter->GetCellData()->SetItem("CenterZ", CenterPoint[2]);
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

    //    p_static_cast_potts_mesh->MapCylinderToPlane();
    //     p_static_cast_potts_mesh->CalculateEdgeLenghts();
    //     p_static_cast_potts_mesh->CalculateLatticeVolumes();

    if (PottsPopulation->IsPottsSimulationPeriodic())
    {
        // Periodic potts population -- need to move center lines 

        std::vector<c_vector<unsigned, 2> > ElementPairing = PottsPopulation->GetElementPairingVector();
        //
        double UnitStep = mWidth / (2 * (mN_D - 1)); // *1/2 because the top of a triangle is 1/2 between the two bottom nodes
        PRINT_VARIABLE(UnitStep)
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

            double TotalNodesInCell = p_element_1->GetNumNodes() + p_element_2->GetNumNodes();

            double NodesInCell1 = p_element_1->GetNumNodes();

            double NodesInCell2 = p_element_2->GetNumNodes();
            TRACE("pre Shuffled")
            PRINT_2_VARIABLES(NodesInCell1, NodesInCell2)

            if (NodesInCell2 > 2 && NodesInCell1 > 2)
            {

                unsigned MoveDirection = 0;
                if (p_element_1->GetNumNodes() > 0.7 * TotalNodesInCell)
                {
                    // There are too mamy nodes in the left cell, need to move the midline left
                    MoveDirection = 1;
                }
                else if (p_element_2->GetNumNodes() > 0.7 * TotalNodesInCell)
                {
                    // There are too mamy nodes in the right cell, need to move the midline right
                    MoveDirection = -1;
                }

                if (MoveDirection != 0)
                {
                    double success = 0;  double NumberOfTrys = 1;
                    while (success == 0)
                    {
                        if (NumberOfTrys > 6)
                        {
                            if (NumberOfTrys > 10)
                            {
                                success = 1;
                                TRACE("Give up on the midline shuffel and call it good")
                            }
                            PRINT_VARIABLE(NumberOfTrys)
                        }

                        // Need to move the center line one unit
                        XCenter += MoveDirection * UnitStep;

                        if (XCenter < 0)
                        {
                            TRACE("WE HAVE FALLEND OF THE LEFT EDGE")
                            XCenter = mWidth;
                        }
                        else if (XCenter > mWidth)
                        {
                            XCenter = 0;
                            TRACE("WE HAVE FALLEND OF THE RIGHT EDGE")
                        }

                        // Now check where things are on each side
                        double LatticeSitesOnLeft = 0;
                        double LatticeSitesOnRight = 0;
                        for (unsigned node_index = p_element_1->GetNumNodes() - 1; node_index < -1; node_index--) // for some reason this needs to go backwards
                        {
                            c_vector<double, DIM> LatticeLocation = p_element_1->GetNodeLocation(node_index);
                            if (LatticeLocation[0] <= XCenter)
                            {
                                LatticeSitesOnLeft += 1;
                            }
                            else
                            {
                                LatticeSitesOnRight += 1;
                            }
                        }

                        // Now check where things are on each side
                        for (unsigned node_index = p_element_2->GetNumNodes() - 1; node_index < -1; node_index--) // for some reason this needs to go backwards
                        {
                            c_vector<double, DIM> LatticeLocation = p_element_2->GetNodeLocation(node_index);
                            if (LatticeLocation[0] <= XCenter)
                            {
                                LatticeSitesOnLeft += 1;
                            }
                            else
                            {
                                LatticeSitesOnRight += 1;
                            }
                        }

                        double TotalNumberOfNodesAfterSwitch = LatticeSitesOnRight + LatticeSitesOnLeft;
                        //moving left because too much in left

                        if (MoveDirection == 1 && LatticeSitesOnRight > 0.55 * TotalNumberOfNodesAfterSwitch)
                        { //Trying to move right because too much on right
                            success = 1;
                        }
                        else if (MoveDirection == -1 && LatticeSitesOnLeft > 0.55 * TotalNumberOfNodesAfterSwitch)
                        {
                            // Trying to move left because too much on left
                            success = 1;
                        }
                        NumberOfTrys += 1;
                    }
                }

                // Now I  need to loop over each of the nodes in each and if they are on the wrong side of the line they need to be switched --- I will need to remember which is left and right buhhh
                // Assumed i is left and i+1 is right

                std::vector<unsigned> NodesToRemoveFromElement1;

                for (unsigned node_index = p_element_1->GetNumNodes() - 1; node_index < -1; node_index--) // for some reason this needs to go backwards
                {
                    c_vector<double, DIM> LatticeLocation = p_element_1->GetNodeLocation(node_index);
                    // If lattice is within certain distance of center and on the right
                    if (XCenter < LatticeLocation[0] && LatticeLocation[0] < XCenter + 5 * UnitStep)
                    {

                        // Need to think about changing but also should be checking first
                        // TRACE("Lattice needs to be removed from the left cell and added to the right")
                        p_element_2->AddNode(p_element_1->GetNode(node_index)); // Need to record this node/lattice site by its global identifier to add it to the other cell
                        NodesToRemoveFromElement1.push_back(node_index); // record this node/lattice site by local identifier to delete
                    }
                }
                for (std::vector<unsigned>::iterator iter = NodesToRemoveFromElement1.begin(); iter != NodesToRemoveFromElement1.end(); iter++)
                {
                    p_element_1->DeleteNode(*iter);
                }

                std::vector<unsigned> NodesToRemoveFromElement2;
                for (unsigned node_index = p_element_2->GetNumNodes() - 1; node_index < -1; node_index--) // for some reason this needs to go backwards
                {
                    c_vector<double, DIM> LatticeLocation = p_element_2->GetNodeLocation(node_index);

                    // THis is the right cell, so everthing should be greater than the center point.

                    if (XCenter > LatticeLocation[0] && LatticeLocation[0] > XCenter - 5 * UnitStep)
                    {

                        // TRACE("Lattice needs to be removed from the right and added to the left ")
                        p_element_1->AddNode(p_element_2->GetNode(node_index));

                        NodesToRemoveFromElement2.push_back(node_index);
                    }
                }

                for (std::vector<unsigned>::iterator iter = NodesToRemoveFromElement2.begin(); iter != NodesToRemoveFromElement2.end(); iter++)
                {

                    p_element_2->DeleteNode(*iter);
                }

                TRACE("Have Shuffled")
                p_cell_1->GetCellData()->SetItem("Center", XCenter);
                p_cell_2->GetCellData()->SetItem("Center", XCenter);
                NodesInCell1 = p_element_1->GetNumNodes();

                NodesInCell2 = p_element_2->GetNumNodes();
                PRINT_2_VARIABLES(NodesInCell1, NodesInCell2)
                 TRACE("Finished shuffling the mid-line")
            }
        }

    }else
    { // Not periodic 

    // I need the center point of each cell, and the area and perimeter 

       for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
             unsigned Cell_Index= rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

             double Perimeter = p_static_cast_potts_mesh->GetPerimeterOfElement(Cell_Index);
             double Area = p_static_cast_potts_mesh->GetVolumeOfElement(Cell_Index) ; 


            cell_iter->GetCellData()->SetItem("Area", Area);
    		cell_iter->GetCellData()->SetItem("Perimeter", Perimeter);


            PottsElement<DIM>* p_element = PottsPopulation->GetElement(Cell_Index);
            
            // Get center point -- this will let me track the center point 
            c_vector<double, DIM> CenterPoint; 
            for (unsigned node_index = p_element->GetNumNodes() - 1; node_index < -1; node_index--) // for some reason this needs to go backwards
            {
                CenterPoint += p_element->GetNodeLocation(node_index);
            }
             CenterPoint /=p_element->GetNumNodes();

             cell_iter->GetCellData()->SetItem("CenterX", CenterPoint[0]);
             cell_iter->GetCellData()->SetItem("CenterY", CenterPoint[1]);
             cell_iter->GetCellData()->SetItem("CenterZ", CenterPoint[2]);
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
