
/*

This code controles the membrame properties. 

At the start of a simulation, the intial conditions delievered at 20(?)X less than original starting configuration of the vessel 
The system then develops to an equlibirum (that will be the realisitc size of a network). After this point the memnbrane properties will 
be variable.

For the inital impletmentation of this code, the membrane properties in the center will be changed to decrease the vessel size to 2X the intial 
condition, with some smoothing on the vessel properties just outside the region. 

After this, the membrane constants will be dependent on the number of cells in the region and the denoted mutation state. 

*/

#include "MembranePropertiesModifier.hpp"
#include <algorithm>
#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"

#include <math.h>
#include "Debug.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MembranePropertiesModifier<ELEMENT_DIM, SPACE_DIM>::MembranePropertiesModifier()
        : AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MembranePropertiesModifier<ELEMENT_DIM, SPACE_DIM>::~MembranePropertiesModifier()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MembranePropertiesModifier<ELEMENT_DIM, SPACE_DIM>::SetMembranePropeties(std::map<double, c_vector<long double, 4> > GrowthMaps, double Strength, bool Hetrogeneous, double StepSize)
{
    mGrowthMaps = GrowthMaps;
    // TRACE("New membrane properties")
    // PRINT_2_VARIABLES(mGrowthMaps[10][3], GrowthMaps[10][3]);

    mOn = 1;
    mAchievedTargetK = 0;
    mStrength = Strength;
    mHetro = Hetrogeneous;
    mStepSize = StepSize;
    mCounter =50;
    mThreshold =50;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MembranePropertiesModifier<ELEMENT_DIM, SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
    TRACE("In the SetUpSolve");
    // PRINT_2_VARIABLES(ELEMENT_DIM, SPACE_DIM);
    assert(SPACE_DIM == 3); assert(ELEMENT_DIM == 2);
    MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);
    MAKE_PTR(EmptyBasementMatrix, p_Basement); //Mutation to mark nodes one basement
    MAKE_PTR(HasEndothelialCell, p_EC); //Mutation to mark nodes with ECs
    MAKE_PTR(LostEndothelialCell, p_NextToBasement); //Mutation to mark nodes with ECs next to the basement
    c_vector<long double, 3> Node_location;
    double MinZ = 1e-3;
    
    std::set<unsigned> MutantNodeIndices; std::set<unsigned> EdgeMutantNodeIndices;
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        Node_location = rCellPopulation.GetNode(node_index)->rGetLocation();
        if (abs(Node_location[2]) < MinZ) // These are the nodes along the lower edge and need to be marked as mutated
        {
            cell_iter->SetMutationState(p_Basement);
            MutantNodeIndices.insert(node_index);
            mSamplebasementNode = node_index;
        }
        else
        {
            cell_iter->SetMutationState(p_EC);
            mSampleECNode = node_index;
        }
        cell_iter->GetCellData()->SetItem("ShearModulus", mGrowthMaps[mStrength](2));
        cell_iter->GetCellData()->SetItem("AreaDilationModulus", mGrowthMaps[mStrength](1));
        cell_iter->GetCellData()->SetItem("AreaConstant", mGrowthMaps[mStrength](0));
        cell_iter->GetCellData()->SetItem("BendingConstant", mGrowthMaps[mStrength](3));
        // PRINT_3_VARIABLES(mGrowthMaps[10](2), mGrowthMaps[10](1), mGrowthMaps[10](0));
    }
    mMaxBasementZ = -10;
    mMinBasementZ = 0;
    // PRINT_2_VARIABLES( mMaxBasementZ, mHetro);
    if (mHetro)
    {
        // Loop over the mutants and get the nodes along the boundary
        for (std::set<unsigned>::iterator iter = MutantNodeIndices.begin();
             iter != MutantNodeIndices.end();
             ++iter)
        {
            // Need to check if any of the neighbours are an EC, if so
            // Save the EC as a mutant on the boundary cell_iter->SetMutationState(p_NextToBasement);
            mBasementNodes.push_back(*iter);
            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
            std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringLocationIndices(p_cell); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
            for (std::set<unsigned>::iterator N_iter = neighbouring_node_indices.begin();
                 N_iter != neighbouring_node_indices.end();
                 ++N_iter)
            {
                CellPtr N_cell = rCellPopulation.GetCellUsingLocationIndex(*N_iter);
                if (N_cell->GetMutationState()->template IsType<HasEndothelialCell>())
                {
                    // This cell is type EC and next to the basement matrix, which means the
                    // this basement mutant is on the edge :)
                    EdgeMutantNodeIndices.insert(*iter);
                    // Node_location = rCellPopulation.GetNode(*iter)->rGetLocation();

                    Node_location = rCellPopulation.GetNode(*N_iter)->rGetLocation();
                    if (mMinBasementZ > Node_location[2])
                    {
                        mMinBasementZ = Node_location[2];
                    }
                    else if (mMaxBasementZ < Node_location[2])
                    {
                        mMaxBasementZ = Node_location[2];
                    }
                }
            }
        }

        // std::map<unsigned, std::vector<unsigned> > NeighboursBoundary;
        std::vector<unsigned> NeighboursBoundary;
        double MaxZ = 0;
        double Minz = 100;
        CellPtr EcCell;
        // vector<unsigned> NeighbourCollection;
        for (std::set<unsigned>::iterator iter = EdgeMutantNodeIndices.begin();
             iter != EdgeMutantNodeIndices.end();
             ++iter)
        {
            // Just collect the boundaries for maxiumum ease of interpretation later

            // Here we have an edge with a EC region, so we should grab all the ECs in the
            // between the nodies location and  Z away in either direction that is an EC
            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
    
            std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringLocationIndices(p_cell); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);

            for (std::set<unsigned>::iterator N_iter = neighbouring_node_indices.begin();
                 N_iter != neighbouring_node_indices.end();
                 ++N_iter)
            {
                CellPtr N_cell = rCellPopulation.GetCellUsingLocationIndex(*N_iter);
                if (N_cell->GetMutationState()->template IsType<HasEndothelialCell>())
                {
                    // Set this as a boundary
                    // N_cell->SetMutationState(p_NextToBasement);
                    // Save in the collection vector
                    NeighboursBoundary.push_back(*N_iter);
                    // EcCell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                    // EcCell->SetMutationState(p_NextToBasement);
                    std::set<unsigned> neighbouring_node_indices2 = rCellPopulation.GetNeighbouringLocationIndices(N_cell); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
                    for (std::set<unsigned>::iterator N_iter_2 = neighbouring_node_indices2.begin();
                         N_iter_2 != neighbouring_node_indices2.end();
                         ++N_iter_2)
                    {
                        CellPtr N_cell_2 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_2);
                        if (N_cell_2->GetMutationState()->template IsType<HasEndothelialCell>())
                        {

                            NeighboursBoundary.push_back(*N_iter_2);

                            std::set<unsigned> neighbouring_node_indices3 = rCellPopulation.GetNeighbouringLocationIndices(N_cell_2); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);

                            for (std::set<unsigned>::iterator N_iter_3 = neighbouring_node_indices3.begin();
                                 N_iter_3 != neighbouring_node_indices3.end();
                                 ++N_iter_3)
                            {

                                CellPtr N_cell_3 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_3);

                                if (N_cell_3->GetMutationState()->template IsType<HasEndothelialCell>())
                                {
                                    //  TRACE("if endothilial state");
                                    NeighboursBoundary.push_back(*N_iter_3);

                                    std::set<unsigned> neighbouring_node_indices4 = rCellPopulation.GetNeighbouringLocationIndices(N_cell_3); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
                                    for (std::set<unsigned>::iterator N_iter_4 = neighbouring_node_indices4.begin();
                                         N_iter_4 != neighbouring_node_indices4.end();
                                         ++N_iter_4)
                                    {

                                        CellPtr N_cell_4 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_4);
                                        if (N_cell_4->GetMutationState()->template IsType<HasEndothelialCell>())
                                        {
                                            NeighboursBoundary.push_back(*N_iter_4);

                                            // --------------------
                                            std::set<unsigned> neighbouring_node_indices5 = rCellPopulation.GetNeighbouringLocationIndices(N_cell_4); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
                                            for (std::set<unsigned>::iterator N_iter_5 = neighbouring_node_indices5.begin();
                                                 N_iter_5 != neighbouring_node_indices5.end();
                                                 ++N_iter_5)
                                            {

                                                CellPtr N_cell_5 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_5);
                                                if (N_cell_5->GetMutationState()->template IsType<HasEndothelialCell>())
                                                {
                                                    NeighboursBoundary.push_back(*N_iter_5);

                                                    // --------------------
                                                    std::set<unsigned> neighbouring_node_indices6 = rCellPopulation.GetNeighbouringLocationIndices(N_cell_5); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
                                                    for (std::set<unsigned>::iterator N_iter_6 = neighbouring_node_indices6.begin();
                                                         N_iter_6 != neighbouring_node_indices6.end();
                                                         ++N_iter_6)
                                                    {

                                                        CellPtr N_cell_6 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_6);
                                                        if (N_cell_6->GetMutationState()->template IsType<HasEndothelialCell>())
                                                        {
                                                            NeighboursBoundary.push_back(*N_iter_6);

                                                            // -------------------
                                                            // -------------------


                                                            std::set<unsigned> neighbouring_node_indices7 = rCellPopulation.GetNeighbouringLocationIndices(N_cell_6); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
                                                            for (std::set<unsigned>::iterator N_iter_7 = neighbouring_node_indices7.begin();
                                                                N_iter_7 != neighbouring_node_indices7.end();
                                                                ++N_iter_7)
                                                            {

                                                                CellPtr N_cell_7 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_7);
                                                                if (N_cell_7->GetMutationState()->template IsType<HasEndothelialCell>())
                                                                {
                                                                    NeighboursBoundary.push_back(*N_iter_7);

                                                                     std::set<unsigned> neighbouring_node_indices8 = rCellPopulation.GetNeighbouringLocationIndices(N_cell_7); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
                                                                    for (std::set<unsigned>::iterator N_iter_8 = neighbouring_node_indices8.begin();
                                                                        N_iter_8 != neighbouring_node_indices8.end();
                                                                        ++N_iter_8)
                                                                    {

                                                                        CellPtr N_cell_8 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_8);
                                                                        if (N_cell_8->GetMutationState()->template IsType<HasEndothelialCell>())
                                                                        {
                                                                            NeighboursBoundary.push_back(*N_iter_8);


                                                                            //  std::set<unsigned> neighbouring_node_indices9 = rCellPopulation.GetNeighbouringLocationIndices(N_cell_8); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
                                                                            //     for (std::set<unsigned>::iterator N_iter_9 = neighbouring_node_indices9.begin();
                                                                            //         N_iter_9 != neighbouring_node_indices9.end();
                                                                            //         ++N_iter_9)
                                                                            //     {

                                                                            //         CellPtr N_cell_9 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_9);
                                                                            //         if (N_cell_9->GetMutationState()->template IsType<HasEndothelialCell>())
                                                                            //         {
                                                                            //             NeighboursBoundary.push_back(*N_iter_9);


                                                                            //             std::set<unsigned> neighbouring_node_indices10 = rCellPopulation.GetNeighbouringLocationIndices(N_cell_9); // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);  // std::set<unsigned> neighbouring_node_indices = rCellPopulation.GetNeighbouringNodeIndices(node_index);
                                                                            //                     for (std::set<unsigned>::iterator N_iter_10 = neighbouring_node_indices10.begin();
                                                                            //                         N_iter_10 != neighbouring_node_indices10.end();
                                                                            //                         ++N_iter_10)
                                                                            //                     {

                                                                            //                         CellPtr N_cell_10 = rCellPopulation.GetCellUsingLocationIndex(*N_iter_10);
                                                                            //                         if (N_cell_10->GetMutationState()->template IsType<HasEndothelialCell>())
                                                                            //                         {
                                                                            //                             NeighboursBoundary.push_back(*N_iter_10);
                                                                            //                         }
                                                                            //                     }


                                                                                
                                                                            //         }
                                                                            //     }




                                                                        }
                                                                    }





                                                                }
                                                            }


                                                        }
                                                    }
                                                    // -------------------
                                                }
                                            }
                                            // -------------------
                                        }
                                    }
                                }
                            }
                        }
                        // Loop over the neighbours for several levels
                    }
                }
            }
        }

        sort(NeighboursBoundary.begin(), NeighboursBoundary.end());
        NeighboursBoundary.erase(unique(NeighboursBoundary.begin(), NeighboursBoundary.end()), NeighboursBoundary.end());

        for (std::vector<unsigned>::iterator it = NeighboursBoundary.begin(); it != NeighboursBoundary.end(); ++it)
        {
            CellPtr BoundaryCell = rCellPopulation.GetCellUsingLocationIndex(*it);
            BoundaryCell->SetMutationState(p_NextToBasement);
            c_vector<double, SPACE_DIM> location2 = rCellPopulation.GetLocationOfCellCentre(BoundaryCell);

            // Set cell propoerties about how far it is from basement
            double distance = location2[2]; //location1[2] - location2[2];
            if (distance < MinZ)
            {
                MinZ = distance;
            }
            else if (distance > MaxZ)
            {
                MaxZ = distance;
            }
            BoundaryCell->GetCellData()->SetItem("DistanceFromBasement", distance);
        }

        // TRACE("Out of loop");
        // MinZ and MaxZ are the minimum and maximum edges of the edge region
        mMinZ = MinZ;
        mMaxZ = MaxZ;
        // PRINT_2_VARIABLES(mMinZ, mMaxZ);

        mNodesNextToBasement = NeighboursBoundary;
        // mBasementNodes = MutantNodeIndices;
    }
    else
    {
        TRACE("NOw in the else part")
        mMinZ = 0;
        mMaxZ = 0;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MembranePropertiesModifier<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    

    if (mHetro)
    {
        PRINT_VARIABLE(mHetro)
        assert(SPACE_DIM == 3);
        // See if any edges are too long and if so divide them
        double num_cells = rCellPopulation.GetNumRealCells();

        // MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_rCellPopulation = dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(rCellPopulation));
        if (mOn == 1 )
        {
            if (mCounter ==mThreshold)
            {
                UpdateCellData(rCellPopulation);
                mCounter =0;
            }
            else
            {
                mCounter +=1;
            }
        }
    }
    
    
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MembranePropertiesModifier<ELEMENT_DIM, SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
   
    /*
        This needs to be properly implimented 
        -- Temp --
        Initally the membrane properties are set for 20X Growth, 
        When the radius gets to 20X inital, the membrane properties in the center will be changed to decreased

        -- To Impliment -- 
        Need to adapt this to sweep over all nodes and adapt the membrane properties of the ones that have lost or gained a Potts lattice
	*/
    assert(SPACE_DIM == 3); // Currently assumes that SPACE_DIM = 3

    CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(mSamplebasementNode);

    double Step_Kbs = p_cell->GetCellData()->GetItem("ShearModulus") + mStepSize*10;
    double Step_Kba = p_cell->GetCellData()->GetItem("AreaDilationModulus") +mStepSize;// 1e-15;
    double Step_KbA = p_cell->GetCellData()->GetItem("AreaConstant") + mStepSize;// 1e-12;

    double Kbs = std::min((double)Step_Kbs, (double)mGrowthMaps[1.2](2));
    double Kba = std::min((double)Step_Kba, (double)mGrowthMaps[1.2](1));
    double KbA = std::min((double)Step_KbA, (double)mGrowthMaps[1.2](0));
    // double Kbb = std::min((double)Step_Kb, (double)mGrowthMaps[1.2](3));
    // PRINT_4_VARIABLES(Kbs, Kba, KbA, Kbb)
    double Ks = mGrowthMaps[mStrength](2);
    double Ka = mGrowthMaps[mStrength](1);
    double KA = mGrowthMaps[mStrength](0);
    // double KB = mGrowthMaps[mStrength](3);

    // Quadradtic --  Fit a quadratic curve to smooth the membrane parameters between the basement and the ECs
    // K = a(Z-b)^2 +c
    // Left side
        // The basement Z
        double Zb = mMinBasementZ;

        // The endothelial Z mMinZ;
        double b1 = Zb;                 //  mMinZ + Zb;
        double c1S = Kbs;               // Ks - a1S * Zb * Zb;
        double a1S = (Ks - c1S) / ((mMinZ - b1)*(mMinZ - b1));

        //Alpha 
        double c1a = Kba;               //Ka - a1a * Zb * Zb;
        double a1a = (Ka - c1a) / ((mMinZ - b1)*(mMinZ - b1));//           (Kba - Ka) / (mMinZ * mMinZ - Zb * Zb);
        

        //Area  
        double c1A = KbA;               //Ka - a1a * Zb * Zb;
        double a1A = (KA - c1A) / ((mMinZ - b1)*(mMinZ - b1));//           (Kba - Ka) / (mMinZ * mMinZ - Zb * Zb);
    
        // Right side

        
        // The basement Z
        Zb = mMaxBasementZ;

        // The endothelial Z   mMaxZ;

        double b2 = Zb;             //mMaxZ + Zb;
        double c2S = c1S;           // Ks - a1S * Zb * Zb;
        double a2S = (Ks - c2S) / ((mMaxZ - b2)*(mMaxZ - b2));               // (Kbs - Ks) / (mMaxZ * mMaxZ - Zb * Zb);
        

        double c2a = c1a;           // Ks - a1S * Zb * Zb;
        double a2a = (Ka - c2a) / ((mMaxZ - b2)*(mMaxZ - b2));               // (Kbs - Ks) / (mMaxZ * mMaxZ - Zb * Zb);
       
        double c2A = c1A;           // Ks - a1S * Zb * Zb;
        double a2A = (KA - c2A) / ((mMaxZ - b2)*(mMaxZ - b2));               // (Kbs - Ks) / (mMaxZ * mMaxZ - Zb * Zb);
       
    //    double c2B = c1B;           // Ks - a1S * Zb * Zb;
    //    double a2B = (KB - c2B) / ((mMaxZ - b2)*(mMaxZ - b2));               // (Kbs - Ks) / (mMaxZ * mMaxZ - Zb * Zb);
       
        double ShearMod;
        double AlphaMod;
        double AreaMod;
        double BendingMod;

    if (mAchievedTargetK == 0)
    {
        // MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);

        for (std::vector<unsigned>::iterator it = mBasementNodes.begin(); it != mBasementNodes.end(); ++it)
        {

            CellPtr cell_iter = rCellPopulation.GetCellUsingLocationIndex(*it);
            cell_iter->GetCellData()->SetItem("ShearModulus", Kbs);
            cell_iter->GetCellData()->SetItem("AreaDilationModulus", Kba);
            cell_iter->GetCellData()->SetItem("AreaConstant", KbA );
         
            // cell_iter->GetCellData()->SetItem("BendingConstant", Kbb);


        }

        for (std::vector<unsigned>::iterator it = mNodesNextToBasement.begin(); it != mNodesNextToBasement.end(); ++it)
        {

            CellPtr cell_iter = rCellPopulation.GetCellUsingLocationIndex(*it);

            // Change the membrane parameters of the nodes neighbouring the basement mutants -- do this linearly in time and quadratically in space
            c_vector<double, SPACE_DIM> location = rCellPopulation.GetLocationOfCellCentre(cell_iter);
            if (location[2] < 0)
            {
                //K =a(Z -b)^2 +c
                ShearMod = a1S * pow(location[2] - b1, 2) + c1S;
                AlphaMod = a1a * pow(location[2] - b1, 2) + c1a;
                AreaMod = a1A * pow(location[2] - b1, 2) + c1A;
                // BendingMod = a1B * pow(location[2] - b1, 2) + c1B;
            }
            else
            {
                ShearMod = a2S * pow(location[2] - b2, 2) + c2S;
                AlphaMod = a2a * pow(location[2] - b2, 2) + c2a;
                AreaMod = a2A * pow(location[2] - b2, 2) + c2A;
                // BendingMod = a2B * pow(location[2] - b2, 2) + c2B;
            }

            cell_iter->GetCellData()->SetItem("ShearModulus", ShearMod);
            cell_iter->GetCellData()->SetItem("AreaDilationModulus", AlphaMod);
            cell_iter->GetCellData()->SetItem("AreaConstant", AreaMod);
            // cell_iter->GetCellData()->SetItem("BendingConstant", BendingMod);
        }

            if (Step_KbA > mGrowthMaps[1.2](0) && Step_Kba > mGrowthMaps[1.2](1) && Step_Kbs > mGrowthMaps[1.2](2))
            {
                mAchievedTargetK = 1;
                TRACE("Area, shear and dilation hit limit");
                 TRACE("Simulation can stop very soon :) ");
                mOn = 0;
            } //else if (Step_KbA > mGrowthMaps[4](0) && Step_Kba > mGrowthMaps[4](1) && Step_Kbs > mGrowthMaps[4](2))
            // {
            //      mThreshold =150;
            // }

    }
    else if (mAchievedTargetK == 1)
    {
        TRACE("Target achieved, shouldnt get here ")
    }

}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MembranePropertiesModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class MembranePropertiesModifier<1, 1>;
template class MembranePropertiesModifier<1, 2>;
template class MembranePropertiesModifier<2, 2>;
template class MembranePropertiesModifier<1, 3>;
template class MembranePropertiesModifier<2, 3>;
template class MembranePropertiesModifier<3, 3>;
// Serialization for Boost >= 1.36

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MembranePropertiesModifier)

// Left side
// double Zb = mMinBasementZ;
// double m1S = (Ks - Kbs) / (mMinZ - Zb);
// double c1S = Ks - mMinZ * m1S;

// double m1a = (Ka - Kba) / (mMinZ - Zb);
// double c1a = Ka - mMinZ * m1a;

// double m1A = (KA - KbA) / (mMinZ - Zb);
// double c1A = KA - mMinZ * m1A;

// // Right side
// Zb = mMaxBasementZ;//1e-3;

// double m2S = (Ks - Kbs) / (mMaxZ - Zb);
// double c2S = Ks - m2S * mMaxZ;

// double m2a = (Ka - Kba) / (mMaxZ - Zb);
// double c2a = Ka - m2a * mMaxZ;

// double m2A = (KA - KbA) / (mMaxZ - Zb);
// double c2A = KA - m2A * mMaxZ;

// // Quadradtic

// // Left side
// // Zb = -1e-3;
// Zb =  mMinBasementZ;
// m1S = (Ks - Kbs) / (mMinZ*mMinZ - Zb*Zb);
// c1S =  Ks - m1S * mMinZ*mMinZ;

// m1a = (Ka - Kba) / (mMinZ*mMinZ - Zb*Zb);
// c1a =  Ka - m1a * mMinZ*mMinZ;

// m1A = (KA - KbA) / (mMinZ*mMinZ - Zb*Zb);
// c1A =  KA - m1A * mMinZ*mMinZ;

// // Right side
// // Zb = 1e-3;
// Zb = mMaxBasementZ;
// m2S = (Ks - Kbs) / (mMaxZ*mMaxZ - Zb*Zb);
// c2S =  Ks - m2S * mMaxZ*mMaxZ;

// m2a = (Ka - Kba) / (mMaxZ*mMaxZ - Zb*Zb);
// c2a =  Ka - m2a * mMaxZ*mMaxZ;

// m2A = (KA - KbA) / (mMaxZ*mMaxZ - Zb*Zb);
// c2A =  KA - m2A * mMaxZ*mMaxZ;

// Zb = -1e-3;
// Kbs = mGrowthMaps[1.2](2);
// Kba = mGrowthMaps[1.2](2);
// KbA = mGrowthMaps[1.2](2);

// / // Quadradtic y =a(x -b)^2 +c

//     // Left side
//     // Zb = -1e-3;
//     double Zb
//     = mMinBasementZ;
// double dz = Zb - mMinZ;

// double a1S = (Kbs - Ks) / (mMinZ * mMinZ - Zb * Zb); // (-Ks + Kbs)/( pow(mMinZ- Zb,2) -1 );
// double b1 = mMinZ + Zb; //Zb;//
// double c1S = Ks - a1S * Zb * Zb;

// double a1a = (Kba - Ka) / (mMinZ * mMinZ - Zb * Zb); //  (Ka - Kba)/( pow(mMinZ- Zb,2) -1 );
// double c1a = Ka - a1a * Zb * Zb; //Kba - a1a;

// double a1A = (KbA - KA) / (mMinZ * mMinZ - Zb * Zb); //   (KA - KbA)/( pow(mMinZ- Zb,2) -1 );
// double c1A = KA - a1A * Zb * Zb; //KbA - a1A;

// // RIght side
// // Zb = -1e-3;
// Zb = mMaxBasementZ;
// dz = Zb - mMaxZ;

// double a2S = (Kbs - Ks) / (mMaxZ * mMaxZ - Zb * Zb); // (-Ks + Kbs)/( pow(mMaxZ- Zb,2) -1 );
// double b2 = mMaxZ + Zb; //Zb;//
// // c1S = Kbs - a1S;
// double c2S = Ks - a1S * Zb * Zb;

// double a2a = (Kba - Ka) / (mMaxZ * mMaxZ - Zb * Zb); //  (Ka - Kba)/( pow(mMaxZ- Zb,2) -1 );
// double c2a = Ka - a1a * Zb * Zb; //Kba - a1a;

// double a2A = (KbA - KA) / (mMaxZ * mMaxZ - Zb * Zb); //   (KA - KbA)/( pow(mMinZ- Zb,2) -1 );
// double c2A = KA - a1A * Zb * Zb; //KbA - a1A;

// // Quadradtic

// MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);
// for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
//      cell_iter != rCellPopulation.End();
//      ++cell_iter)
// {
//     if (cell_iter->GetMutationState()->template IsType<EmptyBasementMatrix>())
//     {
//         cell_iter->GetCellData()->SetItem("ShearModulus", mGrowthMaps[2](2));
//         cell_iter->GetCellData()->SetItem("AreaDilationModulus", mGrowthMaps[2](1));
//         cell_iter->GetCellData()->SetItem("AreaConstant", mGrowthMaps[2](0));
//     }
//     else if (cell_iter->GetMutationState()->template IsType<LostEndothelialCell>())
//     {
//         // SO now easy K = m1 Z + c1
//         c_vector<double, SPACE_DIM> location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
//         if (location[2] < 0)
//         {
//             ShearMod = m1S * location[2] + c1S;
//             AlphaMod = m1a * location[2] + c1a;
//             AreaMod = m1A * location[2] + c1A;

//             // ShearMod = m1S *pow(location[2]  - c1S,2);
//             // AlphaMod = m1a *pow(location[2]  - c1a,2);
//             // AreaMod  = m1A *pow(location[2]  - c1A,2);
//         }
//         else
//         {
//             ShearMod = m2S * location[2] + c2S;
//             AlphaMod = m2a * location[2] + c2a;
//             AreaMod = m2A * location[2] + c2A;

//             // ShearMod = m2S *pow(location[2]  - c2S,2);
//             // AlphaMod = m2a *pow(location[2]  - c2a,2);
//             // AreaMod  = m2A *pow(location[2]  - c2A,2);
//         }

//         cell_iter->GetCellData()->SetItem("ShearModulus", ShearMod);
//         cell_iter->GetCellData()->SetItem("AreaDilationModulus", AlphaMod);
//         cell_iter->GetCellData()->SetItem("AreaConstant", AreaMod);
//     }
// }
