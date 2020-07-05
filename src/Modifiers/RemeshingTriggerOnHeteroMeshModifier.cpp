
/*

    This code controles the membrame properties and can trigger remeshing 

    -- I needed to put the remeshing trigger in this modifier (it was originally in a separate and simple modifier), but when it is a hetro mesh, I need to 
    update the cell properties, which I need to know where the hetero boundaries are. Simplest and easiest way to do this 

    At the start of a simulation, the intial conditions delievered at 20(?)X less than original starting configuration of the vessel 
    The system then develops to an equlibirum (that will be the realisitc size of a network). After this point the memnbrane properties will 
    be variable.

    For the inital impletmentation of this code, the membrane properties in the center will be changed to decrease the vessel size to 2X the intial 
    condition, with some smoothing on the vessel properties just outside the region. 

    After this, the membrane constants will be dependent on the number of cells in the region and the denoted mutation state. 

*/

#include "RemeshingTriggerOnHeteroMeshModifier.hpp"
#include <algorithm>
#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"

#include <math.h>
#include "Debug.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::RemeshingTriggerOnHeteroMeshModifier()
        : AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::~RemeshingTriggerOnHeteroMeshModifier()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetMembranePropeties(std::map<double, c_vector<long double, 4> > GrowthMaps, double Strength, bool Hetrogeneous, double StepSize, double SetupSolve)
{
    mGrowthMaps = GrowthMaps;
    mStrength = Strength;
    mHetro = Hetrogeneous;
    mStepSize = StepSize;
    mOn =  Hetrogeneous;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetMembraneStrength(double Strength)
{
    mStrength = Strength;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetThreshold(double Threshold)
{
    mThreshold = Threshold;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::Boundaries(c_vector<double, 3> UpperPlaneNormal, c_vector<double, 3> UpperPlanePoint, c_vector<double, 3> LowerPlaneNormal, c_vector<double, 3> LowerPlanePoint)
{
    // Set the boundary planes for this hetro region, set an upper and a lower bound.
    std::vector<  c_vector<double, 3> > CurrentBoundary;
    CurrentBoundary.push_back(UpperPlaneNormal);
    CurrentBoundary.push_back(UpperPlanePoint);

    CurrentBoundary.push_back(LowerPlaneNormal);
    CurrentBoundary.push_back(LowerPlanePoint);

    mBoundaries.push_back(CurrentBoundary);
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetRemeshingInterval(int RemeshingInterval)
{
    mRemeshing = 1;
    mRemeshingInterval = RemeshingInterval;
}




template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
        assert(SPACE_DIM ==3);
        MAKE_PTR(EmptyBasementMatrix, p_Basement);
        MAKE_PTR(HasEndothelialCell, p_EC);
        c_vector<double, 3> Node_location;

        for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node_location = rCellPopulation.GetNode(node_index)->rGetLocation();

            cell_iter->SetMutationState(p_EC);

            // Here we need to determine where exatly they 
            // Need to loop over the given boundaries, and see which planes these are defined within
            // std::vector<std::vector<  c_vector<double, 3> > > mBoundaries;
            // PRINT_VARIABLE(mBoundaries.size())
            if (mHetro)
            {
                for (unsigned i = 0; i<mBoundaries.size(); i++)
                {
                    std::vector<  c_vector<double, 3> > BasementRegion =  mBoundaries[i];
                    c_vector<double, 3> UpperPlane = BasementRegion[0];
                    c_vector<double, 3> UpperPoint = BasementRegion[1];

                    c_vector<double, 3> LowerPlane = BasementRegion[2];
                    c_vector<double, 3> LowerPoint = BasementRegion[3];

                    // Vector connecting the node to upper plane
                    c_vector<double, 3> NodeToUpperPlane = Node_location - UpperPoint;
                    c_vector<double, 3> NodeToLowerPlane = Node_location - LowerPoint;

                    c_vector<double, 3> a1 = NodeToLowerPlane;
                    c_vector<double, 3> b1 = - Node_location + LowerPlane;

                    c_vector<double, 3> a2 = NodeToUpperPlane;
                    c_vector<double, 3> b2 = - Node_location + UpperPoint;

                    double a1_1 = abs(inner_prod(NodeToLowerPlane,b1 )/norm_2(b1));
                    double a2_1 = abs(inner_prod(NodeToUpperPlane,b2 )/norm_2(b2));


                    double DotToUpperPlane = inner_prod(NodeToUpperPlane,UpperPlane );
                    double DotToLowerPlane = inner_prod(NodeToLowerPlane,LowerPlane );
                    double radius = 0.1;
                    // PRINT_2_VARIABLES(DotToUpperPlane,DotToLowerPlane)
                    if (DotToLowerPlane > 0 && DotToUpperPlane > 0)
                    {
                        if ( abs(norm_2(Node_location -UpperPoint)) <radius ||  abs(norm_2(Node_location -LowerPlane)) <radius  )
                        {
                            // if (abs(Node_location[1])>0.1   )
                            cell_iter->SetMutationState(p_Basement);
                            mBasementNodes.push_back(node_index);
                            mSamplebasementNode = node_index;
                            // cell_iter->GetCellData()->SetItem("Mut", 1);
                            break;
                        }
                        
                    }
                }
            }
           
            cell_iter->GetCellData()->SetItem("ShearModulus", mGrowthMaps[mStrength](2));
            cell_iter->GetCellData()->SetItem("AreaDilationModulus", mGrowthMaps[mStrength](1));
            cell_iter->GetCellData()->SetItem("AreaConstant", mGrowthMaps[mStrength](0));
            cell_iter->GetCellData()->SetItem("BendingConstant", mGrowthMaps[mStrength](3));
        }
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetBendingForce(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, double BendingConstant)
{
    std::set<unsigned> MutantNodeIndices;
    std::set<unsigned> EdgeMutantNodeIndices;
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        cell_iter->GetCellData()->SetItem("BendingConstant", BendingConstant);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
//    TRACE("UpdateAtEndOfTimeStep")
     HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation = static_cast<HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);
    if(mRemeshing)
    {
         if (mExecute>=mRemeshingInterval)
        {
            pCellPopulation->ExecuteHistoryDependentRemeshing();
            UpdateCellData(rCellPopulation);
            mExecute = 0;

        } 
        mExecute +=1;
        
    }
    // PRINT_VARIABLE(mHetro)

    if (mHetro)
    {
        StepChange(rCellPopulation);
        // assert(SPACE_DIM == 3);
        // // See if any edges are too long and if so divide them
        // // double num_cells = rCellPopulation.GetNumRealCells();

        // if (mHetro == 1)
        // {
        //     if (mCounter == mThreshold)
        //     {
        //         StepChange(rCellPopulation);
        //         mCounter = 0;
        //     }
        //     else
        //     {
        //         mCounter += 1;
        //     }
        // }
    }
    // TRACE("Update finished")
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    mBasementNodes.clear();

        TRACE("in UpdateCellData")
        assert(SPACE_DIM ==3 );
        MAKE_PTR(EmptyBasementMatrix, p_Basement);
        MAKE_PTR(HasEndothelialCell, p_EC);
        c_vector<long double, 3> Node_location;
        
        
        for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node_location = rCellPopulation.GetNode(node_index)->rGetLocation();
            bool HaveSetMembraneParameters=0; 

            cell_iter->SetMutationState(p_EC);

            // Here we need to determine where exatly they 
            // Need to loop over the given boundaries, and see which planes these are defined within
            // std::vector<std::vector<  c_vector<double, SPACE_DIM> > > mBoundaries;
            // PRINT_VARIABLE(mBoundaries.size())
            // PRINT_VARIABLE(mHetro)
            if (mHetro)
            {
                for (unsigned i = 0; i<mBoundaries.size(); i++)
                {
                    std::vector<  c_vector<double, 3> > BasementRegion =  mBoundaries[i];
                    c_vector<double, 3> UpperPlane = BasementRegion[0];
                    c_vector<double, 3> UpperPoint = BasementRegion[1];

                    c_vector<double, 3> LowerPlane = BasementRegion[2];
                    c_vector<double, 3> LowerPoint = BasementRegion[3];

                    // Vector connecting the node to upper plane
                    c_vector<double, 3> NodeToUpperPlane = Node_location - UpperPoint;
                    c_vector<double, 3> NodeToLowerPlane = Node_location - LowerPoint;

                    double DotToUpperPlane = inner_prod(NodeToUpperPlane,UpperPlane );
                    
                    double DotToLowerPlane = inner_prod(NodeToLowerPlane,LowerPlane );

                    // PRINT_2_VARIABLES(DotToUpperPlane,DotToLowerPlane)
                    if ( DotToLowerPlane > 0 && DotToUpperPlane > 0)
                    {
                        cell_iter->SetMutationState(p_Basement);
                        
                        mBasementNodes.push_back(node_index);
                        mSamplebasementNode = node_index;
                        // cell_iter->GetCellData()->SetItem("Mut", 1);

                        cell_iter->GetCellData()->SetItem("ShearModulus",mKbs);
                        cell_iter->GetCellData()->SetItem("AreaDilationModulus", mKba);
                        cell_iter->GetCellData()->SetItem("AreaConstant", mKbA);
                        cell_iter->GetCellData()->SetItem("BendingConstant",0);
                        HaveSetMembraneParameters = 1;
                        break;
                    }
                }
            }

            if( HaveSetMembraneParameters ==0)
            {
                cell_iter->GetCellData()->SetItem("ShearModulus", mGrowthMaps[mStrength](2));
                cell_iter->GetCellData()->SetItem("AreaDilationModulus", mGrowthMaps[mStrength](1));
                cell_iter->GetCellData()->SetItem("AreaConstant", mGrowthMaps[mStrength](0));
                cell_iter->GetCellData()->SetItem("BendingConstant", mGrowthMaps[mStrength](3));
            }
        }   
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::StepChange(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{

//     /*
//         This needs to be properly implimented 
//         -- Temp --
//         Initally the membrane properties are set for 20X Growth, 
//         When the radius gets to 20X inital, the membrane properties in the center will be changed to decreased

//         -- To Impliment -- 
//         Need to adapt this to sweep over all nodes and adapt the membrane properties of the ones that have lost or gained a Potts lattice
// 	*/
//     assert(SPACE_DIM == 3); // Currently assumes that SPACE_DIM = 3
//         // TRACE("StepChange")
//     CellPtr p_Sample_Basement_cell = rCellPopulation.GetCellUsingLocationIndex(mSamplebasementNode);

//     double Step_Kbs = p_Sample_Basement_cell->GetCellData()->GetItem("ShearModulus") + mStepSize * 10;
//     double Step_Kba = p_Sample_Basement_cell->GetCellData()->GetItem("AreaDilationModulus") + mStepSize; // 1e-15;
//     double Step_KbA = p_Sample_Basement_cell->GetCellData()->GetItem("AreaConstant") + mStepSize; // 1e-12;

//     double Kbs = std::min((double)Step_Kbs, (double)mGrowthMaps[1](2));
//     double Kba = std::min((double)Step_Kba, (double)mGrowthMaps[1](1));
//     double KbA = std::min((double)Step_KbA, (double)mGrowthMaps[1](0));

//     if (mAchievedTargetK == 0)
//     {
//         for (std::vector<unsigned>::iterator it = mBasementNodes.begin(); it != mBasementNodes.end(); ++it)
//         {
//             CellPtr cell_iter = rCellPopulation.GetCellUsingLocationIndex(*it);
//             cell_iter->GetCellData()->SetItem("ShearModulus", Kbs);
//             cell_iter->GetCellData()->SetItem("AreaDilationModulus", Kba);
//             cell_iter->GetCellData()->SetItem("AreaConstant", KbA);
//             // cell_iter->GetCellData()->SetItem("BendingConstant", Kbb);
//         }
//     }
//     else if (mAchievedTargetK == 1)
//     {
//         TRACE("Target achieved, shouldnt get here ")
//     }



//   TRACE("StepChange")
    /*
        This needs to be properly implimented 
        -- Temp --
        Initally the membrane properties are set for 20X Growth, 
        When the radius gets to 20X inital, the membrane properties in the center will be changed to decreased

        -- To Impliment -- 
        Need to adapt this to sweep over all nodes and adapt the membrane properties of the ones that have lost or gained a Potts lattice
	*/


    CellPtr p_Sample_Basement_cell = rCellPopulation.GetCellUsingLocationIndex(mSamplebasementNode);
    mKbs = (double)mGrowthMaps[1](2);
    mKba = (double)mGrowthMaps[1](1);
    mKbA = (double)mGrowthMaps[1](0);

        for (std::vector<unsigned>::iterator it = mBasementNodes.begin(); it != mBasementNodes.end(); ++it)
        {
            CellPtr cell_iter = rCellPopulation.GetCellUsingLocationIndex(*it);
            cell_iter->GetCellData()->SetItem("ShearModulus", mKbs);
            cell_iter->GetCellData()->SetItem("AreaDilationModulus", mKba);
            cell_iter->GetCellData()->SetItem("AreaConstant", mKbA);
            //cell_iter->GetCellData()->SetItem("BendingConstant", 0);
        }
    mOn  = 0;
}






template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::UpdateCellData_HillStep(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{

    /*
        This needs to be properly implimented  -- THis gives Hill step 
        -- Temp --
        Initally the membrane properties are set for 20X Growth, 
        When the radius gets to 20X inital, the membrane properties in the center will be changed to decreased

        -- To Impliment -- 
        Need to adapt this to sweep over all nodes and adapt the membrane properties of the ones that have lost or gained a Potts lattice
	*/
    assert(SPACE_DIM == 3); // Currently assumes that SPACE_DIM = 3

    CellPtr p_Sample_Basement_cell = rCellPopulation.GetCellUsingLocationIndex(mSamplebasementNode);

    double Step_Kbs = p_Sample_Basement_cell->GetCellData()->GetItem("ShearModulus") + mStepSize * 10;
    double Step_Kba = p_Sample_Basement_cell->GetCellData()->GetItem("AreaDilationModulus") + mStepSize; // 1e-15;
    double Step_KbA = p_Sample_Basement_cell->GetCellData()->GetItem("AreaConstant") + mStepSize; // 1e-12;

    double Kbs = std::min((double)Step_Kbs, (double)mGrowthMaps[1.2](2));
    double Kba = std::min((double)Step_Kba, (double)mGrowthMaps[1.2](1));
    double KbA = std::min((double)Step_KbA, (double)mGrowthMaps[1.2](0));

    // double Kbb = std::min((double)Step_Kb, (double)mGrowthMaps[1.2](3));
    // PRINT_4_VARIABLES(Kbs, Kba, KbA, Kbb)
    double Ks = mGrowthMaps[mStrength](2);
    double Ka = mGrowthMaps[mStrength](1);
    double KA = mGrowthMaps[mStrength](0);
    // double KB = mGrowthMaps[mStrength](3);

    // H(x) = A/(1+e^-x)+c
    // H(x) = (Km-Ke)(1+e^-X)/(1+e^-x)+Ke

    double X = mMaxZ;
    double As = (Kbs - Ks);
    double Aa = (Kba - Ka);
    double AA = (KbA - KA);

    double Bs = -2 / (3 * X) * log((Kbs - Ks) / (0.01 * Ks));
    double BA = -2 / (3 * X) * log((KbA - KA) / (0.01 * KA));
    double Ba = -2 / (3 * X) * log((Kba - Ka) / (0.01 * Ka));

    double ShearMod;
    double AlphaMod;
    double AreaMod;
    // double BendingMod;

    if (mAchievedTargetK == 0)
    {
        for (std::vector<unsigned>::iterator it = mBasementNodes.begin(); it != mBasementNodes.end(); ++it)
        {
            CellPtr cell_iter = rCellPopulation.GetCellUsingLocationIndex(*it);
            cell_iter->GetCellData()->SetItem("ShearModulus", Kbs);
            cell_iter->GetCellData()->SetItem("AreaDilationModulus", Kba);
            cell_iter->GetCellData()->SetItem("AreaConstant", KbA);
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
                ShearMod = As / (1 + exp(Bs * location[2])) + Ks;
                AlphaMod = Aa / (1 + exp(BA * location[2])) + Ka;
                AreaMod = AA / (1 + exp(Ba * location[2])) + KA;
                // BendingMod = a1B * pow(location[2] - b1, 2) + c1B;
            }
            else
            {
                ShearMod = As / (1 + exp(-Bs * location[2])) + Ks;
                AlphaMod = Aa / (1 + exp(-BA * location[2])) + Ka;
                AreaMod = AA / (1 + exp(-Ba * location[2])) + KA;
                // BendingMod = a2B * pow(location[2] - b2, 2) + c2B;
            }
            // PRINT_3_VARIABLES(ShearMod,AlphaMod, location[2] )
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
            mHetro = 0;
        }
    }
    else if (mAchievedTargetK == 1)
    {
        TRACE("Target achieved, shouldnt get here ")
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class RemeshingTriggerOnHeteroMeshModifier<1, 1>;
template class RemeshingTriggerOnHeteroMeshModifier<1, 2>;
template class RemeshingTriggerOnHeteroMeshModifier<2, 2>;
template class RemeshingTriggerOnHeteroMeshModifier<1, 3>;
template class RemeshingTriggerOnHeteroMeshModifier<2, 3>;
template class RemeshingTriggerOnHeteroMeshModifier<3, 3>;
// Serialization for Boost >= 1.36

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(RemeshingTriggerOnHeteroMeshModifier)



