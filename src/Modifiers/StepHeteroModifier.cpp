#include "StepHeteroModifier.hpp"
#include <algorithm>
#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"
#include <cxxtest/TestSuite.h>
#include <math.h>
#include "Debug.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
StepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::StepHeteroModifier()
        : AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>()
{
             //AreaConstant           AreaDilationModulus        ShearModulus    
        mGrowthMaps =  { {1, Create_c_vector(pow(10, -7), pow(10, -8.4), pow(10, -8), 1e-11) },
        // {0, Create_c_vector(pow(10, -7), pow(10, -6), pow(10, -5), 1e-10)}
        {0, Create_c_vector(pow(10, -7), pow(10, -6), pow(10, -5), 1e-11)}
    };
    // mGrowthMaps =  { {1, Create_c_vector(pow(10, -6), pow(10, -6), pow(10, -5), 1e-14) },
    //     {0, Create_c_vector(pow(10, -6), pow(10, -5.5), pow(10, -4), 1e-14)}
    // };
    PRINT_VECTOR(mGrowthMaps[0])
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
StepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::~StepHeteroModifier()
{
}


    
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{

    if (mSetUpSolve ==1)
    {
        
      TRACE("SETUPSOLVE")
      UpdateCellData(rCellPopulation);

      for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
        {
                cell_iter->GetCellData()->SetItem("CollpasingRegion", 0);
        }
    mSetUpSolve = 0;

    }
        

      

    
      
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    assert(ELEMENT_DIM ==2 &&  SPACE_DIM == 3);
    HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation = static_cast<HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);
    if (mHetro) // Need to make sure that the stiffer regions get stiffer at the right step size
    {
        //  TRACE("hetero")
        if (mCounter ==mMaxCounter) // TODO this needs setting with a member variable
        {
            StepChange(rCellPopulation);
            mCounter =0;
        }else 
        {
            mCounter+=1;
        }
    }


    double NumberOfIterations = 100;
    if (mSlowIncreaseInMembraneStrength)/// Membrane parameters need to slowly increase :) 
    {
        if ( mSteps < NumberOfIterations +1)
        {
            mStepSize = NumberOfIterations; 
            if (mCounter ==100) // TODO this needs setting with a member variable
            {
                SlowIncreaseInMembraneParameters(rCellPopulation);
                mCounter =0;
                mSteps +=1;
                if (mSteps > 20 && mSteps <31)
                { //TRACE("Hit 0")
                    mSteps +=2;}
                else if(mSteps > 30)
                 {mSteps +=10;
                //  TRACE("Hit A")
                 }
                //  else if(mSteps > 80)
                //  {mSteps +=20;
                // //  TRACE("Hit B")
                //  }
            }else 
            {
                mCounter+=1;
            }
        }

        if (mSteps >= NumberOfIterations)
        {
            mSlowIncreaseInMembraneStrength = 0;
            UpdateCellData(rCellPopulation); // This will set all the membrane constants to their final values :) 
        }
                

    }


}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
        assert(ELEMENT_DIM ==2 &&  SPACE_DIM == 3);
        mBasementNodes.clear();
        // mDistanceToEndothelialRegion.clear();

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
            cell_iter->GetCellData()->SetItem("MembraneState", 0);
            // Need to loop over the given boundaries, and see which planes these are defined within --> std::vector<std::vector<  c_vector<double, 3> > > mBoundaries;
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
    
                    double radius = 1;//0.005; // XXX TODO Radius threshold needs fixing

                    if (DotToLowerPlane >= 0 && DotToUpperPlane >= 0)
                    {
                        if (norm_2(NodeToUpperPlane) <radius ||  norm_2(NodeToLowerPlane)<radius  )
                        {
                            cell_iter->GetCellData()->SetItem("MembraneState", 10);
                            cell_iter->SetMutationState(p_Basement);
                            mBasementNodes.push_back(node_index);
                            mSamplebasementNode = node_index;
                                                                                    // DistanceToUpperPlane,DistanceToLowerPlane
                            // mDistanceToEndothelialRegion[node_index] = Create_c_vector(norm_2(NodeToUpperPlane),norm_2(NodeToLowerPlane));
                            break;
                        }
                        
                    }
                }
            }
           
            if (mSlowIncreaseInMembraneStrength ==1) // Incase there needs to be a slow increase :S 
            {
                cell_iter->GetCellData()->SetItem("ShearModulus", mStartingParameterForSlowIncrease);
                cell_iter->GetCellData()->SetItem("AreaDilationModulus", mStartingParameterForSlowIncrease);
                cell_iter->GetCellData()->SetItem("AreaConstant", mStartingParameterForSlowIncrease);
                cell_iter->GetCellData()->SetItem("BendingConstant", 0);
                 
            }else // Set now
            {    
                cell_iter->GetCellData()->SetItem("ShearModulus", mGrowthMaps[mStrength](2));
                cell_iter->GetCellData()->SetItem("AreaDilationModulus", mGrowthMaps[mStrength](1));
                cell_iter->GetCellData()->SetItem("AreaConstant", mGrowthMaps[mStrength](0));
                cell_iter->GetCellData()->SetItem("BendingConstant", mGrowthMaps[mStrength](3));
            }
             
        }
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetUpdateFrequency(double MaxCounter)
{
    mMaxCounter = MaxCounter ;
    mCounter = MaxCounter-1;
}




template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetStartingParameterForSlowIncrease(double StartingParameterForSlowIncrease)
{
    mStartingParameterForSlowIncrease =StartingParameterForSlowIncrease;
}

/////------------------------------------------------------------------

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetMembranePropeties(std::map<double, c_vector<long double, 4> > GrowthMaps, double Strength)
{
    mGrowthMaps = GrowthMaps;
    mStrength = Strength;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetMembraneStrength(double Strength)
{
    mStrength = Strength;
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetStepSize(double StepSize)
{
    mStepSize = StepSize;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::Boundaries(c_vector<double, 3> UpperPlaneNormal, c_vector<double, 3> UpperPlanePoint, c_vector<double, 3> LowerPlaneNormal, c_vector<double, 3> LowerPlanePoint)
{

    // Upper plane is defined as the one upstream and the lower is downstream
    // Set the boundary planes for this hetro region, set an upper and a lower bound.
    std::vector<  c_vector<double, 3> > CurrentBoundary;
    CurrentBoundary.push_back(UpperPlaneNormal);
    CurrentBoundary.push_back(UpperPlanePoint);

    CurrentBoundary.push_back(LowerPlaneNormal);
    CurrentBoundary.push_back(LowerPlanePoint);

    mBoundaries.push_back(CurrentBoundary);

    mHetro = 1;
    mSetUpSolve =1;
    TRACE("GOt here")
    
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetSlowIncreaseInMembraneStrength(bool SlowIncreaseInMembraneStrength, double TimeStepSize)
{
    mSlowIncreaseInMembraneStrength = SlowIncreaseInMembraneStrength;
    mCounter =0;
    mSteps = 1;
}
  


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::StepChange(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    assert(ELEMENT_DIM ==2 &&  SPACE_DIM == 3);
    /*
        Step Change -- linear change temporally 
	*/
 
    CellPtr p_Sample_Basement_cell = rCellPopulation.GetCellUsingLocationIndex(mSamplebasementNode);
    double Step_Kbs = p_Sample_Basement_cell->GetCellData()->GetItem("ShearModulus") + mStepSize;
    double Step_Kba = p_Sample_Basement_cell->GetCellData()->GetItem("AreaDilationModulus") + mStepSize; // 1e-15;
    double Step_KbA = p_Sample_Basement_cell->GetCellData()->GetItem("AreaConstant") + mStepSize; // 1e-12;
 
    mStepSize *=1.05;

    for (unsigned i = 0; i<mBoundaries.size(); i++)
    {   
        TRACE("Doing A step Change")
        double K_ShearMod         =  mGrowthMaps[0](2);
        double K_AreaDilationMod  =  mGrowthMaps[0](1);
        double K_AreaMod          =  mGrowthMaps[0](0);
              

        K_ShearMod = std::min((double)Step_Kbs, K_ShearMod);
        K_AreaDilationMod = std::min((double)Step_Kba, K_AreaDilationMod);
        K_AreaMod = std::min((double)Step_KbA, K_AreaMod);

        for (std::vector<unsigned>::iterator it = mBasementNodes.begin(); it != mBasementNodes.end(); ++it)
        {
            CellPtr cell_iter = rCellPopulation.GetCellUsingLocationIndex(*it);
           
            cell_iter->GetCellData()->SetItem("ShearModulus", K_ShearMod);
            cell_iter->GetCellData()->SetItem("AreaDilationModulus", K_AreaDilationMod);
            cell_iter->GetCellData()->SetItem("AreaConstant", K_AreaMod);
        }
    }

}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SlowIncreaseInMembraneParameters(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    /*
    Linear increase in membrane properties over time, This is Nothing complicated, this is to stop putting on a huge force at the start of an archieved simulation 
	*/
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
        cell_iter != rCellPopulation.End();
        ++cell_iter)
    {
        cell_iter->GetCellData()->SetItem("AreaConstant", mSteps * mGrowthMaps[mStrength](0)/mStepSize);
        cell_iter->GetCellData()->SetItem("AreaDilationModulus",mSteps *  mGrowthMaps[mStrength](1)/mStepSize);
        cell_iter->GetCellData()->SetItem("ShearModulus", mSteps * mGrowthMaps[mStrength](2)/mStepSize);
        cell_iter->GetCellData()->SetItem("BendingConstant", mSteps * mGrowthMaps[mStrength](3)/mStepSize);

   }
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetBendingForce(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, double BendingConstant)
{
    assert(ELEMENT_DIM ==2 &&  SPACE_DIM == 3);
    // std::set<unsigned> MutantNodeIndices;
    // std::set<unsigned> EdgeMutantNodeIndices;
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        cell_iter->GetCellData()->SetItem("BendingConstant", BendingConstant);
    }
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    assert(ELEMENT_DIM ==2 &&  SPACE_DIM == 3);
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
// template class StepHeteroModifier<1, 1>;
// template class StepHeteroModifier<1, 2>;
// template class StepHeteroModifier<2, 2>;
template class StepHeteroModifier<2, 2>;
template class StepHeteroModifier<2, 3>;
template class StepHeteroModifier<3, 3>;
// Serialization for Boost >= 1.36

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(StepHeteroModifier)





