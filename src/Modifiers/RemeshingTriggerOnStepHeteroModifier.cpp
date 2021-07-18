
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

#include "RemeshingTriggerOnStepHeteroModifier.hpp"
#include <algorithm>
#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"
#include <cxxtest/TestSuite.h>
#include <math.h>
#include "Debug.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::RemeshingTriggerOnStepHeteroModifier()
        : AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>()
{
                   //AreaConstant           AreaDilationModulus        ShearModulus    
    mGrowthMaps =  { {1, Create_c_vector(pow(10, -7), pow(10, -8.4), pow(10, -8), 1e-9) },
                    {0.5, Create_c_vector(pow(10, -7), pow(10, -8), pow(10, -8),  1e-7) },
                    // {0, Create_c_vector(pow(10, -7), pow(10, -6), pow(10, -5), 1e-10)}
                      {0, Create_c_vector(pow(10, -7), pow(10, -6), pow(10, -5), 1e-8)}
                     };  
    // mGrowthMaps =  { {1, Create_c_vector(pow(10, -6), pow(10, -6), pow(10, -5), 1e-14) },
    //     {0, Create_c_vector(pow(10, -6), pow(10, -5.5), pow(10, -4), 1e-14)}
    // };
    PRINT_VECTOR(mGrowthMaps[0])
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::~RemeshingTriggerOnStepHeteroModifier()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{

    if (mSetUpSolve ==1)
    {
      TRACE("SETUPSOLVE")
      UpdateCellData(rCellPopulation);
      rCellPopulation.SetBinningRegions();

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
void RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
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


    // Here the remeshing happens
    if(mRemeshing)
    {
        if (mRemeshingTrigger == "AspectRatio")
        {
            // Check all the aspect ratios  and trigger remeshing if needed
            std::vector<double> AspectRatioVector = pCellPopulation->MinimumElementAspectRatio();

            std::vector<double> quartilesNeeded = { 0.25, 0.5, 0.75};
            std::vector<double> quartiles = pCellPopulation->Quantile(AspectRatioVector, quartilesNeeded);
            double MinimumAspectRatio = *std::min_element( std::begin(AspectRatioVector), std::end(AspectRatioVector) );

            typename std::vector<double>::iterator Iterator = quartiles.begin();
            double Quartile1 = *Iterator;

            if (mStepsSinceLastRemesh>1000) 
            {
                
                // if (  MinimumAspectRatio <0.2 || Quartile1 < 0.5)
                if (  MinimumAspectRatio <0.3 || Quartile1 < 0.5)
                {
                    TRACE("REMESHING")
                    PRINT_2_VARIABLES(MinimumAspectRatio,Quartile1 )
                    // TRACE(" AR is too smalle, remeshing now :) ")
                    pCellPopulation->ExecuteHistoryDependentRemeshing();
                    UpdateCellData(rCellPopulation);
                    if (mHetro) // Need to make sure that the stiffer regions get stiffer at the right step size
                    {
                        TRACE("Need to update the membrane strenght for the new mesh, this next method has not yet been written ")
                        SetMembraneStrenghtOnNewMesh(rCellPopulation);
                    }
                    mStepsSinceLastRemesh =1;
                 }
                 

            }else
            {
                mStepsSinceLastRemesh +=1;
            }
        
                
        }else
        {
            if (mExecute>=mRemeshingInterval)
            {
                pCellPopulation->ExecuteHistoryDependentRemeshing();
                UpdateCellData(rCellPopulation);
                if (mHetro) // Need to make sure that the stiffer regions get stiffer at the right step size
                {
                    TRACE("Need to update the membrane strenght for the new mesh, this next method has not yet been written ")
                    SetMembraneStrenghtOnNewMesh(rCellPopulation);
                }
                mExecute = 0;
            } 
            mExecute +=1;
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
void RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
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
                cell_iter->GetCellData()->SetItem("BendingConstant", mGrowthMaps[mStrength](3));
                 
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
void RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetUpdateFrequency(double MaxCounter)
{
    mMaxCounter = MaxCounter ;
    mCounter = MaxCounter-1;
}




template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetStartingParameterForSlowIncrease(double StartingParameterForSlowIncrease)
{
    mStartingParameterForSlowIncrease =StartingParameterForSlowIncrease;
}

/////------------------------------------------------------------------

// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetMembranePropeties(std::map<double, c_vector<long double, 4> > GrowthMaps, double Strength, bool Hetrogeneous, double StepSize, double SetupSolve)
// {
//     mGrowthMaps = GrowthMaps;
//     mStrength = Strength;
//     mHetro = Hetrogeneous;
//     // mStepSize = StepSize;
//     mOn =  Hetrogeneous;
//     mSetUpSolve = SetupSolve;
// }

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetMembranePropeties(std::map<double, c_vector<long double, 4> > GrowthMaps, double Strength)
{
    mGrowthMaps = GrowthMaps;
    mStrength = Strength;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetMembraneStrength(double Strength)
{
    mStrength = Strength;
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetThreshold(double Threshold)
{
    // Interval between increasing membrane stiffness 
    mThreshold = Threshold;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::Boundaries(c_vector<double, 3> UpperPlaneNormal, c_vector<double, 3> UpperPlanePoint, c_vector<double, 3> LowerPlaneNormal, c_vector<double, 3> LowerPlanePoint)
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
void RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetRemeshingInterval(int RemeshingInterval)
{
    mRemeshing = 1;
    mRemeshingInterval = RemeshingInterval;
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetMembraneStrenghtOnNewMesh(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{

    TRACE("Jess needs to do this!!!!");
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetSlowIncreaseInMembraneStrength(bool SlowIncreaseInMembraneStrength, double TimeStepSize)
{
    mSlowIncreaseInMembraneStrength = SlowIncreaseInMembraneStrength;
    // mTimeStepSize = 1e-6;//TimeStepSize;
    mCounter =0;
    mSteps = 1;
}
  


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::StepChange(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
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
void RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SlowIncreaseInMembraneParameters(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
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

        // PRINT_2_VARIABLES(mSteps * mGrowthMaps[mStrength](0)/mStepSize, mGrowthMaps[mStrength](0))
        // PRINT_2_VARIABLES(mSteps * mGrowthMaps[mStrength](1)/mStepSize,mSteps * mGrowthMaps[mStrength](2)/mStepSize) 
    }
}





// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// c_vector<c_vector<double, 3>, 2> RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::PlateauDistributionFuction(double Length)
// {
//     // The spatial smoothing in the changing Variables will spatially follow a PlateauDistributionFuction (temproally linear)
//     // I will need to do this for each boundary, so variables in each function will need to be saved in a vector of some sorts, but I will deal with that later 
//     // Unforunitly this will need to be done for each variable. the a value sould be the same, but k will change dependening on which variable, which I will be more clear about now 

//     // ks will be same for each function, which is good, a will vary on the width of the gap 

//      /*  M(x) = k/(1+(bx)^2a) +Me*/ 

//     double Basement_ShearMod = mGrowthMaps[mBasementMembraneStrength](2);
//     double Basement_AreaDilationMod = mGrowthMaps[mBasementMembraneStrength](1);
//     double Basement_AreaMod = mGrowthMaps[mBasementMembraneStrength](0);
    
//     double Endo_ShearMod = mGrowthMaps[mStrength](2);
//     double Endo_AreaDilationMod = mGrowthMaps[mStrength](1);
//     double Endo_AreaMod = mGrowthMaps[mStrength](0);

//     double K_ShearMod = Basement_ShearMod- Endo_ShearMod;
//     double K_AreaDilationMod = Basement_AreaDilationMod- Endo_AreaDilationMod;
//     double K_AreaMod = Basement_AreaMod - Endo_AreaMod;

//     // double a = 8; // The larger a is, the steaper the increase is 
//     double b = 2/(Length) *std::log(1/0.01 -1)/std::log(ma*2);// Update this to be 2/L --- it is the same thing, just rounded for the thesis :) 

//     c_vector<c_vector<double, 3>, 2> SpatialFuncitonCoefficents;
//     SpatialFuncitonCoefficents[0]=Create_c_vector(K_ShearMod, K_AreaDilationMod, K_AreaMod);
//     SpatialFuncitonCoefficents[1]=Create_c_vector(ma,b,0);
//     // PRINT_VARIABLE(ma,b,Length)

//     /*  M(x) = k/(1+x^2a) */
//     return SpatialFuncitonCoefficents;
// }





template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetBasementMembraneStrength(double Strength)
{
  
    mBasementMembraneStrength = Strength;
}

 


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::exec(const char* cmd)
{
    std::array<char, 128> buffer;
    std::string result;
    int return_code = -1;
    auto pclose_wrapper = [&return_code](FILE* cmd) { return_code = pclose(cmd); };
    { // scope is important, have to make sure the ptr goes out of scope first
        const std::unique_ptr<FILE, decltype(pclose_wrapper)> pipe(popen(cmd, "r"), pclose_wrapper);
        if (pipe)
        {
            while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
            {
                result += buffer.data();
            }
        }
    }
     std::string::size_type sz;     // alias of size_t
     double Answer = std::stod(result,&sz);
    return Answer;
}




template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetBendingForce(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, double BendingConstant)
{
    assert(ELEMENT_DIM ==2 &&  SPACE_DIM == 3);
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
void RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::SetRemeshingTrigger(std::string RemeshingTrigger)
{

    if (RemeshingTrigger == "AspectRatio" ||RemeshingTrigger == "aspectratio" || RemeshingTrigger == "aspectratio" || RemeshingTrigger == "AR" || RemeshingTrigger == "ar")
    {
        mRemeshingTrigger = "AspectRatio";
    }
    else if (RemeshingTrigger == "time"||RemeshingTrigger == "Time"||RemeshingTrigger == "T" || RemeshingTrigger == "t" )
    {
        mRemeshingTrigger = "Time";
    }
    else 
    {
        TRACE("Remeshing trigger must be aspect ratio or time. Please enter correct option")
        TS_ASSERT(1==0);
    }

      
}





template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnStepHeteroModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    assert(ELEMENT_DIM ==2 &&  SPACE_DIM == 3);
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
// template class RemeshingTriggerOnStepHeteroModifier<1, 1>;
// template class RemeshingTriggerOnStepHeteroModifier<1, 2>;
// template class RemeshingTriggerOnStepHeteroModifier<2, 2>;
template class RemeshingTriggerOnStepHeteroModifier<2, 2>;
template class RemeshingTriggerOnStepHeteroModifier<2, 3>;
template class RemeshingTriggerOnStepHeteroModifier<3, 3>;
// Serialization for Boost >= 1.36

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(RemeshingTriggerOnStepHeteroModifier)





