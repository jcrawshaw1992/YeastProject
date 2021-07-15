
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
#include <cxxtest/TestSuite.h>
#include <math.h>
#include "Debug.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::RemeshingTriggerOnHeteroMeshModifier()
        : AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>()
{
    // //AreaConstant           AreaDilationModulus        ShearModulus    
    //     mGrowthMaps =  {  { 10, Create_c_vector(pow(10, -7), pow(10, -10), pow(10, -8), 1e-14) }, 
    //     { 8, Create_c_vector(pow(10,-7), pow(10,-8),pow(10, -8.5) , 1e-14 ) },
    //     { 6, Create_c_vector(pow(10,-6), pow(10,-10),pow(10, -7.5) , 1e-14 )},
    //     {5, Create_c_vector(pow(10, -7), pow(10, -7.5), pow(10, -8), 1e-14) },
    //     {4, Create_c_vector(pow(10, -7), pow(10, -7.5), pow(10, -7.5), 1e-14) },
    //     {3, Create_c_vector(pow(10, -6), pow(10, -10), pow(10, -6.5), 1e-14) },
    //     {2, Create_c_vector(pow(10, -5.500), pow(10, -7.5), pow(10, -6.5), 1e-14) },
    //     {1, Create_c_vector(pow(10, -5.5000), pow(10, -5.5000), pow(10, -5.5000), 1e-14)}
    // };

                                            //AreaConstant AreaDilation    ShearModulus    
        mGrowthMaps =  { {1, Create_c_vector(pow(10, -6), pow(10, -6), pow(10, -6), 1e-11) }, // Trying to Collapse
        // mGrowthMaps =  { {1, Create_c_vector(pow(10, -7.5), pow(10, -9), pow(10, -8.1), 1e-11) },// INitial state
        // mGrowthMaps =  { {1, Create_c_vector(pow(10, -8), pow(10, -9), pow(10, -9), 1e-11) },// Allowed deformation
        //  mGrowthMaps =  { {1, Create_c_vector(pow(10, -6), pow(10, -7), pow(10, -7.5), 1e-11) },
        {0, Create_c_vector(pow(10, -7), pow(10, -6), pow(10, -5), 1e-11)}// {0, Create_c_vector(pow(10, -7), pow(10, -6), pow(10, -5), 1e-10)}
         };
    PRINT_VECTOR(mGrowthMaps[0])

}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::~RemeshingTriggerOnHeteroMeshModifier()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory)
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


    }
        

      

    
      
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetRemeshingTrigger(std::string RemeshingTrigger)
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
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
  
    assert(ELEMENT_DIM ==2 &&  SPACE_DIM == 3);
    HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation = static_cast<HistoryDepMeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);
    // Here the remeshing happens
    if(mRemeshing)
    {
        // if (mRemeshingTrigger == "AspectRatio")
        // {
        //     // Check all the aspect ratios  and trigger remeshing if needed
        //     std::vector<double> AspectRatioVector = pCellPopulation->MinimumElementAspectRatio();

        //     std::vector<double> quartilesNeeded = { 0.25, 0.5, 0.75};
        //     std::vector<double> quartiles = pCellPopulation->Quantile(AspectRatioVector, quartilesNeeded);
        //     double MinimumAspectRatio = *std::min_element( std::begin(AspectRatioVector), std::end(AspectRatioVector) );

        //     typename std::vector<double>::iterator Iterator = quartiles.begin();
        //     double Quartile1 = *Iterator;

        //     if (mStepsSinceLastRemesh>200) 
        //     {
                
        //         // if (  MinimumAspectRatio <0.2 || Quartile1 < 0.5)
        //         if (  MinimumAspectRatio <0.2 || Quartile1 < 0.5)
        //         {
                    
        //             // TRACE(" AR is too smalle, remeshing now :) ")
        //             pCellPopulation->ExecuteHistoryDependentRemeshing();
        //             UpdateCellData(rCellPopulation);
        //             if (mHetro) // Need to make sure that the stiffer regions get stiffer at the right step size
        //             {
        //                 TRACE("Need to update the membrane strenght for the new mesh, this next method has not yet been written ")
        //                 SetMembraneStrenghtOnNewMesh(rCellPopulation);
        //             }
        //             mStepsSinceLastRemesh =1;
        //          }
                 

        //     }else
        //     {
        //         mStepsSinceLastRemesh +=1;
        //     }
        
                
        // }else
        // {
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
        //  }
    }

    // TRACE("IS hetero ?")
    // Update if hetero 
    if (mHetro) // Need to make sure that the stiffer regions get stiffer at the right step size
    {
        //  TRACE("hetero")
        // TRACE("ABout to change")
        // PRINT_VARIABLE(mCounter)
        // PRINT_2_VARIABLES(mCounter,mMaxCounter )
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


    // I want all forces on the nodes 
    // TS_ASSERT_DELTA(rCellPopulation.GetNode(1)->rGetAppliedForce()[0], -analytical_force_magnitude, 1e-4);
            
      


}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
        assert(ELEMENT_DIM ==2 &&  SPACE_DIM == 3);
        mBasementNodes.clear();
        mDistanceToEndothelialRegion.clear();

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
                            mDistanceToEndothelialRegion[node_index] = Create_c_vector(norm_2(NodeToUpperPlane),norm_2(NodeToLowerPlane));
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
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetPlateauParameters(double a, double B)
{
    ma = a;
    mB = B;
    PRINT_2_VARIABLES(ma,mB)
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetUpdateFrequency(double MaxCounter)
{
    mMaxCounter = MaxCounter ;
    mCounter = MaxCounter-1;
}




template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetStartingParameterForSlowIncrease(double StartingParameterForSlowIncrease)
{
    mStartingParameterForSlowIncrease =StartingParameterForSlowIncrease;
}

/////------------------------------------------------------------------

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetMembranePropeties(std::map<double, c_vector<long double, 4> > GrowthMaps, double Strength, bool Hetrogeneous, double StepSize, double SetupSolve)
{
    mGrowthMaps = GrowthMaps;
    mStrength = Strength;
    mHetro = Hetrogeneous;
    // mStepSize = StepSize;
    mOn =  Hetrogeneous;
    mSetUpSolve = SetupSolve;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetMembraneStrength(double Strength)
{
    //    mGrowthMaps =  {  { 10, Create_c_vector(pow(10, -8), pow(10, -8), pow(10, -10), 1e-14) },
    //     { 8, Create_c_vector(pow(10,-7.5), pow(10,-8.5),pow(10, -10) , 1e-14 ) },
    //     { 6, Create_c_vector(pow(10,-7.5), pow(10,-7.5),pow(10, -6.5) , 1e-14 )},
    //     {5, Create_c_vector(pow(10, -7), pow(10, -10), pow(10, -10), 1e-14) },
    //     {4, Create_c_vector(pow(10, -7), pow(10, -7.5), pow(10, -10), 1e-14) },
    //     {3, Create_c_vector(pow(10, -6.5), pow(10, -10), pow(10, -10), 1e-14) }, 
    //     {2, Create_c_vector(pow(10, -6.5), pow(10, -6.5), pow(10, -10), 1e-14) },
    //     {1, Create_c_vector(pow(10, -5.5000), pow(10, -5.5000), pow(10, -5.5000), 1e-14)}
    // };
    mStrength = Strength;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetMembraneParameters(double AreaParameter, double DilationParameter, double DeformationParamter, double BendingParameter)
{
    mGrowthMaps[0] =  Create_c_vector( AreaParameter, DilationParameter, DeformationParamter,BendingParameter);
    mStrength = 0;
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetThreshold(double Threshold)
{
    // Interval between increasing membrane stiffness 
    mThreshold = Threshold;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::Boundaries(c_vector<double, 3> UpperPlaneNormal, c_vector<double, 3> UpperPlanePoint, c_vector<double, 3> LowerPlaneNormal, c_vector<double, 3> LowerPlanePoint)
{

    // Upper plane is defined as the one upstream and the lower is downstream
    // Set the boundary planes for this hetro region, set an upper and a lower bound.
    std::vector<  c_vector<double, 3> > CurrentBoundary;
    CurrentBoundary.push_back(UpperPlaneNormal);
    CurrentBoundary.push_back(UpperPlanePoint);

    CurrentBoundary.push_back(LowerPlaneNormal);
    CurrentBoundary.push_back(LowerPlanePoint);

    mBoundaries.push_back(CurrentBoundary);

    double Length = norm_2(UpperPlanePoint -LowerPlanePoint );
    /*  M(x) = k/(1+x^2a) */
    mMembraneFuctionSpatialConstants.push_back(PlateauDistributionFuction(Length));
    PRINT_VARIABLE(Length)


    mHetro = 1;
    // mStepSize = 1e-7;//1e-10;
    mOn =  1;
    TRACE("GOt here")
    
}


// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::BoundariesEdge(c_vector<double, 3> UpperPlaneNormal, c_vector<double, 3> UpperPlanePoint, c_vector<double, 3> LowerPlaneNormal, c_vector<double, 3> LowerPlanePoint)
// {

//     // Upper plane is defined as the one upstream and the lower is downstream
//     // Set the boundary planes for this hetro region, set an upper and a lower bound.
//     std::vector<  c_vector<double, 3> > CurrentBoundary;
//     CurrentBoundary.push_back(UpperPlaneNormal);
//     CurrentBoundary.push_back(UpperPlanePoint);

//     CurrentBoundary.push_back(LowerPlaneNormal);
//     CurrentBoundary.push_back(LowerPlanePoint);

//     mBoundaries.push_back(CurrentBoundary);

//     double Length = norm_2(UpperPlanePoint -LowerPlanePoint );
//     /*  M(x) = k/(1+x^2a) */
//     mMembraneFuctionSpatialConstants.push_back(PlateauDistributionFuction(Length));
//     mHetro = 1;
//     mOn =  1;
//     TRACE("GOt here")
    
// }




template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetRemeshingInterval(int RemeshingInterval)
{
    mRemeshing = 1;
    mRemeshingInterval = RemeshingInterval;
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetMembraneStrenghtOnNewMesh(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{

    TRACE("Jess needs to do this!!!!");
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetSlowIncreaseInMembraneStrength(bool SlowIncreaseInMembraneStrength, double TimeStepSize)
{
    mSlowIncreaseInMembraneStrength = SlowIncreaseInMembraneStrength;
    // mTimeStepSize = 1e-6;//TimeStepSize;
    mCounter =0;
    mSteps = 1;
}
  


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::StepChange(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    assert(ELEMENT_DIM ==2 &&  SPACE_DIM == 3);
    /*
        Tempral slowly increasing membrane properties --  increasing k, see next point

        Spatial change follows this funciton  M(x) = k/(1+(bx)^2a)+c, where x is distance along the line connecting the two planes (x=0 is the)

        SpatialFuncitonCoefficents[0]=Create_c_vector(K_ShearMod, K_AreaDilationMod, K_AreaMod);
        SpatialFuncitonCoefficents[1]=Create_c_vector(a,b,0);

        Each point in the basement region can be projected to its nearest point on the connecting midline between the planes. To do this, paramatise this line as 
        X = X0 + t(X1-X0), then subtract the given node from this equation, differentiate with respect to t to find the minimum distance, solve for t and put his back into 
        the original equation. Mut also make the midpoint of this line 0, need to figure out how to do that ... later 
         So closest point on the line is 
         P = X1 -(X2-X1)[(X1-Node).(X2-X1)/|X2-X1|^2]

         X1 and X2 are the center points of the planes. 

         For future im going to need to record which boundary region each basement node is in 

	*/
    double c_bs =  mGrowthMaps[mStrength](2);
    double c_ba =  mGrowthMaps[mStrength](1);
    double c_bA =  mGrowthMaps[mStrength](0);
    PRINT_VARIABLE(mStepSize)
    CellPtr p_Sample_Basement_cell = rCellPopulation.GetCellUsingLocationIndex(mSamplebasementNode);
    double Step_Kbs = p_Sample_Basement_cell->GetCellData()->GetItem("ShearModulus") + mStepSize;
    double Step_Kba = p_Sample_Basement_cell->GetCellData()->GetItem("AreaDilationModulus") + mStepSize; // 1e-15;
    double Step_KbA = p_Sample_Basement_cell->GetCellData()->GetItem("AreaConstant") + mStepSize; // 1e-12;
    // if (mStepSize <1e-4)
    // {
    mStepSize *=1.05;
    // }

    // double a  = mMembraneFuctionSpatialConstants[0][1][0];
    double b  = mMembraneFuctionSpatialConstants[0][1][1];

    // PRINT_2_VARIABLES(a,b)
    PRINT_2_VARIABLES(ma,b)



    for (unsigned i = 0; i<mBoundaries.size(); i++)
    {   
        TRACE("Doing A step Change")
        c_vector<double, 3> X1 = mBoundaries[i][1]; // UpperPoint -> Upstream 
        c_vector<double, 3> X2 = mBoundaries[i][3]; // LowerPoint -> Downstreamc_vector<double, 3>
        c_vector<double, 3> MidPoint = (X1+X2)/2;

        // Just for shear 
        double K_ShearMod         = mMembraneFuctionSpatialConstants[i][0][0]; // This needs temproal considerations with the K 
        double K_AreaDilationMod  = mMembraneFuctionSpatialConstants[i][0][1];
        double K_AreaMod          = mMembraneFuctionSpatialConstants[i][0][2];


        K_ShearMod = std::min((double)Step_Kbs, K_ShearMod);
        K_AreaDilationMod = std::min((double)Step_Kba, K_AreaDilationMod);
        K_AreaMod = std::min((double)Step_KbA, K_AreaMod);

        for (std::vector<unsigned>::iterator it = mBasementNodes.begin(); it != mBasementNodes.end(); ++it)
        {
            CellPtr cell_iter = rCellPopulation.GetCellUsingLocationIndex(*it);
            c_vector<double, 3> Node = rCellPopulation.GetLocationOfCellCentre(cell_iter);
            c_vector<double, 3> P = X1-(X2-X1)*(inner_prod( X1-Node,X2-X1)/inner_prod( X2-X1, X2-X1)); // So now I have a point along the line 

            double X =norm_2(MidPoint - P); // Distance from midpoint, I dont really think the signage matters... symmetry
            
            double ShearModulus = K_ShearMod/(1+pow(b*X,2*ma)) +c_bs ; 
            double AreaDilationModulus = K_AreaDilationMod/(1+pow(b*X,2*ma)) + c_ba ; 
            double AreaModulus = K_AreaMod/(1+pow(b*X,2*ma)) +c_bA; 

            cell_iter->GetCellData()->SetItem("ShearModulus", ShearModulus);
            cell_iter->GetCellData()->SetItem("AreaDilationModulus", AreaDilationModulus);
            cell_iter->GetCellData()->SetItem("AreaConstant", AreaModulus);
            cell_iter->GetCellData()->SetItem("CollpasingRegion", 1);
        }
    }
    

    //Need to put in something to stop increasing membrane stiffness after its too stiff
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SlowIncreaseInMembraneParameters(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
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





template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<c_vector<double, 3>, 2> RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::PlateauDistributionFuction(double Length)
{
    // The spatial smoothing in the changing Variables will spatially follow a PlateauDistributionFuction (temproally linear)
    // I will need to do this for each boundary, so variables in each function will need to be saved in a vector of some sorts, but I will deal with that later 
    // Unforunitly this will need to be done for each variable. the a value sould be the same, but k will change dependening on which variable, which I will be more clear about now 

    // ks will be same for each function, which is good, a will vary on the width of the gap 

     /*  M(x) = k/(1+(bx)^2a) +Me*/ 

    double Basement_ShearMod = mGrowthMaps[mBasementMembraneStrength](2);
    double Basement_AreaDilationMod = mGrowthMaps[mBasementMembraneStrength](1);
    double Basement_AreaMod = mGrowthMaps[mBasementMembraneStrength](0);
    
    double Endo_ShearMod = mGrowthMaps[mStrength](2);
    double Endo_AreaDilationMod = mGrowthMaps[mStrength](1);
    double Endo_AreaMod = mGrowthMaps[mStrength](0);

    double K_ShearMod = Basement_ShearMod- Endo_ShearMod;
    double K_AreaDilationMod = Basement_AreaDilationMod- Endo_AreaDilationMod;
    double K_AreaMod = Basement_AreaMod - Endo_AreaMod;

    // double a = 8; // The larger a is, the steaper the increase is 
    double b = 2/(Length) *std::log(1/0.01 -1)/std::log(ma*2);// Update this to be 2/L --- it is the same thing, just rounded for the thesis :) 

    c_vector<c_vector<double, 3>, 2> SpatialFuncitonCoefficents;
    SpatialFuncitonCoefficents[0]=Create_c_vector(K_ShearMod, K_AreaDilationMod, K_AreaMod);
    SpatialFuncitonCoefficents[1]=Create_c_vector(ma,b,0);
    // PRINT_VARIABLE(ma,b,Length)

    /*  M(x) = k/(1+x^2a) */
    return SpatialFuncitonCoefficents;
}





template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetBasementMembraneStrength(double Strength)
{
  
    mBasementMembraneStrength = Strength;
}

 


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::exec(const char* cmd)
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
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::SetBendingForce(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, double BendingConstant)
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
void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    assert(ELEMENT_DIM ==2 &&  SPACE_DIM == 3);
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
// template class RemeshingTriggerOnHeteroMeshModifier<1, 1>;
// template class RemeshingTriggerOnHeteroMeshModifier<1, 2>;
// template class RemeshingTriggerOnHeteroMeshModifier<2, 2>;
template class RemeshingTriggerOnHeteroMeshModifier<2, 2>;
template class RemeshingTriggerOnHeteroMeshModifier<2, 3>;
template class RemeshingTriggerOnHeteroMeshModifier<3, 3>;
// Serialization for Boost >= 1.36

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(RemeshingTriggerOnHeteroMeshModifier)






// template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
// void RemeshingTriggerOnHeteroMeshModifier<ELEMENT_DIM, SPACE_DIM>::UpdateCellData_HillStep(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
// {

//     /*
//         This needs to be properly implimented  -- THis gives Hill step 
//         -- Temp --
//         Initally the membrane properties are set for 20X Growth, 
//         When the radius gets to 20X inital, the membrane properties in the center will be changed to decreased

//         -- To Impliment -- 
//         Need to adapt this to sweep over all nodes and adapt the membrane properties of the ones that have lost or gained a Potts lattice
// 	*/
//     assert(SPACE_DIM == 3); // Currently assumes that SPACE_DIM = 3

//     CellPtr p_Sample_Basement_cell = rCellPopulation.GetCellUsingLocationIndex(mSamplebasementNode);

//     double Step_Kbs = p_Sample_Basement_cell->GetCellData()->GetItem("ShearModulus") + mStepSize * 10;
//     double Step_Kba = p_Sample_Basement_cell->GetCellData()->GetItem("AreaDilationModulus") + mStepSize; // 1e-15;
//     double Step_KbA = p_Sample_Basement_cell->GetCellData()->GetItem("AreaConstant") + mStepSize; // 1e-12;

//     double Kbs = std::min((double)Step_Kbs, (double)mGrowthMaps[1.2](2));
//     double Kba = std::min((double)Step_Kba, (double)mGrowthMaps[1.2](1));
//     double KbA = std::min((double)Step_KbA, (double)mGrowthMaps[1.2](0));

//     // double Kbb = std::min((double)Step_Kb, (double)mGrowthMaps[1.2](3));
//     // PRINT_4_VARIABLES(Kbs, Kba, KbA, Kbb)
//     double Ks = mGrowthMaps[mStrength](2);
//     double Ka = mGrowthMaps[mStrength](1);
//     double KA = mGrowthMaps[mStrength](0);
//     // double KB = mGrowthMaps[mStrength](3);

//     // H(x) = A/(1+e^-x)+c
//     // H(x) = (Km-Ke)(1+e^-X)/(1+e^-x)+Ke

//     double X = mMaxZ;
//     double As = (Kbs - Ks);
//     double Aa = (Kba - Ka);
//     double AA = (KbA - KA);

//     double Bs = -2 / (3 * X) * log((Kbs - Ks) / (0.01 * Ks));
//     double BA = -2 / (3 * X) * log((KbA - KA) / (0.01 * KA));
//     double Ba = -2 / (3 * X) * log((Kba - Ka) / (0.01 * Ka));

//     double ShearMod;
//     double AlphaMod;
//     double AreaMod;
//     // double BendingMod;

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
//         for (std::vector<unsigned>::iterator it = mNodesNextToBasement.begin(); it != mNodesNextToBasement.end(); ++it)
//         {

//             CellPtr cell_iter = rCellPopulation.GetCellUsingLocationIndex(*it);
//             // Change the membrane parameters of the nodes neighbouring the basement mutants -- do this linearly in time and quadratically in space
//             c_vector<double, SPACE_DIM> location = rCellPopulation.GetLocationOfCellCentre(cell_iter);

//             if (location[2] < 0)
//             {
//                 //K =a(Z -b)^2 +c
//                 ShearMod = As / (1 + exp(Bs * location[2])) + Ks;
//                 AlphaMod = Aa / (1 + exp(BA * location[2])) + Ka;
//                 AreaMod = AA / (1 + exp(Ba * location[2])) + KA;
//                 // BendingMod = a1B * pow(location[2] - b1, 2) + c1B;
//             }
//             else
//             {
//                 ShearMod = As / (1 + exp(-Bs * location[2])) + Ks;
//                 AlphaMod = Aa / (1 + exp(-BA * location[2])) + Ka;
//                 AreaMod = AA / (1 + exp(-Ba * location[2])) + KA;
//                 // BendingMod = a2B * pow(location[2] - b2, 2) + c2B;
//             }
//             // PRINT_3_VARIABLES(ShearMod,AlphaMod, location[2] )
//             cell_iter->GetCellData()->SetItem("ShearModulus", ShearMod);
//             cell_iter->GetCellData()->SetItem("AreaDilationModulus", AlphaMod);
//             cell_iter->GetCellData()->SetItem("AreaConstant", AreaMod);
//             // cell_iter->GetCellData()->SetItem("BendingConstant", BendingMod);
//         }
//         if (Step_KbA > mGrowthMaps[1.2](0) && Step_Kba > mGrowthMaps[1.2](1) && Step_Kbs > mGrowthMaps[1.2](2))
//         {
//             mAchievedTargetK = 1;
//             TRACE("Area, shear and dilation hit limit");
//             TRACE("Simulation can stop very soon :) ");
//             mHetro = 0;
//         }
//     }
//     else if (mAchievedTargetK == 1)
//     {
//         TRACE("Target achieved, shouldnt get here ")
//     }
// }
