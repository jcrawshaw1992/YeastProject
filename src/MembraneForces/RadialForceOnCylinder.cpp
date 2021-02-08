#include "RadialForceOnCylinder.hpp"

#include "UblasCustomFunctions.hpp"
#include "MathsFunctions.hpp"
#include "HistoryDepMeshBasedCellPopulation.hpp"


RadialForceOnCylinder::RadialForceOnCylinder()
        : AbstractForce<2, 3>()
{
}


void RadialForceOnCylinder::SetPressure(double Pressure)
{
      mStrength =  Pressure;
}



void RadialForceOnCylinder::SetRadiusThreshold(double RadialThreshold)
{
    mRadialThreshold = RadialThreshold;
    mGrowthThreshold =1;
}



void RadialForceOnCylinder::AddForceContribution(AbstractCellPopulation<2, 3>& rCellPopulation)
{
     for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
        {
                unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                Node<3>* p_node = rCellPopulation.GetNode(node_index);
                c_vector<double, 3> Position = p_node->rGetLocation();
                Position[2]=0;
                
                double Radius = norm_2(Position);
                Position /=Radius;

                c_vector<double, 3> Force = mStrength *Position; // / norm_2(cell_location);
            
                rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(Force);
                cell_iter->GetCellData()->SetItem("RadialForceOnCylinder", norm_2(Force) );
                cell_iter->GetCellData()->SetItem("RadialForceOnCylinder_X", Force[0] );
                cell_iter->GetCellData()->SetItem("RadialForceOnCylinder_Y", Force[1] );
                cell_iter->GetCellData()->SetItem("RadialForceOnCylinder_Z", Force[2] );
               
        }
}
 



void RadialForceOnCylinder::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<2, 3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(RadialForceOnCylinder)
 