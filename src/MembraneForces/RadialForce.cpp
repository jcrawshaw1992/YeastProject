#include "RadialForce.hpp"

#include "UblasCustomFunctions.hpp"
#include "MathsFunctions.hpp"
#include "HistoryDepMeshBasedCellPopulation.hpp"


RadialForce::RadialForce()
        : AbstractForce<2, 3>()
{
}


void RadialForce::SetPressure(double Pressure)
{
      mStrength =  Pressure;
}



void RadialForce::SetRadiusThreshold(double RadialThreshold)
{
    mRadialThreshold = RadialThreshold;
    mGrowthThreshold =1;
}



void RadialForce::AddForceContribution(AbstractCellPopulation<2, 3>& rCellPopulation)
{
        HistoryDepMeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<HistoryDepMeshBasedCellPopulation<2, 3>*>(&rCellPopulation);        
        std::map<unsigned, c_vector<double, 3> > ForceMap;
        for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
        {
                unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                Node<3>* p_node = rCellPopulation.GetNode(node_index);
                c_vector<double, 3> Position = p_node->rGetLocation();
                
                double Radius = norm_2(Position);
                Position /=Radius;

                if (Radius < mRadialThreshold && mGrowthThreshold ==1)
                {
                    c_vector<double, 3> Force = mStrength *Position; // / norm_2(cell_location);
                
                    rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(Force);
                    cell_iter->GetCellData()->SetItem("RadialForce", norm_2(Force) );
                    cell_iter->GetCellData()->SetItem("RadialForce_X", Force[0] );
                    cell_iter->GetCellData()->SetItem("RadialForce_Y", Force[1] );
                    cell_iter->GetCellData()->SetItem("RadialForce_Z", Force[2] );
                }else if (Radius > mRadialThreshold && mGrowthThreshold ==1)
                {
                    c_vector<double, 3> Force = -mStrength *Position; // / norm_2(cell_location);
                
                    rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(Force);
                    cell_iter->GetCellData()->SetItem("RadialForce", norm_2(Force) );
                    cell_iter->GetCellData()->SetItem("RadialForce_X", Force[0] );
                    cell_iter->GetCellData()->SetItem("RadialForce_Y", Force[1] );
                    cell_iter->GetCellData()->SetItem("RadialForce_Z", Force[2] );
                }
            
                
        }
}
 



void RadialForce::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<2, 3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(RadialForce)
 