#include "PressureForce.hpp"
#include <cmath>


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PressureForce<ELEMENT_DIM, SPACE_DIM>::PressureForce()
    : AbstractForce<ELEMENT_DIM, SPACE_DIM>()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PressureForce<ELEMENT_DIM, SPACE_DIM>::~PressureForce()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PressureForce<ELEMENT_DIM, SPACE_DIM>::AddForceContribution(std::vector<c_vector<double, SPACE_DIM> >& rForces,
                                                 AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{

    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();

    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);



        c_vector<double, SPACE_DIM> location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        c_vector<double, SPACE_DIM> radial_location = location;
        radial_location[SPACE_DIM-1] = 0;

        double height = location[SPACE_DIM-1];
        double radial_distance = norm_2(radial_location);

        double magnitude = 2.0*(tanh(height+0.5-(current_time-7.5)) - tanh(height-0.5-(current_time-7.5)));

        cell_iter->GetCellData()->SetItem("pressure", magnitude);

        if (radial_distance > 1e-4)
        {
            rForces[node_index] += magnitude*radial_location/radial_distance;
        }
        
//        // Ozzys "Pulsing Ball" force
//        if (distance > 1e-4)
//        {
//            rForces[node_index] += (sin(current_time*0.2*M_PI)+1.0) * 0.2 * location/distance;
//        }
        
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PressureForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // No parameters to include

    // Call method on direct parent class
    AbstractForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////
template class PressureForce<2,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(PressureForce)
