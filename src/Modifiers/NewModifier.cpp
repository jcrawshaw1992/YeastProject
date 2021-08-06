#include "NewModifier.hpp"
#include <algorithm>
#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"
#include <cxxtest/TestSuite.h>
#include <math.h>
#include "Debug.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
NewModifier<ELEMENT_DIM, SPACE_DIM>::NewModifier()
        : AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM>()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
NewModifier<ELEMENT_DIM, SPACE_DIM>::~NewModifier()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NewModifier<ELEMENT_DIM, SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
  
    TRACE("SetupSolve")
    
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NewModifier<ELEMENT_DIM, SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
  
  TRACE("UpdateAtEndOfTimeStep")

}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NewModifier<ELEMENT_DIM, SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    TRACE("UpdateCellData")
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NewModifier<ELEMENT_DIM,SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}


// Explicit instantiation
template class NewModifier<1,1>;
template class NewModifier<1,2>;
template class NewModifier<2,2>;
template class NewModifier<1,3>;
template class NewModifier<2,3>;
template class NewModifier<3,3>;
// Serialization for Boost >= 1.36

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(NewModifier)

