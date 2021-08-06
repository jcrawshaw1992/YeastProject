#include "NewModifier2.hpp"
#include "MeshBasedCellPopulation.hpp"

#include "Debug.hpp"
#include <math.h>


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
NewModifier2<ELEMENT_DIM,SPACE_DIM>::NewModifier2()
    : AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>()
    {
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
NewModifier2<ELEMENT_DIM,SPACE_DIM>::~NewModifier2()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NewModifier2<ELEMENT_DIM,SPACE_DIM>::SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{


}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NewModifier2<ELEMENT_DIM,SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{


}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NewModifier2<ELEMENT_DIM,SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
	
    
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NewModifier2<ELEMENT_DIM,SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}


// Explicit instantiation
template class NewModifier2<1,1>;
template class NewModifier2<1,2>;
template class NewModifier2<2,2>;
template class NewModifier2<1,3>;
template class NewModifier2<2,3>;
template class NewModifier2<3,3>;
// Serialization for Boost >= 1.36

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(NewModifier2)

