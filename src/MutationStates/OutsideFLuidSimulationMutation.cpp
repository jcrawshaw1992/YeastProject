
#include "OutsideFLuidSimulationMutation.hpp"

OutsideFLuidSimulationMutation::OutsideFLuidSimulationMutation()
    : AbstractCellMutationState(2)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(OutsideFLuidSimulationMutation)
