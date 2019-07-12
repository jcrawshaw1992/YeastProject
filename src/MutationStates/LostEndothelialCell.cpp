
#include "LostEndothelialCell.hpp"

LostEndothelialCell::LostEndothelialCell()
    : AbstractCellMutationState(2)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(LostEndothelialCell)
