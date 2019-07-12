
#include "HasEndothelialCell.hpp"

HasEndothelialCell::HasEndothelialCell()
    : AbstractCellMutationState(4)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(HasEndothelialCell)
