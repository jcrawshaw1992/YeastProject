
#include "EmptyBasementMatrix.hpp"

EmptyBasementMatrix::EmptyBasementMatrix()
    : AbstractCellMutationState(3)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(EmptyBasementMatrix)
