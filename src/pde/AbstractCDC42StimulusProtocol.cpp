#include "AbstractCDC42StimulusProtocol.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::AbstractCDC42StimulusProtocol(const TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh, double stimulusDuration)
    : mStimulusDuration(stimulusDuration)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::~AbstractCDC42StimulusProtocol()
{
}

// Explicit instantiation
template class AbstractCDC42StimulusProtocol<2, 2>;
