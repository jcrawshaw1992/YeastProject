#ifndef ABSTRACTCDC42STIMULUSPROTOCOL_HPP
#define ABSTRACTCDC42STIMULUSPROTOCOL_HPP

#include "TetrahedralMesh.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractCDC42StimulusProtocol
{

public:
    AbstractCDC42StimulusProtocol(const TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh, double stimulusDuration);

    virtual ~AbstractCDC42StimulusProtocol();

    virtual double GetCDC42ActivationRate(const Node<SPACE_DIM>* pNode, double t) const = 0;

protected:
    double mStimulusDuration;
};

#endif
