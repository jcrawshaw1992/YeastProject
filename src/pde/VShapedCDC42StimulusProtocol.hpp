#ifndef VSHAPEDCDC42STIMULUSPROTOCOL_HPP
#define VSHAPEDCDC42STIMULUSPROTOCOL_HPP

#include "GradientCDC42StimulusProtocol.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VShapedCDC42StimulusProtocol : public GradientCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>
{
public:
    VShapedCDC42StimulusProtocol(const TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh, double stimulusDuration);

    ~VShapedCDC42StimulusProtocol();

    double GetCDC42ActivationRate(const Node<SPACE_DIM>* pNode, double t) const;
};

#endif
