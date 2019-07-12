#ifndef GRADIENTCDC42STIMULUSPROTOCOL_HPP
#define GRADIENTCDC42STIMULUSPROTOCOL_HPP

#include "AbstractCDC42StimulusProtocol.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class GradientCDC42StimulusProtocol : public AbstractCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>
{
public:
    GradientCDC42StimulusProtocol(const TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh, double stimulusDuration);

    ~GradientCDC42StimulusProtocol();

    virtual double GetCDC42ActivationRate(const Node<SPACE_DIM>* pNode, double t) const;

    double GetCDC42StimulusGradient() const;

protected:
    double mCDC42StimulusGradient;

    double mMeshXMin;

    double mMeshXMax;

    double mMeshXCentre;

private:
    void ComputeMinMaxCoordinates(const TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh);

};

#endif
