#ifndef RANDOMACTIVATIONCDC42STIMULUSPROTOCOL_HPP
#define RANDOMACTIVATIONCDC42STIMULUSPROTOCOL_HPP

#include "AbstractCDC42StimulusProtocol.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class RandomActivationCDC42StimulusProtocol : public AbstractCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>
{
public:
    RandomActivationCDC42StimulusProtocol(const TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh, double stimulusDuration);

    virtual ~RandomActivationCDC42StimulusProtocol();

    double GetCDC42ActivationRate(const Node<SPACE_DIM>* pNode, double t) const;

    const std::vector<ChastePoint<SPACE_DIM> >& GetNoiseCentre() const;

private:
    void ComputeNoiseCentres(const TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh);

    double Normal_Noise(double x_AXIS, double y_AXIS, double st_deviation_noise, double t) const;

    std::vector<ChastePoint<SPACE_DIM> > mNoiseCentre;

    static const int NOISE_POINTS_NUMBER = 7;
};

#endif
