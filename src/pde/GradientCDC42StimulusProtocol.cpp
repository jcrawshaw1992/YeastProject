#include "GradientCDC42StimulusProtocol.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
GradientCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::GradientCDC42StimulusProtocol(const TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh, double stimulusDuration)
    : AbstractCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>(rMesh, stimulusDuration), mCDC42StimulusGradient(0.05e6)
{
    ComputeMinMaxCoordinates(rMesh);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
GradientCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::~GradientCDC42StimulusProtocol()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GradientCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::GetCDC42ActivationRate(const Node<SPACE_DIM>* pNode, double t) const
{
    if (t < this->mStimulusDuration)
    {
        const ChastePoint<SPACE_DIM>& rX = pNode->GetPoint();

        double dist_x_min = (rX[0] - mMeshXMin);
        double dist_x_max = (mMeshXMax - rX[0]);
        assert(dist_x_min >= 0);
        assert(dist_x_max >= 0);

        // Standard polarization protocol in Maree 2012
        return dist_x_min * mCDC42StimulusGradient;
    }
    else
    {
        return 0.0;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double GradientCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::GetCDC42StimulusGradient() const
{
    return mCDC42StimulusGradient;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GradientCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::ComputeMinMaxCoordinates(const TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
{
    assert(SPACE_DIM == 2);

    double Xmin = DBL_MAX;
    double Xmax = -DBL_MAX;

    for (unsigned i=0; i<rMesh.GetNumNodes(); i++)
      {
         if (rMesh.GetNode(i)->rGetLocation()[0] < Xmin)
         {
            Xmin = rMesh.GetNode(i)->rGetLocation()[0];
         }
         if (rMesh.GetNode(i)->rGetLocation()[0] > Xmax)
         {
                Xmax = rMesh.GetNode(i)->rGetLocation()[0];
         }
      }

    mMeshXMin = Xmin;
    mMeshXMax = Xmax;
    mMeshXCentre = (mMeshXMax + mMeshXMin) / 2.0;
}

// Explicit instantiation
template class GradientCDC42StimulusProtocol<2, 2>;
