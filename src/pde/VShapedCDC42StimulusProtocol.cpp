#include "VShapedCDC42StimulusProtocol.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VShapedCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::VShapedCDC42StimulusProtocol(const TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh, double stimulusDuration)
    : GradientCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>(rMesh, stimulusDuration)
{
    this->mCDC42StimulusGradient = 0.025e6;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VShapedCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::~VShapedCDC42StimulusProtocol()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VShapedCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::GetCDC42ActivationRate(const Node<SPACE_DIM>* pNode, double t) const
{
    if (t < this->mStimulusDuration)
    {
        const ChastePoint<SPACE_DIM>& rX = pNode->GetPoint();
        double abs_dist_x_centre = fabs(this->mMeshXCentre - rX[0]);
        bool is_cell_left_half = (rX[0] < this->mMeshXCentre);

        // V-shaped polarisation protocol in Maree 2012 with a bias towards resolving with positive slope
        if (is_cell_left_half)
        {
            return abs_dist_x_centre * this->mCDC42StimulusGradient * 0.8;
        }
        else
        {
            return abs_dist_x_centre * this->mCDC42StimulusGradient;

        }
    }
    else
    {
        return 0.0;
    }
}

// Explicit instantiation
template class VShapedCDC42StimulusProtocol<2, 2>;
