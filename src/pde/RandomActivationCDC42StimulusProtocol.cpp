#include "RandomActivationCDC42StimulusProtocol.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RandomActivationCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::RandomActivationCDC42StimulusProtocol(const TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh, double stimulusDuration)
    : AbstractCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>(rMesh, stimulusDuration),
      mNoiseCentre(NOISE_POINTS_NUMBER) // This line specifies that mNoiseCentre has NOISE_POINTS_NUMBER entries
{
    ComputeNoiseCentres(rMesh);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
RandomActivationCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::~RandomActivationCDC42StimulusProtocol()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double RandomActivationCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::GetCDC42ActivationRate(const Node<SPACE_DIM>* pNode, double t) const
{
    if (t < this->mStimulusDuration)
    {
        const ChastePoint<SPACE_DIM>& rX = pNode->GetPoint();
        double st_deviation_noise = 0.5e-6;
        double time_fraction =  fmod(t, 1.0);

        if (time_fraction < 0.5)
        {
            return Normal_Noise(rX[0], rX[1], st_deviation_noise, t) * 0.6195;  // NORMAL distribution 21% CDC42 IC 2.95  [previously 0.075e-6]
        }
    }

    return 0.0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double RandomActivationCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::Normal_Noise (double x_AXIS, double y_AXIS, double st_deviation_noise, double t) const
{
    double minVal = 0.0; //(1.0/(st_deviation_noise*sqrt(2.0*M_PI))) * exp(-pow((noise_centre),2)/(2.0*pow(st_deviation_noise, 2))); // Considering that the minimum value of the curve tends to 0.0

    double maxVal = 1.0/(st_deviation_noise*sqrt(2.0*M_PI)); // Maximum value is considered when rU[X] matches noise_centre
    double noise_intensity;
    double noise_intensity_map = 0.0;

    for (int index=0; index<NOISE_POINTS_NUMBER; index++)
    {
        noise_intensity_map += exp(-pow((x_AXIS-mNoiseCentre[index][0]),2)/(2.0*pow(st_deviation_noise, 2))) * exp(-pow((y_AXIS-mNoiseCentre[index][1]),2)/(2.0*pow(st_deviation_noise, 2))) ;
    }

    noise_intensity = (1.0/(st_deviation_noise*sqrt(2.0*M_PI))) * noise_intensity_map;
    cout << "noise_centre_intensity " << noise_intensity << endl;

    double normalized=( noise_intensity - minVal) / (maxVal - minVal);
    cout << "noise_centre_intensity_normalized " << normalized << endl;
    return  normalized;

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void RandomActivationCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::ComputeNoiseCentres(const TetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
{
    double random_node;
    double max_nodes= rMesh.GetNumNodes();
    clock_t clock_time = clock();
    srand (clock_time);
    for (int index=0; index<NOISE_POINTS_NUMBER; index++)
    {
        random_node = rand() % rMesh.GetNumNodes();
        cout << "random nodes "<<random_node << "max_nodes" << max_nodes << endl;
        mNoiseCentre[index] = rMesh.GetNode(random_node)->GetPoint();
        cout << "noise_centre_x " <<mNoiseCentre[index][0]<< "noise_centre_y " << mNoiseCentre[index][1] << endl;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const std::vector<ChastePoint<SPACE_DIM> >& RandomActivationCDC42StimulusProtocol<ELEMENT_DIM, SPACE_DIM>::GetNoiseCentre() const
{
    return mNoiseCentre;
}

// Explicit instantiation
template class RandomActivationCDC42StimulusProtocol<2, 2>;
