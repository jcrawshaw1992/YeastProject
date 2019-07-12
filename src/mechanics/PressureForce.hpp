#ifndef PRESSUREFORCE_HPP_
#define PRESSUREFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"

// Needed here to avoid serialization errors (on Boost<1.37)
#include "CellLabel.hpp"

/**
 * A pressure force class to model spherically radiating pressure.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class PressureForce  : public AbstractForce<ELEMENT_DIM, SPACE_DIM>
{
friend class TestForces;

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     * Constructor.
     */
    PressureForce();

    /**
     * Destructor.
     */
    ~PressureForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rForces reference to vector of forces on nodes
     * @param rCellPopulation reference to the cell population
     *
     */
    void AddForceContribution(std::vector<c_vector<double, SPACE_DIM> >& rForces,
                              AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(PressureForce)

#endif /*PRESSUREFORCE_HPP_*/
