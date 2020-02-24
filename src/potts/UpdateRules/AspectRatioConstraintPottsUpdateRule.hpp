/*


    Wrapped Aspect ratio code 

    There needs to be some relationship between the cell elongation and the flow -- no flow no elongation, high flow, high elongation 
*/

#ifndef ASPECTRATIOCONSTRAINTUPDATERULE_HPP_
#define ASPECTRATIOCONSTRAINTUPDATERULE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractWrappedPottsUpdateRule.hpp"
#include "WrappedPottsBasedCellPopulation.hpp"

// Needed here to avoid serialization errors (on Boost<1.37)
#include "CellLabel.hpp"

/**
 * A volume constraint update rule class for use in Potts based simulations.
 *
 * Note this currently assumes cells don't grow, i.e the target volume is constant
 * for each cell over time.
 */
template<unsigned DIM>
class AspectRatioConstraintPottsUpdateRule : public AbstractWrappedPottsUpdateRule<DIM>
{
friend class TestPottsUpdateRules;

private:

    /**
     * Cell aspect ratio energy parameter.
     * Default set in the constructor.
     * \todo provide units
     */
    double mAspectRatioEnergyParameter;
    double mOrientationParameter;

    /**
     * Non-dimensional target aspect ratio of a mature (fully-grown) cell,
     * given in number of lattice sites.
     * Default set in the constructor.
     */
    double mTargetAspectRatio;

    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractWrappedPottsUpdateRule<DIM> >(*this);
        archive & mAspectRatioEnergyParameter;
        archive & mOrientationParameter;
        archive & mTargetAspectRatio;
    }

public:

    /**
     * Constructor.
     */
    AspectRatioConstraintPottsUpdateRule();

    /**
     * Destructor.
     */
    ~AspectRatioConstraintPottsUpdateRule();

    /**
     * Overridden EvaluateHamiltonianContribution() method
     *
     * Uses sum_elements alpha (V_i - V_i^T)^2.
     *
     * @param currentNodeIndex The index of the current node/lattice site
     * @param targetNodeIndex The index of the target node/lattice site
     * @param rCellPopulation The cell population
     *
     * @return The difference in the Hamiltonian with the configuration of the target node
     * having the same spin as the current node with the current configuration. i.e H_1-H_0
     */
    double EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                           unsigned targetNodeIndex,
                                           WrappedPottsBasedCellPopulation<DIM>& rCellPopulation);

    /**
     * @return mAspectRatioEnergyParameter
     */
    double GetAspectRatioEnergyParameter();

    /**
     * Set mAspectRatioEnergyParameter.
     *
     * @param aspectRatioEnergyParameter the new value of mAspectRatioEnergyParameter
     */
    void SetAspectRatioEnergyParameter(double aspectRatioEnergyParameter);


    double GetOrientationParameter();
    void SetOrientationParameter(double OrientationParameter);




    /**
     * @return mTargetAspectRatio
     */
    double GetTargetAspectRatio() const;

    /**
     * Set mTargetAspectRatio.
     *
     * @param targetAspectRatio the new value of mTargetAspectRatio
     */
    void SetTargetAspectRatio(double targetAspectRatio);

    /**
     * Overridden OutputUpdateRuleParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputUpdateRuleParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AspectRatioConstraintPottsUpdateRule)

#endif /*ASPECTRATIOCONSTRAINTUPDATERULE_HPP_*/
