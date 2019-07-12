
#ifndef ArbitraryCurvatureConstraintPottsUpdateRule_HPP_
#define ArbitraryCurvatureConstraintPottsUpdateRule_HPP_


#include "PottsArbitrarySurfaceIn3DMesh.hpp"


#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractPottsUpdateRule.hpp"
#include "PottsBasedCellPopulation.hpp"

// Needed here to avoid serialization errors (on Boost<1.37)
#include "CellLabel.hpp"

/**
 * An area surface constraint update rule class for use in Potts based simulations.
 *
 * Note this currently assumes cells don't grow, i.e the target surface area is constant
 * for each cell over time.
 */
template<unsigned DIM>
class ArbitraryCurvatureConstraintPottsUpdateRule : public AbstractPottsUpdateRule<DIM>
{
friend class TestPottsUpdateRules;

private:

    /**
     * Cell deformation energy parameter.
     * Set to the default value 0.5 in the constructor.
     * \todo provide units
     */
    double mCurvatureEnergyParameter;

    /**
     * Non-dimensional target volume of a mature (fully-grown) cell,
     * given in number of lattice sites.
     * Set to the default value 16 in the constructor.
     */
    double mTargetCurvature;

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

        archive & boost::serialization::base_object<AbstractPottsUpdateRule<DIM> >(*this);

        archive & mCurvatureEnergyParameter;
        archive & mTargetCurvature;

    }

public:

    /**
     * Constructor.
     */
    ArbitraryCurvatureConstraintPottsUpdateRule();

    /**
     * Destructor.
     */
    ~ArbitraryCurvatureConstraintPottsUpdateRule();

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
                                           PottsBasedCellPopulation<DIM>& rCellPopulation);


    /*
     * Loops over the enternal edges and sums the lenght to give the perimeter of the element
     */

    // double GetCurvatureOfElement(unsigned pottsElementIndex, PottsElement<DIM>* pCurrentElement );
   


    /**
     * Set mDeformationEnergyParameter.
     *
     * @param deformationEnergyParameter the new value of mDeformationEnergyParameter
     */
    void SetCurvatureEnergyParameter(double CurvatureEnergyParameter);

    
    
    /**
     * @return mMatureCellTargetCurvature
     */
    double GetTargetCurvature() const;

    /**
     * Set mMatureCellTargetCurvature.
     *
     * @param matureCellTargetCurvature the new value of mMatureCellTargetCurvature
     */
    void SetTargetCurvature(double TargetCurvature);


   /**
     * @return mDeformationEnergyParameter
     */
    double GetCurvatureEnergyParameter();


    /**
     * Overridden OutputUpdateRuleParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputUpdateRuleParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ArbitraryCurvatureConstraintPottsUpdateRule)

#endif /*ARBITRARYCurvatureCONSTRAINTUPDATERULE_HPP_*/
