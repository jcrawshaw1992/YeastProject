/*
 * RhoBasedContractionUpdateRule.hpp
 *
 *  Created on: 7 May 2015
 *      Author: mobernabeu
 */

#ifndef RHOBASEDCONTRACTIONUPDATERULE_HPP_
#define RHOBASEDCONTRACTIONUPDATERULE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractPottsUpdateRule.hpp"
#include "PottsBasedCellPopulation.hpp"


/**
 * A Rho-level based constraint update rule class for use in Potts based simulations.
 *
 * the target Rho-level or threshold is constant
 * for each cell over time.
 */

template<unsigned DIM>
class RhoBasedContractionUpdateRule : public AbstractPottsUpdateRule<DIM>
{
friend class TestPottsUpdateRules;

private:

    /**
     * Cell deformation energy parameter.
     * Set to the default value 0.5 in the constructor.
     * \todo provide units
     */
    double mDeformationEnergyParameter;

    /**
     * Non-dimensional target volume of a mature (fully-grown) cell,
     * given in number of lattice sites.
     * Set to the default value 16 in the constructor.
     */
    double mMatureCellTargetRhoContractionForces;

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
        archive & mDeformationEnergyParameter;
        archive & mMatureCellTargetRhoContractionForces;
    }


public:

    RhoBasedContractionUpdateRule();

    ~RhoBasedContractionUpdateRule();

    double GetRhoValueForNode(unsigned nodeIndex, unsigned pottsElementIndex, PottsBasedCellPopulation<DIM>& rCellPopulation);

    /**
     * Calculate the contribution to the Hamiltonian.
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

    /**
     * @return mDeformationEnergyParameter
     */
    double GetDeformationEnergyParameter();

    /**
     * Set mDeformationEnergyParameter.
     *
     * @param deformationEnergyParameter the new value of mDeformationEnergyParameter
     */
    void SetDeformationEnergyParameter(double deformationEnergyParameter);

    /**
     * @return mMatureCellTargetRhoDependentContractionForces
     */
    double GetMatureCellTargetRhoContractionForces() const;

    /**
     * Set mMatureCellTargetRhoContractionForces
     *
     * @param matureCellTargetRhoContractionForces the new value of mMatureCellTargetRhoContractionForces
     */
    void SetMatureCellTargetRhoContractionForces(double matureCellTargetRhoContractionForces);

    /**
     * Output update rule parameters to file.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile a file stream
     */
    void OutputUpdateRuleParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RhoBasedContractionUpdateRule)

#endif /* RHOBASEDCONTRACTIONUPDATERULE_HPP */
