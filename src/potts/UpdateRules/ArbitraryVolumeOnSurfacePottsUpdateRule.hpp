/*

Potts volume constrant on a Wrapped potts population 
*/

#ifndef ArbitraryVolumeOnSurfacePottsUpdateRule_HPP_
#define ArbitraryVolumeOnSurfacePottsUpdateRule_HPP_

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
class ArbitraryVolumeOnSurfacePottsUpdateRule : public AbstractWrappedPottsUpdateRule<DIM>
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
    double mMatureCellTargetVolume;

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
        archive & mDeformationEnergyParameter;
        archive & mMatureCellTargetVolume;
    }

public:

    /**
     * Constructor.
     */
    ArbitraryVolumeOnSurfacePottsUpdateRule();

    /**
     * Destructor.
     */
    ~ArbitraryVolumeOnSurfacePottsUpdateRule();

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



    // void CalculateLatticeVolumes(MutableMesh<2, DIM>& rMesh);

    double EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                           unsigned targetNodeIndex,
                                           WrappedPottsBasedCellPopulation<DIM>& rCellPopulation);


    // double GetVolumeOfElement(unsigned pottsElementIndex, PottsElement<DIM>* p_potts_element);

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
     * @return mMatureCellTargetVolume
     */
    double GetMatureCellTargetVolume() const;

    /**
     * Set mMatureCellTargetVolume.
     *
     * @param matureCellTargetVolume the new value of mMatureCellTargetVolume
     */
    void SetMatureCellTargetVolume(double matureCellTargetVolume);

    /**
     * Overridden OutputUpdateRuleParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputUpdateRuleParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ArbitraryVolumeOnSurfacePottsUpdateRule)

#endif /*ArbitraryVolumeOnSurfacePottsUpdateRule_HPP_*/
