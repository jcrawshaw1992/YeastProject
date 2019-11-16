/*

ArbitraryWrappedAdhesionPottsUpdateRule.hpp 

By Jess

*/

#ifndef ArbitraryWrappedAdhesionPottsUpdateRule_HPP_
#define ArbitraryWrappedAdhesionPottsUpdateRule_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractWrappedPottsUpdateRule.hpp"
#include "WrappedPottsBasedCellPopulation.hpp"

// Needed here to avoid serialization errors (on Boost<1.37)
#include "CellLabel.hpp"

/**
 * An adhesion update rule for use in cell-based simulations
 * using the cellular Potts model.
 */
template<unsigned DIM>
class ArbitraryWrappedAdhesionPottsUpdateRule : public AbstractWrappedPottsUpdateRule<DIM>
{
friend class TestPottsUpdateRules;

private:

    /**
     * Cell-cell adhesion energy parameter.
     */
    double mCellCellAdhesionEnergyParameter;

    /**
     * Cell-boundary adhesion energy parameter.
     */
    double mCellBoundaryAdhesionEnergyParameter;

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
        archive & mCellCellAdhesionEnergyParameter;
        archive & mCellBoundaryAdhesionEnergyParameter;
    }

public:

    /**
     * Constructor.
     */
    ArbitraryWrappedAdhesionPottsUpdateRule();

    /**
     * Destructor.
     */
    virtual ~ArbitraryWrappedAdhesionPottsUpdateRule();

    /**
     * Overridden EvaluateHamiltonianContribution() method
     *
     * Uses  sum_adjacentsites (1-delta(spin(i),spin(j))) gamma(spin(i),spin(j))
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
     * Method to calculate the specific interaction between 2 cells can be overridden in
     * child classes to  implement differential adhesion .etc.
     *
     * @param pCellA pointer to the 1st cell
     * @param pCellB pointer to the 2nd cell
     *
     * @return The cell cell interaction adhesion energy between the two cells
     */
    virtual double GetCellCellAdhesionEnergy(CellPtr pCellA, CellPtr pCellB);

    /**
     * Method to calculate the specific interaction between cell and medium can be overridden in
     * child classes to  implement differential adhesion .etc.
     *
     * @param pCell pointer to the cell
     *
     * @return Cell boundary interaction adhesion energy for the cell
     */
    virtual double GetCellBoundaryAdhesionEnergy(CellPtr pCell);

    /**
      * @return mCellCellAdhesionEnergyParameter
      */
    double GetCellCellAdhesionEnergyParameter();

    /**
     * @return mCellBoundaryAdhesionEnergyParameter
     */
    double GetCellBoundaryAdhesionEnergyParameter();

    /**
     * Set mCellCellAdhesionEnergyParameter.
     *
     * @param cellCellAdhesionEnergyEnergyParameter the new value of mCellCellAdhesionEnergyParameter
     */
    void SetCellCellAdhesionEnergyParameter(double cellCellAdhesionEnergyEnergyParameter);

    /**
     * Set mCellBoundaryAdhesionEnergyParameter.
     *
     * @param cellBoundaryAdhesionEnergyParameter the new value of mCellBoundaryAdhesionEnergyParameter
     */
    void SetCellBoundaryAdhesionEnergyParameter(double cellBoundaryAdhesionEnergyParameter);

    /**
     * Overridden OutputUpdateRuleParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputUpdateRuleParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ArbitraryWrappedAdhesionPottsUpdateRule)

#endif /*ARBITRARYADHESIONUPDATERULE_HPP_*/
