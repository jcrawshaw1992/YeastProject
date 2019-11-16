/*

Abstract update rule for the wrapped Potts population

By Jess

*/

#ifndef ABSTRACTWRAPPEDPOTTSUPDATERULE_HPP_
#define ABSTRACTWRAPPEDPOTTSUPDATERULE_HPP_

// #include "AbstractUpdateRule.hpp"
#include "AbstractUpdateRule.hpp"
#include "WrappedPottsBasedCellPopulation.hpp"


template<unsigned DIM>
class WrappedPottsBasedCellPopulation; // Circular definition <- This is neccessary for some reason, ask James 

/**
 * An abstract Potts update rule class, for use in cellular Potts model simulations.
 */
template<unsigned DIM>
class AbstractWrappedPottsUpdateRule : public AbstractUpdateRule<DIM>
{
    /** Needed for serialization. */
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
        archive & boost::serialization::base_object<AbstractUpdateRule<DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    AbstractWrappedPottsUpdateRule();

    /**
     * Destructor.
     */
    virtual ~AbstractWrappedPottsUpdateRule();

    /**
     * Calculate the contribution to the Hamiltonian.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param currentNodeIndex The index of the current node/lattice site
     * @param targetNodeIndex The index of the target node/lattice site
     * @param rCellPopulation The cell population
     *
     * @return The difference in the Hamiltonian with the configuration of the target node
     * having the same spin as the current node with the current configuration. i.e H_1-H_0
     */
    virtual double EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                   unsigned targetNodeIndex,
                                                   WrappedPottsBasedCellPopulation<DIM>& rCellPopulation)=0;

    /**
     * Overridden OutputUpdateRuleParameters() method.
     *
     * @param rParamsFile a file stream
     */
    virtual void OutputUpdateRuleParameters(out_stream& rParamsFile);
};

#endif /*AbstractUpdateRule_HPP_*/
