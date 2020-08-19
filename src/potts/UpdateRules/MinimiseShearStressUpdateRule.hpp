/*

Copyright (c) 2005-2014, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef MinimiseShearStressUpdateRule_HPP_
#define MinimiseShearStressUpdateRule_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "PottsBasedCellPopulation.hpp"

#include "PottsArbitrarySurfaceIn3DMesh.hpp"

#include "WrappedPottsBasedCellPopulation.hpp"

#include "AbstractWrappedPottsUpdateRule.hpp"

// Needed here to avoid serialization errors (on Boost<1.37)
#include "CellLabel.hpp"

/**
 * A mechanotaxis update rule class for use in Potts based simulations.
 *
 * Delta H is computed based on information about the traction experienced by a
 * given lattice site stored in PottsArbitrarySurfaceIn3DMesh, hence the rule
 * is not suitable for simulations with PottsMesh.
 *
 */
template<unsigned DIM>
class MinimiseShearStressUpdateRule : public AbstractWrappedPottsUpdateRule<DIM>
{


private:

    /**
     * Parameter to control amount of mechanotaxis
     * Default set in the constructor.
     * \todo provide units
     */
    double mShearMinimisationCorrelationParameter;


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
        archive & mShearMinimisationCorrelationParameter;
    }

public:

    /**
     * Constructor.
     */
    MinimiseShearStressUpdateRule();

    /**
     * Destructor.
     */
    ~MinimiseShearStressUpdateRule();

    /**
     * Overridden EvaluateHamiltonianContribution() method
     *
     * Favours moves along the direction of the shear stress vector.
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
     * @return mTractionCorrelationParameter
     */
    double GetShearMinimisationCorrelationParameter();

    /**
     * Set mTractionCorrelationParameter.
     *
     * @param tractionCorrelationParameter the new value of mTractionCorrelationParameter
     */
    void SetShearMinimisationCorrelationParameter(double ShearMinimisationCorrelationParameter);

    /**
     * Overridden OutputUpdateRuleParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputUpdateRuleParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MinimiseShearStressUpdateRule<3>)

#endif /*MECHANOTAXISPOTTSUPDATERULE_HPP_*/