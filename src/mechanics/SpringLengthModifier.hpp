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

#ifndef SPRINGLENGTHMODIFIER_HPP_
#define SPRINGLENGTHMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"
#include "PetscTools.hpp"
#include "XdrFileReader.hpp"
#include <map>

/**
 * A modifier class which at each simulation time step
 *
 * To be used in conjunction with the AppliedForce Class
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class SpringLengthModifier : public AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>
{
private:

    std::map<std::pair<unsigned, unsigned>, unsigned> mDivisionsApplied;

    /*
     * Store the original rest length
     */
	std::map<std::pair<unsigned, unsigned>, double> mOriginalRestLength;

    /*
     * Store the target rest length proportion.
     */
    std::map<std::pair<unsigned, unsigned>, double> mTargetRestLengthProportion;

	/*
	 * When to shrink perpendicular springs
	 */
	double mShearThreshold;

	 /*
	  * Ammount to reduce length of perpendicular springs by each timestep
	  */
	double mReductionFactor;

	/*
	 * Max number of times to shrink springs.
	 * If to large can cause elements to become degenerate.
	 */
	unsigned mMaxDivisions;

	bool mFirstRun;

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
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM> >(*this);
        archive & mDivisionsApplied;
        archive & mOriginalRestLength;
        archive & mTargetRestLengthProportion;
        archive & mShearThreshold;
        archive & mReductionFactor;
        archive & mMaxDivisions;
        archive & mFirstRun;
    }

public:

    /**
	 * Default constructor.
	 */
    SpringLengthModifier();

    /**
     * Destructor.
     */
    virtual ~SpringLengthModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory);

	double GetShearThreshold();

	void SetShearThreshold(double shearThreshold);

    double GetReductionFactor();

    void SetReductionFactor(double reductionFactor);

	double GetMaxDivisions();

	void SetMaxDivisions(double maxDivisions);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SpringLengthModifier);

#endif /*SPRINGLENGTHMODIFIER_HPP_*/
