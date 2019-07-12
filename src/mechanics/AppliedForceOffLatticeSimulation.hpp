/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef APPLIEDFORCEOFFLATTICESIMULATION_HPP_
#define APPLIEDFORCEOFFLATTICESIMULATION_HPP_

#include <map>
#include "ChasteSerialization.hpp"
#include "OffLatticeSimulation.hpp"
#include "PetscTools.hpp"
#include "XdrFileReader.hpp"

/**
 * Subclass of OffLatticeSimulation in which the tractions applied to the cells are
 * loaded from a file and passed to a CellData structure to be used with the AppliedForce Class.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AppliedForceOffLatticeSimulation : public OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>
{
private :

	/**
	 * 	Store the loaded traction data; this is the applied traction at a given point.
	 */
	std::vector < c_vector<double,3> > mAppliedTractions;

	/**
	 * 	Store the loaded tangential component of the traction data; this is the applied tangential traction at a given point.
	 */
	std::vector < c_vector<double,3> > mAppliedTangentTractions;

	/**
	 * 	Store the loaded traction data; this is the location that the traction is defined at.
	 */
	std::vector < c_vector<double,3> > mAppliedPosition;

	bool mResetTractionsOnCells;

	std::string mTractionFile;
    
        double mEdgeDivisionThreshold;

        std::map<std::pair<unsigned, unsigned>, unsigned> mDivisionsApplied;

    double mShearThreshold;

    unsigned mMaxDivisions;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM> >(*this);
        archive & mResetTractionsOnCells;
        archive & mTractionFile;
        archive & mEdgeDivisionThreshold;
        archive & mDivisionsApplied;
        archive & mShearThreshold;
        archive & mMaxDivisions;
    }

    /**
     * Overridden SetupSolve() method. Calls UpdateCellData().
     */
    void SetupSolve();

    /**
     * Overridden UpdateAtEndOfTimeStep() method. Calls UpdateCellData().
     */
    void UpdateAtEndOfTimeStep();

    /**
     * Load the applied tractions and store these in the CellData.
     */
    void UpdateCellData();

    /**
     * Return the force applied to a point in space.
     *
     * @param location the location in space where the force is to be calculated
     */
    c_vector<double, SPACE_DIM> GetAppliedForce(c_vector<double, SPACE_DIM> location);

    void LoadTractionFromFile();

public:

    /**
     * Default constructor.
     *
     * @param rCellPopulation A cell population facade class (contains a mesh and cells)
     * @param deleteCellPopulationInDestructor Whether to delete the cell population on destruction to
     *     free up memory (defaults to false)
     * @param edgeDivisionThreshold dividion threshold.
     * @param initialiseCells whether to initialise cells (defaults to true, set to false when loading from an archive)
     * @param reloadTractions whether to load tractions from file (defaults to true, set to false when loading from an archive)
     */
     AppliedForceOffLatticeSimulation(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                   std::string tractionFile,
                                   double edgeDivisionThreshold,
                                   double shearThreshold,
                                   bool deleteCellPopulationInDestructor=false,
                                   bool initialiseCells=true);

     /**
      * Destructor.
      */
    ~AppliedForceOffLatticeSimulation();

    /**
     * @return mResetTractionsOnCells
     */
    bool GetResetTractionsOnCells();

    /**
     * Set mResetTractionsOnCells.
     *
     * @param resetTractionsOnCells the new value of mResetTractionsOnCells
     */
    void SetResetTractionsOnCells(bool resetTractionsOnCells, std::string tractionFile="");

};


// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AppliedForceOffLatticeSimulation)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct an AppliedForceOffLatticeSimulation.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const AppliedForceOffLatticeSimulation<ELEMENT_DIM,SPACE_DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise an AppliedForceOffLatticeSimulation.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, AppliedForceOffLatticeSimulation<ELEMENT_DIM,SPACE_DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance, last two variables set extra
    // member variables to be deleted as they are loaded from archive and to not initialise sells.
    ::new(t)AppliedForceOffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>(*p_cell_population, "", true, false);
}
}
} // namespace

#endif /*APPLIEDFORCEOFFLATTICESIMULATION_HPP_*/
