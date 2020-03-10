/*

Copyright (c) 2005-2013, University of Oxford.
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

#ifndef GTPASEPDEMULTIPLEMESHESMODIFIER_HPP_
#define GTPASEPDEMULTIPLEMESHESMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"
#include "TetrahedralSubsetMesh.hpp"
#include "AbstractCDC42StimulusProtocol.hpp"


/**
 * A modifier class in which an Parabolic PDE is solved and the results are stored in CellData.
 */
template<unsigned DIM>
class GTPasePDEMultipleMeshesModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
       archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
        // Note that archiving of mSolution is handled by the methods save/load_construct_data
        archive & mpFeMesh;
        archive & mDependentVariableName;
        archive & mOutputDirectory;
    }

    /** The solution to the PDE problem at the current timestep
     */
    std::vector<Vec> mSolution; //TODO NEED TO ARCHIVE THIS see PdeandBoundaryCondidtion!!!

    /** Pointer to the mesh to solve the pde on **/
    std::vector<TetrahedralSubsetMesh<DIM,DIM>*> mpFeMesh;

    std::vector<TetrahedralSubsetMesh<DIM,DIM>*> mPreviousFeMesh;

    /** For use in PDEs where we know what the quantity for which we are solving is called,
    * e.g. oxygen concentration.
    */
    std::string mDependentVariableName;

    /** Store the output directory name
    */
    std::string mOutputDirectory;

    /** Keeps track of which files were written for each cell in order to create .pvd files at the end of the simulation*/
    std::map<std::string, std::vector<std::string> > mFilesCollection;

public:

    /**
     * Constructor.
     *
     * @param pde the pde to solve
     *
     */
    GTPasePDEMultipleMeshesModifier();

    /**
     * Destructor.
     */
    virtual ~GTPasePDEMultipleMeshesModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specifies what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden UpdateAtEndOfOutputTimeStep() method.
     *
     * Specifies what to do in the simulation at the end of each output timestep.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);


    /**
     * Overridden SetupSolve() method.
     *
     * Specifies what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper Method to Update the mCellPdeElementMap
     *
     * This method should be called before sending the element map to a PDE class
     * to ensure map is up to date.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellPdeElementMap(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Helper method to Setup a FE Mesh based on the cell population
     *
     * @param rCellPopulation reference to the cell population
     */
    void SetupFeMesh(AbstractCellPopulation<DIM,DIM>& rCellPopulation);


    /**
     * Helper method to copy the CellData to the PDE solution
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateSolutionVector(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    bool HasCellChangedShape(unsigned cellNumber);

    Vec InterpolateSolutionBetweenMeshes(TetrahedralSubsetMesh<DIM,DIM>* newMesh, TetrahedralSubsetMesh<DIM,DIM>* oldMesh, Vec previousSolution);

    void GetNodeSetIntersectionAndDifference(TetrahedralSubsetMesh<DIM,DIM>* newMesh, TetrahedralSubsetMesh<DIM,DIM>* oldMesh, std::vector<unsigned>& commonNodes, std::vector<unsigned>& newNodes);

    const std::map<unsigned, unsigned> ComputeGlobalToLocalMap(TetrahedralSubsetMesh<DIM,DIM>* mesh);

    /**
     * Helper method to copy the PDE solution to CellData
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);

    Vec AverageNodeWiseSolutionAtElements(const Vec nodeWiseSolution, TetrahedralSubsetMesh<DIM, DIM>* mesh) const;

    void EnforceMassConservation(TetrahedralSubsetMesh<DIM,DIM>* previousMesh, TetrahedralSubsetMesh<DIM,DIM>* newMesh, const Vec previousSolution, Vec newSolution) const;

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(GTPasePDEMultipleMeshesModifier)

#endif /*GTPASEPDEMULTIPLEMESHESMODIFIER_HPP_*/
