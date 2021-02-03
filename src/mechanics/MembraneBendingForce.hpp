/*

Copyright (c) 2005-2015, University of Oxford.
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

#ifndef MembraneBendingForce_HPP_
#define MembraneBendingForce_HPP_

#include "AbstractForce.hpp"
#include "MeshBasedCellPopulation.hpp"

/**
 * Membrane Stiffness Force!
 */
class MembraneBendingForce : public AbstractForce<2, 3>
{
private:

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
        archive & boost::serialization::base_object<AbstractForce<2,3> >(*this);
        archive & mMembraneStiffness;
        // archive & mOriginalAngles;
    }

protected:

    double mMembraneStiffness;

    // std::map<std::pair<unsigned, unsigned>, double> mOriginalAngles;

    // bool CalculateElementNormalsInital(MutableMesh<2, 3>& rMesh, std::pair<Node<3>*, Node<3>*> edge,
    //                              std::pair<c_vector<double, 3>, c_vector<double, 3> >& nonUnitNormals,
    //                              std::pair<Node<3>*,  Node<3>*>& otherNodes,
    //                             MeshBasedCellPopulation<2, 3>* p_cell_population );

    bool CalculateElementNormals(MutableMesh<2, 3>& rMesh, std::pair<Node<3>*, Node<3>*> edge,
                                 std::pair<c_vector<double, 3>, c_vector<double, 3> >& nonUnitNormals,
                                 std::pair<Node<3>*,  Node<3>*>& otherNodes);

public:

    /**
     * Constructor.
     */
    MembraneBendingForce();

    void SetMembraneStiffness(double membraneStiffnes);

    // double GetMembraneStiffness() const;

    // double GetOriginalAngle(std::pair<Node<3>*, Node<3>*> edge);

    // double SetOriginalAngle(std::pair<Node<3>*, Node<3>*> edge, double angle);

    // void SetupInitialMembrane(MutableMesh<2,3>& rMesh, AbstractCellPopulation<2, 3>& rCellPopulation);

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<2,3>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};


// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MembraneBendingForce)

#endif /*MembraneBendingForce_HPP_*/
