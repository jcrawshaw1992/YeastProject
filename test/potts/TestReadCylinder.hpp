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
#ifndef TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_
#define TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"

#include <cmath>
#include <cstdio>
#include <ctime>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellsGenerator.hpp"
#include "OnLatticeSimulation.hpp"
#include "OffLatticeSimulation.hpp"
// #include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "WrappedPottsBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"



#include "AbstractCellBasedTestSuite.hpp"
#include "PottsMeshGenerator.hpp"


#include "Warnings.hpp"
#include "LogFile.hpp"
#include "SmartPointers.hpp"

#include "PottsMeshFromMutableMeshGenerator.hpp"
#include "VtkMeshReader.hpp"
#include "RandomNumberGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"


#include "Debug.hpp"

#include "PottsArbitraryMeshFromMutableMeshGenerator.hpp"

#include "FixedG1GenerationalCellCycleModel.hpp"


#include "PottsArbitraryMeshFromMutableMeshGenerator.hpp"


#include "VtkMeshReader.hpp"
#include "MutableMesh.hpp"


// Needed for NodesOnlyMesh.
#include "PetscSetupAndFinalize.hpp"



class TestPottsOnSurface : public AbstractCellBasedTestSuite
{
public:

    void TestReadCylindricalMesh() throw(Exception)
     {
        VtkMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/SimpleTestCylinder/config.vtu");
    	// VtkMeshReader<2,3> mesh_reader("projects/VascularRemodelling/test/data/bifurcation_cut/config.vtu");
		MutableMesh<2,3> mutable_mesh;
		mutable_mesh.ConstructFromMeshReader(mesh_reader);
     }


};

#endif /*TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_*/
