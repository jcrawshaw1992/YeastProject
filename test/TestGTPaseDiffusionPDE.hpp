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

#ifndef TESTGTPASEDIFFUSIONPDE_HPP_
#define TESTGTPASEDIFFUSIONPDE_HPP_

#include <cxxtest/TestSuite.h>
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "GTPasePDESystemSolver.hpp"
#include "GradientCDC42StimulusProtocol.hpp"
#include "FlowBasedCDC42StimulusProtocol.hpp"

class TestGTPaseDiffusionPDE : public CxxTest::TestSuite
{
public:
    void TestNewGTPasePDESystemSolver() throw (Exception)
    {
        // Define a 2D simulation
        static const unsigned DIM = 2;

        // Circular mesh
        TrianglesMeshReader<DIM,DIM> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<DIM,DIM> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        // Original disc is centred in (0,0) with R=1. Scale it to be roughly 300um^2 as in Potts simulation
        mesh.Translate(1, 1);
        mesh.Scale(9.772e-6, 9.772e-6);

        // CDC42 stimulus (any class inheriting from AbstractCDC42StimulusProtocol)
        GradientCDC42StimulusProtocol<DIM,DIM> stimulus_protocol(mesh, 10.0);
        //FlowBasedCDC42StimulusProtocol<DIM,DIM> stimulus_protocol(mesh, 10.0, Create_c_vector(-1, 0));

        // PDE system solver
        GTPasePDESystemSolver<DIM, DIM> solver(&mesh, &stimulus_protocol);

        // Set simulation start and end time, and time step
        solver.SetTimes(0.0, 480.0); // 8 minutes as in Fig.2 of Maree 2012
        solver.SetTimeStep(0.01);

        // Specify an output directory and filename prefix for our results file:
        solver.SetOutputDirectoryAndPrefix("GTPasePDESystemSolver","results");

        /*
         * Compilar y ejecutar: scons b=GccOpt_ndebug_4 ts=projects/VascularRemodelling/test/TestGTPaseDiffusionPDE.hpp
         *
         * When an output directory has been specified, the solver writes output in HDF5 format. To
         * convert this to another output format, we call the relevant method. Here, we convert
         * the output to VTK format.
         *
         * In the output directory run:
         *
         *   ~/Documents/SIMULATION/workspace/Chaste/python/utils/AddVtuTimeAnnotations.py results.vtu results_timeseries.vtu
         */
        solver.SetOutputToVtk(true);

        /*
         * For big problems, you may want to specify how often to write the data, telling the
         * solver to output results to file every n-th timestep.
         */
        solver.SetPrintingTimestepMultiple(10);

        /* We are now ready to solve the system. */
        solver.Solve();

        HeartEventHandler::Headings();
        HeartEventHandler::Report();
    }

};

#endif // TESTWRITINGPDESOLVERSTUTORIAL_HPP_
