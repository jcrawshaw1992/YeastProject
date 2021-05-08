/*

Copyright (c) 2005-2020, University of Oxford.
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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */

#ifndef TESTCREATINGANDUSINGANEWSRNMODELTUTORIAL_HPP_
#define TESTCREATINGANDUSINGANEWSRNMODELTUTORIAL_HPP_

/*
 * = An example showing how to create a new subcellular reaction network (SRN) model and use it in a cell-based simulation. =
 *
 * == Introduction ==
 *
 * In the previous cell-based Chaste tutorials, we used existing cell-cycle and SRN models to define how cells
 * proliferate and update and subcellular model. In this tutorial, we show how to create a new SRN model class, and how this
 * can be used in a cell-based simulation.
 *
 * == Including header files ==
 *
 * We begin by including the necessary header files. */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/* The next header includes the Boost shared_ptr smart pointer, and defines some useful
 * macros to save typing when using it. */
#include "SmartPointers.hpp"
/* The next header includes the NEVER_REACHED macro, used in one of the methods below. */
#include "Exception.hpp"

/* The next header defines a base class for ode-based SRN models.
 * Our new SRN model will inherit from this abstract class. */
#include "AbstractOdeSrnModel.hpp"

/* These headers specify the methods to solve the ODE system. */
#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"

/* This header specifies the ODE solvers. */
#include "CellCycleModelOdeSolver.hpp"

/* The following headers are needed for checkpointing. */
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

/* The remaining header files define classes that will be used in the cell-based
 * simulation test. We have encountered each of these header files in previous cell-based Chaste
 * tutorials. */
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "OffLatticeSimulation.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"


//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"
#include "Debug.hpp"


/*
 * This completes the code for {{{MySrnModel}}}. Note that usually this code would
 * be separated out into a separate declaration in a .hpp file and definition in a .cpp file.
 *
 * === The Tests ===
 *
 * We now define the test class, which inherits from {{{AbstractCellBasedTestSuite}}}.
 */
class TestCreatingAndUsingANewSrnModelTutorial : public AbstractCellBasedTestSuite
{
public:


    void TestOffLatticeSimulationWithMySrnModel()
    {
        TRACE("It works");
        /* We use the honeycomb vertex mesh generator to create a vertex mesh.
         */
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        /* Next, we create some cells. First, define the cells vector. */
        std::vector<CellPtr> cells;
        /* We must create a shared_ptr to a {{{CellMutationState}}} with which to bestow the cells.
         * We make use of the macro MAKE_PTR to do this: the first argument is the class and
         * the second argument is the name of the shared_ptr. */
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        /* Then we loop over the nodes. */
        // for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        // {
        //     /* For each node we create a cell with our SRN model and simple Stochastic cell cycle model. */
        //     // UniformG1GenerationalCellCycleModel* p_cell_cycle_model = new UniformG1GenerationalCellCycleModel();
        //     // MySrnModel* p_srn_model = new MySrnModel;


        //     /* We choose to initialise the concentrations to random levels in each cell. */
        //     std::vector<double> initial_conditions;
        //     initial_conditions.push_back(1.0-2.0*RandomNumberGenerator::Instance()->ranf());
        //     initial_conditions.push_back(1.0-2.0*RandomNumberGenerator::Instance()->ranf());
        //     p_srn_model->SetInitialConditions(initial_conditions);

        //     CellPtr p_cell(new Cell(p_state, p_cell_cycle_model, p_srn_model));
        //     p_cell->SetCellProliferativeType(p_stem_type);


        //     /* Now, we define a random birth time, chosen from [-T,0], where
        //      * T is the typical cell cycle duration
        //      */
        //     double birth_time = - RandomNumberGenerator::Instance()->ranf() * p_cell_cycle_model->GetAverageStemCellCycleTime();
        //     /* We then set the birth time and push the cell back into the vector of cells. */
        //     p_cell->SetBirthTime(birth_time);
        //     cells.push_back(p_cell);
        // }

        /* Now that we have defined the mesh and cells, we can define the cell population, forces, areas modifier, and simulation
         * in the same way as the other tutorials. */
        // VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // OffLatticeSimulation<2> simulator(cell_population);
        // simulator.SetOutputDirectory("TestOffLatticeSimulationWithMySrnModel");
        // simulator.SetEndTime(10.0);
        // simulator.SetSamplingTimestepMultiple(50);

        // MAKE_PTR(NagaiHondaForce<2>, p_force);
        // simulator.AddForce(p_force);

        // // MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        // // simulator.AddSimulationModifier(p_growth_modifier);

        // /* Finally to run the simulation, we call {{{Solve()}}}. */
        // simulator.Solve();
    }
    /*
     * To visualize the results, use Paraview. See the UserTutorials/VisualizingWithParaview tutorial for more information
     *
     * Load the file {{{/tmp/$USER/testoutput/TestOffLatticeSimulationWithMySrnModel/results_from_time_0/results.pvd}}},
     * and color by {{{x}}}.
     */
};

#endif /*TESTCREATINGANDUSINGANEWSRNMODELTUTORIAL_HPP_*/
