
#ifndef TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_CORRECTEDDRAG_HPP_
#define TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_CORRECTEDDRAG_HPP_

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <cstdio>
#include <ctime>
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "CellsGenerator.hpp"
#include "Debug.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"

#include "OffLatticeSimulation.hpp"

#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"

#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"

#include "MutableMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "projects/VascularRemodelling/src/MembraneForces/MembraneShearForce.hpp"

#include "projects/VascularRemodelling/src/MembraneForces/MembraneStiffnessForce.hpp"
 #include "MembraneSurfaceForce.hpp"
  #include "EdgeCorrectionForce.hpp"

#include "AppliedForce.hpp"
#include "AppliedForceModifier.hpp"
#include "ConstantPressure.hpp"
#include "PeriodicRectangleMeshGenerator.hpp"

#include "Honeycomb3DMeshGenerator.hpp"






static const double M_TIME_FOR_SIMULATION = 0.01;//40; //50
static const double M_SAMPLING_TIME_STEP = 100; //50
static const double M_TIME_STEP = 0.002;


class SurfaceAreaForceOnCylinder : public AbstractCellBasedTestSuite
{

private:
    double mLastStartTime;
    //double mEndTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime) / (CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:
    void TestAreaFroceDragCorrectedEqui() throw(Exception)
    {
        double scale = 1e-3; 
        
        double N_D = 20;
        double N_Z = 30;
        double Width = 20; //5e-3
        double Length = 30; //30e-3

        // PeriodicRectangleMeshGenerator generator(N_D, N_Z,  Width, Length);
        // MutableMesh<2, 3>* p_mesh = generator.GetMesh();
        // std::vector<unsigned > BoundaryVector = generator.GetBoundaryVector();

        Honeycomb3DMeshGenerator generator(N_D, N_Z,Width , Length);
        MutableMesh<2, 3>* p_mesh = generator.GetMesh();


        std::vector<unsigned> BoundaryVector = generator.GetBoundaryVector();
        

        

        std::stringstream out;

        std::string output_directory = "MakeARectangleMesh/NotPeriodic";

        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        MAKE_PTR(BetaCateninOneHitCellMutationState, p_state); //Mutation to mark nodes

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Create a cell population
        MeshBasedCellPopulation<2, 3> cell_population(*p_mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(0.15); //(M_TIME_FOR_SIMULATION);
        simulator.SetDt(0.1);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetUpdateCellPopulationRule(false); // No remeshing.

        std::vector<unsigned>::iterator Bound_iter = BoundaryVector.begin();

         for (typename MeshBasedCellPopulation<2,3>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned Cell_Index= cell_population.GetLocationIndexUsingCell(*cell_iter);
            cell_iter->GetCellData()->SetItem("CellNumber", Cell_Index);
            cell_iter->GetCellData()->SetItem("BoundaryIndex", *Bound_iter);
            std::advance(Bound_iter, 1);
        }
        simulator.Solve();

        // To reset before looping: this is usually done by the SetUp and TearDown methods
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
    }
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/
