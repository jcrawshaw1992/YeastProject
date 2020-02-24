#ifndef TestBendingForceParameterSweep_HPP_
#define TestBendingForceParameterSweep_HPP_

#include <cmath>
#include <cstdio>
#include <ctime>
#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "CellsGenerator.hpp"
#include "Debug.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "projects/VascularRemodelling/src/FixedRegionBoundaryCondition.hpp"

#include "MutableMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "MembraneForcesBasic.hpp"
#include "MembraneForcesBasicCylinder.hpp"
#include "MembranePropertiesSecModifier.hpp"

// #include "RadialForce.hpp"

#include "CellMutationStatesWriter.hpp"
#include "EmptyBasementMatrix.hpp"
#include "HasEndothelialCell.hpp"
#include "LostEndothelialCell.hpp"
#include "OutwardsPressure.hpp"
#include "VtkMeshReader.hpp"

#include "BoundariesModifier.hpp"

#include "projects/VascularRemodelling/src/MembraneForces/MembraneStiffnessForce.hpp"

static const double M_TIME_FOR_SIMULATION = 100; //40; //50
static const double M_SAMPLING_TIME_STEP = 1000; //50
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
    void TestHetroPlexus() throw(Exception)
    {
        std::string mesh_file = "/Users/jcrawshaw/docker-polnet-master/Plexus_Course/Plexus_2.vtu";

        // This data file is in mm??
       
        // std::string output_directory = "CollpasePlexusToInitialConfig/ClosedMesh_AreaSquared/";
        std::string output_directory = "ExampleOfMeshCoarseningWithGrowth/";

        std::string load_dir = "ExampleOfMeshCoarseningWithGrowth";
        double start_time = 360;
        double simulation_time = 3*360;

        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(load_dir, start_time);
        TRACE("After Load");

        p_simulator->SetEndTime(start_time + simulation_time);

        p_simulator->Solve();

        TRACE("Before Save");
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        TRACE("After Save");

        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(start_time);
    }
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/