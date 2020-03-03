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
#include "projects/VascularRemodelling/src/FixedRegionBoundaryCondition.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"

#include "MutableMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "MembraneForcesBasicCylinder.hpp"
#include "MembraneForcesBasic.hpp"
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
        double scale = 0.2;//1e-3; // 
        // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
        // in um will be too large and break chaste without carefull playing with or a tiny time step

        std::string mesh_file = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/ShrunkPlexusWithLongInlets/ShrunkPlexusWithLongInlets/config.vtu";

        // This data file is in mm??
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> p_mesh;
        p_mesh.ConstructFromMeshReader(mesh_reader);
        p_mesh.Scale(1.0 * scale, 1.0 * scale, 1.0 * scale); // so distances are back in original scale
        
        std::string output_directory = "ScalledMesh/";

        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        MeshBasedCellPopulation<2, 3> cell_population(p_mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        // cell_population.AddCellWriter<CellMutationStatesWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(0.021);
        simulator.SetDt(0.02); // 0.005
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetUpdateCellPopulationRule(false); // No remeshing.

     
        simulator.Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);

    }


 
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/