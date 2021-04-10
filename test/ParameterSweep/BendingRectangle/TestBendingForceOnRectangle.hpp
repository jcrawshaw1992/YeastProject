#ifndef TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_CORRECTEDDRAG_HPP_
#define TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_CORRECTEDDRAG_HPP_
//TestCylindricalGrowthDeformableMembrane
#include <cxxtest/TestSuite.h>

#include "HoneycombMeshGenerator.hpp"
#include <cmath>
#include <cstdio>
#include <ctime>
#include "CellsGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"

#include "OffLatticeSimulation.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "CellIdWriter.hpp"
#include "CellLabel.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CommandLineArguments.hpp"
#include "Debug.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"
#include "VtkMeshWriter.hpp"

#include "MutableMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"

#include "Debug.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "MembraneStiffnessForceForTrianglePair.hpp"
#include "Honeycomb3DMeshGeneratorBentRectangle.hpp"

static const double M_TIME_FOR_SIMULATION = 0.5; //50
static const double M_SAMPLING_TIME_STEP =1; //50
static const double M_TIME_STEP = 0.0002; //50

static const double Pressure = 1.0666e2;
static const double Ka = 20;
static const double Ks = 10;

static const unsigned N_D = 10;
static const unsigned N_Z = 10;

class TestingAreaForce_IncreasingMeshRefinment : public AbstractCellBasedTestSuite
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
    void TestBeindingForceOnRectangle() throw(Exception)
    {
        TRACE("Bending");
        std::string output_directory = "TestBendingForce/";

        double Length = 1.5e-3;
        // HoneycombMeshGenerator generator(20, 30, 1, Length);
        // MutableMesh<2,2>* p_mesh = generator.GetMesh();

        Honeycomb3DMeshGeneratorBentRectangle generator(20, 30, 1e-3, 1e-3);
         MutableMesh<2,3>* p_mesh = generator.GetMesh();
    

        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
      
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Create a cell population
        MeshBasedCellPopulation<2, 3> cell_population(*p_mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);

        /* Radial location writer   */
        // cell_population.AddCellWriter<RadialLocationWriter>();

        cell_population.CalculateRestLengths();

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);
        simulator.SetDt(M_TIME_STEP);
        simulator.SetSamplingTimestepMultiple(M_SAMPLING_TIME_STEP);
        simulator.SetUpdateCellPopulationRule(false); // No remeshing.

        //  -----------------------------
        //  Bending Force
        //  ----------------------------

        double membrane_constant = 1e-10;
        boost::shared_ptr<MembraneStiffnessForceForTrianglePair> p_membrane_force(new MembraneStiffnessForceForTrianglePair());
        p_membrane_force->SetupInitialMembrane(*p_mesh);
        p_membrane_force->SetMembraneStiffness(membrane_constant);
        simulator.AddForce(p_membrane_force);


        // Need to also a hyperelastic force :S 

              // // //Create a plane boundary to represent the inlet and pass them to the simulation
        // // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, MaxZ * unit_vector<double>(3, 2), unit_vector<double>(3, 2), 10));
        // // simulator.AddCellPopulationBoundaryCondition(p_condition_1);

        // // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, MinZ * unit_vector<double>(3, 2), unit_vector<double>(3, 2), 10));
        // // simulator.AddCellPopulationBoundaryCondition(p_condition_2);

        simulator.Solve();

        // To reset before looping: this is usually done by the SetUp and TearDown methods
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
    }
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/
