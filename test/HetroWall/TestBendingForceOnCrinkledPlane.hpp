#ifndef TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_CORRECTEDDRAG_HPP_
#define TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_CORRECTEDDRAG_HPP_
//TestCylindricalGrowthDeformableMembrane
#include <cxxtest/TestSuite.h>

// #include "HoneycombMeshGenerator.hpp"

#include <cmath>
#include <cstdio>
#include <ctime>
#include "CellsGenerator.hpp"
// #include "CylindricalHoneycombMeshGenerator.hpp"
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

#include "Debug.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "MembraneStiffnessForce.hpp"
#include "MembraneSurfaceForce.hpp"
#include "MembraneShearForce.hpp"


#include "/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/src/Honeycomb3DMeshGenerator.hpp"

static const double M_TIME_FOR_SIMULATION = 400; //50
static const double M_SAMPLING_TIME_STEP =10000; //50
static const double M_TIME_STEP = 0.0002; //50

static const double Pressure = 1.0666e2;
static const double Ka = 20;
static const double Ks = 10;

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

 void TestBeindingForceOnGenericRectangleStretched() throw(Exception)
    {
        TRACE("Bending");

          unsigned N_D[1] = {40 };//,160};
       
        for (unsigned N_D_index = 0; N_D_index < 1; N_D_index++)
        {
                double N_Z = 40;
                 std::stringstream out;
                out << N_D[N_D_index] << "_" << N_Z;
                std::string mesh_size = out.str();
                std::string output_directory = "TestBendingForce/" + mesh_size;

                Honeycomb3DMeshGenerator generator( N_D[N_D_index], N_Z, 20,20);
                MutableMesh<2, 3>* mesh = generator.GetMesh();

                // Create cells
                MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
            
                std::vector<CellPtr> cells;
                CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
                cells_generator.GenerateBasicRandom(cells, mesh->GetNumNodes(), p_differentiated_type);

                // Create a cell population
                MeshBasedCellPopulation<2, 3> cell_population(*mesh, cells);
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

                double membrane_constant = 10;
                boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
                p_membrane_force->SetupInitialMembrane(*mesh, cell_population);
                p_membrane_force->SetMembraneStiffness(membrane_constant);
                simulator.AddForce(p_membrane_force);


                // //    -----------------------------
                // //    Surface Area Force
                // //    ----------------------------

                // boost::shared_ptr<MembraneSurfaceForce> p_surface_force(new MembraneSurfaceForce());
                // p_surface_force->SetupInitialAreas(cell_population);
                // p_surface_force->SetMembraneStiffness(1);
                // simulator.AddForce(p_surface_force);

                simulator.Solve();

                // To reset before looping: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
            }
    }
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/
