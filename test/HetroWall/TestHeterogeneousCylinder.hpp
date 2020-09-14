#ifndef TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_CORRECTEDDRAG_HPP_
#define TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_CORRECTEDDRAG_HPP_
//TestCylindricalGrowthDeformableMembrane
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

#include "projects/Jess/src/DeformableModel/MembraneShearForce.hpp"
#include "projects/Jess/src/DeformableModel/MembraneStiffnessForce.hpp"
#include "projects/Jess/src/DeformableModel/MembraneSurfaceForce.hpp"

#include "AppliedForce.hpp"
#include "AppliedForceModifier.hpp"

#include "AppliedPressureModifier.hpp"

static const double M_TIME_FOR_SIMULATION = 4; //50
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
        unsigned N_D = 20;
        unsigned N_Z = N_D * 1.5;
        double Length = 12e-3; //12e-3;
        double trans = -6e-3;
        double MaxZ = Length + trans;
        double MinZ = trans;

        // MutableMesh<2, 3>* p_mesh = p_mesh_base;
        Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, 1.5e-3, Length);
        // Honeycomb3DCylinderMeshGenerator generator(N_D[N_D_index], N_Z, 1.5e-3, Length);
        MutableMesh<2, 3>* p_mesh = generator.GetMesh();
        p_mesh->Translate(trans * unit_vector<double>(3, 2));

        std::stringstream out;
        std::string output_directory = "TestHeterogeneousCylinder"; // + Parameters + "/";

        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        MAKE_PTR(BetaCateninOneHitCellMutationState, p_state); //Mutation to mark nodes

        MAKE_PTR(WildTypeCellMutationState, p_WildTypeState); //Mutation to mark nodes


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
        simulator.SetEndTime(10); //(M_TIME_FOR_SIMULATION);
        simulator.SetDt(0.002);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetUpdateCellPopulationRule(false); // No remeshing.


        /*
        -----------------------------
         Bending Force
        ----------------------------
        */
        double membrane_constant = 2 * 1e-15;
        boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
        p_membrane_force->SetupInitialMembrane(*p_mesh);
        // p_membrane_force->SetMembraneStiffness(membrane_constant);
        simulator.AddForce(p_membrane_force);

        /*
        -----------------------------
        Shearing Force 
        ----------------------------
        */

        double ElasticShearModulus = 0.5 * 1e-4; // n/m 0.5 * 1e-7;
        double AreaDilationModulus = 5e-9;
        bool HetroMembrane = 1;
        double ElasticShearModulusB = ElasticShearModulus /10;
        double AreaDilationModulusB = AreaDilationModulus/1000;

        boost::shared_ptr<MembraneShearForce> p_shear_force(new MembraneShearForce());
        p_shear_force->SetupMembraneConfiguration(cell_population);
        // p_shear_force->SetAreaDilationModulus(AreaDilationModulus);
        // p_shear_force->SetElasticShearModulus(ElasticShearModulus);
        // p_shear_force->SetHetrogeneousMembrane(HetroMembrane, ElasticShearModulusB,  AreaDilationModulusB);
        simulator.AddForce(p_shear_force);

        /*
        -----------------------------
        Surface Area Force
        ----------------------------
        */

        double Area_constant = 5 * 1e-6;
        
        boost::shared_ptr<MembraneSurfaceForce> p_surface_force(new MembraneSurfaceForce());
        p_surface_force->SetupInitialAreas(cell_population);
        // p_surface_force->SetMembraneStiffness(Area_constant);
        simulator.AddForce(p_surface_force);

        /*       
        -----------------------------
        Tractionforce 
        ----------------------------
        */

        std::string working_directory = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/DeformingCylindricalMembraneWithFluidFlow/";

        std::string PressureFile = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/DeformingCylindricalMembraneWithFluidFlow/results/Extracted/surface-pressure.xtr";

        // Create an Applied Force modifier to couple to Flow
        std::string traction_file = working_directory + "results/Extracted/surface-tractions.xtr";
        boost::shared_ptr<AppliedForceModifier<2, 3> > p_force_modifier(new AppliedForceModifier<2, 3>());
        p_force_modifier->SetResetTractionsOnCells(true, traction_file);
        p_force_modifier->SetMembraneConstants(ElasticShearModulus , AreaDilationModulus, Area_constant, membrane_constant );
        simulator.AddSimulationModifier(p_force_modifier);

        boost::shared_ptr<AppliedForce<2, 3> > p_pressure_force(new AppliedForce<2, 3>());
        simulator.AddForce(p_pressure_force);

        //Create a plane boundary to represent the inlet and pass them to the simulation
        c_vector<long double, 3> Boundary1 = Create_c_vector(0, 0, -0.0058);
        c_vector<long double, 3> Normal1 = Create_c_vector(0, 0, 1);

        c_vector<long double, 3> Boundary2 = Create_c_vector(0, 0, 0.0059);
        c_vector<long double, 3> Normal2 = Create_c_vector(0, 0, -1);

        // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, MaxZ * unit_vector<double>(3, 2), unit_vector<double>(3, 2), 10));

        boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary1, Normal1, 10));
        simulator.AddCellPopulationBoundaryCondition(p_condition_1);

        boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary2, Normal2, 10));

        // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, MinZ * unit_vector<double>(3, 2), unit_vector<double>(3, 2), 10));
        simulator.AddCellPopulationBoundaryCondition(p_condition_2);

        simulator.Solve();

        // To reset before looping: this is usually done by the SetUp and TearDown methods
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
    }
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/
