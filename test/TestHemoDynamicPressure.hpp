#ifndef TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_CORRECTEDDRAG_HPP_
#define TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_CORRECTEDDRAG_HPP_
//TestCylindricalGrowthDeformableMembrane
#include <cxxtest/TestSuite.h>
#include "HemodynamicPressure.hpp"
#include "HemodynamicTraction.hpp"

#include <cmath>
#include <cstdio>
#include <ctime>
#include "CellsGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"

#include "AppliedForceModifier.hpp"

#include "PressureForce.hpp"

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
#include "FixedRegionBoundaryCondition.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"

#include "VoronoiVertexMeshGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"
#include "VtkMeshWriter.hpp"
#include "AppliedForce.hpp"

// #include "BetaCateninOneHitCellMutationState.hpp"
// // #include "MembraneSurfaceForce.hpp"



#include "MutableMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"


#include "Debug.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "UblasCustomFunctions.hpp"

// #include "PressureForce.hpp"
// #include "AppliedForceModifier.hpp"

static const double M_TIME_FOR_SIMULATION = 0.01; //50
static const double M_SAMPLING_TIME_STEP = 100; //50


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
    void TestSurfaceAreaForceCylinderGrowth() throw(Exception)
    {
            unsigned N_D = 10;
            unsigned N_Z = N_D * 1.5;
            double Length = 12e-3;
            double trans = 0;
 
             Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, 1.5e-3, Length);
             MutableMesh<2, 3>* p_mesh = generator.GetMesh();
             p_mesh->Translate(trans * unit_vector<double>(3, 2));
 
            std::string output_directory = "BuildINgHemeLBTractionForce1/";

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
            simulator.SetDt(0.002);
            simulator.SetSamplingTimestepMultiple(M_SAMPLING_TIME_STEP);
            simulator.SetUpdateCellPopulationRule(false); // No remeshing.

         
        /*       
        -----------------------------
        Radial pressure force 
        // ----------------------------
        // */
            // MAKE_PTR_ARGS(RadialForce, p_radial_force, (Pressure));
            // simulator.AddForce(p_radial_force);

        // Create an Applied Force modifier to couple to Flow
         
            // // boost::shared_ptr<HemodynamicPressure> p_force_modifier(new HemodynamicPressure());
            // MAKE_PTR_ARGS(HemodynamicPressure, p_force_modifier);
            // //     // simulator.AddSimulationModifier(p_force_modifier);


            // boost::shared_ptr<HemodynamicPressure> p_pressure(new HemodynamicPressure());
            //  std::string PressureFile = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/RegularCylinder/results/Extracted/surface-pressure.xtr"; 
            //  p_pressure->LoadPressureFromFile( PressureFile);
            //  simulator.AddForce(p_pressure);

            //  boost::shared_ptr<HemodynamicTraction> p_Traction(new HemodynamicTraction());
             std::string traction_file = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/RegularCylinder/results/Extracted/surface-tractions.xtr";
            //  p_Traction->LoadTractionFromFile(traction_file);
            //  simulator.AddForce(p_Traction);

        //     std::string traction_file = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/RegularCylinder/results/Extracted/surface-tractions.xtr"; 
         
        // boost::shared_ptr<AppliedForceModifier<2,3> > p_force_modifier(new AppliedForceModifier<2,3>());
        // // p_force_modifier->SetResetTractionsOnCells(false,"");
        // // // p_force_modifier->SetEdgeDivisionThreshold(edge_division_threshold);
        // // p_force_modifier->SetResetTractionsOnCells(true, traction_file);
        // simulator.AddSimulationModifier(p_force_modifier);

        // boost::shared_ptr<PressureForce> p_pressure_force(new PressureForce());
        // simulator.AddForce(p_pressure_force);

        boost::shared_ptr<AppliedForceModifier<2,3> > p_force_modifier(new AppliedForceModifier<2,3>());
        p_force_modifier->SetResetTractionsOnCells(true, traction_file);
        simulator.AddSimulationModifier(p_force_modifier);

        boost::shared_ptr<AppliedForce<2,3>> p_pressure_force(new AppliedForce<2,3>());
        simulator.AddForce(p_pressure_force);

        boost::shared_ptr<AppliedForce<2,3>> p_Npressure_force(new AppliedForce<2,3>());
        // simulator.AddForce(p_pressure_force);



            //Create a plane boundary to represent the inlet and pass them to the simulation

            boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition_1(new FixedRegionBoundaryCondition<2,3>(&cell_population,5.0e-3*unit_vector<double>(3,2),unit_vector<double>(3,2),10));
            simulator.AddCellPopulationBoundaryCondition(p_condition_1);

            boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition_2(new FixedRegionBoundaryCondition<2,3>(&cell_population,5.0e-3*unit_vector<double>(3,2),-unit_vector<double>(3,2),10));
            simulator.AddCellPopulationBoundaryCondition(p_condition_2);

            simulator.Solve();

            // To reset before looping: this is usually done by the SetUp and TearDown methods
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
    }
};



/*
            -----------------------------
            Bending Force
            ----------------------------
*/

// double membrane_constant = 1e-12;
// boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
// p_membrane_force->SetupInitialMembrane(*p_mesh);
// p_membrane_force->SetMembraneStiffness(membrane_constant);
// simulator.AddForce(p_membrane_force);

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/

