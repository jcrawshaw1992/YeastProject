
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

#include "MembraneShearForce.hpp"

// #include "EdgeCorrectionForce.hpp"
#include "MembraneSurfaceForce.hpp"

// #include "AppliedForce.hpp"
// #include "AppliedForceModifier.hpp"
#include "ConstantPressure.hpp"

static const double M_TIME_FOR_SIMULATION = 0.01; //40; //50
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
    void TestBendingForceParameterSweep2() throw(Exception)
    {
        double scale = 1e3;
        unsigned N_D = 80;
        unsigned N_Z = N_D * 1.5;
        

        double Length = 20e-6 * scale; //12e-3; //12e-3
        double Radius = 5e-6 * scale / 10;
        double trans = -Length / 2;

        double KA[9] = {pow(10, -10), pow(10, -9.5), pow(10, -9), pow(10, -8.5), pow(10, -8), pow(10, -7.5), pow(10, -7), pow(10, -6.5), pow(10, -6) }; //,160};
        double KAlpha[9] = {pow(10, -10), pow(10, -9.5), pow(10, -9), pow(10, -8.5), pow(10, -8), pow(10, -7.5), pow(10, -7), pow(10, -6.5), pow(10, -6) }; //,160};
        // double Ks[1] = { pow(10, -8)};//, pow(10, -6) };
        double NumberOfSims = 9*9;
        // PRINT_VARIABLE(NumberOfSims );
        for (unsigned KA_index = 0; KA_index < 9; KA_index++)
        {
            
            for (unsigned KAlpha_index = 0; KAlpha_index < 9; KAlpha_index++)
            {
                
                    // PRINT_3_VARIABLES(KA_index,KAlpha_index, Ks_index);
                    double ElasticShearModulus = pow(10, -12);// Ks[Ks_index]; // 104.4e-05;
                    double AreaDilationModulus = KAlpha[KAlpha_index]; // //7.0795e-06  (80% defomation ); // 3.1623e-05 (20% defomation ) ;//0.9e-4;
                    double Area_constant = KA[KA_index]; //6.3096e-06; (80% defomation );// 3.5481e-05;// 0.9e-4;

                    // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
                    // in um will be too large and break chaste without carefull playing with or a tiny time step

                    Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
                    MutableMesh<2, 3>* p_mesh = generator.GetMesh();
                    
                    p_mesh->Translate(trans * unit_vector<double>(3, 2));

                    std::stringstream out;
                    out <<  KAlpha_index << "_"<< KA_index;
                    std::string mesh_size = out.str();
                    std::string output_directory = "MembraneParameterSweep_Ks12/" + mesh_size;


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
                    // simulator.SetEndTime(60);
                    // if (Ks_index ==0)
                    // {
                    //     simulator.SetEndTime(35);
                    // }else if (Ks_index ==1)
                    // {
                         simulator.SetEndTime(50);
                    // }
                    
                 
        
                    simulator.SetDt(0.01); // 0.005
                    simulator.SetSamplingTimestepMultiple(1000);
                    simulator.SetUpdateCellPopulationRule(false); // No remeshing.

        /*
        -----------------------------
        Constant Pressure ballance 
        ----------------------------
        */

                    double P_blood = 0.0021; // Pa ==   1.6004e-05 mmHg
                    double P_tissue = 0.0015; // Pa == 1.1000e-05 mmHg
                    double TransmuralPressure = 6.6715e-04;//P_blood - P_tissue;

                    // double pressure = (133.322 * 16 + 144 / 2) * 1e-12 * scale * scale; // 0.00002;//1.0666e4; // to match 80mmhg
                    MAKE_PTR_ARGS(ConstantPressure, p_ConstantPressure, (TransmuralPressure));
                    simulator.AddForce(p_ConstantPressure);

                    /*
        -----------------------------
        Shearing Force
        ----------------------------
        */

                    //    % 80% DEFORMATION

                    // double ElasticShearModulus = 1e-6;// 104.4e-05;
                    // double AreaDilationModulus =  GrowthMaps[30](1);// //7.0795e-06  (80% defomation ); // 3.1623e-05 (20% defomation ) ;//0.9e-4;
                    // double membrane_constant = 0 ;//0.75e-13;
                    // double Area_constant = GrowthMaps[30](0); //6.3096e-06; (80% defomation );// 3.5481e-05;// 0.9e-4;

                    double Scalling = 1;

                    boost::shared_ptr<MembraneShearForce> p_shear_force(new MembraneShearForce());
                    p_shear_force->SetScallingShear(Scalling);

                    p_shear_force->SetAreaDilationModulus(AreaDilationModulus);

                    p_shear_force->SetElasticShearModulus(ElasticShearModulus);
                    p_shear_force->SetupMembraneConfiguration(cell_population);

                    simulator.AddForce(p_shear_force);

                    /*
        -----------------------------
        Surface Area Force
        ----------------------------
        */

                    // double Area_constant = 0 ;//1 * 1e-11;//  0.5 * 1e-7;

                    boost::shared_ptr<MembraneSurfaceForce> p_surface_force(new MembraneSurfaceForce());
                    p_surface_force->SetScallingArea(Scalling);
                    p_surface_force->SetMembraneStiffness(Area_constant);
                    p_surface_force->SetupInitialAreas(cell_population);

                    simulator.AddForce(p_surface_force);

                    /*
         -----------------------------
          Bending Force
         ----------------------------
//        */
                    // // double membrane_constant = 0 ;//1e-11;
                    // boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
                    //  p_membrane_force->SetScallingBending(Scalling);

                    // p_membrane_force->SetMembraneStiffness(membrane_constant);
                    // p_membrane_force->SetupInitialMembrane(*p_mesh, cell_population);

                    // simulator.AddForce(p_membrane_force);

                    // boost::shared_ptr<EdgeCorrectionForce> p_EdgeCorrectionForce(new EdgeCorrectionForce());
                    // p_EdgeCorrectionForce->SetMeshType(1, N_D, N_Z );
                    // simulator.AddForce(p_EdgeCorrectionForce);

                    // //         /*
                    // //         -----------------------------
                    // //         Boundaries
                    // //         ----------------------------
                    //         */

                    //Create a plane boundary to represent the inlet and pass them to the simulation
                    c_vector<long double, 3> Boundary1 = Create_c_vector(0, 0, -10e-6 * scale);
                    c_vector<long double, 3> Normal1 = Create_c_vector(0, 0, 1);

                    c_vector<long double, 3> Boundary2 = Create_c_vector(0, 0, 10e-6 * scale);
                    c_vector<long double, 3> Normal2 = Create_c_vector(0, 0, -1);

                    boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary1, Normal1, 10));
                    simulator.AddCellPopulationBoundaryCondition(p_condition_1);

                    boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary2, Normal2, 10));
                    simulator.AddCellPopulationBoundaryCondition(p_condition_2);

                    simulator.Solve();

                    // To reset before looping: this is usually done by the SetUp and TearDown methods
                    SimulationTime::Instance()->Destroy();
                    SimulationTime::Instance()->SetStartTime(0.0);
            }
        }
    }
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/

