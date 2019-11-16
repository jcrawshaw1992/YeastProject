/*
This code sweeps along the combinations of step sizes to see what is stable
*/

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
#include "MeshBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"

#include "MutableMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "MembraneForcesBasic.hpp"
#include "MembranePropertiesModifier.hpp"

#include "CellMutationStatesWriter.hpp"
#include "EmptyBasementMatrix.hpp"
#include "HasEndothelialCell.hpp"
#include "LostEndothelialCell.hpp"
#include "OutwardsPressure.hpp"

// #include "projects/VascularRemodelling/src/mechanics/MembraneStiffnessForce.hpp"

#include "projects/VascularRemodelling/src/MembraneForces/MembraneStiffnessForce.hpp"

static const double M_TIME_FOR_SIMULATION = 100; //40; //50
static const double M_SAMPLING_TIME_STEP = 100; //50
static const double M_TIME_STEP = 0.002;
// static const double Sampeling;

static const double start_time = 60;
static const std::string load_dir = "DevelopingCollapsingModifier/SetUpArchiving/";


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
    void TestCollapseCode1() throw(Exception)
    {


        // for (int kB = -170; kB < 129; --kB)
        for (int kB = -160; kB <181; kB-=5)
            {
               
                double Dt = 2e-6;
                double Step = 1e-9;
                double duration = 10;
                std::stringstream out;
                out << "KB_" << kB ;
                std::string Parameters = out.str();

                std::string output_directory = "SweepingBendingForceSmallBendingParameters/"+ Parameters;

                OffLatticeSimulation<2,3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2,3>, 3 >::Load(load_dir , start_time);
                p_simulator->SetEndTime(start_time + duration);
                p_simulator->SetOutputDirectory(output_directory);
                p_simulator->SetDt(Dt); // 0.005
                p_simulator->SetSamplingTimestepMultiple(1000);

                double BendingConst =0;
                std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
                GrowthMaps[10] = Create_c_vector(pow(10, -6.9), pow(10, -8.2459), pow(10, -9), 0);// KA, Kalpha  Ks  Kb
                GrowthMaps[1.2] = Create_c_vector(pow(10, -6.2), pow(10, -5.8360), pow(10, -7), BendingConst);



                /*
                -----------------------------
                Alter the membrane properties modifier; 
                ----------------------------
                */
                boost::shared_ptr<MembranePropertiesModifier<2,3> > p_Membrane_modifier(new MembranePropertiesModifier<2,3>());
                p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 10, 1, Step);
                p_simulator->AddSimulationModifier(p_Membrane_modifier);
          

                /*
                -----------------------------
                Bending Force
                ----------------------------
                */

                double scale = 1e3;
                unsigned N_D = 60;
                unsigned N_Z = N_D * 2;
                double CollapseFactor =10;

                double Length = 120e-6 * scale; //12e-3; //12e-3
                double Radius = 5e-6 * scale/CollapseFactor;
                double trans = -Length / 2;

                Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
                MutableMesh<2, 3>* p_Inital_mesh = generator.GetMesh();
                
                p_Inital_mesh->Translate(trans * unit_vector<double>(3, 2));
                double M_BendingConstant = pow(10,kB/10);
                PRINT_2_VARIABLES(kB, M_BendingConstant)
            
                boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
                p_membrane_force->SetupInitialMembrane( *p_Inital_mesh , p_simulator->rGetCellPopulation() );
                p_membrane_force->SetMembraneStiffness(M_BendingConstant,N_D, N_Z );
                p_simulator->AddForce(p_membrane_force);
                p_simulator->Solve();
                TRACE("Before Save");
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
                TRACE("After Save");

                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(60.0);
            }

    }








};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/





// std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
//                 //         KA,          Kalpha           Ks                                 Kb
// GrowthMaps[10] = Create_c_vector(pow(10, -6.9), pow(10, -8.2459), pow(10, -9), 0);
// // GrowthMaps[8] = Create_c_vector(pow(10, -6.9), pow(10, -8.0160), pow(10, -9), BendingConst);
// // GrowthMaps[6] = Create_c_vector(pow(10, -6.9), pow(10, -7.7300), pow(10, -9), BendingConst);
// // GrowthMaps[5] = Create_c_vector(pow(10, -6.9341), pow(10, -7.7), pow(10, -8), BendingConst);
// // GrowthMaps[4] = Create_c_vector(pow(10, -6.9), pow(10, -7.4224), pow(10, 8), BendingConst);
// // GrowthMaps[2] = Create_c_vector(pow(10, -6.8), pow(10, -6.8124), pow(10, -7), BendingConst);
// // GrowthMaps[1.5] = Create_c_vector(pow(10, -6.5), pow(10, -6.3491), pow(10, -7), BendingConst);
// GrowthMaps[1.2] = Create_c_vector(pow(10, -6.2), pow(10, -5.8360), pow(10, -7), BendingConst);
