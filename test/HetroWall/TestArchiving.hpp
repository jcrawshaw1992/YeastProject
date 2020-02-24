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

#include "MembraneForcesBasicCylinder.hpp"
#include "MembranePropertiesModifier.hpp"

#include "CellMutationStatesWriter.hpp"
#include "CellStiffnessWriter.hpp"
#include "EmptyBasementMatrix.hpp"
#include "HasEndothelialCell.hpp"
#include "LostEndothelialCell.hpp"
#include "OutwardsPressure.hpp"

#include "projects/VascularRemodelling/src/MembraneForces/MembraneStiffnessForce.hpp"

static const double M_TIME_FOR_SIMULATION = 60; //40; //50
static const double M_SAMPLING_TIME_STEP = 300; //3000; //50
static const double M_TIME_STEP = 0.01;//  0.001;// 0.0001;
static const double BendingConst = 8e-10;

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
    void TestAllForces() throw(Exception)
    {

        TRACE("Before Load");
        std::string load_dir = "DevelopingCollapsingCylinder/SetUpArchiving/";
        double start_time = 60;

    for (int kB = -160; kB <= -60; kB+=5)
     {
        std::stringstream out;
        out << "KB_" << kB ;
        std::string Parameters = out.str();

        std::string output_directory = "TESTBENDINGFORCE/" +Parameters +"/";
        
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(load_dir, start_time);
        TRACE("After Load");

        p_simulator->SetEndTime(start_time + M_TIME_FOR_SIMULATION);
        p_simulator->SetOutputDirectory(output_directory);
        p_simulator->SetDt(M_TIME_STEP); // M_TIME_STEP
        p_simulator->SetSamplingTimestepMultiple(M_SAMPLING_TIME_STEP);

        /*
        -----------------------------
        Alter the membrane properties modifier; 
        ----------------------------
        */

        std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
            //         KA,          Kalpha           Ks                                 Kb
        GrowthMaps[10] = Create_c_vector(pow(10, -6.9), pow(10, -8.2459), pow(10, -9), 0);
        GrowthMaps[8] = Create_c_vector(pow(10, -6.9), pow(10, -8.0160), pow(10, -9), BendingConst);
        GrowthMaps[6] = Create_c_vector(pow(10, -6.9), pow(10, -7.7300), pow(10, -9), BendingConst);
        GrowthMaps[5] = Create_c_vector(pow(10, -6.9341), pow(10, -7.7), pow(10, -8), BendingConst);
        GrowthMaps[4] = Create_c_vector(pow(10, -6.9), pow(10, -7.4224), pow(10, 8), BendingConst);
        GrowthMaps[2] = Create_c_vector(pow(10, -6.8), pow(10, -6.8124), pow(10, -7), BendingConst);
        GrowthMaps[1.5] = Create_c_vector(pow(10, -6.5), pow(10, -6.3491), pow(10, -7), BendingConst);
        GrowthMaps[1.2] = Create_c_vector(pow(10, -6.2), pow(10, -5.8360), pow(10, -7), BendingConst);

        // SHould be trying to get this now

        boost::shared_ptr<MembranePropertiesModifier<2, 3> > p_Membrane_modifier(new MembranePropertiesModifier<2, 3>());
        p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 10, 1, 1e-9);
        p_simulator->AddSimulationModifier(p_Membrane_modifier);

        /*
        -----------------------------
        Bending Force
        ----------------------------
        */

        double scale = 1e3;
        unsigned N_D = 30;
        unsigned N_Z = N_D * 3;
        double CollapseFactor = 10;

        double Length = 50e-6 * scale; //12e-3; //12e-3
        double Radius = 5e-6 * scale / CollapseFactor;
        double trans = -Length / 2;

        Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
        MutableMesh<2, 3>* p_Inital_mesh = generator.GetMesh();

        p_Inital_mesh->Translate(trans * unit_vector<double>(3, 2));
        double M_BendingConstant = pow(10, kB/10);

        boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
        p_membrane_force->SetupInitialMembrane(*p_Inital_mesh, p_simulator->rGetCellPopulation());
        p_membrane_force->SetMembraneStiffness(M_BendingConstant, N_D, N_Z);
        p_simulator->AddForce(p_membrane_force);

        p_simulator->rGetCellPopulation().AddCellWriter<CellStiffnessWriter>();

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