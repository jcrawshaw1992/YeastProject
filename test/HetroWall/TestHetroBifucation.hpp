
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


#include "MutableMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "VtkMeshReader.hpp"
#include "MembraneForcesBasic.hpp"
#include "MembraneHetroModifier.hpp"

// #include "ConstantPressure.hpp"
#include "EmptyBasementMatrix.hpp"
#include "HasEndothelialCell.hpp"
#include "LostEndothelialCell.hpp"
#include "CellMutationStatesWriter.hpp"

#include "OutwardsPressure.hpp"

#include "projects/VascularRemodelling/src/MembraneForces/MembraneStiffnessForce.hpp"

static const double M_TIME_FOR_SIMULATION = 100; //40; //50
static const double M_SAMPLING_TIME_STEP = 100; //50
static const double M_TIME_STEP = 0.0002;


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
    void TestCollapseCode() throw(Exception)
    {
        double scale = 1e3;
       

            // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
            // in um will be too large and break chaste without carefull playing with or a tiny time step

        // Birfucation
        std::string mesh_file = "projects/VascularRemodelling/test/data/bifurcation_cut/Scalled/configChaste.vtu";

        // This data file is in mm??
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> p_mesh;
        p_mesh.ConstructFromMeshReader(mesh_reader);
        // scale = 1e-3; // so distances are in m
        // p_mesh.Scale(1.0 * scale, 1.0 * scale, 1.0 * scale); // so distances are back in original scal

            std::string output_directory = "DevelopingCollapsingCylinderLongRough/HetroBirfucation/";

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
            simulator.SetEndTime(M_TIME_FOR_SIMULATION);
            simulator.SetDt(M_TIME_STEP); // 0.005
            simulator.SetSamplingTimestepMultiple(M_SAMPLING_TIME_STEP);
            simulator.SetUpdateCellPopulationRule(false); // No remeshing.


        /*
        -----------------------------
        Constant Pressure ballance 
        ----------------------------
        */

        double P_blood = 0.0021;  // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.0015; // Pa == 1.1000e-05 mmHg

        double TransmuralPressure = 6.6715e-4 ;//P_blood - P_tissue;

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure(TransmuralPressure);
        simulator.AddForce(p_ForceOut);



        /*
        -----------------------------
        MembraneProperties Modifier
        ----------------------------
        */
       double bending =1e-8;

        std::map<double, c_vector<long double, 4> > GrowthMaps;  // From matlab sweep results 
                                    //         KA,          Kalpha           Ks
        GrowthMaps[10] = Create_c_vector(pow(10,-6.9), pow(10,-8.2459), pow(10, -9), bending );
        GrowthMaps[8] = Create_c_vector(pow(10,-6.9), pow(10,-8.0160),  pow(10, -9) , bending );
        GrowthMaps[6] = Create_c_vector(pow(10,-6.9), pow(10,-7.7300),  pow(10, -9) , bending );


        // GrowthMaps[5] = Create_c_vector(pow(10,-6.9341)*1e-3, pow(10,-7.7)*1e-2,pow(10, -8) *1e-2, bending );
        
         GrowthMaps[5] = Create_c_vector(pow(10,-6.9341), pow(10,-7.7),pow(10, -8), bending );
       GrowthMaps[4] = Create_c_vector(pow(10,-6.9), pow(10,-7.4224),pow(10, 8), bending  );

        GrowthMaps[2] = Create_c_vector(pow(10,-6.8), pow(10,-6.8124),pow(10, -7) , bending );
        GrowthMaps[1.5] = Create_c_vector(pow(10,-6.5), pow(10,-6.3491),pow(10, -7) , bending );
        // GrowthMaps[1.2] =  Create_c_vector(pow(10,-6.2)*1e-3, pow(10,-5.8360)*1e-2,pow(10, -7)*1e-2, bending  );
GrowthMaps[1.2] =  Create_c_vector(pow(10,-6.2), pow(10,-5.8360),pow(10, -7), bending  );

        boost::shared_ptr<MembraneHetroModifier<2,3> > p_Membrane_modifier(new MembraneHetroModifier<2,3>());
        p_Membrane_modifier->SetMembranePropeties(GrowthMaps,2,1);
        p_Membrane_modifier->SetupSolve(cell_population, output_directory);
        simulator.AddSimulationModifier(p_Membrane_modifier);

        // /*
        // -----------------------------
        // SMembrane forces
        // ----------------------------
        // */

        boost::shared_ptr<MembraneForcesBasic> p_shear_force(new MembraneForcesBasic());
        p_shear_force->SetupMembraneConfiguration(cell_population);
        simulator.AddForce(p_shear_force);
    
         /*
        -----------------------------
        Bending Force
        ----------------------------
        */
        // // double membrane_constant = 0 ;//1e-11;
        // boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());

        // p_membrane_force->SetMembraneStiffness(bending);
        // p_membrane_force->SetupInitialMembrane(p_mesh);

        // simulator.AddForce(p_membrane_force);

        // boost::shared_ptr<EdgeCorrectionForce> p_EdgeCorrectionForce(new EdgeCorrectionForce());
        // p_EdgeCorrectionForce->SetMeshType(1, N_D, N_Z );
        // simulator.AddForce(p_EdgeCorrectionForce);

         /*
        -----------------------------
        Boundary conditions
        ----------------------------
        */

        // //Create a plane boundary to represent the inlet and pass them to the simulation

        c_vector<long double, 3> Boundary1 = Create_c_vector( 0.05605755520275098,  0.15575785260332847, 0.1866408328303961) ;//    c_vector<long double, 3> Boundary1 = Create_c_vector(12.48824184355049, 34.75302061558864, 41.78949821113195);
         c_vector<long double, 3> Normal1 = -Create_c_vector( -0.09536117814727178, 0.6421813385077393, 0.7605980371883513 )  ;// c_vector<long double, 3> Normal1 = Create_c_vector(-0.05607774749413225, 0.762765339938692, 0.6442393362906335);
        double Radius1 = 0.008988;// 1.7;

        c_vector<long double, 3> Boundary2 = Create_c_vector(0.05004181587627815, 0.21343413271111694,0.22811353187212016);// 12.597373380655702, 48.382440094438316, 42.984851419357064);
        c_vector<long double, 3> Normal2 = -Create_c_vector(-0.06301117186316568, -0.9774912340593714,-0.20134666512638721) ;      //-0.04847413662454751, -0.989768366942236, -0.13419701143842153);
        double Radius2 = 0.006160;// .2;

        c_vector<long double, 3> Boundary3 = Create_c_vector(0.05758695640342623, 0.2171917430630741,0.19389525030000318);       //-0.04710924565104161, -0.9898343148598956, -0.13419667693363677);
        c_vector<long double, 3> Normal3 = -Create_c_vector(0.009110299095533111, -0.9961614103846121,-0.08706001901521743);        //-0.04710924565104161, -0.9898343148598956, -0.13419667693363677);
        double Radius3 = 0.00614345;  //1.3;

        boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary1, Normal1, Radius1));
        simulator.AddCellPopulationBoundaryCondition(p_condition_1);

        boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary2, Normal2, Radius2));
        simulator.AddCellPopulationBoundaryCondition(p_condition_2);
        boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_3(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary3, Normal3, Radius3));
        simulator.AddCellPopulationBoundaryCondition(p_condition_3);


        simulator.Solve();
        TRACE("SolveFinished")
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);
        TRACE("SAVED")

        // To reset before looping: this is usually done by the SetUp and TearDown methods
        // SimulationTime::Instance()->Destroy();
        // SimulationTime::Instance()->SetStartTime(0.0);

        
   
    }
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/

