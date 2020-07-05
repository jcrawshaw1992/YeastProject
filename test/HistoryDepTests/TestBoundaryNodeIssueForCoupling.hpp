
#ifndef TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_CORRECTEDDRAG_HPP_
#define TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_CORRECTEDDRAG_HPP_

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

#include "projects/VascularRemodelling/src/MembraneForces/MembraneShearForce.hpp"

#include "projects/VascularRemodelling/src/MembraneForces/MembraneStiffnessForce.hpp"
 #include "MembraneSurfaceForce.hpp"
  #include "EdgeCorrectionForce.hpp"

#include "AppliedForce.hpp"
#include "AppliedForceModifier.hpp"
#include "ConstantPressure.hpp"



static const double M_TIME_FOR_SIMULATION = 0.01;//40; //50
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
        unsigned N_D = 100;
        unsigned N_Z = N_D*0.5;
        double scale = 1e3;
        double Length = 60e-6 * scale;//12e-3; //12e-3;
        double trans = -Length/2;
        double MaxZ = Length + trans;
        double MinZ = trans;

        double Radius =  5e-6 * scale;

        Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, 1e-3, Length);
        MutableMesh<2, 3>* p_mesh = generator.GetMesh();
        p_mesh->Translate(trans * unit_vector<double>(3, 2));


        

        std::stringstream out;

        std::string output_directory = "Mutations";//Shrinking"; // + Parameters + "/";

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

        c_vector<long double, 3> Node_location;
        for (AbstractCellPopulation<2, 3>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
            {
                unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
                Node_location = cell_population.GetNode(node_index)->rGetLocation();
                if (Node_location[2] == MinZ) // These are the nodes along the lower edge and need to be marked as mutated
                {
                    cell_iter->SetMutationState(p_state);
                }
                else if (Node_location[2] == MaxZ) // These are the nodes along the upper edge and need to be marked as mutated
                {
                    cell_iter->SetMutationState(p_state);
                }
            }



        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(10); //(M_TIME_FOR_SIMULATION);
        simulator.SetDt(0.0001);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetUpdateCellPopulationRule(false); // No remeshing.

        /*       
        -----------------------------
        Tractionforce 
        ----------------------------
        */

        std::string working_directory = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/DeformingCylinder/";

         // Create an Applied Force modifier to couple to Flow
        std::string traction_file = working_directory + "results/Extracted/surface-tractions.xtr";
        boost::shared_ptr<AppliedForceModifier<2, 3> > p_force_modifier(new AppliedForceModifier<2, 3>());

        p_force_modifier->SetResetTractionsOnCells(true, traction_file);
		p_force_modifier-s>SetupVessel(cell_population, output_directory);
        simulator.AddSimulationModifier(p_force_modifier);
        
        p_force_modifier->SetMembraneConstants(ElasticShearModulus , AreaDilationModulus, Area_constant, membrane_constant );
        
            double pressure = (133.322 * 16 + 144/2)*1e-12 * scale*scale;// 0.00002;//1.0666e4; // to match 80mmhg
            MAKE_PTR_ARGS(ConstantPressure, p_ConstantPressure, (pressure));
            simulator.AddForce(p_ConstantPressure);


//         /*
//         -----------------------------
//         Shearing Force 
//         ----------------------------
//         */

        double ElasticShearModulus = 4.4e-05;
        double  AreaDilationModulus = 0.9e-4; 
        double  membrane_constant =0.75e-13 ; 
        double Area_constant = 0.9e-4;

        double Scalling =1; 
        

        boost::shared_ptr<MembraneShearForce> p_shear_force(new MembraneShearForce());
        p_shear_force->SetScallingShear(Scalling);
        
        p_shear_force->SetAreaDilationModulus(AreaDilationModulus);
        
        p_shear_force->SetElasticShearModulus(ElasticShearModulus);
        p_shear_force->SetupMembraneConfiguration(cell_population);
        
        simulator.AddForce(p_shear_force);


//         /*
//         -----------------------------
//         Surface Area Force
//         ----------------------------
 //       */


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
        c_vector<long double, 3> Boundary1 = Create_c_vector(0, 0, -30e-6 * scale);
        c_vector<long double, 3> Normal1 = Create_c_vector(0, 0, 1);

        c_vector<long double, 3> Boundary2 = Create_c_vector(0, 0, 30e-6 * scale);
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
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/
