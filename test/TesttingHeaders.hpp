
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

// #include "BetaCateninOneHitCellMutationState.hpp"
// #include "WildTypeCellMutationState.hpp"

#include "MutableMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "MembranePropertiesModifier.hpp"
#include "MembraneForcesBasic.hpp"

// #include "ConstantPressure.hpp"
#include "EmptyBasementMatrix.hpp"
 #include "HasEndothelialCell.hpp"
#include "LostEndothelialCell.hpp"

  #include "CellMutationStatesWriter.hpp"
// //  #include "CellMutationStatesCountWriter.hpp"

//  // Cell population writers
// // #include "CellMutationStatesCountWriter.hpp"


// #include "CellIdWriter.hpp"

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
        unsigned N_D = 40;
        unsigned N_Z = N_D * 1.5;
        double CollapseFactor =10;

        double Length = 20e-6 * scale; //12e-3; //12e-3
        double Radius = 5e-6 * scale/CollapseFactor;
        double trans = -Length / 2;

            // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
            // in um will be too large and break chaste without carefull playing with or a tiny time step

            Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
            MutableMesh<2, 3>* p_mesh = generator.GetMesh();
            
            p_mesh->Translate(trans * unit_vector<double>(3, 2));


            std::string output_directory = "DevelopingCollapsingModifier/Third/";


            // Create cells
             MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

            // Create a cell population
            MeshBasedCellPopulation<2, 3> cell_population(*p_mesh, cells);
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.SetOutputMeshInVtk(true);

            // Here are the problem classes 
        //    cell_population.AddCellWriter<CellIdWriter>();
            // cell_population.AddCellWriter<CellMutationStatesWriter>();
            // AddPopulationWriter<CellMutationStatesWriter>();
            cell_population.AddCellWriter<CellMutationStatesWriter>();
// Cell based write

            // Set up cell-based simulation
            OffLatticeSimulation<2, 3> simulator(cell_population);
            simulator.SetOutputDirectory(output_directory);
            simulator.SetEndTime(100);
            simulator.SetDt(0.005); // 0.005
            simulator.SetSamplingTimestepMultiple(5);
            simulator.SetUpdateCellPopulationRule(false); // No remeshing.

// cellmutationstatewriter

// addtopopwriter
            // Set up the mutants 
            MAKE_PTR(EmptyBasementMatrix, p_Basement); //Mutation to mark nodes one basement 
            MAKE_PTR(HasEndothelialCell, p_EC); //Mutation to mark nodes with ECs
            MAKE_PTR(LostEndothelialCell, p_NextToBasement); //Mutation to mark nodes with ECs next to the basement 
            c_vector<long double, 3> Node_location;
            double MinZ = 2e-3;
         for (AbstractCellPopulation<2, 3>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
            {
                unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
                Node_location = cell_population.GetNode(node_index)->rGetLocation();
                if (abs(Node_location[2]) < MinZ) // These are the nodes along the lower edge and need to be marked as mutated
                {
                    cell_iter->SetMutationState(p_Basement);
                    TRACE("BING");
                }
                else 
                {
                    cell_iter->SetMutationState(p_EC);
                }
            }

        /* In order to visualize labelled cells we need to use the following command.*/
        // cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        /*
        -----------------------------
        Constant Pressure ballance 
        ----------------------------
        */

        // double P_blood = 0.0021; // Pa ==   1.6004e-05 mmHg
        // double P_tissue = 0.0015; // Pa == 1.1000e-05 mmHg

        // double TransmuralPressure = 6.6715e-04;//P_blood - P_tissue;
        // MAKE_PTR_ARGS(RadialForce, p_radial_force, (TransmuralPressure));
        // simulator.AddForce(p_radial_force);

        /*
        -----------------------------
        MembraneProperties Modifier
        ----------------------------
        */


        // std::map<double, c_vector<long double, 3> > GrowthMaps;  // From matlab sweep results 
        //                             //         KA,          Kalpha           Ks
        // GrowthMaps[10] = Create_c_vector(pow(10,-6.9), pow(10,-8.2459),pow(10, -9) );
        // GrowthMaps[8] = Create_c_vector(pow(10,-6.9), pow(10,-8.0160),pow(10, -9) );
        // GrowthMaps[6] = Create_c_vector(pow(10,-6.9), pow(10,-7.7300),pow(10, -9) );


        // GrowthMaps[5] = Create_c_vector(pow(10,-6.9341), pow(10,-7.7),pow(10, -8) );
        // GrowthMaps[4] = Create_c_vector(pow(10,-6.9), pow(10,-7.4224),pow(10, 8) );

        // GrowthMaps[2] = Create_c_vector(pow(10,-6.8), pow(10,-6.8124),pow(10, -7) );
        // GrowthMaps[1.5] = Create_c_vector(pow(10,-6.5), pow(10,-6.3491),pow(10, -7) );
        // GrowthMaps[1.2] =  Create_c_vector(pow(10,-6.2), pow(10,-5.8360),pow(10, -7) );

        // 1.2   -6.2000   -5.8360   -7
        // 1.5   -6.5000   -6.3491   -7
        // 2     -6.8000   -6.8124   -7
        
        // 4   -6.9000   -7.4224      -8
        // 5    -6.9341   -7.7000      -8

        // 6   -6.9000   -7.7300      -9
        // 8   -6.9000   -8.0160      -9
        // 10  -6.9000   -8.2459      -9


    
        // boost::shared_ptr<MembranePropertiesModifier<2,3> > p_Membrane_modifier(new MembranePropertiesModifier<2,3>());
        // p_Membrane_modifier->SetMembranePropeties(GrowthMaps);
        // simulator.AddSimulationModifier(p_Membrane_modifier);

        /*
        -----------------------------
        SMembrane forces
        ----------------------------
        */

        // double ElasticShearModulus = 1e-10;// 104.4e-05;
        // double AreaDilationModulus =  GrowthMaps[CollapseFactor](1);// //7.0795e-06  (80% defomation ); // 3.1623e-05 (20% defomation ) ;//0.9e-4;
        // double membrane_constant = 0 ;//0.75e-13;
        // double Area_constant = GrowthMaps[CollapseFactor](0); //6.3096e-06; (80% defomation );// 3.5481e-05;// 0.9e-4;

        // boost::shared_ptr<MembraneForcesBasic> p_shear_force(new MembraneForcesBasic());
        // p_shear_force->SetupMembraneConfiguration(cell_population);
        // simulator.AddForce(p_shear_force);
    
         /*
        -----------------------------
        Bending Force
        ----------------------------
        */
        // // double membrane_constant = 0 ;//1e-11;
        // boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
        //  p_membrane_force->SetScallingBending(Scalling);

        // p_membrane_force->SetMembraneStiffness(membrane_constant);
        // p_membrane_force->SetupInitialMembrane(*p_mesh, cell_population);

        // simulator.AddForce(p_membrane_force);

        // boost::shared_ptr<EdgeCorrectionForce> p_EdgeCorrectionForce(new EdgeCorrectionForce());
        // p_EdgeCorrectionForce->SetMeshType(1, N_D, N_Z );
        // simulator.AddForce(p_EdgeCorrectionForce);

        /*
        -----------------------------
        Boundaries
        ----------------------------
        */

        // //Create a plane boundary to represent the inlet and pass them to the simulation
        // c_vector<long double, 3> Boundary1 = Create_c_vector(0, 0, -10e-6 * scale);
        // c_vector<long double, 3> Normal1 = Create_c_vector(0, 0, 1);

        // c_vector<long double, 3> Boundary2 = Create_c_vector(0, 0, 10e-6 * scale);
        // c_vector<long double, 3> Normal2 = Create_c_vector(0, 0, -1);

        // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary1, Normal1, 10));
        // simulator.AddCellPopulationBoundaryCondition(p_condition_1);

        // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary2, Normal2, 10));
        // simulator.AddCellPopulationBoundaryCondition(p_condition_2);

        // simulator.Solve();

   
    }
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/

/*       
        -----------------------------
        Tractionforce 
        ----------------------------
        */

// std::string working_directory = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/DeformingCylinder/";

//  // Create an Applied Force modifier to couple to Flow
// std::string traction_file = working_directory + "results/Extracted/surface-tractions.xtr";
// boost::shared_ptr<AppliedForceModifier<2, 3> > p_force_modifier(new AppliedForceModifier<2, 3>());

// p_force_modifier->SetResetTractionsOnCells(true, traction_file);
// p_force_modifier->SetupVessel(cell_population, output_directory);
// simulator.AddSimulationModifier(p_force_modifier);

// p_force_modifier->SetMembraneConstants(ElasticShearModulus , AreaDilationModulus, Area_constant, membrane_constant );

//  c_vector<long double, 3> Node_location;
//         for (AbstractCellPopulation<2, 3>::Iterator cell_iter = cell_population.Begin();
//                  cell_iter != cell_population.End();
//                  ++cell_iter)
//             {
//                 unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
//                 Node_location = cell_population.GetNode(node_index)->rGetLocation();
//                 if (Node_location[2] == MinZ) // These are the nodes along the lower edge and need to be marked as mutated
//                 {
//                     cell_iter->SetMutationState(p_state);
//                 }
//                 else if (Node_location[2] == MaxZ) // These are the nodes along the upper edge and need to be marked as mutated
//                 {
//                     cell_iter->SetMutationState(p_state);
//                 }
//             }

// std::map<unsigned, c_vector<double, 2> > GrowthMaps;
//                                 //   KA,     Kalpha
// GrowthMaps[30] = Create_c_vector(pow(10,-9.7), pow(10,-9.7454));
// GrowthMaps[20] = Create_c_vector(pow(10,-9.3), pow(10,-9.3983));
// GrowthMaps[10] = Create_c_vector(pow(10,-8.8), pow(10,-8.8349));
// GrowthMaps[8] = Create_c_vector(pow(10,-8.5), pow(10,-8.5504));
// GrowthMaps[6] = Create_c_vector(pow(10,-8.1), pow(10,-8.1939));
// GrowthMaps[5] = Create_c_vector(pow(10,-7.9), pow(10,-7.9032));
// GrowthMaps[4] = Create_c_vector(pow(10,-7.6), pow(10,-7.6017));
// GrowthMaps[2] = Create_c_vector(pow(10,-6.6), pow(10,-6.603));
