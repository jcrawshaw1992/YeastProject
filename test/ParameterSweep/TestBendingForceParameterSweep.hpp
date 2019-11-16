
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

// #include "MembraneShearForce.hpp"
#include "MembraneForcesBasic.hpp"

// #include "EdgeCorrectionForce.hpp"
// #include "MembraneSurfaceForce.hpp"

// #include "AppliedForce.hpp"
// #include "AppliedForceModifier.hpp"
#include "ConstantPressure.hpp"

static const double M_TIME_FOR_SIMULATION = 100; //40; //50
static const double M_SAMPLING_TIME_STEP = 100; //50
static const double M_TIME_STEP = 0.002;


class RadialForce : public AbstractForce<2, 3>
{
private:
    double mStrength;

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractForce<2, 3> >(*this);
        archive& mStrength;
    }

public:
    RadialForce(double strength = 1.0)
            : AbstractForce<2, 3>(),
              mStrength(strength)
    {
        assert(mStrength > 0.0);
    }

    void AddForceContribution(AbstractCellPopulation<2, 3>& rCellPopulation)
    {
        MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);
        

       // Calculate midpoint
        c_vector<double, 3> centroid = zero_vector<double>(3);
        for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            centroid += rCellPopulation.GetNode(node_index)->rGetLocation();
           // Area += rCellPopulation.GetVolumeOfCell(*cell_iter);
        }
        centroid /= rCellPopulation.GetNumRealCells();
        
        
        for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {

            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<3>* p_node = rCellPopulation.GetNode(node_index);
            c_vector<double, 3> cell_location = p_node->rGetLocation() - centroid;
            cell_location(2) = 0.0;
            cell_location /= norm_2(cell_location);

            c_vector<long double, 3> Normal = zero_vector<long double>(3);
            long double Area = 0;
            std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();
            assert(containing_elements.size() > 0);
            for (std::set<unsigned>::iterator iter = containing_elements.begin();
                 iter != containing_elements.end();
                 ++iter)
            {
                Node<3>* pNode0 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(0));
                Node<3>* pNode1= p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(1));
                Node<3>* pNode2 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(2));

                c_vector<long double, 3> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
                c_vector<long double, 3> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3

                c_vector<long double, 3> normalVector = VectorProduct(vector_12, vector_13);
                Area+= norm_2(normalVector)/6;

            }
            
             c_vector<long double, 3> force = mStrength * cell_location; // / norm_2(cell_location);
         //   PRINT_VECTOR(force);

            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(force);

            cell_iter->GetCellData()->SetItem("applied_force_x", force[0]);
            cell_iter->GetCellData()->SetItem("applied_force_y", force[1]);
            cell_iter->GetCellData()->SetItem("applied_force_z", force[2]);
            cell_iter->GetCellData()->SetItem("applied_force_mag", norm_2(force));

            // cell_iter->GetCellData()->SetItem("norm_x", normal[0]);
            // cell_iter->GetCellData()->SetItem("norm_y", normal[1]);
            // cell_iter->GetCellData()->SetItem("norm_z", normal[2]);
        }
    }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<Strength>" << mStrength << "</Strength>\n";
        AbstractForce<2, 3>::OutputForceParameters(rParamsFile);
    }
};





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
    void TestBendingForceParameterSweep() throw(Exception)
    {
        double scale = 1e3;
        unsigned N_D = 70;
        unsigned N_Z = N_D * 1;
        

        double Length = 5e-6 * scale; //12e-3; //12e-3
        double Radius = 5e-6 * scale /10;
        double trans = -Length / 2;

        std::map<unsigned, c_vector<long double, 2> > GrowthMaps;
                                        //   KA,     Kalpha
        GrowthMaps[30] = Create_c_vector(pow(10,-9), pow(10, -8.993866625781013 ));
        GrowthMaps[20] = Create_c_vector(pow(10,-9.3), pow(10,-9.3983));
        GrowthMaps[10] = Create_c_vector(pow(10,-8.8), pow(10,-8.8349));
        GrowthMaps[8] = Create_c_vector(pow(10,-8.5), pow(10,-8.5504));
        GrowthMaps[6] = Create_c_vector(pow(10,-8.1), pow(10,-8.1939));
        GrowthMaps[5] = Create_c_vector(pow(10,-7.9), pow(10,-7.9032));
        GrowthMaps[4] = Create_c_vector(pow(10,-7.6), pow(10,-7.6017));
        GrowthMaps[2] = Create_c_vector(pow(10,-6.6), pow(10,-6.603));

        // -9.5000   -9.5572
        // -8.800000000000002  -8.873376863585321

                    // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
                    // in um will be too large and break chaste without carefull playing with or a tiny time step

                    Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
                    MutableMesh<2, 3>* p_mesh = generator.GetMesh();
                    
                    p_mesh->Translate(trans * unit_vector<double>(3, 2));


                    std::string output_directory = "TestingParameterSweep/";


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
                    simulator.SetEndTime(50);
                    simulator.SetDt(0.0005); // 0.005
                    simulator.SetSamplingTimestepMultiple(500);
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
                    // MAKE_PTR_ARGS(ConstantPressure, p_ConstantPressure, (TransmuralPressure));
                    // simulator.AddForce(p_ConstantPressure);

                     MAKE_PTR_ARGS(RadialForce, p_radial_force, (TransmuralPressure));
            simulator.AddForce(p_radial_force);


                    /*
        -----------------------------
        SMembrane forces
        ----------------------------
        */

                    //    % 80% DEFORMATION

                    double ElasticShearModulus = 1e-10;// 104.4e-05;
                    double AreaDilationModulus =  GrowthMaps[30](1);// //7.0795e-06  (80% defomation ); // 3.1623e-05 (20% defomation ) ;//0.9e-4;
                    double membrane_constant = 0 ;//0.75e-13;
                    double Area_constant = GrowthMaps[30](0); //6.3096e-06; (80% defomation );// 3.5481e-05;// 0.9e-4;

                PRINT_2_VARIABLES(Area_constant ,AreaDilationModulus);


              


                    boost::shared_ptr<MembraneForcesBasic> p_shear_force(new MembraneForcesBasic());
                    p_shear_force->SetAreaDilationModulus(AreaDilationModulus);
                    p_shear_force->SetElasticShearModulus(ElasticShearModulus);
                    p_shear_force->SetMembraneStiffness(Area_constant);
                    p_shear_force->SetupMembraneConfiguration(cell_population);
                    simulator.AddForce(p_shear_force);

    
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
                    c_vector<long double, 3> Boundary1 = Create_c_vector(0, 0, -2.5e-6 * scale);
                    c_vector<long double, 3> Normal1 = Create_c_vector(0, 0, 1);

                    c_vector<long double, 3> Boundary2 = Create_c_vector(0, 0, 2.5e-6 * scale);
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
