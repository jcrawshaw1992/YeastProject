
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

#include "MembranePropertiesModifier.hpp"
#include "MembraneForcesBasic.hpp"

// #include "ConstantPressure.hpp"
#include "EmptyBasementMatrix.hpp"
#include "HasEndothelialCell.hpp"
#include "LostEndothelialCell.hpp"
#include "CellMutationStatesWriter.hpp"

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

            c_vector<double, 3> Normal = zero_vector<double>(3);
            double Area = 0;
            std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();
            assert(containing_elements.size() > 0);
            for (std::set<unsigned>::iterator iter = containing_elements.begin();
                 iter != containing_elements.end();
                 ++iter)
            {
                Node<3>* pNode0 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(0));
                Node<3>* pNode1= p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(1));
                Node<3>* pNode2 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(2));

                c_vector<double, 3> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
                c_vector<double, 3> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3

                c_vector<double, 3> normalVector = VectorProduct(vector_12, vector_13);
                Area+= norm_2(normalVector)/6;

            }
            
             c_vector<double, 3> force = mStrength * cell_location; // / norm_2(cell_location);
         //   PRINT_VECTOR(force);

            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(force);

            cell_iter->GetCellData()->SetItem("applied_force_x", force[0]);
            cell_iter->GetCellData()->SetItem("applied_force_y", force[1]);
            cell_iter->GetCellData()->SetItem("applied_force_z", force[2]);
            cell_iter->GetCellData()->SetItem("applied_force_mag", norm_2(force));
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
    void TestCollapseCode() throw(Exception)
    {
        double scale = 1e3;
        unsigned N_D = 60;
        unsigned N_Z = N_D * 2;
        double CollapseFactor =10;

        double Length = 120e-6 * scale; //12e-3; //12e-3
        double Radius = 5e-6 * scale/CollapseFactor;
        double trans = -Length / 2;

            // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
            // in um will be too large and break chaste without carefull playing with or a tiny time step

            Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
            MutableMesh<2, 3>* p_mesh = generator.GetMesh();
            
            p_mesh->Translate(trans * unit_vector<double>(3, 2));


            std::string output_directory = "DevelopingCollapsingCylinderLong/SetUpArchiving/";


            // Create cells
             MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

            // Create a cell population
            MeshBasedCellPopulation<2, 3> cell_population(*p_mesh, cells);
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.SetOutputMeshInVtk(true);
            cell_population.AddCellWriter<CellMutationStatesWriter>();

            // Set up cell-based simulation
            OffLatticeSimulation<2, 3> simulator(cell_population);
            simulator.SetOutputDirectory(output_directory);
            simulator.SetEndTime(1);
            simulator.SetDt(0.01); // 0.005
            simulator.SetSamplingTimestepMultiple(100);
            simulator.SetUpdateCellPopulationRule(false); // No remeshing.

        /*
        -----------------------------
        Constant Pressure ballance 
        ----------------------------
        */

        double P_blood = 0.0021; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.0015; // Pa == 1.1000e-05 mmHg

        double TransmuralPressure = 6.6715e-04;//P_blood - P_tissue;
        MAKE_PTR_ARGS(RadialForce, p_radial_force, (TransmuralPressure));
        simulator.AddForce(p_radial_force);

        /*
        -----------------------------
        MembraneProperties Modifier
        ----------------------------
        */


        std::map<double, c_vector<long double, 4> > GrowthMaps;  // From matlab sweep results 
                                    //         KA,          Kalpha           Ks
        GrowthMaps[10] = Create_c_vector(pow(10,-6.9), pow(10,-8.2459),pow(10, -9) , 0);
        GrowthMaps[8] = Create_c_vector(pow(10,-6.9), pow(10,-8.0160),pow(10, -9) , 0);
        GrowthMaps[6] = Create_c_vector(pow(10,-6.9), pow(10,-7.7300),pow(10, -9) , 0);


        GrowthMaps[5] = Create_c_vector(pow(10,-6.9341), pow(10,-7.7),pow(10, -8), 0 );
        GrowthMaps[4] = Create_c_vector(pow(10,-6.9), pow(10,-7.4224),pow(10, 8), 0 );

        GrowthMaps[2] = Create_c_vector(pow(10,-6.8), pow(10,-6.8124),pow(10, -7) , 0);
        GrowthMaps[1.5] = Create_c_vector(pow(10,-6.5), pow(10,-6.3491),pow(10, -7) , 0);
        GrowthMaps[1.2] =  Create_c_vector(pow(10,-6.2), pow(10,-5.8360),pow(10, -7) , 0);
    
        boost::shared_ptr<MembranePropertiesModifier<2,3> > p_Membrane_modifier(new MembranePropertiesModifier<2,3>());
        p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 10, 1, 1);
        simulator.AddSimulationModifier(p_Membrane_modifier);

        /*
        -----------------------------
        SMembrane forces
        ----------------------------
        */

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

        //Create a plane boundary to represent the inlet and pass them to the simulation
        c_vector<double, 3> Boundary1 = Create_c_vector(0, 0, -60e-6 * scale);
        c_vector<double, 3> Normal1 = Create_c_vector(0, 0, 1);

        c_vector<double, 3> Boundary2 = Create_c_vector(0, 0, 60e-6 * scale);
        c_vector<double, 3> Normal2 = Create_c_vector(0, 0, -1);

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

