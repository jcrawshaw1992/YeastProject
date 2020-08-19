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
#include "MembranePropertiesSecModifier.hpp"

// #include "RadialForce.hpp"

#include "CellMutationStatesWriter.hpp"
#include "EmptyBasementMatrix.hpp"
#include "HasEndothelialCell.hpp"
#include "LostEndothelialCell.hpp"
#include "OutwardsPressure.hpp"

#include "projects/VascularRemodelling/src/MembraneForces/MembraneStiffnessForce.hpp"


static const double M_TIME_FOR_SIMULATION = 100; //40; //50
static const double M_SAMPLING_TIME_STEP = 1000; //50
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
    void TestSetUpArchiving() throw(Exception)
    {
        double scale = 1e3;
        unsigned N_D = 30;
        unsigned N_Z = N_D * 1.5;
        double CollapseFactor = 5;

        double Length = 20e-6 * scale; //12e-3; //12e-3
        double Radius = 4e-6 * scale / CollapseFactor;
        double trans = -Length / 2;

         for (int i = -16; i < -15; i++)
        {


            // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
            // in um will be too large and break chaste without carefull playing with or a tiny time step

            Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
            MutableMesh<2, 3>* p_mesh = generator.GetMesh();

            p_mesh->Translate(trans * unit_vector<double>(3, 2));


            std::stringstream out;
            out << "Iter_" << i;
            std::string Iteration = out.str();
            std::string output_directory = "BendingForceOnCollapsingCylinder/"+Iteration;

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
            simulator.SetEndTime(10);
            simulator.SetDt(0.01); // 0.005
            simulator.SetSamplingTimestepMultiple(100);
            simulator.SetUpdateCellPopulationRule(false); // No remeshing.

              for (AbstractCellPopulation<2, 3>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
                cell_iter != simulator.rGetCellPopulation().End();
                ++cell_iter)
            {

                cell_iter->GetCellData()->SetItem("Boundary", 0);
            }

            /*
            -----------------------------
            Constant Pressure ballance 
            ----------------------------
            */

            double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
            double P_tissue = 0.001466542; // Pa == 1.1000e-05 mmHg

            double TransmuralPressure =  P_blood - P_tissue;
            boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
            p_ForceOut->SetPressure(TransmuralPressure);
            simulator.AddForce(p_ForceOut);

            /*
            -----------------------------
            MembraneProperties Modifier
            ----------------------------
            */

            double BendingConst =pow(10, i);
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

            boost::shared_ptr<MembranePropertiesSecModifier<2, 3> > p_Membrane_modifier(new MembranePropertiesSecModifier<2, 3>());
            p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 5, 0,1e-10, 1); 
            p_Membrane_modifier->SetupSolve(cell_population, output_directory);
            simulator.AddSimulationModifier(p_Membrane_modifier);

            /*
            -----------------------------
            SMembrane forces
            ----------------------------
            */
            boost::shared_ptr<MembraneForcesBasicCylinder> p_shear_force(new MembraneForcesBasicCylinder());
            p_shear_force->SetupMembraneConfiguration(cell_population);
            simulator.AddForce(p_shear_force);

            /*
            -----------------------------
            Bending forces
            ----------------------------
            */

            boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
            p_membrane_force->SetupInitialMembrane(*p_mesh, simulator.rGetCellPopulation());
            p_membrane_force->SetMembraneStiffness(BendingConst,N_D, N_Z );
            simulator.AddForce(p_membrane_force);

            /*
            -----------------------------
            Boundaries
            ----------------------------
            */

            //Create a plane boundary to represent the inlet and pass them to the simulation
            c_vector<double, 3> Boundary1 = Create_c_vector(0, 0, -10e-6 * scale);
            c_vector<double, 3> Normal1 = Create_c_vector(0, 0, 1);

            c_vector<double, 3> Boundary2 = Create_c_vector(0, 0, 10e-6 * scale);
            c_vector<double, 3> Normal2 = Create_c_vector(0, 0, -1);

            boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary1, Normal1, 10));
            simulator.AddCellPopulationBoundaryCondition(p_condition_1);
            boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary2, Normal2, 10));
            simulator.AddCellPopulationBoundaryCondition(p_condition_2);

            simulator.Solve();

            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);

            TRACE("After Growth To Equi");
            
            simulator.SetEndTime(10+ 100);
            // simulator.SetOutputDirectory(output_directory);
            simulator.SetDt(0.001);
            simulator.SetSamplingTimestepMultiple(200);

            std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = simulator.GetSimulationModifiers()->begin();
            boost::shared_ptr<MembranePropertiesSecModifier<2,3> > force_modifier = boost::static_pointer_cast<MembranePropertiesSecModifier<2, 3> >(*iter);

            
            force_modifier->SetBendingForce(simulator.rGetCellPopulation() , BendingConst);

            // c_vector<double, 3> UpperPlaneNormal = Create_c_vector(0,0,1);
            // c_vector<double, 3> UpperPlanePoint = Create_c_vector(0,0,2);

            // c_vector<double, 3> LowerPlanePoint = Create_c_vector(0,0,0);
            // c_vector<double, 3> LowerPlaneNormal = Create_c_vector(0,0,1);
            
            // force_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);
            // force_modifier->SetBendingForce(cell_population, BendingConst);
            // force_modifier->SetMembranePropeties(GrowthMaps, 5, 1,1e-9, 1);
            // force_modifier->SetupSolve(cell_population, output_directory); 




            simulator.Solve();

            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);
        }

    }
    //  void offTestSweepBendingParameter() throw(Exception)
    // {
    //     std::string load_dir = "DevelopingHetroCylinder/SetChange/SetUpArchiving"; 
    //     double start_time = 40;
    //     double simulation_time = 10;

    //      std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
    //             //         KA,          Kalpha           Ks                                 Kb

    //     for (int i = -16; i < -11; i++)
    //     {
    //         PRINT_VARIABLE(i)
    //         std::stringstream out;
    //         out << "Iter_" << i;
    //         std::string Iteration = out.str();

    //         std::string output_directory = "DevelopingHetroCylinder/SweepingBendingForces/"+Iteration;
    //         std::cout<< load_dir <<endl;
            
    //         OffLatticeSimulation<2,3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2,3>, 3 >::Load(load_dir , start_time);
    //         TRACE("After Load");
            
    //         p_simulator->SetEndTime(start_time + simulation_time);
    //         p_simulator->SetOutputDirectory(output_directory);
    //         p_simulator->SetDt(0.0005);
    //         p_simulator->SetSamplingTimestepMultiple(1);



    //         for (AbstractCellPopulation<2, 3>::Iterator cell_iter = p_simulator->rGetCellPopulation().Begin();
    //             cell_iter != p_simulator->rGetCellPopulation().End();
    //             ++cell_iter)
    //         {

    //             cell_iter->GetCellData()->SetItem("Boundary", 0);
    //         }

    //         double BendingConst =pow(10, i);

    //         GrowthMaps[5] = Create_c_vector(pow(10, -6.9341), pow(10, -7.7), pow(10, -8), BendingConst);
    //         GrowthMaps[1.2] = Create_c_vector(pow(10, -6.2), pow(10, -5.8360), pow(10, -7), BendingConst);

    //         std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
    //         boost::shared_ptr<MembranePropertiesSecModifier<2,3> > p_force_modifier = boost::static_pointer_cast<MembranePropertiesSecModifier<2, 3> >(*iter);
    //         p_force_modifier->SetMembranePropeties(GrowthMaps, 5, 1,1e-9, 1);  

    //         c_vector<double, 3> UpperPlaneNormal = Create_c_vector(0,0,1);
    //         c_vector<double, 3> UpperPlanePoint = Create_c_vector(0,0,0);

    //         c_vector<double, 3> LowerPlanePoint = Create_c_vector(0,0,0);
    //         c_vector<double, 3> LowerPlaneNormal = Create_c_vector(0,0,1);

        
    //     p_force_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);
    //     p_force_modifier->SetBendingForce(cell_population, BendingConst);

    //     simulator.AddSimulationModifier(p_force_modifier);





            
    //         p_force_modifier->SetBendingForce(p_simulator->rGetCellPopulation() , BendingConst);

    //         p_simulator->Solve();

    //         TRACE("Before Save");
    //         CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
    //         TRACE("After Save");

    //         SimulationTime::Instance()->Destroy();
    //         SimulationTime::Instance()->SetStartTime(start_time);
    //     }
    // }

};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/