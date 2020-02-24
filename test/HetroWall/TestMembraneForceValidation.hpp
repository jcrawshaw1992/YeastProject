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
#include "projects/VascularRemodelling/src/FixedRegionBoundaryCondition.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"

#include "MutableMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "MembraneForcesBasicCylinder.hpp"
#include "MembraneForcesBasic.hpp"
#include "MembranePropertiesSecModifier.hpp"

// #include "RadialForce.hpp"

#include "CellMutationStatesWriter.hpp"
#include "EmptyBasementMatrix.hpp"
#include "HasEndothelialCell.hpp"
#include "LostEndothelialCell.hpp"
#include "OutwardsPressure.hpp"
#include "VtkMeshReader.hpp"

#include "BoundariesModifier.hpp"

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



     void TestForcesOnCylinder_Stretched() throw(Exception)
    {
        unsigned Refinment[4] = {20,40,60,80};//,160};
    
        for (unsigned N_D_index = 0; N_D_index <4; N_D_index++)
        {
            unsigned N_D = Refinment[N_D_index];
            unsigned N_Z = Refinment[N_D_index]*0.75;//

            double Length = 20*1e-3;
            double Radius = 5*1e-3;

            Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
            MutableMesh<2, 3>* mesh = generator.GetMesh();
        
            std::stringstream out;
            out << N_D;
            std::string mesh_size = out.str();
            std::string output_directory = "TestForceValidation/Stretched/" + mesh_size;

            // Create cells
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells,mesh->GetNumNodes(), p_differentiated_type);

            // Create a cell population
            MeshBasedCellPopulation<2, 3> cell_population(*mesh, cells);
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.SetOutputMeshInVtk(true);
            // cell_population.AddCellWriter<CellMutationStatesWriter>();

            // Set up cell-based simulation
            OffLatticeSimulation<2, 3> simulator(cell_population);
            simulator.SetOutputDirectory(output_directory);
            simulator.SetEndTime(2000);
            simulator.SetDt(0.02); // 0.005
            simulator.SetSamplingTimestepMultiple(1000);
            simulator.SetUpdateCellPopulationRule(false); // No remeshing.

            /*
            -----------------------------
            Constant Pressure ballance 
            ----------------------------
            */
            double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
            double P_tissue = 0.001466542; // Pa == 1.1000e-05 mmHg

            double TransmuralPressure =  P_blood - P_tissue;

            boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
            p_ForceOut->SetPressure(TransmuralPressure*1e3);
            simulator.AddForce(p_ForceOut);

            /*
            -----------------------------
            MembraneProperties Modifier
            ----------------------------
            */

            double BendingConst = 1e-11;
            std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
                    //         KA,          Kalpha           Ks                                 Kb
            GrowthMaps[10] = Create_c_vector(pow(10, -6.9), pow(10, -8.2459), pow(10, -9), 0);
            GrowthMaps[5] = Create_c_vector(pow(10, -6.9341), pow(10, -7.7), pow(10, -8), BendingConst);
            GrowthMaps[1.5] = Create_c_vector(pow(10, -6.5), pow(10, -6.3491), pow(10, -7), BendingConst);
            GrowthMaps[1.2] = Create_c_vector(pow(10, -6.2)*1e3, pow(10, -5.8360)*1e2, pow(10, -7)*1e2, BendingConst*1e1);
            GrowthMaps[2] = Create_c_vector(pow(10, -6.8), pow(10, -6.8124), pow(10, -7), BendingConst);


            boost::shared_ptr<MembranePropertiesSecModifier<2, 3> > p_Membrane_modifier(new MembranePropertiesSecModifier<2, 3>());
            p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 2, 0,10, 1); 
            p_Membrane_modifier->SetupSolve(cell_population,output_directory);

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
            Boundaries
            ----------------------------
            */

            //Create a plane boundary to represent the inlet and pass them to the simulation
                c_vector<double, 3> Boundary1 = Create_c_vector(0, 0, 0);
                c_vector<double, 3> Normal1 = Create_c_vector(0, 0, 1);

                c_vector<double, 3> Boundary2 = Create_c_vector(0, 0,12*1e-3);
                c_vector<double, 3> Normal2 = Create_c_vector(0, 0, -1);

                boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary1, Normal1, 2));
                simulator.AddCellPopulationBoundaryCondition(p_condition_1);

                boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary2, Normal2, 2));
                simulator.AddCellPopulationBoundaryCondition(p_condition_2);


                simulator.Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);

                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
            }
    }
        


     void TestForcesOnCylinder_Equi() throw(Exception)
    {
        unsigned Refinment[4] = {20,40,60,80};//,160};
    
        for (unsigned N_D_index = 0; N_D_index <4; N_D_index++)
        {
            unsigned N_D = Refinment[N_D_index];
            unsigned N_Z = Refinment[N_D_index]*1.5;//

            double Length = 20*1e-3;
            double Radius = 5*1e-3;

            Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
            MutableMesh<2, 3>* mesh = generator.GetMesh();
        
            std::stringstream out;
            out << N_D;
            std::string mesh_size = out.str();
            std::string output_directory = "TestForceValidation/Equi/" + mesh_size;

            // Create cells
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells,mesh->GetNumNodes(), p_differentiated_type);

            // Create a cell population
            MeshBasedCellPopulation<2, 3> cell_population(*mesh, cells);
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.SetOutputMeshInVtk(true);
            // cell_population.AddCellWriter<CellMutationStatesWriter>();

            // Set up cell-based simulation
            OffLatticeSimulation<2, 3> simulator(cell_population);
            simulator.SetOutputDirectory(output_directory);
            simulator.SetEndTime(2000);
            simulator.SetDt(0.02); // 0.005
            simulator.SetSamplingTimestepMultiple(1000);
            simulator.SetUpdateCellPopulationRule(false); // No remeshing.

            /*
            -----------------------------
            Constant Pressure ballance 
            ----------------------------
            */
            double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
            double P_tissue = 0.001466542; // Pa == 1.1000e-05 mmHg

            double TransmuralPressure =  P_blood - P_tissue;

            boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
            p_ForceOut->SetPressure(TransmuralPressure*1e3);
            simulator.AddForce(p_ForceOut);

            /*
            -----------------------------
            MembraneProperties Modifier
            ----------------------------
            */

            double BendingConst = 1e-11;
            std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
                    //         KA,          Kalpha           Ks                                 Kb
            GrowthMaps[10] = Create_c_vector(pow(10, -6.9), pow(10, -8.2459), pow(10, -9), 0);
            GrowthMaps[5] = Create_c_vector(pow(10, -6.9341), pow(10, -7.7), pow(10, -8), BendingConst);
            GrowthMaps[1.5] = Create_c_vector(pow(10, -6.5), pow(10, -6.3491), pow(10, -7), BendingConst);
            GrowthMaps[1.2] = Create_c_vector(pow(10, -6.2)*1e3, pow(10, -5.8360)*1e2, pow(10, -7)*1e2, BendingConst*1e1);
            GrowthMaps[2] = Create_c_vector(pow(10, -6.8), pow(10, -6.8124), pow(10, -7), BendingConst);


            boost::shared_ptr<MembranePropertiesSecModifier<2, 3> > p_Membrane_modifier(new MembranePropertiesSecModifier<2, 3>());
            p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 2, 0,10, 1); 
            p_Membrane_modifier->SetupSolve(cell_population,output_directory);

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
            Boundaries
            ----------------------------
            */

            //Create a plane boundary to represent the inlet and pass them to the simulation
                c_vector<double, 3> Boundary1 = Create_c_vector(0, 0, 0);
                c_vector<double, 3> Normal1 = Create_c_vector(0, 0, 1);

                c_vector<double, 3> Boundary2 = Create_c_vector(0, 0,12*1e-3);
                c_vector<double, 3> Normal2 = Create_c_vector(0, 0, -1);

                boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary1, Normal1, 2));
                simulator.AddCellPopulationBoundaryCondition(p_condition_1);

                boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary2, Normal2, 2));
                simulator.AddCellPopulationBoundaryCondition(p_condition_2);


                simulator.Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);

                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
            }
    }
     
 
     void TestForcesOnCylinder_Squashed() throw(Exception)
    {
        unsigned Refinment[4] = {20,40,60,80};//,160};
    
        for (unsigned N_D_index = 0; N_D_index <4; N_D_index++)
        {
            unsigned N_D = Refinment[N_D_index];
            unsigned N_Z = Refinment[N_D_index]*3;//

            double Length = 20*1e-3;
            double Radius = 5*1e-3;

            Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
            MutableMesh<2, 3>* mesh = generator.GetMesh();
        
            std::stringstream out;
            out << N_D;
            std::string mesh_size = out.str();
            std::string output_directory = "TestForceValidation/Squashed/" + mesh_size;

            // Create cells
            MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells,mesh->GetNumNodes(), p_differentiated_type);

            // Create a cell population
            MeshBasedCellPopulation<2, 3> cell_population(*mesh, cells);
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.SetOutputMeshInVtk(true);
            // cell_population.AddCellWriter<CellMutationStatesWriter>();

            // Set up cell-based simulation
            OffLatticeSimulation<2, 3> simulator(cell_population);
            simulator.SetOutputDirectory(output_directory);
            simulator.SetEndTime(2000);
            simulator.SetDt(0.02); // 0.005
            simulator.SetSamplingTimestepMultiple(1000);
            simulator.SetUpdateCellPopulationRule(false); // No remeshing.

            /*
            -----------------------------
            Constant Pressure ballance 
            ----------------------------
            */
            double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
            double P_tissue = 0.001466542; // Pa == 1.1000e-05 mmHg

            double TransmuralPressure =  P_blood - P_tissue;

            boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
            p_ForceOut->SetPressure(TransmuralPressure*1e3);
            simulator.AddForce(p_ForceOut);

            /*
            -----------------------------
            MembraneProperties Modifier
            ----------------------------
            */

            double BendingConst = 1e-11;
            std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
                    //         KA,          Kalpha           Ks                                 Kb
            GrowthMaps[10] = Create_c_vector(pow(10, -6.9), pow(10, -8.2459), pow(10, -9), 0);
            GrowthMaps[5] = Create_c_vector(pow(10, -6.9341), pow(10, -7.7), pow(10, -8), BendingConst);
            GrowthMaps[1.5] = Create_c_vector(pow(10, -6.5), pow(10, -6.3491), pow(10, -7), BendingConst);
            GrowthMaps[1.2] = Create_c_vector(pow(10, -6.2)*1e3, pow(10, -5.8360)*1e2, pow(10, -7)*1e2, BendingConst*1e1);
            GrowthMaps[2] = Create_c_vector(pow(10, -6.8), pow(10, -6.8124), pow(10, -7), BendingConst);


            boost::shared_ptr<MembranePropertiesSecModifier<2, 3> > p_Membrane_modifier(new MembranePropertiesSecModifier<2, 3>());
            p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 2, 0,10, 1); 
            p_Membrane_modifier->SetupSolve(cell_population,output_directory);

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
            Boundaries
            ----------------------------
            */

            //Create a plane boundary to represent the inlet and pass them to the simulation
                c_vector<double, 3> Boundary1 = Create_c_vector(0, 0, 0);
                c_vector<double, 3> Normal1 = Create_c_vector(0, 0, 1);

                c_vector<double, 3> Boundary2 = Create_c_vector(0, 0,12*1e-3);
                c_vector<double, 3> Normal2 = Create_c_vector(0, 0, -1);

                boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary1, Normal1, 2));
                simulator.AddCellPopulationBoundaryCondition(p_condition_1);

                boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary2, Normal2, 2));
                simulator.AddCellPopulationBoundaryCondition(p_condition_2);


                simulator.Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);

                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
            }
    }
     
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/