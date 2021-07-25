#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "PetscSetupAndFinalize.hpp"

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"

#include "Debug.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "SmartPointers.hpp"
#include "VtkMeshReader.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"

#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "HistoryDepMutableMesh.hpp"

// #include "AppliedForceModifier.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "HemeLBForce.hpp"
#include "MembraneDeformationForce.hpp"
#include "MembraneDeformationForceOnCylinder.hpp"
#include "OutwardsPressureWithBreaks.hpp"
#include "OutwardsPressure.hpp"
#include "RemeshingTriggerOnStepHeteroModifier.hpp"
#include "StepHeteroModifier.hpp"
#include "MembraneBendingForce.hpp"

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:
    
    
    void TestFSICylinder_Hetero() throw(Exception)
    {
        std::string output_dir = "FSICylinder/Fine/Hetero4";
        std::string Archieve = "FSICylinder/Fine";
    
        double SamplingTimestepMultiple = 2000;
        double EndTime = 45;
        double scale = 1e3;
        double Length = 50e-6 * scale;
        double Radius = 0.5e-6 * scale;
        double dt = 0.001;
        double FSIIterations = 2000;


        // Load and fix any settings in the simulator
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieve, EndTime);

        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

        /* Remove the constant pressure force   */
        p_simulator->RemoveAllForces();
        p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);

        /*
        -----------------------------
        Add the HemeLB Force
        ----------------------------
        */


        c_vector<double, 3> PlaneNormal1 = Create_c_vector(0, 0, 1);
        c_vector<double, 3> Point1 = Create_c_vector(0, 0, 0.0002e-6 * scale);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(0, 0, -1);
        c_vector<double, 3> Point2 = Create_c_vector(0, 0, Length - 0.0002e-6 * scale);

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = (0.002133152 - 0.001466542) * 1.001; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later
        double OutletPressure = (0.002133152 - 0.001466542) * (0.999);

        boost::shared_ptr<HemeLBForce<2, 3> > p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(EndTime);
        p_ForceOut->SetFluidSolidIterations(FSIIterations);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir+"/HemeLBForce/", p_simulator->rGetCellPopulation());
        p_simulator->AddForce(p_ForceOut);

        // boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        // p_ForceOut->SetPressure(P_blood - P_tissue);
        // // p_ForceOut->SetRadiusThreshold(10 * Radius);
        // p_simulator->AddForce(p_ForceOut);


        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_simulator->AddForce(p_shear_force);


        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_simulator->AddForce(p_membrane_force);
 
        /* 
        -----------------------------
        Edit  RemeshingTriggerOnStepHeteroModifier
        ----------------------------
         */

        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2,3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);

        
        std::map<double, c_vector<long double, 4> >  GrowthMaps =  { {1, Create_c_vector(pow(10, -8), pow(10, -9), pow(10, -8.05), 1e-9) },
                                                                     {0,  Create_c_vector(pow(10, -5), pow(10, -4), pow(10, -5), 1e-9)}    };
 

        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
        // Upstream 
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0,0,20e-6* scale);
        c_vector<double, 3> UpperPlaneNormal = Create_c_vector(0,0,1);
        // Down stream
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0,0,30e-6 * scale);
        c_vector<double, 3> LowerPlaneNormal = Create_c_vector(0,0,-1);
        p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);
        p_Mesh_modifier->SetUpdateFrequency(10/dt);
        p_Mesh_modifier->SetmSetUpSolve(1);


        for (int j =1; j<=40; j++)
        {

            
            for (int i =1; i<=20; i++)
            { 
                // static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
                EndTime +=5;
                p_simulator->SetEndTime(EndTime);

                p_simulator->Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
            }

            dt/=5 ;  SamplingTimestepMultiple*= 5; FSIIterations*=5;
 
            p_ForceOut->SetFluidSolidIterations(FSIIterations);
            p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
            p_simulator->SetDt(dt);
            p_Mesh_modifier->SetUpdateFrequency(10/dt);


        }


       
    }

  
  
     void offTestFSICylinder_GrowToEqui() throw(Exception)
    {

        double EndTime = 30;
        double scale = 1e3;
        double Length = 50e-6 * scale;
        double Radius = 0.5e-6 * scale;

        unsigned N_D = 80;
        unsigned N_Z = 100;

        std::string output_dir = "FSICylinder/Fine/";

        // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
        // in um will be too large and break chaste without carefull playing with or a tiny time step

        Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
        MutableMesh<2, 3>* p_mesh = generator.GetMesh();
        HistoryDepMutableMesh<2, 3>* mesh = static_cast<HistoryDepMutableMesh<2, 3>*>(p_mesh);

        // Create the cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh->GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(*mesh, cells);
        cell_population.SetChasteOutputDirectory(output_dir, 0);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        // Set population to output all data to results files
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetDt(0.05);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(EndTime);

        /*
        -----------------------------
        RemeshingTriggerOnStepHeteroModifier
        ----------------------------
        */
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnStepHeteroModifier<2, 3>());
    
        std::map<double, c_vector<long double, 4> > GrowthMaps =  { {4, Create_c_vector(pow(10, -6), pow(10, -7), pow(10, -7), 1e-11) }, // Trying to Collapse
                         {3, Create_c_vector(pow(10, -8), pow(10, -9), pow(10, -8.05), 1e-11) },// Initial state 
                         {2, Create_c_vector(pow(10, -8), pow(10, -9), pow(10, -9), 1e-11) },// Allowed deformation
                         {1, Create_c_vector(pow(10, -6), pow(10, -7), pow(10, -7.5), 1e-11) },
                         {0, Create_c_vector(pow(10, -7), pow(10, -6), pow(10, -15), 1e-11)}// {0, Create_c_vector(pow(10, -7), pow(10, -6), pow(10, -5), 1e-10)}
                        };

        p_Mesh_modifier->SetMembranePropeties(GrowthMaps,3);
        // p_Mesh_modifier->SetMembraneStrength(3);
        simulator.AddSimulationModifier(p_Mesh_modifier);


        /*
        -----------------------------
        HemeLB Force
        ----------------------------
        */
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(0, 0, 1);
        c_vector<double, 3> Point1 = Create_c_vector(0, 0, 0.1e-6 * scale);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(0, 0, -1);
        c_vector<double, 3> Point2 = Create_c_vector(0, 0, Length - 0.1e-6 * scale);

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */

        boost::shared_ptr<OutwardsPressureWithBreaks> p_ForceOut(new OutwardsPressureWithBreaks());
        p_ForceOut->SetPressure(P_blood - P_tissue);
        p_ForceOut->SetRadiusThreshold(10 * Radius);
        simulator.AddForce(p_ForceOut);

        /*
        -----------------------------
        Boundary conditions
        ----------------------------
        */
        std::vector<c_vector<double, 3> > boundary_plane_points;
        std::vector<c_vector<double, 3> > boundary_plane_normals;

        boundary_plane_points.push_back(Point1);
        boundary_plane_normals.push_back(PlaneNormal1);

        boundary_plane_points.push_back(Point2);
        boundary_plane_normals.push_back(PlaneNormal2);

        for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
            boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(&cell_population, boundary_plane_points[boundary_id], boundary_plane_normals[boundary_id], 0.5));
            simulator.AddCellPopulationBoundaryCondition(p_condition);
        }

        simulator.Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);
    }

    
    
    void offTestFSICylinder_AssertEqui() throw(Exception)
    {
        std::string output_dir = "FSICylinder/Fine/";
    

        double EndTime = 30;
        double scale = 1e3;
        double Length = 50e-6 * scale;
        double Radius = 0.5e-6 * scale;


        // Load and fix any settings in the simulator
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(output_dir, EndTime);

        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

        /* Remove the constant pressure force   */
        // p_simulator->RemoveForce(0); // TRACE("RemoveForce will only work with the edit I made in OffLatticeSimulation.cpp line 69" )
        p_simulator->RemoveAllForces();
        p_simulator->SetEndTime(EndTime + 15);
        p_simulator->SetSamplingTimestepMultiple(500);
        p_simulator->SetDt(0.005);
        p_simulator->SetOutputDirectory(output_dir);

        /*
        -----------------------------
        Add the HemeLB Force
        ----------------------------
        */



        c_vector<double, 3> PlaneNormal1 = Create_c_vector(0, 0, 1);
        c_vector<double, 3> Point1 = Create_c_vector(0, 0, 0.0002e-6 * scale);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(0, 0, -1);
        c_vector<double, 3> Point2 = Create_c_vector(0, 0, Length - 0.0002e-6 * scale);

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = (0.002133152 - 0.001466542) * 1.001; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later
        double OutletPressure = (0.002133152 - 0.001466542) * (0.999);

        // boost::shared_ptr<HemeLBForce<2, 3> > p_ForceOut(new HemeLBForce<2, 3>());
        // p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        // p_ForceOut->Inlets(PlaneNormal2, Point2, OutletPressure, "Outlet");
        // p_ForceOut->SetStartTime(EndTime);
        // p_ForceOut->SetFluidSolidIterations(1000);
        // p_ForceOut->SetUpHemeLBConfiguration(output_dir+"Hetero/HemeLBForce/", p_simulator->rGetCellPopulation());
        // p_simulator->AddForce(p_ForceOut);



        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure(P_blood - P_tissue);
        // p_ForceOut->SetRadiusThreshold(10 * Radius);
        p_simulator->AddForce(p_ForceOut);



        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_simulator->AddForce(p_shear_force);


        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_simulator->AddForce(p_membrane_force);
 
        /* 
        -----------------------------
        Edit  RemeshingTriggerOnStepHeteroModifier
        ----------------------------
         */

        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2,3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);


       std::map<double, c_vector<long double, 4> >  GrowthMaps =  { {1, Create_c_vector(pow(10, -8), pow(10, -9), pow(10, -8.05), 1e-14) },
                                                                    {0,  Create_c_vector(pow(10, -6), pow(10, -7), pow(10, -7), 1e-9)}    };


        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);

        p_simulator->Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
    }

  
    
 


};

#endif /*TESTRELAXATION_HPP_*/
