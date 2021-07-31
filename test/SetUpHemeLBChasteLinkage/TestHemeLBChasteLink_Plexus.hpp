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
    
    
    void offTestFSICylinder_HeteroContinued() throw(Exception)
    {
        std::string output_dir = "FSICylinder/Medium/Hetro5";
        std::string Archieve = "FSICylinder/Medium/Hetro4";
    
        double SamplingTimestepMultiple = 100;
        double EndTime = 30;
        double scale = 1e3;
        double Length = 50e-6 * scale;
        double Radius = 0.5e-6 * scale;
        double dt = 0.0005;
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

        
        std::map<double, c_vector<long double, 4> >  GrowthMaps =  { {1, Create_c_vector(pow(10, -8), pow(10, -9), pow(10, -8.05), 1e-10 ) },
                                                                     {0,  Create_c_vector(pow(10, -5), pow(10, -4), pow(10, -5), 1e-10 )}    };
 

        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);

        p_Mesh_modifier->SetStepSize(pow(10, -8));

        // Upstream
        // p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);
        p_Mesh_modifier->SetUpdateFrequency(2/dt);
        p_Mesh_modifier->SetmSetUpSolve(0);

       
        for (int j =1; j<=40; j++)
        {

            
            for (int i =1; i<=10; i++)
            { 
                // static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
                EndTime +=5;
                p_simulator->SetEndTime(EndTime);

                p_simulator->Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
            }

            dt/=2 ;  SamplingTimestepMultiple*= 2; 
            FSIIterations*=2;
 
            p_ForceOut->SetFluidSolidIterations(FSIIterations);
            p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
            p_simulator->SetDt(dt);
            p_Mesh_modifier->SetUpdateFrequency(2/dt);


        }


       
    }

  

    void offTestFSICylinder_Hetero() throw(Exception)
    {
        std::string output_dir = "FSICylinder/Medium/Hetro4";
        std::string Archieve = "FSICylinder/Medium";
    
        double SamplingTimestepMultiple = 100;
        double EndTime = 20;
        double scale = 1e3;
        double Length = 50e-6 * scale;
        double Radius = 0.5e-6 * scale;
        double dt = 0.005;
        double FSIIterations = 200;


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

        
        std::map<double, c_vector<long double, 4> >  GrowthMaps =  { {1, Create_c_vector(pow(10, -8), pow(10, -9), pow(10, -8.05), 1e-10 ) },
                                                                     {0,  Create_c_vector(pow(10, -5), pow(10, -4), pow(10, -5), 1e-10 )}    };
 

        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);

        p_Mesh_modifier->SetStepSize(pow(10, -8));

        // Upstream 
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0,0,20e-6* scale);
        c_vector<double, 3> UpperPlaneNormal = Create_c_vector(0,0,1);
        // Down stream
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0,0,30e-6 * scale);
        c_vector<double, 3> LowerPlaneNormal = Create_c_vector(0,0,-1);
        p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);
        p_Mesh_modifier->SetUpdateFrequency(2/dt);
        p_Mesh_modifier->SetmSetUpSolve(1);

       
        for (int j =1; j<=40; j++)
        {

            
            for (int i =1; i<=10; i++)
            { 
                // static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
                EndTime +=5;
                p_simulator->SetEndTime(EndTime);

                p_simulator->Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
            }

            dt/=2 ;  SamplingTimestepMultiple*= 2; 
            FSIIterations*=2;
 
            p_ForceOut->SetFluidSolidIterations(FSIIterations);
            p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
            p_simulator->SetDt(dt);
            p_Mesh_modifier->SetUpdateFrequency(2/dt);


        }


       
    }

  
    void TestContinuingHomoArchieve() throw(Exception)
   {
        std::string Archieved = "DeformingPlexus/BendingForceOnly/MoreCorse3/";
        std::string output_dir = "DeformingPlexus/BendingForceOnly/Hetero_FSI/";

        double SamplingStep = 25;
        double dt = 0.0001*5;
        double RemeshingTime = 50;
        double EdgeLength =0.0004; // 0.0003;
 

        double EndTime = 4.3;
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);

        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
        // p_simulator->RemoveAllForces();

        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */
        // double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        // double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        // boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        // p_ForceOut->SetPressure(2 * (P_blood - P_tissue) / 3);
        // p_simulator->AddForce(p_ForceOut);

        /*
        -----------------------------
        Update membrane properties
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 



        for (int j = 0; j < 10; j++)
        {
            for (int i = 0; i <= 2; i++)
            {
                static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);
                EndTime += 0.5;
                p_simulator->SetEndTime(EndTime);
                p_simulator->Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
            }
            // (p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(0.0003);
            // dt /= 2;
            // SamplingStep *= 2;
            // RemeshingTime /= 2;
            // p_simulator->SetSamplingTimestepMultiple(SamplingStep);
            // p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 
            // p_simulator->SetDt(dt);
        }
    }

    



};

#endif /*TESTRELAXATION_HPP_*/
