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
// #include "VtkMeshReader.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"

#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "HistoryDepMutableMesh.hpp"


#include "FixedRegionBoundaryCondition.hpp"
#include "MembraneDeformationForceOnCylinder.hpp"
#include "OutwardsPressure.hpp"
#include "StepHeteroModifier.hpp"
#include "MembraneBendingForce.hpp"

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:


    void offTestRunningArchieve2() throw(Exception)
    {
        std::string Archieved = "StepChangeHetroCylinder";
        double EndTime = 13;//9;
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved,EndTime);

        double dt= 0.01;
        double NewEndTime = 1;
     
        double SamplingTimestepMultiple = 100;
        std::string output_dir = "StepChangeHetroCylinder/CollapseWithBending/LargeBending/";
        
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

        p_simulator->RemoveAllForces();
        p_simulator->SetEndTime(EndTime + NewEndTime);
        p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
    

        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForceOnCylinder> p_shear_force(new MembraneDeformationForceOnCylinder());
        p_simulator->AddForce(p_shear_force);


        /*
        -----------------------------
        Outwards force
        ----------------------------
        */


        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        boost::shared_ptr<OutwardsPressure> p_force(new OutwardsPressure());
        p_force->SetPressure(P_blood - P_tissue);
        // p_simulator->AddForce(p_force);


        p_simulator->Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);


    }



    void TestRunningArchieve2() throw(Exception)
    {

        double EndTime = 36;
        std::string output_dir = "StepChangeHetroCylinder/CollapseWithBending/SmallBending/";
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(output_dir, EndTime);

        double dt= 0.001;
        double NewEndTime = EndTime;
        
        
        double SamplingTimestepMultiple = 500;

        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

        p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);

        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<StepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<StepHeteroModifier<2, 3> >(*iter);
        p_Mesh_modifier->SetUpdateFrequency(0.5/dt);

          for (int i =1; i<=50; i++)
        { 
            static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
            NewEndTime +=1;
            p_simulator->SetEndTime(NewEndTime);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }

        dt/= 10;SamplingTimestepMultiple *= 10;
        p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
        p_Mesh_modifier->SetUpdateFrequency(0.5/dt);

        for (int i =1; i<=40; i++)
        { 
            static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
            NewEndTime +=4;
            p_simulator->SetEndTime(NewEndTime);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }



    }




    // void TestIntroduceHetro() throw(Exception)
    // {
    //     std::string Archieved = "StepChangeHetroCylinder";
    //     double EndTime = 5;//9;
    //     OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved,EndTime);

    //     c_vector<double, 3> UpperPlaneNormal =  Create_c_vector(0, 0, 1);
    //     c_vector<double, 3> UpperPlanePoint = Create_c_vector(0, 0, 15e-6 * 1e3);
    //     c_vector<double, 3> LowerPlaneNormal =  Create_c_vector(0, 0, -1);
    //     c_vector<double, 3> LowerPlanePoint = Create_c_vector(0, 0,25e-6 * 1e3);

    //     double dt= 0.001;
    //     double NewEndTime = EndTime;
    
    //     double SamplingTimestepMultiple = 100;
    
    //     std::string output_dir = "StepChangeHetroCylinder/CollapseWithBending/LargeBending/";
        
    //     /* Update the ouput directory for the population  */
    //     static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
    
    //     p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
    //     p_simulator->SetDt(dt);
    //     p_simulator->SetOutputDirectory(output_dir);

    //     /* 
    //     -----------------------------
    //     Update membrane properties
    //     ----------------------------
    //     */
    //     std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
    //     boost::shared_ptr<StepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<StepHeteroModifier<2, 3> >(*iter);
    //     p_Mesh_modifier->Boundaries(UpperPlaneNormal,UpperPlanePoint,  LowerPlaneNormal,LowerPlanePoint );
    //     // std::map<double, c_vector<long double, 4> >  GrowthMaps =  { {1, Create_c_vector(pow(10, -7), pow(10, -8.4), pow(10, -8), 1e-14) },
    //     // {0, Create_c_vector(pow(10, -7), pow(10, -8.4), pow(10, -5), 1e-14)}    };
    //     // p_Mesh_modifier->SetMembranePropeties( GrowthMaps, 1);

    //     p_Mesh_modifier->SetUpdateFrequency(0.5/dt);


    //     boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
    //     p_membrane_force->SetMembraneStiffness(pow(10, -9));
    //     p_simulator->AddForce(p_membrane_force);


    //     for (int i =1; i<=5; i++)
    //     { 
    //         static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
    //         NewEndTime +=2;
    //         p_simulator->SetEndTime(NewEndTime);

    //         p_simulator->Solve();
    //         CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
    //     }

    //     dt= 0.001;SamplingTimestepMultiple = 500;
    
    //     p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
    //     p_simulator->SetDt(dt);
    //     p_simulator->SetOutputDirectory(output_dir);
    //     p_Mesh_modifier->SetUpdateFrequency(0.5/dt);

    //     for (int i =1; i<=40; i++)
    //     { 
    //         static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
    //         NewEndTime +=4;
    //         p_simulator->SetEndTime(NewEndTime);

    //         p_simulator->Solve();
    //         CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
    //     }




    // }

    //     void offTestIntroduceHetro2() throw(Exception)
    // {
    //     std::string Archieved = "StepChangeHetroCylinder/Collapse/";
    //     double EndTime = 65;
    //     OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved,EndTime);


    //     double dt= 0.0001;
    //     double NewEndTime = EndTime;
    
    //     double SamplingTimestepMultiple = 5000;
    
    //     std::string output_dir = "StepChangeHetroCylinder/Collapse/";
        
    //     /* Update the ouput directory for the population  */
    //     static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
    
    //     p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
    //     p_simulator->SetDt(dt);
    //     p_simulator->SetOutputDirectory(output_dir);

    //     /* 
    //     -----------------------------
    //     Update membrane properties
    //     ----------------------------
    //     */
    //     std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
    //     boost::shared_ptr<StepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<StepHeteroModifier<2, 3> >(*iter);
    //     p_Mesh_modifier->SetUpdateFrequency(0.5/dt);


    //     for (int i =1; i<=40; i++)
    //     { 
    //         static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
    //         NewEndTime +=4;
    //         p_simulator->SetEndTime(NewEndTime);

    //         p_simulator->Solve();
    //         CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
    //     }


    //     dt= 0.00001;SamplingTimestepMultiple = 5000;
    
    //     p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
    //     p_simulator->SetDt(dt);
    //     p_simulator->SetOutputDirectory(output_dir);
    //     p_Mesh_modifier->SetUpdateFrequency(0.5/dt);


    //     for (int i =1; i<=40; i++)
    //     { 
    //         static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
    //         NewEndTime +=4;
    //         p_simulator->SetEndTime(NewEndTime);

    //         p_simulator->Solve();
    //         CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
    //     }



    // }


    //   void offTestIntroduceHetro() throw(Exception)
    // {
    //     std::string Archieved = "StepChangeHetroCylinder/Collapse/";
    //     double EndTime = 85;
    //     OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved,EndTime);


    //     double dt= 0.0004;
    //     double NewEndTime = EndTime;
    
    //     double SamplingTimestepMultiple = 10;
    
    //     std::string output_dir = "StepChangeHetroCylinder/Collapse/";
        
    //     /* Update the ouput directory for the population  */
    //     static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
    
    //     p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
    //     p_simulator->SetDt(dt);
    //     p_simulator->SetOutputDirectory(output_dir);

    //     /* 
    //     -----------------------------
    //     Update membrane properties
    //     ----------------------------
    //     */
    //     std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
    //     boost::shared_ptr<StepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<StepHeteroModifier<2, 3> >(*iter);
    //     p_Mesh_modifier->SetUpdateFrequency(50/dt);


    //     for (int i =1; i<=40; i++)
    //     { 
    //         static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
    //         NewEndTime +=4;
    //         p_simulator->SetEndTime(NewEndTime);

    //         p_simulator->Solve();
    //         CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
    //     }


    //     dt= 0.00001;SamplingTimestepMultiple = 5000;
    
    //     p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
    //     p_simulator->SetDt(dt);
    //     p_simulator->SetOutputDirectory(output_dir);
    //     p_Mesh_modifier->SetUpdateFrequency(0.5/dt);


    //     for (int i =1; i<=40; i++)
    //     { 
    //         static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
    //         NewEndTime +=4;
    //         p_simulator->SetEndTime(NewEndTime);

    //         p_simulator->Solve();
    //         CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
    //     }



    // }





};

#endif /*TESTRELAXATION_HPP_*/