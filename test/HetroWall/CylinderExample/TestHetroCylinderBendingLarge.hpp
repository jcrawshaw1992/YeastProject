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
#include "OutwardsPressureWithBreaks.hpp"
// #include "StepHeteroModifier.hpp"
#include "StepHeteroModifier.hpp"
#include "MembraneBendingForce.hpp"

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:

    void offTestRunningArchieve2() throw(Exception)
    {
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load("StepChangeHetroCylinder", 13);


        double dt= 0.001;
        double NewEndTime = 13;
        double EndTime = 13;
        
        double SamplingTimestepMultiple = 100;
        std::string output_dir = "StepChangeHetroCylinder/CollapseWithBending/LargeBending/";


        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

        // p_simulator->SetEndTime(EndTime + NewEndTime);
        p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
    
     

        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_membrane_force->SetMembraneStiffness(pow(10, 8));
        p_simulator->AddForce(p_membrane_force);

        /* 
        -----------------------------
        Update membrane properties
        ----------------------------
        */

        c_vector<double, 3> UpperPlaneNormal =  Create_c_vector(0, 0, 1);
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0, 0, 15e-6 * 1e3);
        c_vector<double, 3> LowerPlaneNormal =  Create_c_vector(0, 0, -1);
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0, 0,25e-6 * 1e3);

        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<StepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<StepHeteroModifier<2, 3> >(*iter);
        p_Mesh_modifier->Boundaries(UpperPlaneNormal,UpperPlanePoint,  LowerPlaneNormal,LowerPlanePoint );
        p_Mesh_modifier->SetUpdateFrequency(0.5/dt);



          for (int i =1; i<=15; i++)
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

    
        

        p_simulator->Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);


    }



    void TestRunningArchieve2() throw(Exception)
    {

        std::string output_dir = "StepChangeHetroCylinder/CollapseWithBending/LargeBending/";
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(output_dir, 28);


        double dt= 0.001;
        double NewEndTime = 28;
        double EndTime = 28;
        
        double SamplingTimestepMultiple = 500;

        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

        // p_simulator->SetEndTime(EndTime + NewEndTime);
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



    // void TestRunningArchieve3() throw(Exception)
    // {

    //     // std::string Archieved = ;
    //     OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load("StepChangeHetroCylinder", 7);
    //     // Load and fix any settings in the simulator

   
    //     double DilationParameter=8.3;
    //     double AreaParameter=7.3;
    //     double DeformationParamter=8.1;

     
    //     double NewEndTime = 1;
    //     double EndTime = 7;
        
    //     double SamplingTimestepMultiple = 100;
    //     std::string output_dir = "StepChangeHetroCylinder/";
        
    //     /* Update the ouput directory for the population  */
    //     static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
    //     static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

    //     p_simulator->SetEndTime(EndTime + NewEndTime);
    //     p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
    //     p_simulator->SetOutputDirectory(output_dir);
    
    //     /* 
    //     -----------------------------
    //     Update membrane properties
    //     ----------------------------
    //     */
    //     std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
    //     boost::shared_ptr<StepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<StepHeteroModifier<2, 3> >(*iter);
    //     std::map<double, c_vector<long double, 4> > GrowthMaps;
    //     GrowthMaps[1] = Create_c_vector(pow(10, -AreaParameter), pow(10, -DilationParameter), pow(10, -DeformationParamter), 0);
    //     GrowthMaps[0] = Create_c_vector(pow(10, -6), pow(10, -6), pow(10, -6), 0);
    //     // Strength,hetro,stepsize, setupsolve
    //     // GrowthMaps, Strength, Hetrogeneous,  StepSize,SetupSolve
    //     p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 0, 1e-18, 1);
    //     p_simulator->Solve();
    //     CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);


    // }


    //  void TestIntroduceHetro1() throw(Exception)
    // {
    //     std::string Archieved = "StepChangeHetroCylinder/";
    //     OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, 8);

    //     c_vector<double, 3> UpperPlaneNormal =  Create_c_vector(0, 0, 1);
    //     // c_vector<double, 3> UpperPlanePoint = Create_c_vector(12e-6 * scale, 0, 0);
    //     c_vector<double, 3> UpperPlanePoint = Create_c_vector(0, 0, 15e-6 * 1e3);
    //     c_vector<double, 3> LowerPlaneNormal =  Create_c_vector(0, 0, -1);
    //     c_vector<double, 3> LowerPlanePoint = Create_c_vector(0, 0,25e-6 * 1e3);

    //     double DilationParameter=8.4;
    //     double AreaParameter=7;
    //     double DeformationParamter=8;

     
    //     double dt= 0.001;
    //     double NewEndTime = 100;
    //     double EndTime = 8;
        
    //     double SamplingTimestepMultiple = 100;
    //     std::string output_dir = "StepChangeHetroCylinder/B_2/";
        
    //     /* Update the ouput directory for the population  */
    //     static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
    //     static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

    //     p_simulator->SetEndTime(EndTime + NewEndTime);
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
    //     std::map<double, c_vector<long double, 4> > GrowthMaps;
    //     GrowthMaps[1] = Create_c_vector(pow(10, -AreaParameter), pow(10, -DilationParameter), pow(10, -DeformationParamter), 0);
    //     GrowthMaps[0] = Create_c_vector(pow(10, -6), pow(10, -6), pow(10, -5), 0);
    //     p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 1, 1e-10, 1);
    //     p_Mesh_modifier->Boundaries(UpperPlaneNormal,UpperPlanePoint,  LowerPlaneNormal,LowerPlanePoint );
    //     p_Mesh_modifier->SetBasementMembraneStrength(0);
    //     p_Mesh_modifier->SetPlateauParameters(8, 2);
    //     p_Mesh_modifier->SetUpdateFrequency(100);

    //     p_simulator->Solve();
    //     CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);

    // }


    //  void offTestIntroduceHetro2() throw(Exception)
    // {
    //     std::string Archieved = "StepChangeHetroCylinder/";
    //     OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, 8);

    //     c_vector<double, 3> UpperPlaneNormal =  Create_c_vector(0, 0, 1);
    //     // c_vector<double, 3> UpperPlanePoint = Create_c_vector(12e-6 * scale, 0, 0);
    //     c_vector<double, 3> UpperPlanePoint = Create_c_vector(0, 0, 15e-6 * 1e3);
    //     c_vector<double, 3> LowerPlaneNormal =  Create_c_vector(0, 0, -1);
    //     c_vector<double, 3> LowerPlanePoint = Create_c_vector(0, 0,25e-6 * 1e3);

    //     double DilationParameter=8.4;
    //     double AreaParameter=7;
    //     double DeformationParamter=8;

     
    //     double dt= 0.001;
    //     double NewEndTime = 100;
    //     double EndTime = 8;
        
    //     double SamplingTimestepMultiple = 100;
    //     std::string output_dir = "StepChangeHetroCylinder/B_0.5/";
        
    //     /* Update the ouput directory for the population  */
    //     static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
    //     static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

    //     p_simulator->SetEndTime(EndTime + NewEndTime);
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
    //     std::map<double, c_vector<long double, 4> > GrowthMaps;
    //     GrowthMaps[1] = Create_c_vector(pow(10, -AreaParameter), pow(10, -DilationParameter), pow(10, -DeformationParamter), 0);
    //     GrowthMaps[0] = Create_c_vector(pow(10, -6), pow(10, -6), pow(10, -5), 0);
    //     p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 1, 1e-10, 1);
    //     p_Mesh_modifier->Boundaries(UpperPlaneNormal,UpperPlanePoint,  LowerPlaneNormal,LowerPlanePoint );
    //     p_Mesh_modifier->SetBasementMembraneStrength(0);
    //     p_Mesh_modifier->SetPlateauParameters(8, 0.5);
    //     p_Mesh_modifier->SetUpdateFrequency(100);

    //     p_simulator->Solve();
    //     CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);

    // }



    //  void offTestIntroduceHetro3() throw(Exception)
    // {
    //     std::string Archieved = "StepChangeHetroCylinder/";
    //     OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, 8);

    //     c_vector<double, 3> UpperPlaneNormal =  Create_c_vector(0, 0, 1);
    //     // c_vector<double, 3> UpperPlanePoint = Create_c_vector(12e-6 * scale, 0, 0);
    //     c_vector<double, 3> UpperPlanePoint = Create_c_vector(0, 0, 15e-6 * 1e3);
    //     c_vector<double, 3> LowerPlaneNormal =  Create_c_vector(0, 0, -1);
    //     c_vector<double, 3> LowerPlanePoint = Create_c_vector(0, 0,25e-6 * 1e3);

    //     double DilationParameter=8.4;
    //     double AreaParameter=7;
    //     double DeformationParamter=8;

     
    //     double dt= 0.001;
    //     double NewEndTime = 100;
    //     double EndTime = 8;
        
    //     double SamplingTimestepMultiple = 100;
    //     std::string output_dir = "StepChangeHetroCylinder/B_4/";
        
    //     /* Update the ouput directory for the population  */
    //     static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
    //     static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

    //     p_simulator->SetEndTime(EndTime + NewEndTime);
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
    //     std::map<double, c_vector<long double, 4> > GrowthMaps;
    //     GrowthMaps[1] = Create_c_vector(pow(10, -AreaParameter), pow(10, -DilationParameter), pow(10, -DeformationParamter), 0);
    //     GrowthMaps[0] = Create_c_vector(pow(10, -6), pow(10, -6), pow(10, -5), 0);
    //     p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 1, 1e-10, 1);
    //     p_Mesh_modifier->Boundaries(UpperPlaneNormal,UpperPlanePoint,  LowerPlaneNormal,LowerPlanePoint );
    //     p_Mesh_modifier->SetBasementMembraneStrength(0);
    //     p_Mesh_modifier->SetPlateauParameters(8, 4);
    //     p_Mesh_modifier->SetUpdateFrequency(100);

    //     p_simulator->Solve();
    //     CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);

    // }




};

#endif /*TESTRELAXATION_HPP_*/
