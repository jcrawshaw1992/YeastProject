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
#include "OutwardsPressureWithBreaks.hpp"
#include "RemeshingTriggerOnHeteroMeshModifier.hpp"

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:
     void TestIntroduceHetro3() throw(Exception)
    {
        std::string Archieved = "HetroCylinder/";
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, 8);

        c_vector<double, 3> UpperPlaneNormal =  Create_c_vector(0, 0, 1);
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0, 0, 15e-6 * 1e3);
        c_vector<double, 3> LowerPlaneNormal =  Create_c_vector(0, 0, -1);
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0, 0,25e-6 * 1e3);

        double DilationParameter=8.4;
        double AreaParameter=7;
        double DeformationParamter=8;

        double dt= 0.001;
        double NewEndTime = 1;
        double EndTime = 8;
        
        double SamplingTimestepMultiple = 100;
    
        std::string output_dir = "HetroCylinder/B_0.0002/";
        
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

        p_simulator->SetEndTime(EndTime + NewEndTime);
        p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);

        /* 
        -----------------------------
        Update membrane properties
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnHeteroMeshModifier<2, 3> >(*iter);
        std::map<double, c_vector<long double, 4> > GrowthMaps;
        GrowthMaps[1] = Create_c_vector(pow(10, -AreaParameter), pow(10, -DilationParameter), pow(10, -DeformationParamter), 0);
        GrowthMaps[0] = Create_c_vector(pow(10, -7), pow(10, -7), pow(10, -5), 0);
        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 1, 1e-10, 1);
        p_Mesh_modifier->Boundaries(UpperPlaneNormal,UpperPlanePoint,  LowerPlaneNormal,LowerPlanePoint );
        p_Mesh_modifier->SetBasementMembraneStrength(0);
        p_Mesh_modifier->SetPlateauParameters(8,0.0002);
        p_Mesh_modifier->SetUpdateFrequency(50);

        p_simulator->Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);


        // ========================

        dt= 0.0001;
        SamplingTimestepMultiple = 1000;
        NewEndTime = EndTime + NewEndTime;
        p_Mesh_modifier->SetUpdateFrequency(500);
        p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        p_simulator->SetDt(dt);

        for (int i =1; i<=3; i++)
        { 
            static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
            NewEndTime +=1;
            p_simulator->SetEndTime(NewEndTime);
            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }


        dt= 0.00001;
        SamplingTimestepMultiple = 10000;
        p_Mesh_modifier->SetUpdateFrequency(5000);
        p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        p_simulator->SetDt(dt);

        for (int i =1; i<=4; i++)
        { 
            static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
            NewEndTime +=1;
            p_simulator->SetEndTime(NewEndTime);
            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }

        dt= 0.000001;
        SamplingTimestepMultiple = 100000;
        p_Mesh_modifier->SetUpdateFrequency(50000);
        p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        p_simulator->SetDt(dt);

        for (int i =1; i<=5; i++)
        { 
            static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
            NewEndTime +=1;
            p_simulator->SetEndTime(NewEndTime);
            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }

    }







};

#endif /*TESTRELAXATION_HPP_*/
