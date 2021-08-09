#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "PetscSetupAndFinalize.hpp"

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"

#include "Debug.hpp"
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
#include "MembraneBendingForce.hpp"
#include "MembraneDeformationForce.hpp"
#include "StaticOutwardsPressure.hpp"
#include "OutwardsPressureWithBreaks.hpp"
#include "OutwardsPressure.hpp"
#include "MembraneBendingForce0TargetAngle.hpp"


#include "RemeshingTriggerOnStepHeteroModifier.hpp"

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:


 void TestRemesh() throw(Exception)
   {


         TRACE("Jess is good")
        double EndTime = 11;
        double scale = 0.05;
        double SamplingStep = 25;
        double dt = 0.0001;
        double RemeshingTime = 50000;
        double EdgeLength = 0.0005/2;//(2e-6 * scale);

        
        std::string output_dir = "DeformingHoneyComb/RemeshingStep/";
        std::string Archieved =  "DeformingHoneyComb/FlatForce3";


        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
 
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetPrintRemeshedIC(1);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetInitialAnlgesAcrossMembrane();
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).ExecuteHistoryDependentRemeshing();
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingIterations(10);

        p_simulator->RemoveAllForces();
        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
        p_simulator->RemoveAllCellPopulationBoundaryConditions();

        EndTime +=1;
        p_simulator->SetEndTime(EndTime);

        p_simulator->Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
    }


  
};

#endif /*TESTRELAXATION_HPP_*/
