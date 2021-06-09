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
     void TestIntroduceHetro() throw(Exception)
    {
        std::string Archieved = "HetroCylinder/";
        double EndTime = 8;
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved,EndTime);

        c_vector<double, 3> UpperPlaneNormal =  Create_c_vector(0, 0, 1);
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0, 0, 15e-6 * 1e3);
        c_vector<double, 3> LowerPlaneNormal =  Create_c_vector(0, 0, -1);
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0, 0,25e-6 * 1e3);

        double DilationParameter=8.4;
        double AreaParameter=7;
        double DeformationParamter=8;

        double dt= 0.005;
        double NewEndTime = EndTime;
    
        double SamplingTimestepMultiple = 100;
    
        std::string output_dir = "HetroCylinder/A_2/";
        
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
    
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
        GrowthMaps[0] = Create_c_vector(pow(10, -6.5), pow(10, -6), pow(10, -5), 0);
        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 1, 1e-10, 1);
        p_Mesh_modifier->Boundaries(UpperPlaneNormal,UpperPlanePoint,  LowerPlaneNormal,LowerPlanePoint );
        p_Mesh_modifier->SetBasementMembraneStrength(0);
        p_Mesh_modifier->SetPlateauParameters(2,2);
        p_Mesh_modifier->SetUpdateFrequency(50);



        std::vector<c_vector<double, 3> > boundary_plane_points;
        std::vector<c_vector<double, 3> > boundary_plane_normals;

        boundary_plane_points.push_back(UpperPlanePoint);
        boundary_plane_normals.push_back(UpperPlaneNormal);

        boundary_plane_points.push_back(LowerPlanePoint);
        boundary_plane_normals.push_back(LowerPlaneNormal);

        for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
            boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(  &(p_simulator->rGetCellPopulation())  , boundary_plane_points[boundary_id], boundary_plane_normals[boundary_id], 0.006));
            p_simulator->AddCellPopulationBoundaryCondition(p_condition);
        }

        for (int i =1; i<=40; i++)
        { 
            static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
            NewEndTime +=2;
            p_simulator->SetEndTime(NewEndTime);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }


    }


 void offTestIntroduceHetroPickedUp() throw(Exception)
    {
        std::string Archieved = "HetroCylinder/A_2/";
        double EndTime = 42;
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved,EndTime);

        double DilationParameter=8.4;
        double AreaParameter=7;
        double DeformationParamter=8;

        double dt= 0.0005; //double dt= 0.0005;
        double NewEndTime = EndTime;
    
        double SamplingTimestepMultiple = 1000;
    
        std::string output_dir = "HetroCylinder/A_2/";
        
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
    
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
        // GrowthMaps[0] = Create_c_vector(pow(10, -7), pow(10, -6), pow(10, -5), 0);
        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1, 1, 1e-10, 0);
        p_Mesh_modifier->SetBasementMembraneStrength(0);
        p_Mesh_modifier->SetUpdateFrequency(500);

        for (int i =1; i<=40; i++)
        { 
            static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
            NewEndTime +=2;
            p_simulator->SetEndTime(NewEndTime);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }


    }







};

#endif /*TESTRELAXATION_HPP_*/
