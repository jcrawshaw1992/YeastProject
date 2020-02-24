
#ifndef TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_CORRECTEDDRAG_HPP_
#define TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_CORRECTEDDRAG_HPP_

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <cstdio>
#include <ctime>
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

#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"

#include "MutableMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "MembraneForcesBasic.hpp"
#include "MembranePropertiesModifier.hpp"

#include "CellMutationStatesWriter.hpp"
#include "EmptyBasementMatrix.hpp"
#include "HasEndothelialCell.hpp"
#include "LostEndothelialCell.hpp"
#include "OutwardsPressure.hpp"

#include <vector>
#include "AppliedForce.hpp"
#include "AppliedForceModifier.hpp"
#include "VtkMeshWriter.hpp"
#include <iostream>
#include <fstream>


static const double M_TIME_FOR_SIMULATION = 0.01; //40; //50
static const double M_SAMPLING_TIME_STEP = 100; //50
static const double M_TIME_STEP = 0.002;


class ScalledPlexus : public AbstractCellBasedTestSuite
{

public:

    void TestAreaFroceDragCorrectedEqui() throw(Exception)
    {

        double mesh_scale = 1;//1e-3;

        // std::string mesh_file = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/Bifurcation/SetUpData/config.vtu";
        // std::string mesh_file = "/Users/jcrawshaw/docker-polnet-master/ScalledMesh/Plexus.vtu";
        std::string mesh_file = "/Users/jcrawshaw/docker-polnet-master/ScalledMesh/Clipped.vtu";
        
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
    
      
        std::stringstream out;
        std::string output_directory = "ScalledPlexus/Third/"; //Shrinking"; // + Parameters + "/";

        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        MAKE_PTR(BetaCateninOneHitCellMutationState, p_state); //Mutation to mark nodes

        // MAKE_PTR(WildTypeCellMutationState, p_WildTypeState); //Mutation to mark nodes


        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        MeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(500); //(M_TIME_FOR_SIMULATION);
        simulator.SetDt(0.02);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetUpdateCellPopulationRule(false); // No remeshing.

        //  -----------------------------
        //  Save original mesh before I shrink it
        // ----------------------------

        VtkMeshWriter<2, 3> mesh_writer(output_directory, "Original", false);
        MutableMesh<2, 3>* p_mesh = &(dynamic_cast<MeshBasedCellPopulation<2, 3>*>(&(simulator.rGetCellPopulation()))->rGetMesh());
        p_mesh->Scale(1.0 / mesh_scale, 1.0 / mesh_scale, 1.0 / mesh_scale); // so distances are back in original scal
        mesh_writer.WriteFilesUsingMesh(*p_mesh);

        // //  -----------------------------
        // //  Parameters
        // //----------------------------

        double ElasticShearModulus = 4.4e-9;
        double AreaDilationModulus = 0.9e-11;
        double membrane_constant = 5;
        // double Area_constant = 0.9e-15;
 
        /*
        -----------------------------
        Constant Pressure ballance 
        ----------------------------
        */
        double P_blood = 0.2133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.1466542; // Pa == 1.1000e-05 mmHg

        double TransmuralPressure =  P_blood - P_tissue;

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure(TransmuralPressure);
        simulator.AddForce(p_ForceOut);



        /*
        -----------------------------
        MembraneProperties Modifier
        ----------------------------
        */

        double BendingConst = 0.001;
        std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
                //         KA,          Kalpha           Ks                                 Kb
        GrowthMaps[10] = Create_c_vector(pow(10, -6.9), pow(10, -8.2459), pow(10, -9), 0);
        GrowthMaps[8] = Create_c_vector(pow(10, -6.9), pow(10, -8.0160), pow(10, -9), BendingConst);
        GrowthMaps[6] = Create_c_vector(pow(10, -6.9), pow(10, -7.7300), pow(10, -9), BendingConst);
        GrowthMaps[5] = Create_c_vector(pow(10, -6.9341), pow(10, -7.7), pow(10, -8), BendingConst);
        GrowthMaps[4] = Create_c_vector(pow(10, -6.9), pow(10, -7.4224), pow(10, 8), BendingConst);

        
        GrowthMaps[2] = Create_c_vector(pow(10, -7.15), pow(10, -9), pow(10, -9), BendingConst);
        GrowthMaps[1.5] = Create_c_vector(pow(10, -6.5), pow(10, -6.3491), pow(10, -7), BendingConst);
        GrowthMaps[1.2] = Create_c_vector(pow(10, -6.2), pow(10, -5.8360), pow(10, -7), BendingConst);


        boost::shared_ptr<MembranePropertiesModifier<2, 3> > p_Membrane_modifier(new MembranePropertiesModifier<2, 3>());
        p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 2, 0,1e-10); 
        p_Membrane_modifier->SetupSolve(cell_population, output_directory);

        // /*
        // -----------------------------
        // SMembrane forces
        // ----------------------------
        // */
        boost::shared_ptr<MembraneForcesBasic> p_shear_force(new MembraneForcesBasic());
        p_shear_force->SetupMembraneConfiguration(cell_population);
        simulator.AddForce(p_shear_force);
       
 


        /*
        -----------------------------
        Boundaries
        ----------------------------
        */

        //Create a plane boundary to represent the inlet and pass them to the simulation

        // c_vector<long double, 3> Boundary1 = Create_c_vector(12.48824184355049, 34.75302061558864, 41.78949821113195);
        // c_vector<long double, 3> Normal1 = Create_c_vector(-0.05607774749413225, 0.762765339938692, 0.6442393362906335);
        // double Radius1 = 1.7;

        // c_vector<long double, 3> Boundary2 = Create_c_vector(12.597373380655702, 48.382440094438316, 42.984851419357064);
        // c_vector<long double, 3> Normal2 = Create_c_vector(-0.04847413662454751, -0.989768366942236, -0.13419701143842153);
        // double Radius2 = 1.2;

        // c_vector<long double, 3> Boundary3 = Create_c_vector(-0.04710924565104161, -0.9898343148598956, -0.13419667693363677);
        // c_vector<long double, 3> Normal3 = Create_c_vector(-0.04710924565104161, -0.9898343148598956, -0.13419667693363677);
        // double Radius3 = 1.3;

        // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary1, Normal1, Radius1));
        // simulator.AddCellPopulationBoundaryCondition(p_condition_1);

        // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary2, Normal2, Radius2));
        // simulator.AddCellPopulationBoundaryCondition(p_condition_2);
        // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_3(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary3, Normal3, Radius3));
        // simulator.AddCellPopulationBoundaryCondition(p_condition_3);

          simulator.Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);

        // To reset before looping: this is usually done by the SetUp and TearDown methods
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
    }
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/
