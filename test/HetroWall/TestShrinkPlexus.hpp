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
    void TestHetroPlexus() throw(Exception)
    {
        double scale = 1e-3; // 
        // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
        // in um will be too large and break chaste without carefull playing with or a tiny time step


        // std::string mesh_file = "projects/VascularRemodelling/test/data/ShrunkPlexus/NewClipped2.vtu";
        std::string mesh_file = "/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/test/data/embryo_plexus/config.vtu";
        // This data file is in mm??
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> p_mesh;
        p_mesh.ConstructFromMeshReader(mesh_reader);
        p_mesh.Scale(1.0 * scale, 1.0 * scale, 1.0 * scale); // so distances are back in original scale
        
        std::string output_directory = "ShrinkingPlexus/LargestBendingForce/";

        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        MeshBasedCellPopulation<2, 3> cell_population(p_mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        // cell_population.AddCellWriter<CellMutationStatesWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(20);
        simulator.SetDt(0.001); // 0.005
        simulator.SetSamplingTimestepMultiple(200);
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
        p_ForceOut->SetPressure(TransmuralPressure);
        simulator.AddForce(p_ForceOut);

        /*
        -----------------------------
        MembraneProperties Modifier
        ----------------------------
        */

        double BendingConst = 1e-10;
        std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
                //         KA,          Kalpha           Ks                                 Kb
        GrowthMaps[10] = Create_c_vector(pow(10, -6.9), pow(10, -8.2459), pow(10, -9), 0);
        GrowthMaps[5] = Create_c_vector(pow(10, -6.9341), pow(10, -7.7), pow(10, -8), BendingConst);
        GrowthMaps[1.5] = Create_c_vector(pow(10, -6.5), pow(10, -6.3491), pow(10, -7), BendingConst);
        GrowthMaps[1.2] = Create_c_vector(pow(10, -6.2)*1e3, pow(10, -5.8360)*1e2, pow(10, -7)*1e2, BendingConst*1e1);
        GrowthMaps[2] = Create_c_vector(pow(10, -6.8), pow(10, -6.8124), pow(10, -7), BendingConst);


        boost::shared_ptr<MembranePropertiesSecModifier<2, 3> > p_Membrane_modifier(new MembranePropertiesSecModifier<2, 3>());
        // p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 2, 1,1e-10, 1); 
        p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 2, 0,10, 1); 

        c_vector<double, 3> UpperPlaneNormal = Create_c_vector(-0.8,-0.5,0.2);
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0,0,0);

        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0.14,0.192,0);
        c_vector<double, 3> LowerPlaneNormal = Create_c_vector(0.8,-0.5,0.1);

        
        p_Membrane_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);
        p_Membrane_modifier->SetBendingForce(cell_population, BendingConst);
        simulator.AddSimulationModifier(p_Membrane_modifier);
    
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
        Bending forces
        ----------------------------
        */

        boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
        p_membrane_force->SetupInitialMembrane(p_mesh, simulator.rGetCellPopulation());
        p_membrane_force->SetMembraneStiffness(BendingConst,30,30 );
        simulator.AddForce(p_membrane_force);

        for (AbstractCellPopulation<2, 3>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_iter->GetCellData()->SetItem("Boundary", 0);
        }

  
        /*
        -----------------------------
        Boundaries
        ----------------------------
        */

        std::vector<c_vector<double,3> > boundary_plane_points;
     
        std::vector<c_vector<double,3> > boundary_plane_normals;

        boundary_plane_points.push_back(Create_c_vector(0.208,0.21638,-0.0010));
        boundary_plane_normals.push_back(Create_c_vector(0.7945,0.60,-0.000055));


        boundary_plane_points.push_back(Create_c_vector(0.1699077963347786, 0.29207 , -0.0006817));
        boundary_plane_normals.push_back(Create_c_vector(0.878339963, 0.4746312, 0.0569565));


        boundary_plane_points.push_back(Create_c_vector(0.108050897, 0.300792494,-0.0023384 ));
        boundary_plane_normals.push_back(Create_c_vector(-0.68914432,0.72444911,-0.0159 ));

        boundary_plane_points.push_back(Create_c_vector( 0.032559564910557455, 0.19958179421238176 , -0.00022190318136226447));
        boundary_plane_normals.push_back(Create_c_vector( -0.9997308436644144, -0.022921836854266218, -0.003581846069338259  ));

        boundary_plane_points.push_back(Create_c_vector(0.0665665,0.1280056,0.00476316));
        boundary_plane_normals.push_back(Create_c_vector(-0.7375899212325719, -0.6727617086929154, -0.057903293625814305));
        // boundary_plane_points.push_back(Create_c_vector(0.0657,0.1322,0.002439 ));
        // boundary_plane_normals.push_back(Create_c_vector(-0.7264,-0.687,-0.0122 ));

        // boundary_plane_points.push_back(Create_c_vector(0.1169868469286675,0.08244644665474296, -0.0020484454135636133 ));
        // boundary_plane_normals.push_back(Create_c_vector(-0.13392910018912516, -0.9907740752468547,-0.02072988039693599  ));

        boundary_plane_points.push_back(Create_c_vector(0.1170383788089726,0.08236980144490846, -0.0025995520652228256));
        boundary_plane_normals.push_back(Create_c_vector(-0.14907205374285917,-0.9886401838485218,-0.019186184426487517));

        boundary_plane_points.push_back(Create_c_vector(0.2531,0.1776,-0.00665 ));
        boundary_plane_normals.push_back(Create_c_vector(0.74382,-0.667,0.0365));

        boundary_plane_points.push_back(Create_c_vector(0.19766687695542723, 0.23294767634236974 , -0.00030336941283469883 ));
        boundary_plane_normals.push_back(Create_c_vector(0.810720950328382, 0.5847206536888926, 0.028866898833778264));

        boundary_plane_points.push_back(Create_c_vector(0.1232435688, 0.08267,-0.00302008 ));
        boundary_plane_normals.push_back(Create_c_vector(-0.12537,-0.9913,-0.0398008));
  

        boundary_plane_points.push_back(Create_c_vector(0.2104411508338335, 0.13041892970979158, 0.005931369791122276));
        boundary_plane_normals.push_back(Create_c_vector(0.7410281143959528,-0.670580874021506,0.03462116507762));


        TRACE("got all the planes")
    
    PRINT_VARIABLE(boundary_plane_normals.size())
        for(unsigned boundary_id = 0; boundary_id < boundary_plane_normals.size(); boundary_id++)
        {
            PRINT_VARIABLE(norm_2(-boundary_plane_normals[boundary_id]))
            PRINT_VECTOR(boundary_plane_points[boundary_id])
           boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population,boundary_plane_points[boundary_id],boundary_plane_normals[boundary_id],0.05));
            simulator.AddCellPopulationBoundaryCondition(p_condition);
        }




        simulator.Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);

    }


 
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/