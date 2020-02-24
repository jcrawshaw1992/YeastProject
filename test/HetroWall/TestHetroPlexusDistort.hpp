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


        std::string mesh_file = "projects/VascularRemodelling/test/data/embryo_plexus/config.vtu";

        // This data file is in mm??
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> p_mesh;
        p_mesh.ConstructFromMeshReader(mesh_reader);
        p_mesh.Scale(1.0 * scale, 1.0 * scale, 1.0 * scale); // so distances are back in original scale
        
        std::string output_directory = "HeterogeneousPlexus/RunningOriginalPlexusDistort2/";

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
        simulator.SetEndTime(200);
        simulator.SetDt(0.002); // 0.005
        simulator.SetSamplingTimestepMultiple(100);
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
        p_ForceOut->SetPressure(-TransmuralPressure);
        simulator.AddForce(p_ForceOut);

        /*
        -----------------------------
        MembraneProperties Modifier
        ----------------------------
        */

        double BendingConst = 1e-17;
        std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
                //         KA,          Kalpha           Ks                                 Kb
        GrowthMaps[10] = Create_c_vector(pow(10, -6.9), pow(10, -8.2459), pow(10, -9), 0);
        GrowthMaps[5] = Create_c_vector(pow(10, -6.9341), pow(10, -7.7), pow(10, -8), BendingConst);
        GrowthMaps[1.5] = Create_c_vector(pow(10, -6.5), pow(10, -6.3491), pow(10, -7), BendingConst);
        GrowthMaps[1.2] = Create_c_vector( 0,0,0, BendingConst) ;
        GrowthMaps[2] = Create_c_vector(pow(10, -6.8), pow(10, -6.8124), pow(10, -7), BendingConst);

//;
        boost::shared_ptr<MembranePropertiesSecModifier<2, 3> > p_Membrane_modifier(new MembranePropertiesSecModifier<2, 3>());
        // p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 2, 1,1e-10, 1); 
        p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 2, 1,10, 1); 

        c_vector<double, 3> UpperPlaneNormal = Create_c_vector(-0.8,-0.5,0.2);
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0,0,0);

        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0.14,0.192,0);
        c_vector<double, 3> LowerPlaneNormal = Create_c_vector(0.8,-0.5,0.1);
        p_Membrane_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);
        p_Membrane_modifier->SetBendingForce(cell_population, BendingConst);

        simulator.AddSimulationModifier(p_Membrane_modifier);


        // /*
        // // -----------------------------
        // // SMembrane forces
        // // ----------------------------
        // // */
        boost::shared_ptr<MembraneForcesBasic> p_shear_force(new MembraneForcesBasic());
        p_shear_force->SetupMembraneConfiguration(cell_population);
        simulator.AddForce(p_shear_force);

        // // /*
        // // -----------------------------
        // // Bending forces
        // // ----------------------------
        // // */

        boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
        p_membrane_force->SetupInitialMembrane(p_mesh, simulator.rGetCellPopulation());
        p_membrane_force->SetMembraneStiffness(BendingConst,30,30 );
        simulator.AddForce(p_membrane_force);

        for (AbstractCellPopulation<2, 3>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {

            cell_iter->GetCellData()->SetItem("Boundary", 0);

                    cell_iter->GetCellData()->SetItem("Mut", 0);
        }

         /*
        -----------------------------
        Boundaries
        ----------------------------
        */

        std::vector<c_vector<double,3> > boundary_plane_points;
        std::vector<c_vector<double,3> > boundary_plane_normals;

        boundary_plane_points.push_back(Create_c_vector(0.2086141261846316, 0.2113633402408825, -0.011480150515788608 ));
        boundary_plane_normals.push_back(Create_c_vector(0.7977338628844153,0.6007614630245445, 0.05202257733059901  ));


        boundary_plane_points.push_back(Create_c_vector(0.16770245865140276, 0.28947813761622154 , 0.0016863697687186426));
        boundary_plane_normals.push_back(Create_c_vector(0.8812648098552295, 0.47238172061006944, -0.015094533589904778 ));


        boundary_plane_points.push_back(Create_c_vector(0.108461040582645, 0.29447132627620465, 0.0017621149095964472 ));
        boundary_plane_normals.push_back(Create_c_vector(-0.6720375227750773, 0.7390518943920804, -0.04656034125535193 ));


        boundary_plane_points.push_back(Create_c_vector(  0.0336,0.2089, -0.00116 ));
        boundary_plane_normals.push_back(Create_c_vector( -0.9978,-0.034,0.0568));

        boundary_plane_points.push_back(Create_c_vector(0.0657,0.1322,0.002439 ));
        boundary_plane_normals.push_back(Create_c_vector(-0.7264,-0.687,-0.0122 ));
        


        boundary_plane_points.push_back(Create_c_vector(0.120955,0.0849,-0.00048 ));
        boundary_plane_normals.push_back(Create_c_vector(-0.13520599,-0.9903,-0.028949 ));

        boundary_plane_points.push_back(Create_c_vector(0.2036,0.12533,-0.00324 ));
        boundary_plane_normals.push_back(Create_c_vector(0.7438,-0.667,0.03332 ));
  
        
        
    

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