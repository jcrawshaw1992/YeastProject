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
    void TestGrowPlexusToEqui() throw(Exception)
    {
        double scale = 1e-3; // 
        // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
        // in um will be too large and break chaste without carefull playing with or a tiny time step


        std::string mesh_file = "projects/VascularRemodelling/test/data/ShrunkPlexus/Plexus2.vtu";
        
        // This data file is in mm??
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> p_mesh;
        p_mesh.ConstructFromMeshReader(mesh_reader);
        p_mesh.Scale(1.0 * scale, 1.0 * scale, 1.0 * scale); // so distances are back in original scale
        
        std::string output_directory = "HeterogeneousPlexus/SetUpArchiving2/";

        // // Need to delete the unacoumpanied nodes
                // std::vector<unsigned> NodesToDelete;
                // double NumberOfNodes = 0;
                // for (typename  MutableMesh<2, 3>::NodeIterator node_iter = p_mesh.GetNodeIteratorBegin();
                //         node_iter != p_mesh.GetNodeIteratorEnd(); ++node_iter)
                //     {
                //         // node_iter->ClearAppliedForce();
                //         std::set<unsigned>& containing_elements = node_iter->rGetContainingElementIndices();
                //         if (containing_elements.size() <1)
                //         { 
                //             NodesToDelete.push_back(node_iter->GetIndex());
                //         }
                //         NumberOfNodes +=1;
                        
                //         // p_mesh.DeleteNode(node_iter->GetIndex());
                //     }
                //     PRINT_2_VARIABLES(NumberOfNodes, NodesToDelete.size())
                //     double counter =0;
                // for (std::vector<unsigned>::iterator it = NodesToDelete.begin(); it != NodesToDelete.end(); ++it)
                //     {
                //         // PRINT_VARIABLE(*it)
                //         // p_mesh.DeleteNode(*it- counter);
                //         p_mesh.DeleteNodePriorToReMesh(*it- counter);
                //         counter +=1;

                //     }
                    

                // std::vector<unsigned> NodesToDelete2;
                //     double NumberOfNodes2= 0;
                // for (typename  MutableMesh<2, 3>::NodeIterator node_iter = p_mesh.GetNodeIteratorBegin();
                //         node_iter != p_mesh.GetNodeIteratorEnd(); ++node_iter)
                //     {
                //         // node_iter->ClearAppliedForce();
                //         std::set<unsigned>& containing_elements = node_iter->rGetContainingElementIndices();
                //         if (containing_elements.size() ==0)
                //         { 
                //             NodesToDelete2.push_back(node_iter->GetIndex());
                //         }
                //         NumberOfNodes2 +=1;
                        
                //         // p_mesh.DeleteNode(node_iter->GetIndex());
                //     }
                //     PRINT_2_VARIABLES(NumberOfNodes2, NodesToDelete2.size())
                //      counter =0;
                // for (std::vector<unsigned>::iterator it = NodesToDelete2.begin(); it != NodesToDelete2.end(); ++it)
                //     {
                //         // PRINT_VARIABLE(*it)
                //         // p_mesh.DeleteNode(*it- counter);
                //         p_mesh.DeleteNodePriorToReMesh(*it- counter);
                //         counter +=1;

                //     }

                //      std::vector<unsigned> NodesToDelete3;
                //     double NumberOfNodes3= 0;
                // for (typename  MutableMesh<2, 3>::NodeIterator node_iter = p_mesh.GetNodeIteratorBegin();
                //         node_iter != p_mesh.GetNodeIteratorEnd(); ++node_iter)
                //     {
                //         // node_iter->ClearAppliedForce();
                //         std::set<unsigned>& containing_elements = node_iter->rGetContainingElementIndices();
                //         if (containing_elements.size() ==0)
                //         { 
                //             NodesToDelete3.push_back(node_iter->GetIndex());
                //         }
                //         NumberOfNodes3 +=1;
                        
                //         // p_mesh.DeleteNode(node_iter->GetIndex());
                //     }
                //     PRINT_2_VARIABLES(NumberOfNodes3, NodesToDelete3.size())





        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        MeshBasedCellPopulation<2, 3> cell_population(p_mesh, cells);

       unsigned num_removed = 0;
       for (AbstractCellPopulation<2, 3>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
        {
    
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
            Node<3>* pNode = cell_population.rGetMesh().GetNode(node_index);

            std::set<unsigned>& containing_elements = pNode->rGetContainingElementIndices();
            if (containing_elements.size()  <1)
            { 
                // TRACE("Mark as dead")
                cell_iter->SetApoptosisTime(0);
                cell_iter->StartApoptosis();
                // PRINT_VARIABLE(cell_iter->IsDead())
            }
        }
        
        cell_population.RemoveDeadCells();
        

        // NodeMap node_map(p_mesh.GetNumAllNodes());
        // p_mesh.ReIndex(node_map);

        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);

        // cell_population.AddCellWriter<CellMutationStatesWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(0.01);
        simulator.SetDt(0.00001); // 0.005
        simulator.SetSamplingTimestepMultiple(1);
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

        // /*
        // -----------------------------
        // MembraneProperties Modifier
        // ----------------------------
        // */

        // double BendingConst = 1e-14;
        // std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
        //         //         KA,          Kalpha           Ks                                 Kb
        // GrowthMaps[10] = Create_c_vector(pow(10, -6.9), pow(10, -8.2459), pow(10, -9), 0);
        // GrowthMaps[5] = Create_c_vector(pow(10, -6.9341), pow(10, -7.7), pow(10, -8), BendingConst);
        // GrowthMaps[1.5] = Create_c_vector(pow(10, -6.5), pow(10, -6.3491), pow(10, -7), BendingConst);
        // GrowthMaps[1.2] = Create_c_vector(pow(10, -6.2), pow(10, -5.8360), pow(10, -7), BendingConst);


        // boost::shared_ptr<MembranePropertiesSecModifier<2, 3> > p_Membrane_modifier(new MembranePropertiesSecModifier<2, 3>());
        // p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 5, 0,1e-10, 1); 

        // c_vector<double, 3> UpperPlaneNormal = Create_c_vector(0.8,0.58,0.058);
        // c_vector<double, 3> UpperPlanePoint = Create_c_vector(580,520,-0.1)*scale;

        // c_vector<double, 3> LowerPlaneNormal = Create_c_vector(0.9,0.45,0.054);
        // c_vector<double, 3> LowerPlanePoint = Create_c_vector(270.7, 597,37.8)*scale;
        // p_Membrane_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);
        // p_Membrane_modifier->SetBendingForce(cell_population, BendingConst);

        // simulator.AddSimulationModifier(p_Membrane_modifier);


        // /*
        // -----------------------------
        // SMembrane forces
        // ----------------------------
        // */
        // boost::shared_ptr<MembraneForcesBasic> p_shear_force(new MembraneForcesBasic());
        // p_shear_force->SetupMembraneConfiguration(cell_population);
        // simulator.AddForce(p_shear_force);

        // /*
        // -----------------------------
        // Bending forces
        // ----------------------------
        // */

        // boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
        // p_membrane_force->SetupInitialMembrane(p_mesh, simulator.rGetCellPopulation());
        // p_membrane_force->SetMembraneStiffness(BendingConst,30,30 );
        // simulator.AddForce(p_membrane_force);

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

        boundary_plane_normals.push_back(Create_c_vector(0.796512814285952,0.5625719395969793, 0.22154040140878425 ));
        boundary_plane_points.push_back(Create_c_vector(0.3709572080903858, 0.6640869252057044, 0.004854026486994977 ));
        
        boundary_plane_points.push_back(Create_c_vector( 0.5256067377437604,0.4620514197517712, -0.005173984009126431));
        boundary_plane_normals.push_back(Create_c_vector( 0.8541629230600764,0.5051514253254935, 0.12340072269273226));
    
        boundary_plane_points.push_back(Create_c_vector( 0.5068979181373336,0.25351647549806283,  0.021576436539306097));
        boundary_plane_normals.push_back(Create_c_vector( 0.4890372370519744, -0.868615353372259, -0.07968656513212781));

        boundary_plane_points.push_back(Create_c_vector(0.3264510482232788, 0.18589376368842103, 0.0007284710196420146 ));
        boundary_plane_normals.push_back(Create_c_vector(0.06079840539314355, -0.9970659090712014,  -0.04650942775039967));

        boundary_plane_points.push_back(Create_c_vector(0.18046257819123496, 0.2475290829034053, 0.006742930992320383  ));
        boundary_plane_normals.push_back(Create_c_vector( -0.7849488907857648, -0.6150729078345657, -0.07443491722390383 ));

        boundary_plane_points.push_back(Create_c_vector(0.11212066476933438, 0.4706722959297537, 0.006937365263584527  ));
        boundary_plane_normals.push_back(Create_c_vector(-0.991200354397287,0.11657903242117686  , 0.0626991757715747 ));

        boundary_plane_points.push_back(Create_c_vector( 0.2990633793409403, 0.6604501670537247, 2.799541929783861e-5  ));
        boundary_plane_normals.push_back(Create_c_vector( 0.713391205627203, -0.6789810757156608, 0.17337152751755425 ));
    

        for(unsigned boundary_id = 0; boundary_id < boundary_plane_normals.size(); boundary_id++)
         {

           boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population,boundary_plane_points[boundary_id],boundary_plane_normals[boundary_id],0.1));
            simulator.AddCellPopulationBoundaryCondition(p_condition);
        }


        simulator.Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);

    }

    void offTestCreatingBendingSimulation() throw(Exception)
    {
        std::string load_dir = "DevelopingHetroCylinder/SetChange/SetUpArchiving/"; 
        std::string output_directory = "DevelopingHetroCylinder/SetChange/HalfCollapse2/";
        double start_time = 40;
        OffLatticeSimulation<2,3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2,3>, 3 >::Load(load_dir , start_time);

        
        p_simulator->SetEndTime(start_time +2);
        p_simulator->SetOutputDirectory(output_directory);
        p_simulator->SetDt(0.001); // M_TIME_STEP
        p_simulator->SetSamplingTimestepMultiple(30);

        double BendingConst =pow(10, -13);
        std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
            //         KA,          Kalpha           Ks                                 Kb
        GrowthMaps[5] = Create_c_vector(pow(10, -6.9341), pow(10, -7.7), pow(10, -8), BendingConst);  GrowthMaps[4] = Create_c_vector(pow(10, -6.9), pow(10, -7.4224), pow(10, 8), BendingConst);  GrowthMaps[2] = Create_c_vector(pow(10, -6.8), pow(10, -6.8124), pow(10, -7), BendingConst);
        GrowthMaps[1.5] = Create_c_vector(pow(10, -6.5), pow(10, -6.3491), pow(10, -7), BendingConst);  GrowthMaps[1.2] = Create_c_vector(pow(10, -6.2), pow(10, -5.8360), pow(10, -7), BendingConst);

        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
	    // assert(boost::dynamic_pointer_cast<MembranePropertiesSecModifier<2,3> >(*iter));
	    boost::shared_ptr<MembranePropertiesSecModifier<2,3> > p_force_modifier = boost::static_pointer_cast<MembranePropertiesSecModifier<2, 3> >(*iter);

		p_force_modifier->SetMembranePropeties(GrowthMaps, 5, 1, 1e-8,1);
        p_force_modifier->SetBendingForce(p_simulator->rGetCellPopulation() , BendingConst);


        p_simulator->Solve();

        TRACE("Before Save");
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        TRACE("After Save");

        SimulationTime::Instance()->Destroy();
   
    }


    void offTestCreatingBendingSimulation3() throw(Exception)
    {
        std::string load_dir = "DevelopingHetroCylinder/SetChange/HalfCollapse2/";
        std::string output_directory = "DevelopingHetroCylinder/SetChange/HalfCollapse2/NoBending/";
        double start_time = 42;
        OffLatticeSimulation<2,3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2,3>, 3 >::Load(load_dir , start_time);

        
        p_simulator->SetEndTime(start_time +10);
        p_simulator->SetOutputDirectory(output_directory);
        p_simulator->SetDt(0.001); // M_TIME_STEP
        p_simulator->SetSamplingTimestepMultiple(30);

        double BendingConst =0;//pow(10, -13);
       
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
	    boost::shared_ptr<MembranePropertiesSecModifier<2,3> > p_force_modifier = boost::static_pointer_cast<MembranePropertiesSecModifier<2, 3> >(*iter);
        p_force_modifier->SetBendingForce(p_simulator->rGetCellPopulation() , BendingConst);
        p_simulator->Solve();

        TRACE("Before Save");
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        TRACE("After Save");

        SimulationTime::Instance()->Destroy();
   
    }

 
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/