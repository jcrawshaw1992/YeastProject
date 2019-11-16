
#ifndef TestPottsModelOn2DSquare_HPP_
#define TestPottsModelOn2DSquare_HPP_

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <cstdio>
#include <ctime>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellsGenerator.hpp"

#include "OnLatticeSimulation.hpp"

#include "FixedG1GenerationalCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"

// Include the relevant cell writers
#include "CellAspectRatioWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellMajorAxisAngleWriter.hpp"

#include "PottsBasedCellPopulation.hpp"

#include "MeshBasedCellPopulation.hpp"

#include "ArbitraryAdhesionPottsUpdateRule.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
// New rule to make the cells wrap
#include "AppliedForceModifier.hpp"
#include "ArbitraryCurvatureConstraintPottsUpdateRule.hpp"

#include "PottsCellPropertiesModifier.hpp"

// #include "TractionDataLoader.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "PottsMeshGenerator.hpp"
// #include "WildTypeCellMutationState.hpp"
#include "AspectRatioConstraintPottsUpdateRule.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "LogFile.hpp"
#include "MechanotaxisPottsUpdateRule.hpp"
#include "PottsMeshFromMutableMeshGenerator.hpp"
#include "SmartPointers.hpp"
#include "Warnings.hpp"
// #include "VtkMeshReader.hpp"
#include "RandomNumberGenerator.hpp"
#include "AbstractForce.hpp"
#include "ArbitrarySurfaceAreaConstraintPottsUpdateRule.hpp"
#include "ArbitraryVolumeConstraintPottsUpdateRule.hpp"
#include "Debug.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "Honeycomb3DMeshGenerator.hpp"

#include "PottsArbitrarySurfaceIn3DMesh.hpp"
#include "PottsMeshFromMutableMeshGeneratorJess.hpp"
#include "UblasCustomFunctions.hpp"

// Needed for NodesOnlyMesh.
#include "PetscSetupAndFinalize.hpp"

#include <iostream>
#include <stdio.h>

// #include "MembraneForcesBasic.hpp"
// #include "MembranePropertiesModifier.hpp"
// #include "projects/VascularRemodelling/src/MembraneForces/MembraneForcesBasic.hpp"

// #include "EmptyBasementMatrix.hpp"
// #include "HasEndothelialCell.hpp"
// #include "LostEndothelialCell.hpp"
// #include "OutwardsPressure.hpp"


class RadialForce : public AbstractForce<2, 3>
{
private:
    double mStrength;

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractForce<2, 3> >(*this);
        archive& mStrength;
    }

public:
    RadialForce(double strength = 1.0)
            : AbstractForce<2, 3>(),
              mStrength(strength)
    {
        assert(mStrength > 0.0);
    }

    void AddForceContribution(AbstractCellPopulation<2, 3>& rCellPopulation)
    {
        // Helper variables
        MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);

        // Calculate midpoint
        c_vector<double, 3> centroid = zero_vector<double>(3);
        for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            centroid += rCellPopulation.GetNode(node_index)->rGetLocation();
        }
        centroid /= rCellPopulation.GetNumRealCells();

        for (AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<3>* p_node = rCellPopulation.GetNode(node_index);
            c_vector<double, 3> cell_location = p_node->rGetLocation() - centroid;
            cell_location(2) = 0.0;
            c_vector<double, 3> force = zero_vector<double>(3);

            //Calculate cell normal (average of element normals)
            c_vector<double, 3> normal = zero_vector<double>(3);

            std::set<unsigned>& containing_elements = p_node->rGetContainingElementIndices();
            assert(containing_elements.size() > 0);
            for (std::set<unsigned>::iterator iter = containing_elements.begin();
                 iter != containing_elements.end();
                 ++iter)
            {

                normal += p_cell_population->rGetMesh().GetElement(*iter)->CalculateNormal();
            }
            normal /= norm_2(normal);

            double cell_area = rCellPopulation.GetVolumeOfCell(*cell_iter);
            // force = mStrength * cell_area * normal; // cell_location / norm_2(cell_location);

            force = -mStrength * normal; // cell_location / norm_2(cell_location);
            cell_iter->GetCellData()->SetItem("area", cell_area);

            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(force);

            cell_iter->GetCellData()->SetItem("applied_force_x", force[0]);
            cell_iter->GetCellData()->SetItem("applied_force_y", force[1]);
            cell_iter->GetCellData()->SetItem("applied_force_z", force[2]);
            cell_iter->GetCellData()->SetItem("applied_force_mag", norm_2(force));

            cell_iter->GetCellData()->SetItem("norm_x", normal[0]);
            cell_iter->GetCellData()->SetItem("norm_y", normal[1]);
            cell_iter->GetCellData()->SetItem("norm_z", normal[2]);
        }
    }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<Strength>" << mStrength << "</Strength>\n";
        AbstractForce<2, 3>::OutputForceParameters(rParamsFile);
    }
};

class TestPottsOnCylinder : public AbstractCellBasedTestSuite
{
public:
    double ComputeMeshTotalArea(TetrahedralMesh<2, 3>& mutable_mesh)
    {
        double total_area = 0.;
        for (TetrahedralMesh<2, 3>::ElementIterator element_iter = mutable_mesh.GetElementIteratorBegin();
             element_iter != mutable_mesh.GetElementIteratorEnd();
             ++element_iter)
        {
            const c_vector<double, 3> AC = element_iter->GetNodeLocation(1) - element_iter->GetNodeLocation(0);
            const c_vector<double, 3> AB = element_iter->GetNodeLocation(2) - element_iter->GetNodeLocation(0);
            total_area += 0.5 * norm_2(VectorProduct(AC, AB));
        }

        PRINT_VARIABLE(total_area);
        return total_area;
    }

    void TestPottsOnCylinderWithDeformation()
    {

        // Regualr square --hexagonal lattie
        // double N_D = 50;
        // double N_Z = 50;
        // Honeycomb3DMeshGenerator generator(N_D, N_Z, 1,1);
        // MutableMesh<2, 3>* mutable_mesh = generator.GetMesh();
        // mutable_mesh->Translate( unit_vector<double>(3, 2));

        // Random cylinder
        //   unsigned N = 581;
        //   std::stringstream out;
        //     out << N;
        //     std::string mesh_size = out.str();
        //     std::string mesh_file = "projects/EMBC2018/test/data/cyl_" + mesh_size + "_nodes.vtu";
        //     // std::string output_directory = "CylinderValidation/Random/" + mesh_size;

        //     // This data file is in mm
        //     VtkMeshReader<2,3> mesh_reader(mesh_file);
        //     MutableMesh<2,3> mutable_mesh;
        //     mutable_mesh.ConstructFromMeshReader(mesh_reader);
        //     double scaling = 1e-3;  // so distances are in m

        // Regular cylinder

        // // scale everything up

        // double scale = 1e1;
        // PRINT_VARIABLE(scale);

        // 			double NumberOfCells =4;
        // 			unsigned N_D = 10;
        // 			unsigned N_Z = 20;//N_D*1.5;
        // 			double Radius = 5e-4 *scale;
        // 			double Length = 30e-4 *scale;//2*96.6e-4 * scale; //12e-3;
        // 			double trans = -Length/2;

        // 			// MutableMesh<2, 3>* p_mesh = p_mesh_base;
        // 			Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
        // 			// Honeycomb3DCylinderMeshGenerator generator(N_D[N_D_index], N_Z, 1.5e-3, Length);
        // 			MutableMesh<2, 3>* mutable_mesh = generator.GetMesh();
        // 			mutable_mesh->Translate(trans * unit_vector<double>(3, 2));

        // Birfucation
        unsigned N = 581;
        std::stringstream out;
        out << N;
        std::string mesh_size = out.str();
        std::string mesh_file = "projects/VascularRemodelling/test/data/bifurcation_cut/Scalled/config.vtu";
        std::string output_directory = "BifucationPottsTest/";

        // This data file is in mm??
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> mutable_mesh;
        mutable_mesh.ConstructFromMeshReader(mesh_reader);
        double scale = 1e-3; // so distances are in m
        mutable_mesh.Scale(1.0 * scale, 1.0 * scale, 1.0 * scale); // so distances are back in original scal

        double Radius = 5e-4 * scale;
        double Length = 30e-4 * scale; //2*96.6e-4 * scale; //12e-3;

        /*
		 * Setup Potts simulation
		 */

        PottsMeshFromMutableMeshGeneratorJess<3> potts_generator(mutable_mesh);
        // PottsMeshFromMutableMeshGeneratorJess<3> potts_generator(* mutable_mesh);
        PottsArbitrarySurfaceIn3DMesh<3>* p_potts_mesh = potts_generator.GetMesh();

        unsigned max_num_cells = 10;
        unsigned num_cells_seeded = 8;
     
     
        double target_volume = ComputeMeshTotalArea(mutable_mesh) / max_num_cells;
        double MeshArea = ComputeMeshTotalArea(mutable_mesh);

        double CellScale =1e-3;
        double AthesticScalling = 5;
		double LongAxis = 10.0/2.0* 3.0*AthesticScalling *CellScale;//* 1e-6 * CellScale;//*scale ;
		double ShortAxis = 5.0/2.0*AthesticScalling*  CellScale;//*1e-6 * CellScale ;//*scale ;
		double CellArea = M_PI * LongAxis  * ShortAxis ;

		double AxisRatio = pow(LongAxis -ShortAxis,2)/ pow(LongAxis +ShortAxis,2);//  (LongAxis -ShortAxis)* (LongAxis -ShortAxis) /(LongAxis +ShortAxis)*(LongAxis +ShortAxis); 
		double CellPerimeter = M_PI * (LongAxis +ShortAxis)* (3 *AxisRatio * 1/(sqrt(-3*AxisRatio+4) +10) +1) ;   
        CellPerimeter/=1.0;//

        PRINT_3_VARIABLES(CellArea, MeshArea, CellPerimeter);

        for (unsigned i = 0; i < num_cells_seeded; i++)
        {
            // Pick a node randomly avoiding boundary nodes
            unsigned num_node;
            do
            {
                num_node = RandomNumberGenerator::Instance()->randMod(p_potts_mesh->GetNumNodes()); // Selects a random node

            } while (p_potts_mesh->GetNode(num_node)->IsBoundaryNode()); // XXX not sure what this does, Im not sure it does anything
            PRINT_VARIABLE(i);

            std::vector<Node<3>*> element_nodes;
            element_nodes.push_back(p_potts_mesh->GetNode(num_node));
            p_potts_mesh->AddElement(new PottsElement<3>(i, element_nodes));
            std::cout << "Seeded at " << num_node << std::endl;
        }

        // Randomly place the cell markers for seeding
        std::vector<CellPtr> potts_cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type); // Set the cell type for the seeds
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> potts_cells_generator; // Set the cell cylce for the seeds
        potts_cells_generator.GenerateBasicRandom(potts_cells, p_potts_mesh->GetNumElements(), p_diff_type); // Assigning the cell type and the cell cycle

        TRACE("Generate Potts cell population");
        // Create cell population linking potts mesh and cells
        PottsBasedCellPopulation<3> potts_population(*p_potts_mesh, potts_cells);

        // Add in all the cell writers
        potts_population.AddCellWriter<CellIdWriter>();
        potts_population.SetNumSweepsPerTimestep(1);
        double Confluence_time =30;
        // Set up Potts simulation
        OnLatticeSimulation<3> potts_simulator(potts_population);
        potts_simulator.SetDt(0.1);
        potts_simulator.SetOutputDirectory("TestPottsOn2dSquare_Potts_Setup");
        potts_simulator.SetSamplingTimestepMultiple(10);
        potts_simulator.SetEndTime(Confluence_time);

        // // Reduce temperature to avoid randomly removing cells before they grow (parameters will need proper fitting)
        potts_population.SetTemperature(1e-9);

        /*
         * Simulate cell seeding from randomly placed cell markers
         */

        MAKE_PTR(ArbitraryVolumeConstraintPottsUpdateRule<3>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(CellArea);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(4e7);
        potts_simulator.AddUpdateRule(p_volume_constraint_update_rule);

        MAKE_PTR(ArbitraryAdhesionPottsUpdateRule<3>, p_adhesion_update_rule);
        p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02 / 1000);
        p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16 / 1000);
        potts_simulator.AddUpdateRule(p_adhesion_update_rule);

        
        MAKE_PTR(ArbitrarySurfaceAreaConstraintPottsUpdateRule<3>, p_area_constraint_update_rule);
        p_area_constraint_update_rule->SetSurfaceAreaEnergyParameter(0.016e2);
        p_area_constraint_update_rule->SetTargetSurfaceArea(CellPerimeter);
        potts_simulator.AddUpdateRule(p_area_constraint_update_rule);

  
        p_potts_mesh->UpdatePottsNodeLocationFromDelaunay();
        p_potts_mesh->GetRadius(Radius);
        potts_simulator.Solve();

      
        potts_simulator.SetStartTime(Confluence_time)

         /*
          Add mechanotaxis to Potts simulation for timesteps to follow
        */

        std::string traction_file = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/Tester/results/Extracted/surface-tractions.xtr";
        p_potts_mesh->TractionDataLoader(traction_file);
        MAKE_PTR(MechanotaxisPottsUpdateRule<3>, p_mechanotaxis_update_rule);
        p_mechanotaxis_update_rule->SetTractionCorrelationParameter(0.016e25); //1e-8);
        potts_simulator.AddUpdateRule(p_mechanotaxis_update_rule);

        /*
        -----------------------------
        Mesh SetUp
        ----------------------------
        */
        // std::vector<CellPtr> cells;
        // CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        // cells_generator.GenerateBasicRandom(cells, mutable_mesh.GetNumNodes(), p_diff_type);

        // // Create a cell population. Links the Mesh and the cells
        // MeshBasedCellPopulation<2, 3> cell_population(mutable_mesh, cells);

        // cell_population.SetWriteVtkAsPoints(true);
        // cell_population.SetOutputMeshInVtk(true);

        // // Set up cell-based simulation
        // SimulationTime::Instance()->Destroy();
        // SimulationTime::Instance()->SetStartTime(0);
        // OffLatticeSimulation<2, 3> simulator(cell_population);
        // simulator.SetOutputDirectory("TestPottsOnVesselGrowth_Growth");
        // simulator.SetDt(0.002);
        // simulator.SetSamplingTimestepMultiple(50);
        // simulator.SetEndTime(1);
        // simulator.SetUpdateCellPopulationRule(false); // No remeshing.

        // // Add a force law for normal stress, i.e. blood pressure only
        // double Tissue_Pressure = 0.001466542; //1.0666e4; // to match 80mmhg
        // MAKE_PTR_ARGS(RadialForce, p_radial_force, (Tissue_Pressure));
        // simulator.AddForce(p_radial_force);

        /*
        -----------------------------
        Shear and Area Force
        ----------------------------
        */

        // std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
        //     //         KA,          Kalpha           Ks
        // GrowthMaps[1.5] = Create_c_vector(pow(10, -6.5), pow(10, -6.3491), pow(10, -7), 0);
        // GrowthMaps[1.2] = Create_c_vector(pow(10, -6.2), pow(10, -5.8360), pow(10, -7), 0);

        // double Strength = 1.2;
        // bool Hetro = 0;

        // /*
        // -----------------------------
        // SMembrane forces
        // ----------------------------
        // */

        // boost::shared_ptr<MembranePropertiesModifier<2, 3> > p_Membrane_modifier(new MembranePropertiesModifier<2, 3>());
        // p_Membrane_modifier->SetMembranePropeties(GrowthMaps, Strength, Hetro);
        // TRACE("Did the SetMembranePropeties")
        // p_Membrane_modifier->SetupSolve(cell_population, output_directory);
        // TRACE("Did the SetupSolve")

        // // boost::shared_ptr<MembraneForcesBasic> p_Membrane_forces(new MembraneForcesBasic());
        // // p_shear_force->SetupMembraneConfiguration(cell_population);
        // // simulator.AddForce(p_shear_force);

        // // /*
        // // -----------------------------
        // // Boundary conditions
        // // ----------------------------
        // // */

        // // // //Create a plane boundary to represent the inlet and pass them to the simulation
        // // c_vector<long double, 3> Boundary1 = Create_c_vector(12.48824184355049, 34.75302061558864, 41.78949821113195);
        // // c_vector<long double, 3> Normal1 = Create_c_vector(-0.05607774749413225, 0.762765339938692, 0.6442393362906335);
        // // double Radius1 = 1.7;

        // // c_vector<long double, 3> Boundary2 = Create_c_vector(12.597373380655702, 48.382440094438316, 42.984851419357064);
        // // c_vector<long double, 3> Normal2 = Create_c_vector(-0.04847413662454751, -0.989768366942236, -0.13419701143842153);
        // // double Radius2 = 1.2;

        // // c_vector<long double, 3> Boundary3 = Create_c_vector(-0.04710924565104161, -0.9898343148598956, -0.13419667693363677);
        // // c_vector<long double, 3> Normal3 = Create_c_vector(-0.04710924565104161, -0.9898343148598956, -0.13419667693363677);
        // // double Radius3 = 1.3;

        // // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary1, Normal1, Radius1));
        // // simulator.AddCellPopulationBoundaryCondition(p_condition_1);

        // // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary2, Normal2, Radius2));
        // // simulator.AddCellPopulationBoundaryCondition(p_condition_2);
        // // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_3(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary3, Normal3, Radius3));
        // // simulator.AddCellPopulationBoundaryCondition(p_condition_3);

        // /*
		//  * We need timestepping for two different simulations
		//  */
        // double potts_timestep = 0.1;
        // double growth_timestep = 0.1;
        // double potts_current_time = 0;
        // double growth_current_time = 0;

        // /*
        // -----------------------------
        // Inital Seeding
        // ----------------------------
        // */


        // double Confluence_time = 10;
        // potts_simulator.SetOutputDirectory("TestPottsOnVesselGrowth_Potts_Setup");
        // potts_simulator.Solve();





        // /// This is the iterations part 
        // potts_simulator.SetOutputDirectory("TestPottsOnVesselGrowth_Potts");
        // potts_simulator.SetSamplingTimestepMultiple(1);
        // unsigned num_iterations = 60;
        // // for (unsigned iter_num = 0; iter_num < num_iterations; iter_num++)
        // // {
        // //     std::cout << "Iteration  " << iter_num << std::endl;
        // //     /*
		// //      * Do one timestep of cell movement
		// //      */

        // //     std::cout << "Potts timestep" << std::endl;

        // //     // Destroy the simulation time that has been left from the setup step or the previous iteration
        //     SimulationTime::Instance()->Destroy();
        //     SimulationTime::Instance()->SetStartTime(potts_current_time);
        //     if (potts_current_time > 0)
        //     {
        //         double dt = simulator.GetDt();
        //         unsigned num_time_steps = (unsigned)(potts_timestep / dt + 0.5);
        //         SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(potts_current_time + potts_timestep, num_time_steps);
        //     }

        //     potts_simulator.SetEndTime(potts_current_time + potts_timestep);
        //     potts_simulator.Solve();
        //     potts_current_time += potts_timestep;

        //     /*
        //      * Do one timestep of vessel mechanics
        //      */
        //     std::cout << "Mechanics timestep" << std::endl;

        //     SimulationTime::Instance()->Destroy();
        //     SimulationTime::Instance()->SetStartTime(growth_current_time);
        //     if (growth_current_time > 0)
        //     {
        //         double dt = simulator.GetDt();
        //         unsigned num_time_steps = (unsigned)(growth_timestep / dt + 0.5);
        //         SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(growth_current_time + growth_timestep, num_time_steps);
        //     }

        //     simulator.SetEndTime(growth_current_time + growth_timestep);

        //     // there is a problem here with the solve???
        //     simulator.Solve();

        //     // // // Update location of nodes in potts mesh from recently computed deformation
        //     p_potts_mesh->UpdatePottsNodeLocationFromDelaunay();

        //     growth_current_time += growth_timestep;
        // }

        // Now need a step where the mesh deforms and we go again
    }
};

#endif /*TESTREPRESENTATIVEPOTTSBASEDONLATTICESIMULATION_HPP_*/







		// * ADDITIONAL UPDATE RULES 
		// */

        // MAKE_PTR(PottsCellPropertiesModifier<3>, p_modifier);
        // potts_simulator.AddSimulationModifier(p_modifier);

        // MAKE_PTR(ArbitraryCurvatureConstraintPottsUpdateRule<3>, p_Curvature_constraint_update_rule);
        // p_Curvature_constraint_update_rule->SetCurvatureEnergyParameter(0.016e18); //(0.016e18); // lambda_p
        // p_Curvature_constraint_update_rule->SetTargetCurvature(20993);// 7993) // p_area_constraint_update_rule->SetTargetSurfaceArea (150e-6); // (150e-6); // P (is 150e-12 in paper but value was found wrong and corrected for 150e-6)
        // potts_simulator.AddUpdateRule(p_Curvature_constraint_update_rule);

        // // Need to have a look at what this is doing, to asserts need to be DIM ==2
        //         MAKE_PTR(AspectRatioConstraintPottsUpdateRule<3>, p_aspect_ratio_update_rule);
        //         p_aspect_ratio_update_rule->SetTargetAspectRatio(4);
        //         p_aspect_ratio_update_rule->SetAspectRatioEnergyParameter(1e-7);
        //         potts_simulator.AddUpdateRule(p_aspect_ratio_update_rule);



        // MAKE_PTR(MechanotaxisPottsUpdateRule<3>, p_mechanotaxis_update_rule);
        // p_mechanotaxis_update_rule->SetTractionCorrelationParameter(0.01e20);//1e-8);
        // potts_simulator.AddUpdateRule(p_mechanotaxis_update_rule);







// Checked that i have a proper mesh

// std::vector<CellPtr> cells;
// CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
// cells_generator.GenerateBasicRandom(cells, mutable_mesh->GetNumNodes(), p_diff_type);

// // Create a cell population. Links the Mesh and the cells
// MeshBasedCellPopulation<2,3> cell_population(* mutable_mesh, cells);

// cell_population.SetWriteVtkAsPoints(true);
// cell_population.SetOutputMeshInVtk(true);
// // cell_population.SetOutputCellProliferativeTypes(true);
// cell_population.CalculateRestLengths();
// // cell_population.SetOutputCellIdData(true);
// //cell_population.SetDampingConstantNormal(0.1);//e5);
// //cell_population.SetDampingConstantMutant(0.1);//e5);

// // Set up cell-based simulation
// SimulationTime::Instance()->Destroy();
// SimulationTime::Instance()->SetStartTime(0);
// TRACE("NewSimulation");
// OffLatticeSimulation<2,3> simulator(cell_population);
// simulator.SetOutputDirectory("TestPottsOnVesselGrowth_Growth");
// simulator.SetDt(0.002);
// simulator.SetSamplingTimestepMultiple(500);
// simulator.SetEndTime(13);
// simulator.SetUpdateCellPopulationRule(false); // No remeshing.
// simulator.Solve();

// /*
//  * Setup growth problem
//  */
// // Create representation of cells required for the mesh based simulation
// TRACE("NewCells");
// std::vector<CellPtr> vessel_cells;
// CellsGenerator<FixedG1GenerationalCellCycleModel, 3> vessel_cells_generator;
// vessel_cells_generator.GenerateBasicRandom(vessel_cells, mutable_mesh->GetNumNodes(), p_diff_type);

// // Create a cell population linking mechanics mesh and cells
// MeshBasedCellPopulation<2,3> vessel_cell_population(*mutable_mesh, vessel_cells);
// vessel_cell_population.SetWriteVtkAsPoints(true);
// vessel_cell_population.SetOutputMeshInVtk(true);

//  // Set up cell-based simulation
// TRACE("NewSimulation");
// OffLatticeSimulation<2,3> vessel_simulator(vessel_cell_population);
// vessel_simulator.SetOutputDirectory("TestPottsOnVesselGrowth_Growth");
// vessel_simulator.SetDt(0.002);
// vessel_simulator.SetSamplingTimestepMultiple(500);
// vessel_simulator.SetUpdateCellPopulationRule(false); // No remeshing.

// This loads tractions into the lattice site objects so it can be used for mechanotaxis.
// @todo the force law for normal stress (a few lines above) should probably be reimplemented reading from file too.
// TractionDataLoader potts_traction_loader("projects/VascularRemodelling/test/data/straight_vessel-surface-tractions.xtr");
// potts_traction_loader.UpdateLatticeSiteData(p_potts_mesh);

/*       
        -----------------------------
        Tractionforce 
        ----------------------------
        */

// // std::string working_directory = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/DeformingCylindricalBasicHetroWallTesting/";
// std::string working_directory = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/FluidFlowInPipe/";
//  std::string output_directory = "AppropriateInitalConditions";//Shrinking"; // + Parameters + "/";

//  // Create an Applied Force modifier to couple to Flow
// std::string traction_file = working_directory + "results/Extracted/surface-tractions.xtr";
// boost::shared_ptr<AppliedForceModifier<2, 3> > p_force_modifier(new AppliedForceModifier<2, 3>());

// p_force_modifier->SetResetTractionsOnCells(true, traction_file);
// p_force_modifier->SetupVessel(vessel_cell_population, output_directory);
// vessel_simulator.AddSimulationModifier(p_force_modifier);
