#ifndef TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_
#define TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_
//TestCylindricalGrowthDeformableMembrane
#include <cxxtest/TestSuite.h>


#include <cstdio>
#include <ctime>
#include <cmath>
#include "OffLatticeSimulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
//#include "GeneralisedLinearSpringForce_CorrectedForArea.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "SmartPointers.hpp"
#include "CellLabel.hpp"
#include "VtkMeshWriter.hpp"
#include "Debug.hpp"
#include "CommandLineArguments.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "CellIdWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "MembraneSurfaceForceTest.hpp"
 //#include "MembraneSurfaceForce.hpp"
// #include "MembraneStiffnessForce.hpp"
// #include "MembraneShearForce.hpp"
#include "NewForce.hpp"

//#include "TriangleVertexMesh.hpp"



#include "Debug.hpp"
#include "PetscSetupAndFinalize.hpp"




#include "VertexMeshReader.hpp"



class RadialForce : public AbstractForce<2,3>
{
private:

    double mStrength;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<2,3> >(*this);
        archive & mStrength;
    }

public:
    RadialForce(double strength=1.0)
        : AbstractForce<2,3>(),
          mStrength(strength)
    {
        assert(mStrength > 0.0);
    }

    void AddForceContribution(AbstractCellPopulation<2,3>& rCellPopulation)
    {
        // Helper variables
        MeshBasedCellPopulation<2,3>* p_cell_population = static_cast<MeshBasedCellPopulation<2,3>*>(&rCellPopulation);


        // Calculate midpoint
        c_vector<double,3> centroid = zero_vector<double>(3);
        for (AbstractCellPopulation<2,3>::Iterator cell_iter = rCellPopulation.Begin();
                     cell_iter != rCellPopulation.End();
                     ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            centroid += rCellPopulation.GetNode(node_index)->rGetLocation();
        }
        centroid /= rCellPopulation.GetNumRealCells();

        // Get the vector from the centroid to the node and add a force in this direction 

        for (AbstractCellPopulation<2,3>::Iterator cell_iter = rCellPopulation.Begin();
                     cell_iter != rCellPopulation.End();
                     ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<3>* p_node = rCellPopulation.GetNode(node_index);
            c_vector<double,3> cell_location = p_node->rGetLocation() - centroid;
            c_vector<double, 3> force = zero_vector<double>(3);

            //Calculate cell normal (average of element normals)
            // c_vector<double,3> normal = zero_vector<double>(3);

            // std::set<unsigned>&  containing_elements = p_node->rGetContainingElementIndices();
            // assert(containing_elements.size()>0);
            // for (std::set<unsigned>::iterator iter = containing_elements.begin();
            //         iter != containing_elements.end();
            //         ++iter)
            // {
            //     // Negative as normals point inwards for these surface meshes
            //     normal += - p_cell_population->rGetMesh().GetElement(*iter)->CalculateNormal();
            // }
            // normal /= norm_2(normal);
            //double RadialLocation  = 1 / norm_2(cell_location);
            // cell_location[0]=1/cell_location[0];
            // cell_location[1]=1/cell_location[1];
            force =  mStrength * (cell_location) ; //mStrength * cell_area * normal;

          //  cell_iter->GetCellData()->SetItem("area", cell_area);
            
            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(force);

            cell_iter->GetCellData()->SetItem("applied_force_x", force[0]);
            cell_iter->GetCellData()->SetItem("applied_force_y", force[1]);
            cell_iter->GetCellData()->SetItem("applied_force_z", force[2]);
            cell_iter->GetCellData()->SetItem("applied_force_mag", norm_2(force));

            // cell_iter->GetCellData()->SetItem("norm_x", normal[0]);
            // cell_iter->GetCellData()->SetItem("norm_y", normal[1]);
            // cell_iter->GetCellData()->SetItem("norm_z", normal[2]);

            cell_iter->GetCellData()->SetItem("radius", norm_2(cell_location));
        }
    }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<Strength>" << mStrength << "</Strength>\n";
        AbstractForce<2,3>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(RadialForce)

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(RadialForce)











class CylinderValidation_CorrectedForDrag : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:

    void TestMembraneForce() throw (Exception)
    {
        double Area_constant = 25;
        unsigned N_D= 20;//{20,40,80};//,160};
        unsigned N_Z = 30;//{15,30,60,120,240};//,,480};
        // double Membrane_Constant_search[4] = {25, 30, 35, 40};
         double AreaDilationModulus_Sweep[4] = {25, 30, 35, 40};
         double ElasticShearModulus_Sweep[4] = {25, 30, 35, 40};

        double AreaDilationModulus = AreaDilationModulus_Sweep[1];
        double ElasticShearModulus = ElasticShearModulus_Sweep[1];
        

         TrianglesMeshReader<2,3> mesh_reader("projects/Jess/test/MembraneForceModel/mesh");
        MutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

            //    Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, 1.5e-3, 12e-3);
            //    MutableMesh<2,3>* p_mesh = generator.GetMesh();
            //    p_mesh->Translate(-6.0e-3*unit_vector<double>(3,2));

                //std::stringstream out;
                //out << membrane_constant;
                //std::string mesh_size = out.str();
               // std::string output_directory = "TestingForce/AreaConst30/AreaConst_" + mesh_size;
            std::string output_directory = "TwoTriangleMesh";



                // Create cells
                MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
                std::vector<CellPtr> cells;
                CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
                //cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);
                cells_generator.GenerateBasicRandom(cells,4, p_differentiated_type);

                // Create a cell population
                MeshBasedCellPopulation<2,3> cell_population(mesh, cells);
                cell_population.SetWriteVtkAsPoints(true);
                cell_population.SetOutputMeshInVtk(true);

                cell_population.AddCellWriter<CellIdWriter>();
                cell_population.AddCellWriter<CellProliferativeTypesWriter>();
                cell_population.AddCellWriter<CellMutationStatesWriter>();

                cell_population.CalculateRestLengths();

                // Set up cell-based simulation
                OffLatticeSimulation<2,3> simulator(cell_population);
                simulator.SetOutputDirectory(output_directory);
                simulator.SetEndTime(30); //50
                simulator.SetDt(0.002);
                simulator.SetSamplingTimestepMultiple(20);
                simulator.SetUpdateCellPopulationRule(false); // No remeshing.





/*
            -----------------------------
            Shearing Force 
            ----------------------------
*/
        //  double AreaDilationModulus = 0.001;
        //  double ElasticShearModulus = 0.001;
        //  boost::shared_ptr<MembraneShearForce> p_shear_force(new MembraneShearForce());
        //  p_shear_force->SetupMembraneConfiguration(cell_population);
        //  p_shear_force->SetAreaDilationModulus(AreaDilationModulus);
        //  p_shear_force->SetElasticShearModulus(ElasticShearModulus);
        //  simulator.AddForce(p_shear_force);




       
/*
            -----------------------------
            Bending Force
            ----------------------------
*/           
   

        // // double membrane_constant = 1e-12;
        // boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
        // p_membrane_force->SetupInitialMembrane(*p_mesh);
        // p_membrane_force->SetMembraneStiffness(membrane_constant);
        // simulator.AddForce(p_membrane_force);
        


 




/*
        
            -----------------------------
            Radial pressure force 
            ----------------------------
*/      
             double Pressure = 0.2; //1.0666e4;//to match 80mmh
        //  boost::shared_ptr<NewForce> p_pressure_force(new NewForce());
        //  p_pressure_force->GetPressure(Pressure);
        //  simulator.AddForce(p_pressure_force);

                MAKE_PTR_ARGS(RadialForce, p_radial_force, (Pressure));
                simulator.AddForce(p_radial_force);
 
/*
            -----------------------------
            Surface Area Force
            ----------------------------
*/           
          
        //double Area_constant = 50;
         boost::shared_ptr<MembraneSurfaceForceTest> p_surface_force(new MembraneSurfaceForceTest());
         p_surface_force->SetupInitialAreas(cell_population);
         p_surface_force->SetMembraneStiffness(Area_constant);
         simulator.AddForce(p_surface_force);




            //  MAKE_PTR(MembraneSurfaceForce, p_surface_force);
            // simulator.AddForce(p_surface_force);


 /*
            -----------------------------
            Spring force, corrected for drag
            ----------------------------
*/  

           // Create a force law and pass it to the simulation
           //boost::shared_ptr<GeneralisedLinearSpringForce_CorrectedForArea<2,3> > p_linearSpring_force(new GeneralisedLinearSpringForce_CorrectedForArea<2,3>());
           // p_linearSpring_force->SetMeinekeSpringStiffness(50.0);
          //Spring_force->SetMeinekeSpringStiffness(50.0);
          //simulator.AddForce(p_linearSpring_force);


               // Create a plane boundary to represent the inlet and pass them to the simulation
                boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition_1(new FixedRegionBoundaryCondition<2,3>(&cell_population,5.0e-3*unit_vector<double>(3,2),unit_vector<double>(3,2),10));
                simulator.AddCellPopulationBoundaryCondition(p_condition_1);

                boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition_2(new FixedRegionBoundaryCondition<2,3>(&cell_population,5.0e-3*unit_vector<double>(3,2),-unit_vector<double>(3,2),10));
                simulator.AddCellPopulationBoundaryCondition(p_condition_2);

                simulator.Solve();

                // To reset before looping: this is usually done by the SetUp and TearDown methods
                SimulationTime::Instance()->Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
    }
 

};



#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/





        // std::vector<Node<3>*> nodes;
        // nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        // nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        // nodes.push_back(new Node<3>(2, false, 0.5, 1.0, 0.0));
        // nodes.push_back(new Node<3>(3, false, 1.5, 1.0, 0.0));

        // // Make two triangular faces out of these nodes
        // std::vector<std::vector<Node<3>*> > nodes_faces(2);

        // unsigned node_indices_face_0[3] = {0, 1, 2};
        // unsigned node_indices_face_1[3] = {1, 2, 3};
        // for (unsigned i=0; i<3; i++)
        // {
        //         nodes_faces[0].push_back(nodes[node_indices_face_0[i]]);
        //         nodes_faces[1].push_back(nodes[node_indices_face_1[i]]);
        // }
        // // Make the faces
        // std::vector<Element<3,3>*> faces;
        // for (unsigned i=0; i<3; i++)
        // {
        //     faces.push_back(new Element<3,3>(i, nodes_faces[i]));
        // }
        // // Make the elements
        // std::vector<Element<3,3>*> faces_element;
        // std::vector<bool> orientations;

        // for (unsigned i=0; i<3; i++)
        // {
        //     faces_element.push_back(faces[i]);
        //     orientations.push_back(true);
        // }
        // // faces_element.push_back(faces[5]);
        // // orientations.push_back(false);

        // std::vector<Element<3,3>*> elements;
        // elements.push_back(new Element<2,3>(0, faces_element, orientations));
    
        // MutableMesh<2,3>* p_mesh = MutableMesh<2,3>MutableMesh(nodes, faces, elements);

