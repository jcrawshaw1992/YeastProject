#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

//  You might like to fix up the boundaries here

#include <cxxtest/TestSuite.h>
#include <cstdio>
#include <ctime>
#include <cmath>

#include "Debug.hpp"

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "HoneycombMeshGenerator.hpp"

#include "CylindricalHoneycombMeshGenerator.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "AppliedForceModifier.hpp"
#include "AppliedForce.hpp"
#include "SpringLengthModifier.hpp"

#include "CellIdWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"

#include "CellsGenerator.hpp"
// #include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
// #include "LinearSpringForce.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "MeshBasedCellPopulation.hpp"

#include "FixedRegionBoundaryCondition.hpp"
#include "SmartPointers.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VtkMeshWriter.hpp"
#include "CommandLineArguments.hpp"
// #include "MembraneStiffnessForce.hpp"
#include "PetscSetupAndFinalize.hpp"

// #include "PressureForce.hpp"
#include "XmlTools.hpp"

#include "UblasCustomFunctions.hpp"


#include "VtkMeshReader.hpp"

#include "MembraneHetroModifier.hpp"

// #include "ConstantPressure.hpp"
#include "EmptyBasementMatrix.hpp"
#include "HasEndothelialCell.hpp"
#include "LostEndothelialCell.hpp"
#include "CellMutationStatesWriter.hpp"
#include "OutwardsPressure.hpp"

#include "MembranePropertiesSecModifier.hpp"
#include "MembraneForcesBasic.hpp"

#include "projects/VascularRemodelling/src/MembraneForces/MembraneStiffnessForce.hpp"

#include "BoundariesModifier.hpp"

#include "LinearSpringForce.hpp"



using namespace xsd::cxx::tree;

class TestRemeshing  : public AbstractCellBasedTestSuite
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

    void ReadHemeLBv3Vector(c_vector<double,3>& vector, xercesc::DOMElement* xmlElement)
    {
        std::stringstream raw_value(X2C(xmlElement->getAttribute(X("value"))));
        char left_par, comma1, comma2, right_par;
        double x, y, z;
        raw_value >> left_par; assert(left_par == '(');
        raw_value >> x;
        raw_value >> comma1; assert(comma1 == ',');
        raw_value >> y;
        raw_value >> comma2; assert(comma2 == ',');
        raw_value >> z;
        raw_value >> right_par; assert(right_par == ')');

        vector[0] = x;
        vector[1] = y;
        vector[2] = z;
      
    }

    void ReadIoletPlanes(std::string hemelbConfigFile,
                         std::vector<c_vector<double,3> >& rBoundaryPlanePoints,
                         std::vector<c_vector<double,3> >& rBoundaryPlaneNormals,
                         std::vector<double>& rBoundaryPlaneRadii)
    {
        
        bool hemelb_config_v3;
        std::vector<xercesc::DOMElement*> iolets;
        try
        {
        	properties<char> props;
            xsd::cxx::xml::dom::auto_ptr<xercesc::DOMDocument> p_doc = XmlTools::XmlTools::ReadXmlFile(hemelbConfigFile, props, false);
            xercesc::DOMElement* p_root_elt = p_doc->getDocumentElement();
            hemelb_config_v3 = (atoi(X2C(p_root_elt->getAttribute(X("version"))).c_str()) == 3);
            iolets = XmlTools::FindElements(p_root_elt, "inlets/inlet");
            std::vector<xercesc::DOMElement*> aux = XmlTools::FindElements(p_root_elt, "outlets/outlet");
            iolets.insert(iolets.end(), aux.begin(), aux.end());
        }
        catch (const exception<char>& e)
        {
            std::cerr << e << std::endl;
            EXCEPTION("XML parsing error in configuration file: " + hemelbConfigFile);
        }

        if (!hemelb_config_v3)
        {
            EXCEPTION("Only HemeLB XML version 3 is supported");
        }

PRINT_VARIABLE(iolets.size());

        for (unsigned iolet_id=0; iolet_id<iolets.size(); iolet_id++)
        {
    
    
            c_vector<double,3> vector;
          
          

            std::vector<xercesc::DOMElement*> points = XmlTools::FindElements(iolets[iolet_id], "position");
          
          
            assert(points.size()==1);
            
            // PRINT_VECTOR(vector);
			ReadHemeLBv3Vector(vector, points[0]);
            rBoundaryPlanePoints.push_back(vector);
    

            std::vector<xercesc::DOMElement*> normals = XmlTools::FindElements(iolets[iolet_id], "normal");
            assert(normals.size()==1);
  
  
			ReadHemeLBv3Vector(vector, normals[0]);
            rBoundaryPlaneNormals.push_back(vector);
    
    

            std::vector<xercesc::DOMElement*> radi = XmlTools::FindElements(iolets[iolet_id], "radius");
            if(radi.size()==1)
            {
          
            	double radius = atof(X2C(radi[0]->getAttribute(X("radius"))).c_str());
            	rBoundaryPlaneRadii.push_back(radius);
        
            }
            else
            {
            	rBoundaryPlaneRadii.push_back(10); // in m
            }
        }
    }

public:

    void TestSetupRunAndSavePipe() throw (Exception)
    {
        
        std::string working_directory = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/ShrunkPlexusWithLongInlets/" ;
        std::string mesh_file = working_directory +"SetUpData/Plexus_MostCourse.vtu" ;//"SetUpData/config.vtu";
		double mesh_scale = 1e-3; 
       
           // Read in the plexus 


        VtkMeshReader<2,3> mesh_reader(mesh_file);
        MutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Scale(mesh_scale,mesh_scale,mesh_scale); // so distances are in mm

        // Honeycomb3DCylinderMeshGenerator generator(30, 40, 0.001, 0.01);
        // MutableMesh<2, 3>* mesh = generator.GetMesh();
        
        // Create the cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        MeshBasedCellPopulation<2,3> cell_population(mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        cell_population.CalculateRestLengths();

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();


        // Set up cell-based simulation
        OffLatticeSimulation<2,3> simulator(cell_population);
        std::string output_dir = "RemeshingComparison/WithoutRemeshing/";
        
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetDt(0.00002e-1);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(500);



        /*
        -----------------------------
       Compressive tissue pressure
        ----------------------------
        */

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.1000e-05 mmHg



        // Create a plane boundary to represent the inlets/outlets and pass them to the simulation


        std::vector<c_vector<double,3> > boundary_plane_points;
        std::vector<c_vector<double,3> > boundary_plane_normals;

        boundary_plane_points.push_back(Create_c_vector(0.19,0.73,0.001));
        boundary_plane_normals.push_back(Create_c_vector(1,0,0));


        boundary_plane_points.push_back(Create_c_vector(0.495,0.35,0.0088));
        boundary_plane_normals.push_back(Create_c_vector(0.709,0.7,-0.052));


        boundary_plane_points.push_back(Create_c_vector(0.755,0.185,0.0007));
        boundary_plane_normals.push_back(Create_c_vector(0.21,0.97,0.022));


// ------------
        boundary_plane_points.push_back(Create_c_vector(1.243, 0.258,-0.007));
        boundary_plane_normals.push_back(Create_c_vector(-0.75,0.665,-0.041));


        boundary_plane_points.push_back(Create_c_vector(1.25,0.868,-0.0106));
        boundary_plane_normals.push_back(Create_c_vector(-0.735,-0.67,-0.021));


        boundary_plane_points.push_back(Create_c_vector(1.06,1.06,0.004));
        boundary_plane_normals.push_back(Create_c_vector(-0.8605,-0.488, 0.144));

        boundary_plane_points.push_back(Create_c_vector(0.6852,1.2,-0.0065));
        boundary_plane_normals.push_back(Create_c_vector(0.685,-0.72,0.103));


        boost::shared_ptr<BoundariesModifier<2, 3> > p_Boundary_modifier(new BoundariesModifier<2, 3>());
        p_Boundary_modifier->CreateBoundaryNodes(cell_population,boundary_plane_normals, boundary_plane_points);
        p_Boundary_modifier->SetupSolve(cell_population,output_dir );
        std::map<unsigned, c_vector<unsigned, 2> > NearestNodesMap = p_Boundary_modifier->GetNeighbouringNodesMap();

        
    
        // // // PRINT_VARIABLE( boundary_plane_points.size());
        for(unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
            c_vector<double, 3> Point = boundary_plane_points[boundary_id];

            boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population,Point,-boundary_plane_normals[boundary_id],0.5));
            simulator.AddCellPopulationBoundaryCondition(p_condition);
            // TRACE("Here")
        }

        //              //Create a plane boundary to represent the inlet and pass them to the simulation
        // c_vector<long double, 3> Boundary1 = Create_c_vector(0, 0, 0.01);
        // c_vector<long double, 3> Normal1 = -Create_c_vector(0, 0, -1);

        // c_vector<long double, 3> Boundary2 = Create_c_vector(0, 0, 0);
        // c_vector<long double, 3> Normal2 = -Create_c_vector(0, 0, 1);
        
        // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary1, Normal1, 1));
        // simulator.AddCellPopulationBoundaryCondition(p_condition_1);

        // boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary2, Normal2, 1));
        // simulator.AddCellPopulationBoundaryCondition(p_condition_2);


        //  Create a force law and pass it to the simulation
        boost::shared_ptr<LinearSpringForce<2,3> > p_linear_force(new LinearSpringForce<2,3>());
        p_linear_force->SetMeinekeSpringStiffness(0.001);
        p_linear_force->SetNearestNeighboursMap(NearestNodesMap);
        simulator.AddForce(p_linear_force);

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure((P_blood-P_tissue)*1e-2);
        p_ForceOut->SetNearestNeighboursMap(NearestNodesMap);
        simulator.AddForce(p_ForceOut);



        

     	simulator.Solve();

        // // Save the set up simulation ready to be executed once flow results are available
        // CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(&simulator);
    
        // // Output the mesh in .vtu format for HemeLB setup tool to pick up (first converted to stl, though).
        // VtkMeshWriter<2,3> mesh_writer(output_dir, "config", false);
        // MutableMesh<2,3>* p_mesh= &(dynamic_cast<MeshBasedCellPopulation<2,3>*>(&(simulator.rGetCellPopulation()))->rGetMesh());
        // p_mesh->Scale(1.0/mesh_scale,1.0/mesh_scale,1.0/mesh_scale); // so distances are back in original scale
        // mesh_writer.WriteFilesUsingMesh(*p_mesh);
    }
};

#endif /*TESTRELAXATION_HPP_*/

