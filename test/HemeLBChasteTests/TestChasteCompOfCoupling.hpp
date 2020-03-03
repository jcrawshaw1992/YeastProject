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
#include "GeneralisedLinearSpringForce.hpp"

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



using namespace xsd::cxx::tree;

class TestSetupFlowInPipe : public AbstractCellBasedTestSuite
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
            	rBoundaryPlaneRadii.push_back(40e-6); // in m
            }
        }
    }

public:

    void TestSetupRunAndSavePipe() throw (Exception)
    {
        
        std::string working_directory = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/ShrunkPlexusWithLongInlets/" ;
		
        std::string mesh_file =  working_directory + "SetUpData/config.vtu";
		std::string mesh_file_0 =  working_directory + "SetUpData/InitalCondition_Plexus.vtu";
        double mesh_scale = 1e-3; 
        std::string hemelb_config_file =  working_directory + "config.xml";
        std::string traction_file =  working_directory + "results/Extracted/surface-tractions.xtr";
        double edge_division_threshold = 1e10;
      
    
        // Read inlet position from the Hemelb input XML file
        std::vector<c_vector<double,3> > boundary_plane_points;
        std::vector<c_vector<double,3> > boundary_plane_normals;
        std::vector<double> boundary_plane_radii;
        // TRACE(hemelb_config_file);
        ReadIoletPlanes(hemelb_config_file, boundary_plane_points, boundary_plane_normals, boundary_plane_radii);
        // TRACE(hemelb_config_file);
        // TRACE(mesh_file);
// 
        VtkMeshReader<2,3> mesh_reader(mesh_file);
        MutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Scale(mesh_scale,mesh_scale,mesh_scale); // so distances are in mm
        
		
		// And do this for the inital condition 
		VtkMeshReader<2,3> mesh_0_reader(mesh_file_0);
        MutableMesh<2,3> mesh_0;
        mesh_0.ConstructFromMeshReader(mesh_0_reader);
        mesh_0.Scale(mesh_scale,mesh_scale,mesh_scale); // so distances are in mm


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



        std::vector<CellPtr> cells2;
        cells_generator.GenerateBasicRandom(cells2, mesh_0.GetNumNodes(), p_differentiated_type);


		 // Create a cell population
        MeshBasedCellPopulation<2,3> cell_0_population(mesh_0, cells2);
		cell_0_population.SetWriteVtkAsPoints(true);
        cell_0_population.SetOutputMeshInVtk(true);
        cell_0_population.CalculateRestLengths();

        // Set population to output all data to results files
        cell_0_population.AddCellWriter<CellMutationStatesWriter>();
        cell_0_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_0_population.AddCellWriter<CellIdWriter>();





        // Set up cell-based simulation
        OffLatticeSimulation<2,3> simulator(cell_population);
        std::string output_dir = "TestNewInitalCondition/";
        
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetDt(0.002);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(10);


     // Create an Applied Force modifier to couple to Flow
        boost::shared_ptr<AppliedForceModifier<2,3> > p_force_modifier(new AppliedForceModifier<2,3>());
        p_force_modifier->SetResetTractionsOnCells(false,"");
        p_force_modifier->SetEdgeDivisionThreshold(edge_division_threshold);
        p_force_modifier->SetResetTractionsOnCells(true, traction_file);
        simulator.AddSimulationModifier(p_force_modifier);
        

        // boost::shared_ptr<AppliedForce<2,3>> p_pressure_force(new AppliedForce<2,3>());
        // simulator.AddForce(p_pressure_force);

        /*
        -----------------------------
       Compressive tissue pressure
        ----------------------------
        */

         double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.1000e-05 mmHg

        // double TransmuralPressure =  P_blood - P_tissue;

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure(-P_tissue);
        simulator.AddForce(p_ForceOut);



    //     /*
    //     -----------------------------
    //     MembraneProperties Modifier
    //     ----------------------------
    //     */

        double BendingConst = 1e-13;
        std::map<double, c_vector<long double, 4> > GrowthMaps;  // From matlab sweep results 
                                    //         KA,          Kalpha           Ks
        GrowthMaps[10] = Create_c_vector(pow(10,-6.9), pow(10,-8.2459),pow(10, -9),BendingConst );
        GrowthMaps[8] = Create_c_vector(pow(10,-6.9), pow(10,-8.0160),pow(10, -9) ,BendingConst);
        GrowthMaps[6] = Create_c_vector(pow(10,-6.9), pow(10,-7.7300),pow(10, -9) ,BendingConst);

        GrowthMaps[5] = Create_c_vector(pow(10,-6.9341), pow(10,-7.7),pow(10, -8) ,BendingConst );
        GrowthMaps[4] = Create_c_vector(pow(10,-6.9), pow(10,-7.4224),pow(10, 8),BendingConst);

        GrowthMaps[2] = Create_c_vector(pow(10,-6.8), pow(10,-6.8124),pow(10, -7) ,BendingConst );
        GrowthMaps[1.5] = Create_c_vector(pow(10,-6.5), pow(10,-6.3491),pow(10, -7) ,BendingConst );
        GrowthMaps[1.2] =  Create_c_vector(pow(10,-6.2), pow(10,-5.8360),pow(10, -7),BendingConst  );

        boost::shared_ptr<MembranePropertiesSecModifier<2, 3> > p_Membrane_modifier(new MembranePropertiesSecModifier<2, 3>());
        p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 1.2, 0,10, 1); 
        p_Membrane_modifier->SetupSolve(cell_population,output_dir);

     

         /*
        -----------------------------
        Bending forces
        ----------------------------
        */

   

        // // // Create a plane boundary to represent the inlets/outlets and pass them to the simulation

        boost::shared_ptr<BoundariesModifier<2, 3> > p_Boundary_modifier(new BoundariesModifier<2, 3>());
        p_Boundary_modifier->CreateBoundaryNodes(cell_population,boundary_plane_normals, boundary_plane_points);
        p_Boundary_modifier->SetupSolve(cell_population,output_dir );
        std::map<unsigned, c_vector<unsigned, 5> > NearestNodesMap = p_Boundary_modifier->GetNeighbouringNodesMap();

   
     


        // PRINT_VARIABLE( boundary_plane_points.size());
        for(unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
            boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population,boundary_plane_points[boundary_id],-boundary_plane_normals[boundary_id],boundary_plane_radii[boundary_id]));
            simulator.AddCellPopulationBoundaryCondition(p_condition);
            TRACE("Here")
        }
         /*
        -----------------------------
        SMembrane forces
        ----------------------------
        //  */
        boost::shared_ptr<MembraneForcesBasic> p_shear_force(new MembraneForcesBasic());
        p_shear_force->SetupMembraneConfiguration(cell_0_population);
        p_shear_force->SetNearestNodesForBoundaryNodes(NearestNodesMap);
        simulator.AddForce(p_shear_force);

        boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
        p_membrane_force->SetupInitialMembrane(mesh_0, cell_0_population);
        p_membrane_force->SetNearestNodesForBoundaryNodesBending(NearestNodesMap);
        p_membrane_force->SetMembraneStiffness(BendingConst,30,30 );
        simulator.AddForce(p_membrane_force);






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

