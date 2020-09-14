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


#include "projects/VascularRemodelling/src/MembraneForces/MembraneShearForce.hpp"
#include "projects/VascularRemodelling/src/MembraneForces/MembraneStiffnessForce.hpp"
 #include "MembraneSurfaceForce.hpp"

 #include "EmptyBasementMatrix.hpp"
 #include "LostEndothelialCell.hpp"
 #include "HasEndothelialCell.hpp"

  #include "EdgeCorrectionForce.hpp"
  #include "BetaCateninOneHitCellMutationState.hpp"



#include "MembraneForcesBasic.hpp"
#include "MembranePropertiesSecModifier.hpp"



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
            PRINT_2_VARIABLES(iolet_id, iolets[iolet_id]);
        }
        for (unsigned iolet_id=0; iolet_id<iolets.size(); iolet_id++)
        {
            PRINT_2_VARIABLES(iolet_id, iolets[iolet_id]);
          
            c_vector<double,3> vector;

            std::vector<xercesc::DOMElement*> points = XmlTools::FindElements(iolets[iolet_id], "position");
       
            assert(points.size()==1);
      
            PRINT_VECTOR(vector);
            PRINT_VARIABLE(points[0]);
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
        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-mesh"));
        std::string mesh_file = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-mesh");
        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-mesh_scale"));
        double mesh_scale = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-mesh_scale").c_str());

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-xml"));
        std::string hemelb_config_file = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-xml");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-div_threshold"));
        double edge_division_threshold = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-div_threshold").c_str());

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-wss_threshold"));
        double wss_threshold = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-wss_threshold").c_str());

       TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-dt"));
        double dt = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-dt").c_str());
        PRINT_VARIABLE(dt);

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-SamplingTimestepMultiple"));
        double SamplingTimestepMultiple = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-SamplingTimestepMultiple").c_str());

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-Scalling"));
        double Scalling = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-Scalling").c_str());

        // Read inlet position from the Hemelb input XML file
        std::vector<c_vector<double,3> > boundary_plane_points;
        std::vector<c_vector<double,3> > boundary_plane_normals;
        std::vector<double> boundary_plane_radii;

        ReadIoletPlanes(hemelb_config_file, boundary_plane_points, boundary_plane_normals, boundary_plane_radii);
  
  
// 
        PRINT_VARIABLE(mesh_file);
        VtkMeshReader<2,3> mesh_reader(mesh_file);
    
        MutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Scale(mesh_scale,mesh_scale,mesh_scale); // so distances are in m

        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        MeshBasedCellPopulation<2,3> cell_population(mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        cell_population.CalculateRestLengths();
        // cell_population.SetDampingConstantNormal(drag_coeficient);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();

        // Iterate over all cells and lable them as HasEndothelialCell
        MAKE_PTR(HasEndothelialCell, p_mutation_state);
        for (AbstractCellPopulation<2, 3>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
            {
                    cell_iter->SetMutationState(p_mutation_state);
                    cell_iter->GetCellData()->SetItem("Boundary", 0);
            }


        // Set up cell-based simulation
        OffLatticeSimulation<2,3> simulator(cell_population);
        std::string output_dir = "./";
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        simulator.SetDt(dt);
        simulator.SetUpdateCellPopulationRule(false);


        // Mark the edge cells as mutant cells 

        double Length = 12e-3; //12e-3;
        double trans = -6e-3;
        double MaxZ = Length + trans;
        double MinZ = trans;
  TRACE("here")

     // Create an Applied Force modifier to couple to Flow
        boost::shared_ptr<AppliedForceModifier<2,3> > p_force_modifier(new AppliedForceModifier<2,3>());
        p_force_modifier->SetResetTractionsOnCells(false,"");
        // p_force_modifier->SetMembraneConstants(ElasticShearModulus , AreaDilationModulus, Area_constant, membrane_constant );
        simulator.AddSimulationModifier(p_force_modifier);

        boost::shared_ptr<AppliedForce<2,3>> p_pressure_force(new AppliedForce<2,3>());
        simulator.AddForce(p_pressure_force);


        double BendingConst = 1e-17;
        std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
                //         KA,          Kalpha           Ks                                 Kb
        GrowthMaps[10] = Create_c_vector(pow(10, -6.9), pow(10, -8.2459), pow(10, -9), 0);
        GrowthMaps[5] = Create_c_vector(pow(10, -6.9341), pow(10, -7.7), pow(10, -8), BendingConst);
        GrowthMaps[1.5] = Create_c_vector(pow(10, -6.5), pow(10, -6.3491), pow(10, -7), BendingConst);
        GrowthMaps[1.2] = Create_c_vector(pow(10, -6.2)*1e3, pow(10, -5.8360)*1e2, pow(10, -7)*1e2, BendingConst*1e1);
        GrowthMaps[2] = Create_c_vector(pow(10, -6.8), pow(10, -6.8124), pow(10, -7), BendingConst);


        boost::shared_ptr<MembranePropertiesSecModifier<2, 3> > p_Membrane_modifier(new MembranePropertiesSecModifier<2, 3>());
        
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
        p_membrane_force->SetupInitialMembrane(mesh, simulator.rGetCellPopulation());
        p_membrane_force->SetMembraneStiffness(BendingConst,30,30 );
        simulator.AddForce(p_membrane_force);

        // Create a plane boundary to represent the inlets/outlets and pass them to the simulationVtkMeshReader
        double inlet_offset = 0.001;
        for(unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
            boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population,boundary_plane_points[boundary_id],-boundary_plane_normals[boundary_id],boundary_plane_radii[boundary_id] +10));
            simulator.AddCellPopulationBoundaryCondition(p_condition);
            // PRINT_2_VARIABLES(-boundary_plane_normals[boundary_id],boundary_plane_radii[boundary_id]);


        }

        // Save the set up simulation ready to be executed once flow results are available
        CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(&simulator);
    
        // Output the mesh in .vtu format for HemeLB setup tool to pick up (first converted to stl, though).
        VtkMeshWriter<2,3> mesh_writer(output_dir, "config", false);   
        MutableMesh<2,3>* p_mesh= &(dynamic_cast<MeshBasedCellPopulation<2,3>*>(&(simulator.rGetCellPopulation()))->rGetMesh());
        p_mesh->Scale(1.0/mesh_scale,1.0/mesh_scale,1.0/mesh_scale); // so distances are back in original scale
        mesh_writer.WriteFilesUsingMesh(*p_mesh);
    }
};

#endif /*TESTRELAXATION_HPP_*/

