
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
#include "Debug.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"

#include "MutableMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "MembraneForcesBasic.hpp"
#include "MembraneHetroModifier.hpp"
#include "VtkMeshReader.hpp"

// #include "ConstantPressure.hpp"
#include "CellMutationStatesWriter.hpp"
#include "EmptyBasementMatrix.hpp"
#include "HasEndothelialCell.hpp"
#include "LostEndothelialCell.hpp"

#include "OutwardsPressure.hpp"

#include "projects/VascularRemodelling/src/MembraneForces/MembraneStiffnessForce.hpp"

#include "CommandLineArguments.hpp"
#include "VtkMeshWriter.hpp"
#include "XmlTools.hpp"

using namespace xsd::cxx::tree;

static const double M_TIME_FOR_SIMULATION = 500; //40; //50
static const double M_SAMPLING_TIME_STEP = 1000000;//200;//100;
static const double M_TIME_STEP =  0.00005;// 0.0001;

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
    void ReadHemeLBv3Vector(c_vector<double, 3>& vector, xercesc::DOMElement* xmlElement)
    {

        std::stringstream raw_value(X2C(xmlElement->getAttribute(X("value"))));

        char left_par, comma1, comma2, right_par;
        double x, y, z;
        raw_value >> left_par;

        raw_value >> x;
        PRINT_VARIABLE(x);
        raw_value >> comma1;
        raw_value >> y;
        raw_value >> comma2;
        raw_value >> z;
        raw_value >> right_par;
        PRINT_4_VARIABLES(left_par, comma1, comma2, right_par)


        vector[0] = x;
        vector[1] = y;
        vector[2] = z;
       
        assert(left_par == '('); assert(comma1 == ','); assert(comma2 == ',');  assert(right_par == ')');
    }

    void ReadIoletPlanes(std::string hemelbConfigFile,
                         std::vector<c_vector<double, 3> >& rBoundaryPlanePoints,
                         std::vector<c_vector<double, 3> >& rBoundaryPlaneNormals,
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
            // std::vector<xercesc::DOMElement*> aux = XmlTools::FindElements(p_root_elt, "outlets/outlet");
            // iolets.insert(iolets.end(), aux.begin(), aux.end());
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
        for (unsigned iolet_id = 0; iolet_id < iolets.size(); iolet_id++)
        {
            c_vector<double, 3> vector;

            std::vector<xercesc::DOMElement*> points = XmlTools::FindElements(iolets[iolet_id], "position");

            assert(points.size() == 1);

            // PRINT_VECTOR(vector);
            PRINT_VARIABLE(points[0]);

            ReadHemeLBv3Vector(vector, points[0]);

            rBoundaryPlanePoints.push_back(vector);

            std::vector<xercesc::DOMElement*> normals = XmlTools::FindElements(iolets[iolet_id], "normal");
            assert(normals.size() == 1);

            ReadHemeLBv3Vector(vector, normals[0]);
            rBoundaryPlaneNormals.push_back(vector);
            std::vector<xercesc::DOMElement*> radi = XmlTools::FindElements(iolets[iolet_id], "radius");
            if (radi.size() == 1)
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
 void TestHomogeneousPlexusGrowth() throw(Exception)
    {
        double scale = 1e3;
        // Plexus
        // std::string mesh_file = "projects/VascularRemodelling/test/data/embryo_plexus/config.vtu";
        
        // std::string mesh_file = "~/docker-polnet-master/GeneratingShrunkMesh/SecondTry/NewScalledMesh.vtu";//ResultantMesh/ShrunkPlexus.vtu";
      
      std::string mesh_file = "projects/VascularRemodelling/test/data/ShrunkPlexus/Plexus.vtu";

        // This data file is in mm??
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> p_mesh;
        p_mesh.ConstructFromMeshReader(mesh_reader);
        scale = 1e-3; // so distances are in m
        p_mesh.Scale(1.0 * scale, 1.0 * scale, 1.0 * scale); // so distances are back in original scal

        std::string output_directory = "CapillaryPlexusDeformation/ShrunkMesh/";

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
        simulator.SetEndTime(100);
        simulator.SetDt(0.0001); //(0.01); // 0.005
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetUpdateCellPopulationRule(false); // No remeshing.

 
        /*
        -----------------------------
        Constant Pressure ballance 
        ----------------------------
        */

        double P_blood = 0.0021; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.0015; // Pa == 1.1000e-05 mmHg

        // For some reason on the plexus the direction of the normals are reversed??
        double TransmuralPressure = -6.6715e-4 * 1e-3; //P_blood - P_tissue;

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure(TransmuralPressure);
        simulator.AddForce(p_ForceOut);

        /*
            -----------------------------
            MembraneProperties Modifier
            ----------------------------
            */
        double bending = 1e-8;

        std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
            //         KA,          Kalpha           Ks

        GrowthMaps[1.2] = Create_c_vector(pow(10, -6.93), pow(10, -7.70), pow(10, -8.00), bending);
        GrowthMaps[1.1] = Create_c_vector(pow(10, -5.9), pow(10, -5.52), pow(10, -7), bending);

        bool Hetrogeneous = 0;
        boost::shared_ptr<MembraneHetroModifier<2, 3> > p_Membrane_modifier(new MembraneHetroModifier<2, 3>());
        p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 1.2, Hetrogeneous);
        p_Membrane_modifier->SetAConstantHetrogenaity(cell_population, output_directory);

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
        Boundary conditions
        ----------------------------
        */

        std::vector<c_vector<double,3> > boundary_plane_points;
        std::vector<c_vector<double,3> > boundary_plane_normals;

        boundary_plane_normals.push_back(Create_c_vector(-0.900883974588,-0.3901023946,-0.190337032807));
        boundary_plane_points.push_back(Create_c_vector(0.000165315,0.000290584647052,4.34131733398e-07));

        boundary_plane_normals.push_back(Create_c_vector(0.682503460834,-0.730650337725,0.0184149377318));
        boundary_plane_points.push_back(Create_c_vector(0.000108234680393,0.00029617095176,1.31970700782e-08));

        boundary_plane_normals.push_back(Create_c_vector(0.999156349137,0.0376419761764,-0.0164216810429));
        boundary_plane_points.push_back(Create_c_vector(3.42362e-05,0.00020895266296,-6.29838561806e-07));

        boundary_plane_normals.push_back(Create_c_vector(-0.81694348867,-0.563757323322,-0.121577204773));
        boundary_plane_points.push_back(Create_c_vector(0.000209418202716,0.000212406847202,-1.55505611998e-06));

        boundary_plane_normals.push_back(Create_c_vector(0.657440300446,0.746256113619,0.104278781332));
        boundary_plane_points.push_back(Create_c_vector(6.6954396938e-05,0.000132966087987,-4.12819038903e-07));

        boundary_plane_normals.push_back(Create_c_vector(0.151533835168,0.988451772619,0.000768117762117));
        boundary_plane_points.push_back(Create_c_vector(0.000119687757385,8.3646048258e-05,-1.5894283281e-06));

        boundary_plane_normals.push_back(Create_c_vector(-0.706473421705,0.70263424701,-0.0848552847671));
        boundary_plane_points.push_back(Create_c_vector(0.000200124062897,0.000124468739077,1.02362047278e-06));

        // for(unsigned boundary_id = 0; boundary_id < 7; boundary_id++)
        // {
        //  TRACE("Add bcs")
        //    boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population,boundary_plane_points[boundary_id]*1e3,-boundary_plane_normals[boundary_id],0.05));
        //     simulator.AddCellPopulationBoundaryCondition(p_condition);
        // }

        simulator.Solve();
        TRACE("SolveFinished")
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);
        TRACE("SAVED")

        // To reset before looping: this is usually done by the SetUp and TearDown methods
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
    }
 

    void offTestHetrogeneousPlexusGrowth() throw(Exception)
    {
        double scale = 1e3;
        // Plexus
        std::string mesh_file = "projects/VascularRemodelling/test/data/embryo_plexus/config.vtu";

        // This data file is in mm??
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> p_mesh;
        p_mesh.ConstructFromMeshReader(mesh_reader);
        scale = 1e-3; // so distances are in m
        p_mesh.Scale(1.0 * scale, 1.0 * scale, 1.0 * scale); // so distances are back in original scal

        std::string output_directory = "CapillaryPlexusDeformation/LowPressure/Hetro/";

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
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);
        simulator.SetDt(M_TIME_STEP); // 0.005
        simulator.SetSamplingTimestepMultiple(M_SAMPLING_TIME_STEP);
        simulator.SetUpdateCellPopulationRule(false); // No remeshing.

 
        /*
        -----------------------------
        Constant Pressure ballance 
        ----------------------------
        */

        double P_blood = 0.0021; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.0015; // Pa == 1.1000e-05 mmHg

        // For some reason on the plexus the direction of the normals are reversed??
        double TransmuralPressure = -6.6715e-4 * 1e-3;; //P_blood - P_tissue;

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure(TransmuralPressure);
        simulator.AddForce(p_ForceOut);

        /*
            -----------------------------
            MembraneProperties Modifier
            ----------------------------
            */
        double bending = 1e-8;

        std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
            //         KA,          Kalpha           Ks

        GrowthMaps[1.2] = Create_c_vector(pow(10, -6.93), pow(10, -7.70), pow(10, -8.00), bending);
        GrowthMaps[1.1] = Create_c_vector(pow(10, -5.9)*1e4, pow(10, -5.52)*1e3, pow(10, -7)*1e2, bending);

        bool Hetrogeneous = 1;
        boost::shared_ptr<MembraneHetroModifier<2, 3> > p_Membrane_modifier(new MembraneHetroModifier<2, 3>());
        p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 1.2, Hetrogeneous);
        p_Membrane_modifier->SetAConstantHetrogenaity(cell_population, output_directory);

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
        Boundary conditions
        ----------------------------
        */

        std::vector<c_vector<double,3> > boundary_plane_points;
        std::vector<c_vector<double,3> > boundary_plane_normals;

        boundary_plane_normals.push_back(Create_c_vector(-0.900883974588,-0.3901023946,-0.190337032807));
        boundary_plane_points.push_back(Create_c_vector(0.000165315,0.000290584647052,4.34131733398e-07));

        boundary_plane_normals.push_back(Create_c_vector(0.682503460834,-0.730650337725,0.0184149377318));
        boundary_plane_points.push_back(Create_c_vector(0.000108234680393,0.00029617095176,1.31970700782e-08));

        boundary_plane_normals.push_back(Create_c_vector(0.999156349137,0.0376419761764,-0.0164216810429));
        boundary_plane_points.push_back(Create_c_vector(3.42362e-05,0.00020895266296,-6.29838561806e-07));

        boundary_plane_normals.push_back(Create_c_vector(-0.81694348867,-0.563757323322,-0.121577204773));
        boundary_plane_points.push_back(Create_c_vector(0.000209418202716,0.000212406847202,-1.55505611998e-06));

        boundary_plane_normals.push_back(Create_c_vector(0.657440300446,0.746256113619,0.104278781332));
        boundary_plane_points.push_back(Create_c_vector(6.6954396938e-05,0.000132966087987,-4.12819038903e-07));

        boundary_plane_normals.push_back(Create_c_vector(0.151533835168,0.988451772619,0.000768117762117));
        boundary_plane_points.push_back(Create_c_vector(0.000119687757385,8.3646048258e-05,-1.5894283281e-06));

        boundary_plane_normals.push_back(Create_c_vector(-0.706473421705,0.70263424701,-0.0848552847671));
        boundary_plane_points.push_back(Create_c_vector(0.000200124062897,0.000124468739077,1.02362047278e-06));

        for(unsigned boundary_id = 0; boundary_id < 7; boundary_id++)
        {
         TRACE("Add bcs")
           boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population,boundary_plane_points[boundary_id]*1e3,-boundary_plane_normals[boundary_id],0.05));
            simulator.AddCellPopulationBoundaryCondition(p_condition);
        }

        simulator.Solve();
        TRACE("SolveFinished")
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);
        TRACE("SAVED")

        // To reset before looping: this is usually done by the SetUp and TearDown methods
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
    }
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/
