#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "UblasCustomFunctions.hpp"
#include <cxxtest/TestSuite.h>
// #include <cstdio>
// #include <ctime>
// #include <cmath>
// #include <vector> 

#include "Debug.hpp"
#include "VtkMeshReader.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "SmartPointers.hpp"
 
#include "AbstractCellBasedTestSuite.hpp"
#include "OffLatticeSimulation.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "HistoryDepMutableMesh.hpp"

// #include "AppliedForceModifier.hpp"
#include "RemeshingTriggerOnHeteroMeshModifier.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "MembraneDeformationForce.hpp"
// #include "OutwardsPressureWithBreaks.hpp"
#include "OutwardsPressure.hpp"
#include "MembraneBendingForce.hpp"
#include "HemeLBForce.hpp"


class TestRemeshing  : public AbstractCellBasedTestSuite
{
public:

  void TestParametersOnBifucation() throw (Exception)
    {


        std::map<double, c_vector<double, 3> >  GrowthMaps = { { 1, Create_c_vector(pow(10, -5.5),pow(10, -6.5),  pow(10, -7.5)) },
                                                            { 2, Create_c_vector(pow(10, -5.5),pow(10, -7.5),  pow(10, -6.5))  },
                                                            { 3, Create_c_vector(pow(10, -5.5),pow(10, -8),  pow(10, -6.5)) },
                                                            { 4, Create_c_vector(pow(10, -5.5),pow(10, -8.5),  pow(10, -6.5))},
                                                            { 5, Create_c_vector(pow(10, -5.5 ),pow(10, -9 ),  pow(10, -6.5 )) },
                                                            { 6, Create_c_vector(pow(10, -5.5 ),pow(10, -6.5 ),  pow(10, -10 )) }, 
                                                            { 7, Create_c_vector(pow(10, -5.5 ),pow(10, -9.5 ),  pow(10, -6.5 ))},
                                                            { 8,  Create_c_vector(pow(10, -5.5 ),pow(10, -10 ),  pow(10, -6.5 ))}};

        std::map<double, double >  BendingMap =   { {1, pow(10, -8) },{2, pow(10, -9) },  {3, pow(10, -10) },   {4, pow(10, -11) }, {5, pow(10, -12) },   {6, pow(10, -13) },   {7, pow(10, -14) }, {8, pow(10, -15) }, {9, pow(10, -16) }};                               

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-ParameterSet"));

        double PSet = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-ParameterSet");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-BendingParameter"));
        double BendingParameter = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-BendingParameter");

        double dt= 0.0005;
         if (CommandLineArguments::Instance()->OptionExists("-dt"))
        {
            dt= CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-dt");
        }
        double EndTime = 100;
        if (CommandLineArguments::Instance()->OptionExists("-EndTime"))
        {
            EndTime = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-EndTime");
        }
        double SamplingTimestepMultiple = 200;
        if (CommandLineArguments::Instance()->OptionExists("-SamplingTimestepMultiple"))
        {
            SamplingTimestepMultiple = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-SamplingTimestepMultiple");
        }
        std::string mesh_file = "/Users/jcrawshaw/Documents/Projects/MeshCollection/SimpleBifucation/mesh.vtu";
        if (CommandLineArguments::Instance()->OptionExists("-MeshFile"))
        {
            mesh_file  = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-MeshFile");
        }

        //////////
        
        double scale = 1e1;//1e-3;  
        double scale2 =  0.15;
        // std::string output_dir = "FSIIdealNetwork/SimpleBifucation_NoRemeshing/";

        std::stringstream out;
        out << "PSet" << PSet << "_Bending_" << BendingParameter;
        std::string ParameterSet = out.str();
        std::string output_dir = "TestForBuldgingInBifucation/"+ParameterSet;

        
      
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Scale(scale*scale2,scale*scale2,scale*scale2);

       // Create the cells 
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a cell population
        HistoryDepMeshBasedCellPopulation<2, 3> cell_population(mesh, cells);
        cell_population.SetChasteOutputDirectory(output_dir, 0);
        cell_population.SetInitialAnlgesAcrossMembrane();
        cell_population.SetRelativePath(output_dir, 0);
        cell_population.SetTargetRemeshingEdgeLength(0.5*1e-2*scale2); //cell_population.SetTargetRemeshingEdgeLength(0.005*scale); 
        cell_population.SetPrintRemeshedIC(1);
        cell_population.SetTargetRemeshingIterations(10);
        cell_population.EdgeLengthVariable(1.001); 
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        cell_population.SetRemeshingSoftwear("CGAL");
        
        // Set population to output all data to results files
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();

        // Binning Functions
        cell_population.SetBinningIntervals(2,2,1);
        // Set up cell-based simulation
        OffLatticeSimulation<2,3> simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        simulator.SetDt(dt);
        simulator.SetUpdateCellPopulationRule(false);
        simulator.SetEndTime(EndTime);

        /*
        -----------------------------
        RemeshingTriggerOnHeteroMeshModifier
        ----------------------------
        */  
        boost::shared_ptr<RemeshingTriggerOnHeteroMeshModifier<2, 3> > p_Mesh_modifier(new RemeshingTriggerOnHeteroMeshModifier<2, 3>());
        p_Mesh_modifier->SetMembraneParameters(GrowthMaps[PSet][0], GrowthMaps[PSet][1], GrowthMaps[PSet][2], BendingMap[BendingParameter]);
        // p_Mesh_modifier->SetRemeshingInterval(200);//
        simulator.AddSimulationModifier(p_Mesh_modifier);

        /*
        -----------------------------
        HemeLB Force
        ----------------------------
        */        

        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0);
        c_vector<double, 3> Point1 = Create_c_vector(0.85,-0.28,0)*scale2;

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(-0.7,-0.7,0);
        c_vector<double, 3> Point2 = Create_c_vector(94,-22,-0.05)*1e-2*scale2;

        c_vector<double, 3> PlaneNormal3 = Create_c_vector(-0.7,0.7,0);
        c_vector<double, 3> Point3 = Create_c_vector(94,-33,-0.05)*1e-2*scale2;


        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure((P_blood - P_tissue));// needs to be negative for server ?? 
        simulator.AddForce(p_ForceOut);

        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        simulator.AddForce(p_shear_force);


       /*
        -----------------------------
        Bending forces
        ----------------------------
        */

        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_membrane_force->SetMembraneStiffness(BendingMap[BendingParameter]);
        simulator.AddForce(p_membrane_force);

        /*
        -----------------------------
        Boundary conditions
        ----------------------------
        */
        std::vector<c_vector<double,3> > boundary_plane_points;
        std::vector<c_vector<double,3> > boundary_plane_normals;

        boundary_plane_points.push_back(Point1);
        boundary_plane_normals.push_back(PlaneNormal1);

        boundary_plane_points.push_back(Point2);
        boundary_plane_normals.push_back(PlaneNormal2);
        
        boundary_plane_points.push_back(Point3);
        boundary_plane_normals.push_back(PlaneNormal3);
        

        for(unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
          boost::shared_ptr<FixedRegionBoundaryCondition<2,3> > p_condition(new FixedRegionBoundaryCondition<2,3>(&cell_population, boundary_plane_points[boundary_id],-boundary_plane_normals[boundary_id],0.1));
          simulator.AddCellPopulationBoundaryCondition(p_condition);
        }
      
     	  simulator.Solve();
        CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(&simulator);
}


};




#endif /*TESTRELAXATION_HPP_*/

