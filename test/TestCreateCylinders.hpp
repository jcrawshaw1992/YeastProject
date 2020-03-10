
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
#include "FixedRegionBoundaryCondition.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"

#include "OffLatticeSimulation.hpp"

#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"

#include "BetaCateninOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"

#include "MutableMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"

// #include "MembraneShearForce.hpp"
#include "MembraneForcesBasic.hpp"

#include "ConstantPressure.hpp"

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

void TestEquationSystem() throw(Exception)
    {

Valuable a1, a2, b1, b2; // init with values

    System sys;
    Variable x,y;
    sys << (x-a1)^2 + (y-b1)^2 - c1; // addin an equation as an equality to 0
    sys << (x-a2)^2 + (y-b2)^2 - c2;

    for(auto& solution : sys.Solve(x))
            std::cout << solution;
    }

void OffTest_CreateBendEllipticalCylinders_WithDiffereingRefinments() throw(Exception)
    {
     

      unsigned Refinment[11] = {20,30,40,50,55,60,65,70,75,80,85};//,160};
       
        for (unsigned N_D_index = 0; N_D_index < 11; N_D_index++)
        {
            std::stringstream out;
            out << Refinment[N_D_index];
            std::string mesh_size = out.str();
        
        unsigned N_D = Refinment[N_D_index];
        unsigned N_Z = Refinment[N_D_index]/1.2;

        double Length = 1.5;
        double Radius = 0.4;

        Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
        MutableMesh<2, 3>* p_mesh = generator.GetMesh();

        std::string output_directory = "CylinderCollection/BunchOfBendEllicpses/" + mesh_size+"/Inital";

        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Create a cell population
        MeshBasedCellPopulation<2, 3> cell_population(*p_mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(0.5);
        simulator.SetDt(0.5); // 0.005
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetUpdateCellPopulationRule(false); // No remeshing.

        simulator.Solve();

            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);


            for (int i=0; i<p_mesh->GetNumNodes(); i++)
            {
                
                c_vector<double,3> InitalLocation =  cell_population.GetNode(i)->rGetLocation();

                double R = sqrt(InitalLocation[0]*InitalLocation[0]  + InitalLocation[1]*InitalLocation[1] );
                double Angle;
                double Scalled_R;
                if (InitalLocation[0] >= 0)
                {
                    Angle = atan(InitalLocation[1]/ InitalLocation[0]);
                    Scalled_R = 1.1*R;
                }
                 else if (InitalLocation[0] < 0 && InitalLocation[1] <= 0)
                {
                  Scalled_R = 1.1*R;
                  Angle = M_PI +atan( InitalLocation[1]/ InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] >= 0)
                {
                   Scalled_R = 1.1*R;
                     Angle = -M_PI +atan( InitalLocation[1]/ InitalLocation[0]);
                }
                
                double X = Scalled_R  * cos(Angle)+ 0.4*InitalLocation[2]*InitalLocation[2] ;//+ InitalLocation[2]-1;
                double Y = Scalled_R*2 * sin(Angle);
        
                c_vector<double,3>  DeformedLocation =Create_c_vector(X,Y, InitalLocation[2]);

                cell_population.GetNode(i)->rGetModifiableLocation() = DeformedLocation;        
            }
            std::string output_directory2 = "CylinderCollection/BunchOfBendEllicpses/" + mesh_size+"/Deformed";
            simulator.SetOutputDirectory(output_directory2);

            simulator.SetEndTime(0.5);
            simulator.SetDt(0.5); // 0.005
            simulator.SetSamplingTimestepMultiple(1);
            simulator.SetUpdateCellPopulationRule(false); // No remeshing.
            simulator.Solve();

            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

        }

    }


 
 
 
 


 void ofoTest_CreateAndDeformCollectionOfCylinders_WithDiffereingRefinments() throw(Exception)
    {
     
      unsigned Refinment[11] = {20,30,40,50,55,60,65,70,75,80,85};//,160};
       
        for (unsigned N_D_index = 0; N_D_index < 11; N_D_index++)
        {
            std::stringstream out;
            out << Refinment[N_D_index];
            std::string mesh_size = out.str();
        
        unsigned N_D = Refinment[N_D_index];
        unsigned N_Z = Refinment[N_D_index]/2;

        double Length = 1.5/2;
        double Radius = 0.5;

        Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
        MutableMesh<2, 3>* p_mesh = generator.GetMesh();

        std::string output_directory = "CylinderCollection/BunchOfCylinder/" + mesh_size+"/Inital";

        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Create a cell population
        MeshBasedCellPopulation<2, 3> cell_population(*p_mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(0.5);
        simulator.SetDt(0.5); // 0.005
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetUpdateCellPopulationRule(false); // No remeshing.

        simulator.Solve();

            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);


            for (int i=0; i<p_mesh->GetNumNodes(); i++)
            {
                
                c_vector<double,3> InitalLocation =  cell_population.GetNode(i)->rGetLocation();

                double R = sqrt(InitalLocation[0]*InitalLocation[0]  + InitalLocation[1]*InitalLocation[1] );
                double Angle;
                double Scalled_R;
                if (InitalLocation[0] >= 0)
                {
                    Angle = atan(InitalLocation[1]/ InitalLocation[0]);
                    Scalled_R = 2*R;
                }
                 else if (InitalLocation[0] < 0 && InitalLocation[1] <= 0)
                {
                  Scalled_R = 2*R;
                  Angle = M_PI +atan( InitalLocation[1]/ InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] >= 0)
                {
                   Scalled_R = 2*R;
                     Angle = -M_PI +atan( InitalLocation[2-1]/ InitalLocation[1-1]);
                }
                
                double X = Scalled_R  * cos(Angle);
                double Y = Scalled_R  * sin(Angle);
        
                c_vector<double,3>  DeformedLocation =Create_c_vector(X,Y, InitalLocation[2]);

                cell_population.GetNode(i)->rGetModifiableLocation() = DeformedLocation;        
            }
            std::string output_directory2 = "CylinderCollection/BunchOfCylinder/" + mesh_size+"/Deformed";
            simulator.SetOutputDirectory(output_directory2);

            simulator.SetEndTime(0.5);
            simulator.SetDt(0.5); // 0.005
            simulator.SetSamplingTimestepMultiple(1);
            simulator.SetUpdateCellPopulationRule(false); // No remeshing.
            simulator.Solve();

            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

        }

    }


 
 
 
 
 
 
 
    void offTestRectLong() throw(Exception)
    {
     
        unsigned N_D = 8;
        unsigned N_Z = 3;

        double L = 4; //12e-3; //12e-3
        double R = 4;

        // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
        // in um will be too large and break chaste without carefull playing with or a tiny time step

        Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, R, L);
        MutableMesh<2, 3>* p_mesh = generator.GetMesh();


        std::string output_directory = "CylinderCollection/Rectangle_Wider";

        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Create a cell population
        MeshBasedCellPopulation<2, 3> cell_population(*p_mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(0.05);
        simulator.SetDt(0.05); // 0.005
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetUpdateCellPopulationRule(false); // No remeshing.

    
        cell_population.GetNode(1)->rGetModifiableLocation()[0] = R/2; cell_population.GetNode(1)->rGetModifiableLocation()[1] = R/2; 
        cell_population.GetNode(3)->rGetModifiableLocation()[0] = -R/2;   cell_population.GetNode(3)->rGetModifiableLocation()[1] = R/2; 
        cell_population.GetNode(5)->rGetModifiableLocation()[0] = -R/2; cell_population.GetNode(5)->rGetModifiableLocation()[1] = -R/2; 
        cell_population.GetNode(7)->rGetModifiableLocation()[0] = R/2;   cell_population.GetNode(7)->rGetModifiableLocation()[1] = -R/2; 
     
        cell_population.GetNode(8)->rGetModifiableLocation() =Create_c_vector(R, 0, L/2); 
        cell_population.GetNode(9)->rGetModifiableLocation() =Create_c_vector(R/2 ,R/2, L/2); 
        cell_population.GetNode(10)->rGetModifiableLocation() =Create_c_vector(0 ,R, L/2); 
        cell_population.GetNode(11)->rGetModifiableLocation() =Create_c_vector(-R/2 ,R/2, L/2); 
        cell_population.GetNode(12)->rGetModifiableLocation() =Create_c_vector(-R ,0, L/2 ); 
        cell_population.GetNode(13)->rGetModifiableLocation() =Create_c_vector(-R/2 ,-R/2 ,L/2 ); 
        cell_population.GetNode(14)->rGetModifiableLocation() =Create_c_vector(0 ,-R, L/2 ); 
        cell_population.GetNode(15)->rGetModifiableLocation() =Create_c_vector(R/2 ,-R/2 ,L/2  ); 


        cell_population.GetNode(16)->rGetModifiableLocation() =Create_c_vector(R, 0, L ); 
        cell_population.GetNode(17)->rGetModifiableLocation() =Create_c_vector(R/2 ,R/2, L ); 
        cell_population.GetNode(18)->rGetModifiableLocation() =Create_c_vector(0 ,R, L ); 
        cell_population.GetNode(19)->rGetModifiableLocation() =Create_c_vector(-R/2 ,R/2 ,L  ); 
        cell_population.GetNode(20)->rGetModifiableLocation() =Create_c_vector(-R ,0, L  ); 
        cell_population.GetNode(21)->rGetModifiableLocation() =Create_c_vector(-R/2 ,-R/2 ,L ); 
        cell_population.GetNode(22)->rGetModifiableLocation() =Create_c_vector(0 ,-R, L ); 
        cell_population.GetNode(23)->rGetModifiableLocation() =Create_c_vector(R/2 ,-R/2 ,L  ); 

        simulator.Solve();

        // To reset before looM_PIng: this is usually done by the SetUp and TearDown methods
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
    }
    void offTestCylinderRect() throw(Exception)
    {
     
        unsigned N_D = 8;
        unsigned N_Z = 3;

        double Length = 2; //12e-3; //12e-3
        double Radius = 0.5;

        // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
        // in um will be too large and break chaste without carefull playing with or a tiny time step

        Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
        MutableMesh<2, 3>* p_mesh = generator.GetMesh();


        std::string output_directory = "CylinderCollection/Rectangle";

        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Create a cell population
        MeshBasedCellPopulation<2, 3> cell_population(*p_mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(0.05);
        simulator.SetDt(0.05); // 0.005
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetUpdateCellPopulationRule(false); // No remeshing.

    
        cell_population.GetNode(1)->rGetModifiableLocation()[0] = 0.25; cell_population.GetNode(1)->rGetModifiableLocation()[1] = 0.25; 
        cell_population.GetNode(3)->rGetModifiableLocation()[0] = -0.25;   cell_population.GetNode(3)->rGetModifiableLocation()[1] = 0.25; 
        cell_population.GetNode(5)->rGetModifiableLocation()[0] = -0.25; cell_population.GetNode(5)->rGetModifiableLocation()[1] = -0.25; 
        cell_population.GetNode(7)->rGetModifiableLocation()[0] = 0.25;   cell_population.GetNode(7)->rGetModifiableLocation()[1] = -0.25; 
     
        cell_population.GetNode(8)->rGetModifiableLocation() =Create_c_vector(0.5, 0, 1 ); 
        cell_population.GetNode(9)->rGetModifiableLocation() =Create_c_vector(0.25 ,0.25, 1); 
        cell_population.GetNode(10)->rGetModifiableLocation() =Create_c_vector(0 ,0.5, 1 ); 
        cell_population.GetNode(11)->rGetModifiableLocation() =Create_c_vector(-0.25 ,0.25 ,1  ); 
        cell_population.GetNode(12)->rGetModifiableLocation() =Create_c_vector(-0.5 ,0, 1  ); 
        cell_population.GetNode(13)->rGetModifiableLocation() =Create_c_vector(-0.25 ,-0.25 ,1  ); 
        cell_population.GetNode(14)->rGetModifiableLocation() =Create_c_vector(0 ,-0.5, 1  ); 
        cell_population.GetNode(15)->rGetModifiableLocation() =Create_c_vector(0.25 ,-0.25 ,1   ); 


        cell_population.GetNode(16)->rGetModifiableLocation() =Create_c_vector(0.5, 0, 2 ); 
        cell_population.GetNode(17)->rGetModifiableLocation() =Create_c_vector(0.25 ,0.25, 2); 
        cell_population.GetNode(18)->rGetModifiableLocation() =Create_c_vector(0 ,0.5, 2 ); 
        cell_population.GetNode(19)->rGetModifiableLocation() =Create_c_vector(-0.25 ,0.25 ,2  ); 
        cell_population.GetNode(20)->rGetModifiableLocation() =Create_c_vector(-0.5 ,0, 2  ); 
        cell_population.GetNode(21)->rGetModifiableLocation() =Create_c_vector(-0.25 ,-0.25 ,2 ); 
        cell_population.GetNode(22)->rGetModifiableLocation() =Create_c_vector(0 ,-0.5, 2  ); 
        cell_population.GetNode(23)->rGetModifiableLocation() =Create_c_vector(0.25 ,-0.25 ,2   ); 





        simulator.Solve();

        // To reset before looping: this is usually done by the SetUp and TearDown methods
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
    } 

};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/







