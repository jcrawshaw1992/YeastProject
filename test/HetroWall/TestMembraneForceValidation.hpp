#ifndef TestBendingForceParameterSweep_HPP_
#define TestBendingForceParameterSweep_HPP_

#include <cmath>
#include <cstdio>
#include <ctime>
#include <cxxtest/TestSuite.h>

#include "RandomNumberGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "CellsGenerator.hpp"
#include "Debug.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"

#include "MutableMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "MembraneForcesBasic.hpp"
#include "MembranePropertiesSecModifier.hpp"

#include "CellMutationStatesWriter.hpp"
#include "VtkMeshReader.hpp"

#include "projects/VascularRemodelling/src/MembraneForces/MembraneStiffnessForce.hpp"
#include "/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/src/Honeycomb3DMeshGenerator.hpp"

#include "FixedRegionBoundaryCondition.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"

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

  
     void TestBendingForceOnCrinckledCylinders() throw(Exception)
    {


    unsigned Refinment[4] = {20,40,60,80};//,160};
    
    for (unsigned N_D_index = 0; N_D_index <4; N_D_index++)
    {
     

        unsigned N_D = Refinment[N_D_index];
        unsigned N_Z = Refinment[N_D_index]*0.75;//

        double Length = 1.5;
        double Radius = 0.5;

        Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
        MutableMesh<2, 3>* mesh = generator.GetMesh();
    
        std::stringstream out;
        out << N_D<< "_" << N_Z;
        std::string mesh_size = out.str();
        std::string output_directory = "TestBendingForce/CrinkledCyliners/" + mesh_size;

        
        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,mesh->GetNumNodes(), p_differentiated_type);

        // Create a cell population
        MeshBasedCellPopulation<2, 3> cell_population(*mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        // cell_population.AddCellWriter<CellMutationStatesWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(5000);
        simulator.SetDt(0.01); // 0.005
        simulator.SetSamplingTimestepMultiple(1000);
        simulator.SetUpdateCellPopulationRule(false); // No remeshing.


        /*
        -----------------------------
        MembraneProperties Modifier
        ----------------------------
        */

        double BendingConst = 3e-7;
        std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
                //         KA,          Kalpha           Ks                                 Kb
        GrowthMaps[10] = Create_c_vector(pow(10, -6.9), 0, 0, BendingConst);
        GrowthMaps[5] = Create_c_vector(pow(10, -6.9341), pow(10, -7.7), pow(10, -8), BendingConst);
        GrowthMaps[1.5] = Create_c_vector(pow(10, -6.5), pow(10, -6.3491), pow(10, -7), BendingConst);
        GrowthMaps[1.2] = Create_c_vector(pow(10, -6.2)*1e3, pow(10, -5.8360)*1e2, pow(10, -7)*1e2, BendingConst*1e1);
        GrowthMaps[2] = Create_c_vector(pow(10, -6.8), pow(10, -6.8124), pow(10, -7), BendingConst);


        boost::shared_ptr<MembranePropertiesSecModifier<2, 3> > p_Membrane_modifier(new MembranePropertiesSecModifier<2, 3>());
        p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 2, 0,10, 1); 

        /*
        -----------------------------
        Bending forces
        ----------------------------
        */

        boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
        p_membrane_force->SetupInitialMembrane(*mesh, simulator.rGetCellPopulation());
        p_membrane_force->SetMembraneStiffness(BendingConst,30,30 );
        simulator.AddForce(p_membrane_force);

          // /*
        // -----------------------------
        // Boundaries
        // ----------------------------
        // */

        //Create a plane boundary to represent the inlet and pass them to the simulation
        c_vector<double, 3> Boundary1 = Create_c_vector(0, 0, 0);
        c_vector<double, 3> Normal1 = Create_c_vector(0, 0, 1);

        c_vector<double, 3> Boundary2 = Create_c_vector(0, 0,1.5);
        c_vector<double, 3> Normal2 = Create_c_vector(0, 0, -1);

        boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_1(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary1, Normal1, 2));
        simulator.AddCellPopulationBoundaryCondition(p_condition_1);

        boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition_2(new FixedRegionBoundaryCondition<2, 3>(&cell_population, Boundary2, Normal2, 2));
        simulator.AddCellPopulationBoundaryCondition(p_condition_2);

        



        // After the inital conditions are set, distort the z component of each node          
        for (int i=0; i<mesh->GetNumNodes(); i++)
            {
                
                c_vector<double,3> InitalLocation =  cell_population.GetNode(i)->rGetLocation();

                double R = sqrt(InitalLocation[0]*InitalLocation[0]  + InitalLocation[1]*InitalLocation[1] );
                double Angle;
                if (InitalLocation[0] >= 0)
                {
                    Angle = atan(InitalLocation[1]/ InitalLocation[0]);
                }
                 else if (InitalLocation[0] < 0 && InitalLocation[1] <= 0)
                {
                  Angle = M_PI +atan( InitalLocation[1]/ InitalLocation[0]);
                }
                else if (InitalLocation[0] < 0 && InitalLocation[1] >= 0)
                {
                     Angle = -M_PI +atan( InitalLocation[1]/ InitalLocation[0]);
                }
                double Pertebation = 0.05*(0.1*RandomNumberGenerator::Instance()->randMod(50)-2.5);
                double X = R  * cos(Angle) + Pertebation ;//+ InitalLocation[2]-1;
                double Y = R * sin(Angle)+Pertebation ;
        
                c_vector<double,3>  DeformedLocation =Create_c_vector(X,Y, InitalLocation[2]);

                cell_population.GetNode(i)->rGetModifiableLocation() = DeformedLocation;        
            }







        simulator.Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);

            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
        }
    }

    void OffTestBendingForceOnBentCylinder() throw(Exception)
    {
     

     
      unsigned Refinment[6] = {20,30,40,50,60,70};//,160};
       
        for (unsigned N_D_index = 0; N_D_index <6; N_D_index++)
        {
     
        std::stringstream out;
        out << Refinment[N_D_index];
        std::string mesh_size = out.str();
        
        unsigned N_D = Refinment[N_D_index];
        unsigned N_Z = Refinment[N_D_index]/1.2;

        double Length = 1.5;
        double Radius = 0.4;

        Honeycomb3DCylinderMeshGenerator generator(N_D, N_Z, Radius, Length);
        MutableMesh<2, 3>* mesh = generator.GetMesh();
        std::string output_directory = "TestBendingForce/BunchOfBentCylinders/" + mesh_size;

        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh->GetNumNodes(), p_differentiated_type);

        // Create a cell population
        MeshBasedCellPopulation<2, 3> cell_population(*mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(100);
        simulator.SetDt(0.02); // 0.005
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetUpdateCellPopulationRule(false); // No remeshing.

            /*
        -----------------------------
        MembraneProperties Modifier
        ----------------------------
        */

        double BendingConst = 0.000000001;
        std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
                //         KA,          Kalpha           Ks                                 Kb
        GrowthMaps[10] = Create_c_vector(pow(10, -6.9), 0, 0, BendingConst);
        GrowthMaps[5] = Create_c_vector(pow(10, -6.9341), pow(10, -7.7), pow(10, -8), BendingConst);
        GrowthMaps[1.5] = Create_c_vector(pow(10, -6.5), pow(10, -6.3491), pow(10, -7), BendingConst);
        GrowthMaps[1.2] = Create_c_vector(pow(10, -6.2)*1e3, pow(10, -5.8360)*1e2, pow(10, -7)*1e2, BendingConst*1e1);
        GrowthMaps[2] = Create_c_vector(pow(10, -6.8), pow(10, -6.8124), pow(10, -7), BendingConst);


        boost::shared_ptr<MembranePropertiesSecModifier<2, 3> > p_Membrane_modifier(new MembranePropertiesSecModifier<2, 3>());
        p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 2, 0,10, 1); 


        // /*
        // -----------------------------
        // SMembrane forces
        // ----------------------------
        // */
        // boost::shared_ptr<MembraneForcesBasic> p_shear_force(new MembraneForcesBasic());
        // p_shear_force->SetupMembraneConfiguration(cell_population);
        // p_shear_force->SetNearestNodesForBoundaryNodes(NearestNodesMap);
        // simulator.AddForce(p_shear_force);

        /*
        -----------------------------
        Bending forces
        ----------------------------
        */

        boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
        p_membrane_force->SetupInitialMembrane(*mesh, simulator.rGetCellPopulation());
        p_membrane_force->SetMembraneStiffness(BendingConst,30,30 );
        simulator.AddForce(p_membrane_force);


       

            for (int i=0; i<mesh->GetNumNodes(); i++)
            {
                
                c_vector<double,3> InitalLocation =  cell_population.GetNode(i)->rGetLocation();
                    
                double X = InitalLocation[0]+ 0.4*InitalLocation[2]*InitalLocation[2] ;//+ InitalLocation[2]-1;
        
                c_vector<double,3>  DeformedLocation =Create_c_vector(X,InitalLocation[1], InitalLocation[2]);

                cell_population.GetNode(i)->rGetModifiableLocation() = DeformedLocation;        
            }
            simulator.Solve();

            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);


        }

    }


 
 
 



    void offTestBendingPlane() throw(Exception)
    {
         double N_Z = 40;
         double N_D = 40;
        std::stringstream out;
        out << N_D<< "_" << N_Z;
        std::string mesh_size = out.str();
        std::string output_directory = "TestBendingForce/" + mesh_size;

        Honeycomb3DMeshGenerator generator( N_D, N_Z, 20,20);
        MutableMesh<2, 3>* mesh = generator.GetMesh(); 
            

        // Create cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,mesh->GetNumNodes(), p_differentiated_type);

        // Create a cell population
        MeshBasedCellPopulation<2, 3> cell_population(*mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);
        cell_population.SetOutputMeshInVtk(true);
        // cell_population.AddCellWriter<CellMutationStatesWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2, 3> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(100);
        simulator.SetDt(0.02); // 0.005
        simulator.SetSamplingTimestepMultiple(400);
        simulator.SetUpdateCellPopulationRule(false); // No remeshing.


        /*
        -----------------------------
        MembraneProperties Modifier
        ----------------------------
        */

        double BendingConst = 0.001;
        std::map<double, c_vector<long double, 4> > GrowthMaps; // From matlab sweep results
                //         KA,          Kalpha           Ks                                 Kb
        GrowthMaps[10] = Create_c_vector(pow(10, -6.9), 0, 0, BendingConst);
        GrowthMaps[5] = Create_c_vector(pow(10, -6.9341), pow(10, -7.7), pow(10, -8), BendingConst);
        GrowthMaps[1.5] = Create_c_vector(pow(10, -6.5), pow(10, -6.3491), pow(10, -7), BendingConst);
        GrowthMaps[1.2] = Create_c_vector(pow(10, -6.2)*1e3, pow(10, -5.8360)*1e2, pow(10, -7)*1e2, BendingConst*1e1);
        GrowthMaps[2] = Create_c_vector(pow(10, -6.8), pow(10, -6.8124), pow(10, -7), BendingConst);


        boost::shared_ptr<MembranePropertiesSecModifier<2, 3> > p_Membrane_modifier(new MembranePropertiesSecModifier<2, 3>());
        // p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 2, 1,1e-10, 1); 
        p_Membrane_modifier->SetMembranePropeties(GrowthMaps, 2, 0,10, 1); 

        c_vector<double, 3> UpperPlaneNormal = Create_c_vector(-0.8,-0.5,0.2);
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0,0,0);

        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0.14,0.192,0);
        c_vector<double, 3> LowerPlaneNormal = Create_c_vector(0.8,-0.5,0.1);

        
        p_Membrane_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);
        p_Membrane_modifier->SetBendingForce(cell_population, BendingConst);
        simulator.AddSimulationModifier(p_Membrane_modifier);


        // /*
        // -----------------------------
        // SMembrane forces
        // ----------------------------
        // */
        // boost::shared_ptr<MembraneForcesBasic> p_shear_force(new MembraneForcesBasic());
        // p_shear_force->SetupMembraneConfiguration(cell_population);
        // p_shear_force->SetNearestNodesForBoundaryNodes(NearestNodesMap);
        // simulator.AddForce(p_shear_force);

        /*
        -----------------------------
        Bending forces
        ----------------------------
        */

        boost::shared_ptr<MembraneStiffnessForce> p_membrane_force(new MembraneStiffnessForce());
        p_membrane_force->SetupInitialMembrane(*mesh, simulator.rGetCellPopulation());
        // p_membrane_force->SetNearestNodesForBoundaryNodesBending(NearestNodesMap);
        p_membrane_force->SetMembraneStiffness(BendingConst,30,30 );
        simulator.AddForce(p_membrane_force);



        // After the inital conditions are set, distort the z component of each node          

        for (int i=0; i< mesh->GetNumNodes(); i++)
        { 
            c_vector<double,3> InitalLocation =  cell_population.GetNode(i)->rGetLocation();
            double Pertebation = 0.1*RandomNumberGenerator::Instance()->randMod(30)-1.5;

            if (std::fmod(InitalLocation[0], 4.0) == 0 || std::fmod(InitalLocation[1], 4.0 )== 0)
            {
                std::set<unsigned> NeighbouringNodeIndices = cell_population.GetNeighbouringNodeIndices(i);
                for (std::set<unsigned>::iterator iter = NeighbouringNodeIndices.begin();
                    iter != NeighbouringNodeIndices.end();
                    ++iter)
                {           
                        cell_population.GetNode(*iter)->rGetModifiableLocation()[2] = 3*Pertebation/4;      
                }
                cell_population.GetNode(i)->rGetModifiableLocation()[2] = Pertebation;       

            }
        }


        simulator.Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(&simulator);

    }


 
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/