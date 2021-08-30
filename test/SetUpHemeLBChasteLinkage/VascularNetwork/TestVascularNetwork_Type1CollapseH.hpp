#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "PetscSetupAndFinalize.hpp"

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"

#include "Debug.hpp"
#include "Honeycomb3DCylinderMeshGenerator.hpp"
#include "SmartPointers.hpp"
#include "VtkMeshReader.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"

#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "HistoryDepMutableMesh.hpp"

// #include "AppliedForceModifier.hpp"
#include "FixedRegionBoundaryCondition.hpp"
#include "EnclosedRegionBoundaryCondition.hpp"
#include "HemeLBForce.hpp"
#include "MembraneDeformationForce.hpp"
#include "MembraneDeformationForceOnCylinder.hpp"
#include "OutwardsPressureWithBreaks.hpp"
#include "OutwardsPressure.hpp"
#include "RemeshingTriggerOnStepHeteroModifier.hpp"
#include "StepHeteroModifier.hpp"
#include "MembraneBendingForce.hpp"

#include "MembraneBendingForceSensitive.hpp"


class TestRemeshing : public AbstractCellBasedTestSuite
{
public:

  void offTestWithHemeLBForce() throw(Exception)
   {

        double AreaParameter = -5;  double DilationParameter = -5.5; double DeformationParamter = -5; double BendingParameter = -5;
        std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), pow(10, BendingParameter)) }, {0,  Create_c_vector(pow(10, -4), pow(10, -4), pow(10, -4),pow(10, BendingParameter))} };

        TRACE("Jess is good")
        double EndTime = 70;
        double scale = 0.0011; 

        double SamplingStep = 25;
        double dt = 0.001;
        double RemeshingTime = 100;
        double FSI_Iterations =50;// 50;
        double EdgeLength =1.1*scale;
        
        /////////////////////////////////////////////////////////////////////////////////////
        //  std::string Archieved = "FSISimulations/VascularNetworkLargerBendingForce/GrowingToEqui/ConstantForceArchiving/";//
         std::string Archieved = "FSISimulations/VascularNetworkLargerBendingForce/HemeLBEqui2/ConstantForceArchiving/";//std::string mesh_file = "/data/vascrem/testoutput/DeformingPlexus/FlatForceFINAL9/results_from_time_3/mesh_50.vtu";
        std::string output_dir = "FSISimulations/VascularNetworkLargerBendingForce/Type1Collapse/H/";
        // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
        // in um will be too large and break chaste without carefull playing with or a tiny time step

        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
 

        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetPrintRemeshedIC(1);

        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
        p_simulator->RemoveAllForces();
        // p_simulator->RemoveAllCellPopulationBoundaryConditions();

        /* 
        -----------------------------
        Edit  RemeshingTriggerOnStepHeteroModifier
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
        p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); // 
        p_Mesh_modifier->SetStepSize(pow(10, -8));
        p_Mesh_modifier->SetmSetUpSolve(1);


        // Thirdcollapse option 
        // Upstream 
        c_vector<double, 3> UpperPlanePoint =  Create_c_vector( 0.148450771419833, 0.17540001846872066,0 ); 
        c_vector<double, 3> UpperPlaneNormal = Create_c_vector(0.6681809321405496, -0.7397789537602159, 0);  
        // Down stream
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0.16467321416038386, 0.1592262610683606,0 ); 
        c_vector<double, 3> LowerPlaneNormal = -Create_c_vector(0.870408755948827, -0.47734167113341636, 0 ); 
  
        p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);
        p_Mesh_modifier->SetRadius(0.015);
        // p_Mesh_modifier->TurnOffRemeshing();
        p_Mesh_modifier->SetUpdateFrequency(1/dt);


        boost::shared_ptr<EnclosedRegionBoundaryCondition<2, 3> > p_condition(new EnclosedRegionBoundaryCondition<2, 3>(&(p_simulator->rGetCellPopulation()) , UpperPlanePoint, UpperPlaneNormal, 0.01)); //0.01));

        p_condition->SetPointOnPlane2( LowerPlanePoint);
        p_condition->SetNormalToPlane2(-LowerPlaneNormal);
        p_simulator->AddCellPopulationBoundaryCondition(p_condition);
       
       /* 
        -----------------------------
        finished editing  RemeshingTriggerOnStepHeteroModifier
        ----------------------------
        */

    
        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_shear_force->SetCollapseType(1);
        p_simulator->AddForce(p_shear_force);

        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_membrane_force->SetMembraneStiffness(pow(10, -7));
        p_membrane_force->SetCollapseType(1);
        p_simulator->AddForce(p_membrane_force);


        /*
        -----------------------------
        Add the HemeLB Force
        ----------------------------
        */
       // Inlet1
        // Inlet1
       c_vector<double, 3> Point1 = Create_c_vector(0.012896438155977795,0.1296022543191808, -0.002236333573976076 );
       c_vector<double, 3> PlaneNormal1 = Create_c_vector(0.9491166468623324, -0.3145145356937025,0 );
        

       c_vector<double, 3> Point2 = Create_c_vector(0.04152531239022141, 0.23542737676454184, -0.0007601220860857502);
       c_vector<double, 3> PlaneNormal2 = Create_c_vector(0.7969315092041309, -0.604064075870372, 0 );
       

       c_vector<double, 3> Point3 = Create_c_vector(0.1651556668051053,0.02512909630785012,0.0011455782250521602);
       c_vector<double, 3> PlaneNormal3 = Create_c_vector(0.0807491790098989,0.9963585613635368,  0.027371285808498066 );

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg

        double InletPressure = P_blood; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later
        double OutletPressure = P_blood * (0.9);

        boost::shared_ptr<HemeLBForce> p_ForceOut(new HemeLBForce());
        p_ForceOut->Inlets(PlaneNormal1, Point1, OutletPressure, "Outlet"); // Issues here 
        p_ForceOut->Inlets(PlaneNormal2, Point2, OutletPressure, "Outlet"); //FIne
        p_ForceOut->Inlets(PlaneNormal3, Point3, InletPressure, "Inlet");// Issues here 
        p_ForceOut->SetStartTime(EndTime);
        p_ForceOut->Network("VascularNetwork");
        p_ForceOut->SetCollapseType(1);
        p_ForceOut->SetFluidSolidIterations(FSI_Iterations);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir+"HemeLBForce/", p_simulator->rGetCellPopulation(),1);
        p_simulator->AddForce(p_ForceOut);


      for (int i =1; i<=50; i++)
        { 
    
            EndTime +=1;
            p_simulator->SetEndTime(EndTime);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }

    }



    // source ~/.bash_profile

  void TestType2Collapse() throw(Exception)
   {

        double AreaParameter = -5;  double DilationParameter = -5.5; double DeformationParamter = -5; double BendingParameter = -5;
        std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), pow(10, BendingParameter)) }, {0,  Create_c_vector(pow(10, -4), pow(10, -4), pow(10, -4),pow(10, BendingParameter))} };

        TRACE("Jess is good")
        double EndTime = 88;
        double scale = 0.0011; 

        double SamplingStep = 1;
        double dt = 0.001;
        double RemeshingTime = 100;
        double FSI_Iterations =50;// 50;
        double EdgeLength =0.9*scale;
        
        /////////////////////////////////////////////////////////////////////////////////////
        //  std::string Archieved = "FSISimulations/VascularNetworkLargerBendingForce/GrowingToEqui/ConstantForceArchiving/";//
        std::string Archieved = "FSISimulations/VascularNetworkLargerBendingForce/Type1Collapse/H";
        std::string output_dir = "FSISimulations/VascularNetworkLargerBendingForce/Type2Collapse/H/";
        // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
        // in um will be too large and break chaste without carefull playing with or a tiny time step

        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
 

        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetPrintRemeshedIC(1);

    // static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).ExecuteHistoryDependentRemeshing();

        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
        p_simulator->RemoveAllForces();
        // p_simulator->RemoveAllCellPopulationBoundaryConditions();

        /* 
        -----------------------------
        Edit  RemeshingTriggerOnStepHeteroModifier
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
        p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); // 
        p_Mesh_modifier->SetStepSize(pow(10, -8));
        p_Mesh_modifier->SetmSetUpSolve(1);


        // Thirdcollapse option 
        // Upstream 
        c_vector<double, 3> UpperPlanePoint =  Create_c_vector( 0.148450771419833, 0.17540001846872066,0 ); 
        c_vector<double, 3> UpperPlaneNormal = Create_c_vector(0.6681809321405496, -0.7397789537602159, 0);  
        // Down stream
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0.16467321416038386, 0.1592262610683606,0 ); 
        c_vector<double, 3> LowerPlaneNormal = -Create_c_vector(0.870408755948827, -0.47734167113341636, 0 ); 
  
        p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);
        p_Mesh_modifier->SetRadius(0.015);
        p_Mesh_modifier->TurnOffRemeshing();
        p_Mesh_modifier->SetUpdateFrequency(1/dt);

        //       // Thirdcollapse option 
        // // Upstream 
        // c_vector<double, 3> UpperPlanePoint =  Create_c_vector( 0.15429402927086544, 0.179490164745615,0 ); 
        // c_vector<double, 3> UpperPlaneNormal = Create_c_vector(0.7229590897477705, -0.6905862771818618,0);  
        // // Down stream
        // c_vector<double, 3> LowerPlanePoint = Create_c_vector( 0.16827330952059658, 0.16613685161594424,0); 
        // c_vector<double, 3> LowerPlaneNormal = -Create_c_vector(0.7229590897477705, -0.6905862771818618,0 ); 


        boost::shared_ptr<EnclosedRegionBoundaryCondition<2, 3> > p_condition(new EnclosedRegionBoundaryCondition<2, 3>(&(p_simulator->rGetCellPopulation()) , UpperPlanePoint, UpperPlaneNormal, 0.01)); //0.01));

        p_condition->SetPointOnPlane2( LowerPlanePoint);
        p_condition->SetNormalToPlane2(-LowerPlaneNormal);
        p_simulator->AddCellPopulationBoundaryCondition(p_condition);
       



    
        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_shear_force->SetCollapseType(2);
        p_simulator->AddForce(p_shear_force);

        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_membrane_force->SetMembraneStiffness(pow(10, -7));
        p_membrane_force->SetCollapseType(2);
        p_simulator->AddForce(p_membrane_force);


        /*
        -----------------------------
        Add the HemeLB Force
        ----------------------------
        */
       // Inlet1
        // Inlet1
       c_vector<double, 3> Point1 = Create_c_vector(0.012896438155977795,0.1296022543191808, -0.002236333573976076 );
       c_vector<double, 3> PlaneNormal1 = Create_c_vector(0.9491166468623324, -0.3145145356937025,0 );
        

       c_vector<double, 3> Point2 = Create_c_vector(0.04152531239022141, 0.23542737676454184, -0.0007601220860857502);
       c_vector<double, 3> PlaneNormal2 = Create_c_vector(0.7969315092041309, -0.604064075870372, 0 );
       

       c_vector<double, 3> Point3 = Create_c_vector(0.1651556668051053,0.02512909630785012,0.0011455782250521602);
       c_vector<double, 3> PlaneNormal3 = Create_c_vector(0.0807491790098989,0.9963585613635368,  0.027371285808498066 );

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg

        double InletPressure = P_blood; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later
        double OutletPressure = P_blood * (0.9);

        boost::shared_ptr<HemeLBForce> p_ForceOut(new HemeLBForce());
        p_ForceOut->Inlets(PlaneNormal1, Point1, OutletPressure, "Outlet"); // Issues here 
        p_ForceOut->Inlets(PlaneNormal2, Point2, OutletPressure, "Outlet"); //FIne
        p_ForceOut->Inlets(PlaneNormal3, Point3, InletPressure, "Inlet");// Issues here 
        p_ForceOut->SetStartTime(EndTime);
        p_ForceOut->Network("VascularNetwork");
        p_ForceOut->SetCollapseType(2);
        p_ForceOut->SetFluidSolidIterations(FSI_Iterations);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir+"HemeLBForce/", p_simulator->rGetCellPopulation(),0);
        p_simulator->AddForce(p_ForceOut);


      for (int i =1; i<=50; i++)
        { 
    
            EndTime +=1;
            p_simulator->SetEndTime(EndTime);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }

    }






  void offTestWithHemeLBForceFiner() throw(Exception)
   {

        double AreaParameter = -5;  double DilationParameter = -5.5; double DeformationParamter = -5; double BendingParameter = -5;
        std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), pow(10, BendingParameter)) }, {0,  Create_c_vector(pow(10, -4), pow(10, -4), pow(10, -4),pow(10, BendingParameter))} };

        
        TRACE("Jess is good")
        double EndTime = 71;
        double scale = 0.0011; 

        double SamplingStep = 25;
        double dt = 0.001;
        double RemeshingTime = 100;
        double FSI_Iterations =100;// 50;
        double EdgeLength =0.6*scale;
        
        /////////////////////////////////////////////////////////////////////////////////////
        std::string output_dir = "FSISimulations/VascularNetworkLargerBendingForce/Type1Collapse/H/MeshCorrection/";
        std::string Archieved = "FSISimulations/VascularNetworkLargerBendingForce/Type1Collapse/H";
        // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
         // this simulaiton is in mm. Have chosen this magnitude because um in m will give me numbers too close to machince presision, and movment
        // in um will be too large and break chaste without carefull playing with or a tiny time step

        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
 

        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetPrintRemeshedIC(1);

    static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).ExecuteHistoryDependentRemeshing();


        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
        p_simulator->RemoveAllForces();
        // p_simulator->RemoveAllCellPopulationBoundaryConditions();

        /* 
        -----------------------------
        Edit  RemeshingTriggerOnStepHeteroModifier
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
        p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); // 
        p_Mesh_modifier->SetStepSize(pow(10, -8));

        p_Mesh_modifier->TurnOffRemeshing();
        p_Mesh_modifier->SetmSetUpSolve(1);


        // Thirdcollapse option 
        // Upstream 
        c_vector<double, 3> UpperPlanePoint =  Create_c_vector( 0.148450771419833, 0.17540001846872066,0 ); 
        c_vector<double, 3> UpperPlaneNormal = Create_c_vector(0.6681809321405496, -0.7397789537602159, 0);  
        // Down stream
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0.16467321416038386, 0.1592262610683606,0 ); 
        c_vector<double, 3> LowerPlaneNormal = -Create_c_vector(0.870408755948827, -0.47734167113341636, 0 ); 
  
        p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);
        p_Mesh_modifier->SetRadius(0.015);
        // p_Mesh_modifier->TurnOffRemeshing();
        p_Mesh_modifier->SetUpdateFrequency(1/dt);


        boost::shared_ptr<EnclosedRegionBoundaryCondition<2, 3> > p_condition(new EnclosedRegionBoundaryCondition<2, 3>(&(p_simulator->rGetCellPopulation()) , UpperPlanePoint, UpperPlaneNormal, 0.01)); //0.01));

        p_condition->SetPointOnPlane2( LowerPlanePoint);
        p_condition->SetNormalToPlane2(-LowerPlaneNormal);
        p_simulator->AddCellPopulationBoundaryCondition(p_condition);
       
       /* 
        -----------------------------
        finished editing  RemeshingTriggerOnStepHeteroModifier
        ----------------------------
        */

    
        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_shear_force->SetCollapseType(1);
        p_simulator->AddForce(p_shear_force);

        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_membrane_force->SetMembraneStiffness(pow(10, -7));
        p_membrane_force->SetCollapseType(1);
        p_simulator->AddForce(p_membrane_force);


        /*
        -----------------------------
        Add the HemeLB Force
        ----------------------------
        */
       // Inlet1
        // Inlet1
       c_vector<double, 3> Point1 = Create_c_vector(0.012896438155977795,0.1296022543191808, -0.002236333573976076 );
       c_vector<double, 3> PlaneNormal1 = Create_c_vector(0.9491166468623324, -0.3145145356937025,0 );
        

       c_vector<double, 3> Point2 = Create_c_vector(0.04152531239022141, 0.23542737676454184, -0.0007601220860857502);
       c_vector<double, 3> PlaneNormal2 = Create_c_vector(0.7969315092041309, -0.604064075870372, 0 );
       

       c_vector<double, 3> Point3 = Create_c_vector(0.1651556668051053,0.02512909630785012,0.0011455782250521602);
       c_vector<double, 3> PlaneNormal3 = Create_c_vector(0.0807491790098989,0.9963585613635368,  0.027371285808498066 );

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg

        double InletPressure = P_blood; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later
        double OutletPressure = P_blood * (0.9);

        boost::shared_ptr<HemeLBForce> p_ForceOut(new HemeLBForce());
        p_ForceOut->Inlets(PlaneNormal1, Point1, OutletPressure, "Outlet"); // Issues here 
        p_ForceOut->Inlets(PlaneNormal2, Point2, OutletPressure, "Outlet"); //FIne
        p_ForceOut->Inlets(PlaneNormal3, Point3, InletPressure, "Inlet");// Issues here 
        p_ForceOut->SetStartTime(EndTime);
        p_ForceOut->Network("VascularNetwork");
        p_ForceOut->SetCollapseType(1);
        p_ForceOut->SetFluidSolidIterations(FSI_Iterations);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir+"HemeLBForce/", p_simulator->rGetCellPopulation(),1);
        p_simulator->AddForce(p_ForceOut);


      for (int i =1; i<=50; i++)
        { 
    
            EndTime +=1;
            p_simulator->SetEndTime(EndTime);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }

    }



};

#endif /*TESTRELAXATION_HPP_*/



