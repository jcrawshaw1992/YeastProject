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


  void TestHemeLBCollapse() throw(Exception)
   {

      
        std::string output_dir = "FSISimulations/Honey/Collapse1_StrongMembraneParameterVariationAdditionalInitialConditionCollapseMoreRemeshing2/ShorterTS2/";

        double AreaParameter = -5;  double DilationParameter = -5.5; double DeformationParamter = -5; double BendingParameter = -7;
        std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), pow(10, BendingParameter)) }, {0,  Create_c_vector(pow(10, -4), pow(10, -4), pow(10, -4),pow(10, BendingParameter))} };
                                
        std::string Archieved =  "FSISimulations/Honey/Collapse1_StrongMembraneParameterVariationAdditionalInitialConditionCollapseMoreRemeshing2/";

        
        double EndTime = 11.8;
        double SamplingStep = 50;
        double dt = 0.0001;
        double RemeshingTime = 500;//50;
        double EdgeLength =0.0005;
        double FSI_Iterations = 500;//50;


        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);

        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetPrintRemeshedIC(1);

        ////////

        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetRelativePath(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetBinningIntervals(15, 10, 1);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).EdgeLengthVariable(1);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingIterations(5);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetRemeshingSoftwear("CGAL");
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetOperatingSystem("server");

        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
        p_simulator->RemoveAllForces();
  

        /* 
        -----------------------------
        Edit  RemeshingTriggerOnStepHeteroModifier
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 
        p_Mesh_modifier->TurnOffRemeshing();   
        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);

        // First collapse option 
        // Upstream 
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0.034963365591332625, -5e-5,0);
        c_vector<double, 3> UpperPlaneNormal = Create_c_vector(1,0,0);
        // Down stream                                                             
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0.04307533991933138, -5e-5,0 );
        c_vector<double, 3> LowerPlaneNormal = -Create_c_vector(1,0,0);

        p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);

        p_Mesh_modifier->SetRadius(0.007);
        p_Mesh_modifier->SetStepSize(pow(10, -8));
        p_Mesh_modifier->SetUpdateFrequency(0.05/dt);
        p_Mesh_modifier->SetmSetUpSolve(1);


        boost::shared_ptr<EnclosedRegionBoundaryCondition<2, 3> > p_condition(new EnclosedRegionBoundaryCondition<2, 3>(&(p_simulator->rGetCellPopulation())  , UpperPlanePoint, UpperPlaneNormal, 0.01)); //0.01));

        p_condition->SetPointOnPlane2( LowerPlanePoint);
        p_condition->SetNormalToPlane2(-LowerPlaneNormal);
        p_simulator->AddCellPopulationBoundaryCondition(p_condition);
       
        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_simulator->AddForce(p_shear_force);

        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_membrane_force->SetMembraneStiffness(pow(10, -7));
        p_simulator->AddForce(p_membrane_force);


        /*
        -----------------------------
        Add the HemeLB Force
        ----------------------------
        */
      
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg

        double InletPressure = P_blood* (1.2); // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later
        double OutletPressure = P_blood * (0.8);

        


        // Inlet1
        c_vector<double, 3> Point1 = Create_c_vector(0.0022124205632440453, -0.01409714283809368,0);
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0 );
        // Inlet2
        c_vector<double, 3> Point2 = Create_c_vector(0.0022124205632440453, 0.0004613958360830862,0);
        c_vector<double, 3> PlaneNormal2 = Create_c_vector(1,0,0 );
        // Outlet1
        c_vector<double, 3> Point3 = Create_c_vector(0.07584026562645844, -0.01409714283809368,0);
        c_vector<double, 3> PlaneNormal3 = Create_c_vector(-1,0,0 );
        // Outlet2
        c_vector<double, 3> Point4 = Create_c_vector(0.07584026562645844, 0.0004613958360830862,0);
        c_vector<double, 3> PlaneNormal4 = Create_c_vector(-1,0,0 );


        boost::shared_ptr<HemeLBForce<2, 3> > p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure,  "Inlet");
        p_ForceOut->Inlets(PlaneNormal3, Point3, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal4, Point4, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(EndTime);
        p_ForceOut->SetFluidSolidIterations(FSI_Iterations);

        p_ForceOut->SetUpHemeLBConfiguration(output_dir+"HemeLBForce/", p_simulator->rGetCellPopulation(),0);
        p_simulator->AddForce(p_ForceOut);


      for (int i =1; i<=1600; i++)
        { 
    
            EndTime +=0.1;
            p_simulator->SetEndTime(EndTime);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }

    }
    /* */


void offTestWithConstantForce() throw(Exception)
   {

      
        std::string output_dir = "FSISimulations/Honey/Collapse1_StrongMembraneParameterVariationAdditionalInitialConditionCollapseMoreRemeshing2/ShorterTS/";

        double AreaParameter = -5;  double DilationParameter = -5.5; double DeformationParamter = -5; double BendingParameter = -7;
        std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), pow(10, BendingParameter)) }, {0,  Create_c_vector(pow(10, -4), pow(10, -4), pow(10, -4),pow(10, BendingParameter))} };
                                
        std::string Archieved = "FSISimulations/Honey/Collapse1_StrongMembraneParameterVariationAdditionalInitialConditionCollapseMoreRemeshing2/";
        
        double EndTime = 11.4;
        double SamplingStep = 100;
        double dt = 0.001;
        double RemeshingTime = 500;//50;
        double EdgeLength =0.0005;
        double FSI_Iterations = 500;//50;


        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
 

        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetPrintRemeshedIC(1);

        ////////

        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetRelativePath(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetBinningIntervals(15, 5, 1);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).EdgeLengthVariable(1.05);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingIterations(5);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetRemeshingSoftwear("CGAL");
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetOperatingSystem("server");

        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
        p_simulator->RemoveAllForces();
  

        /* 
        -----------------------------
        Edit  RemeshingTriggerOnStepHeteroModifier
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 
        p_Mesh_modifier->TurnOffRemeshing();   
        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
        p_Mesh_modifier->SetStepSize(pow(10, -8));

        // First collapse option 
        // Upstream 
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0.034963365591332625, -5e-5,0);
        c_vector<double, 3> UpperPlaneNormal = Create_c_vector(1,0,0);
        // Down stream                                                             
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0.04307533991933138, -5e-5,0 );
        c_vector<double, 3> LowerPlaneNormal = -Create_c_vector(1,0,0);

        p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);

        p_Mesh_modifier->SetRadius(0.007);
        p_Mesh_modifier->SetUpdateFrequency(0.05/dt);
        p_Mesh_modifier->SetmSetUpSolve(1);


        boost::shared_ptr<EnclosedRegionBoundaryCondition<2, 3> > p_condition(new EnclosedRegionBoundaryCondition<2, 3>(&(p_simulator->rGetCellPopulation())  , UpperPlanePoint, UpperPlaneNormal, 0.01)); //0.01));

        p_condition->SetPointOnPlane2( LowerPlanePoint);
        p_condition->SetNormalToPlane2(-LowerPlaneNormal);
        p_simulator->AddCellPopulationBoundaryCondition(p_condition);
       
        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_simulator->AddForce(p_shear_force);

        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_membrane_force->SetMembraneStiffness(pow(10, -7));
        p_simulator->AddForce(p_membrane_force);


        /*
        -----------------------------
        Add the HemeLB Force
        ----------------------------
        */
      
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg

        double InletPressure = P_blood* (1.2); // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later
        double OutletPressure = P_blood * (0.8);

        


        // Inlet1
        c_vector<double, 3> Point1 = Create_c_vector(0.0022124205632440453, -0.01409714283809368,0);
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(1,0,0 );
        // Inlet2
        c_vector<double, 3> Point2 = Create_c_vector(0.0022124205632440453, 0.0004613958360830862,0);
        c_vector<double, 3> PlaneNormal2 = Create_c_vector(1,0,0 );
        // Outlet1
        c_vector<double, 3> Point3 = Create_c_vector(0.07584026562645844, -0.01409714283809368,0);
        c_vector<double, 3> PlaneNormal3 = Create_c_vector(-1,0,0 );
        // Outlet2
        c_vector<double, 3> Point4 = Create_c_vector(0.07584026562645844, 0.0004613958360830862,0);
        c_vector<double, 3> PlaneNormal4 = Create_c_vector(-1,0,0 );


        boost::shared_ptr<HemeLBForce<2, 3> > p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure,  "Inlet");
        p_ForceOut->Inlets(PlaneNormal3, Point3, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal4, Point4, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(EndTime);
        p_ForceOut->SetFluidSolidIterations(FSI_Iterations);

        p_ForceOut->SetUpHemeLBConfiguration(output_dir+"HemeLBForce/", p_simulator->rGetCellPopulation());
        p_simulator->AddForce(p_ForceOut);


      for (int i =1; i<=50; i++)
        { 
    
            EndTime +=0.1;
            p_simulator->SetEndTime(EndTime);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }

    }
    /* */



};

#endif /*TESTRELAXATION_HPP_*/





        /*
        -----------------------------
        Add the HemeLB Force
        ----------------------------
        */
        // Inlet1
    //    c_vector<double, 3> Point1 = Create_c_vector(0.02297977879048401, 0.037075481876763385,  0.0004635334404919405);
    //     c_vector<double, 3> PlaneNormal1 = Create_c_vector(0.999910929986149, 0.0020912833516396026,-0.013181753607835568  );
    //     // Inlet2
    //     c_vector<double, 3> Point2 = Create_c_vector(0.036837793060229,0.055384952301980456,-0.0007597519518940717)  ;
    //     c_vector<double, 3> PlaneNormal2 = Create_c_vector(0.6732805342808156,-0.7380737547778258,  -0.04405059212657047);

    //     // Outlet1 This is the bug bar
    //     c_vector<double, 3> Point3 = Create_c_vector(0.051243051697026636,0.05386889597979771, 0.00016345376729906492 ) ;
    //     c_vector<double, 3> PlaneNormal3 = Create_c_vector( -0.8734008505817445, -0.4862907639633924, 0.026310588875685135  );

    //     // Outlet2
    //     c_vector<double, 3> Point4 = Create_c_vector( 0.05849786867183286,0.04003892834739773, -8.346812953121241e-5) ;
    //     c_vector<double, 3> PlaneNormal4 = Create_c_vector( -0.7771752850914773, -0.6286228836915505, 0.028841746519592842);

    //     // Inlet3
    //     c_vector<double, 3> Point5 = Create_c_vector(0.056296622912376706, 0.020105116941221777, 0.0002243854912119816) ;
    //     c_vector<double, 3> PlaneNormal5 = Create_c_vector(-0.7208815333497773, 0.6929575165100672, -0.011819271867331964);

    //     // Outlet3
    //     c_vector<double, 3> Point6 = Create_c_vector(0.039492636709086544, 0.013468141696478432, -0.00030703284641268327) ;
    //     c_vector<double, 3> PlaneNormal6 = Create_c_vector(0.14577215937619503, 0.9879552859046741, -0.0519117578571096);

    //      // Inlet4
    //     c_vector<double, 3> Point7 = Create_c_vector(0.027954386139121126, 0.02122180981896879, 0.0008357837352478671) ;
    //     c_vector<double, 3> PlaneNormal7 = Create_c_vector(0.6878807670924608, 0.7247343474980099, -0.03975142539484216);


    //     double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg

    //     double InletPressure = P_blood * 1.001; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later
    //     double OutletPressure = P_blood * (0.999);

    //     boost::shared_ptr<HemeLBForce<2, 3> > p_ForceOut(new HemeLBForce<2, 3>());
    //     p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
    //     p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure*1.1, "Inlet");
    //     p_ForceOut->Inlets(PlaneNormal3, Point3, InletPressure*1.2, "Inlet");
    //     p_ForceOut->Inlets(PlaneNormal4, Point4, InletPressure*1, "Inlet");
    //     p_ForceOut->Inlets(PlaneNormal5, Point5, OutletPressure, "Outlet");
    //     p_ForceOut->Inlets(PlaneNormal6, Point6, OutletPressure*0.98, "Outlet");
    //     p_ForceOut->Inlets(PlaneNormal7, Point7, OutletPressure*0.95, "Outlet");
    //     p_ForceOut->SetStartTime(EndTime);
    //     p_ForceOut->SetFluidSolidIterations(FSI_Iterations);
    //     p_ForceOut->SetUpHemeLBConfiguration(output_dir+"HemeLBForce/", simulator.rGetCellPopulation());
    //     simulator.AddForce(p_ForceOut);

       






        // ///////////////////////////////////////////////////////////////////////////////////////

        // std::vector<c_vector<double, 3> > boundary_plane_points1;
        // std::vector<c_vector<double, 3> > boundary_plane_normals1;

        // std::vector<c_vector<double, 3> > boundary_plane_points2;
        // std::vector<c_vector<double, 3> > boundary_plane_normals2;



        // boundary_plane_points1.push_back(Create_c_vector(0.0344839553141627, 0.039166859323769385,  -7.321768764989628e-5   ));
        // boundary_plane_normals1.push_back(Create_c_vector( 0.7474390621015472,  0.662857583536715,  0.04421167710714476  ));
        // /// ----------------------------------------
        // boundary_plane_points2.push_back(Create_c_vector(0.038633725599634766 , 0.042834650295782095, 0.00014677123886970464 ));
        // boundary_plane_normals2.push_back(Create_c_vector( 0.7713091680800281,0.6338501379897751,  -0.057586194578546857  ));
        // /* */
        // boundary_plane_points1.push_back(Create_c_vector(0.03408644744248119, 0.02948610276941591,1.1431840682180057e-5   ));
        // boundary_plane_normals1.push_back(Create_c_vector( -0.13710180943275943, 0.9904392690600551,0 ));
        // /// ----------------------------------------
        // boundary_plane_points2.push_back(Create_c_vector(0.031214753375716997,0.03384381942477147, -0.00020100804806112238  ));
        // boundary_plane_normals2.push_back(Create_c_vector( -0.1365186273602051,  0.9906313098812551, 0.00350317950355177 ));
        // /* */
        // boundary_plane_points1.push_back(Create_c_vector(0.05115507250002135,0.029210202823000585, 0.0013139290500819793  ));
        // boundary_plane_normals1.push_back(Create_c_vector( 0.25675931195415885, 0.9657125663824311, 0.03839133829245353 ));
        // /// ----------------------------------------
        // boundary_plane_points2.push_back(Create_c_vector(0.052237407320433346, 0.033906865036683775, 0.0012929371808254831 ));
        // boundary_plane_normals2.push_back(Create_c_vector( 0.23648919545312502, 0.9712857281725686, 0.026017199738188904 ));
        // /* */
        // boundary_plane_points1.push_back(Create_c_vector(0.04877608409834386,0.03986943093515008, 0.00376937928511186   ));
        // boundary_plane_normals1.push_back(Create_c_vector(-0.7187448933270177, 0.6950578067249691 ,   0.017332732820003426));
        // /// ----------------------------------------
        // boundary_plane_points2.push_back(Create_c_vector(0.04524453117254053, 0.043228289032675365,0.003720997732043512  ));
        // boundary_plane_normals2.push_back(Create_c_vector( -0.6624867475157061, 0.7438111690784874, 0.08863551274830604 ));
        // /* */
        // boundary_plane_points1.push_back(Create_c_vector(0.03663965409412324, 0.020155302012919925, -0.0009673119597267015     ));
        // boundary_plane_normals1.push_back(Create_c_vector( 0.894371990942059,  -0.44368042276915565,  0.05697740139494638 ));
        // /// ----------------------------------------
        // boundary_plane_points2.push_back(Create_c_vector(0.037791247760704444, 0.01958401895618103, -0.0008939478391491307));
        // boundary_plane_normals2.push_back(Create_c_vector( 0.894371990942059,  -0.44368042276915565,  0.05697740139494638 ));
        // /* */
        // boundary_plane_points1.push_back(Create_c_vector(0.04535333451092035,  0.01827253994971488, -0.000826775696576238   ));
        // boundary_plane_normals1.push_back(Create_c_vector( 0.8800927838947658,  0.47304154554359323,  0.040845904642813985   ));
        // /// ----------------------------------------
        // boundary_plane_points2.push_back(Create_c_vector(0.04777378462420049,  0.02435646728406592, -0.0015296906390526728  ));
        // boundary_plane_normals2.push_back(Create_c_vector(0.8999000530240061,0.43609399958642875, -0.001384952021894634    ));


        // /* */
        // boundary_plane_points1.push_back(Create_c_vector(0.03308530881950349, 0.02149761438499426, -0.00018009428504821387));
        // boundary_plane_normals1.push_back(Create_c_vector(0.9261098741471802 ,  -0.376954441509806, -0.015028307726570638 ));
        // /// ----------------------------------------
        // boundary_plane_points2.push_back(Create_c_vector(0.034300824847584066,0.021164260640778476, -0.0002909367534368054 ));
        // boundary_plane_normals2.push_back(Create_c_vector(0.9691730622969827, -0.2311539334429166, -0.08527270590145869 ));

       


        // unsigned counter =0;
        // for (unsigned boundary_id = 0; boundary_id < boundary_plane_points1.size()-1; boundary_id++)
        // {
        //     counter +=1;
        //     boost::shared_ptr<EnclosedRegionBoundaryCondition<2, 3> > p_condition(new EnclosedRegionBoundaryCondition<2, 3>(&(simulator.rGetCellPopulation()) , boundary_plane_points1[boundary_id], boundary_plane_normals1[boundary_id], 0.01)); //0.01));

        //     p_condition->SetPointOnPlane2( boundary_plane_points2[boundary_id]);
        //     p_condition->SetNormalToPlane2(boundary_plane_normals2[boundary_id]);


        //      simulator.AddCellPopulationBoundaryCondition(p_condition);
        // }

        // boost::shared_ptr<EnclosedRegionBoundaryCondition<2, 3> > p_condition1(new EnclosedRegionBoundaryCondition<2, 3>(&(simulator.rGetCellPopulation()) , boundary_plane_points1[counter], boundary_plane_normals1[counter], 0.004));

        // p_condition1->SetPointOnPlane2( boundary_plane_points2[counter]);
        // p_condition1->SetNormalToPlane2(boundary_plane_normals2[counter]);

        // simulator.AddCellPopulationBoundaryCondition(p_condition1);

        ///////////////////////////////////////////////////////////////////////////////////////



