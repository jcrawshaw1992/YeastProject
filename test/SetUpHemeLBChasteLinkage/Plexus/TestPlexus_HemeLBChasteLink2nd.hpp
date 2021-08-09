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

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:



  void TestContinuingHomoArchieve2collapse() throw(Exception)
   {
        double DilationParameter = -7;
        double AreaParameter = -7.9;
        double DeformationParamter = -6.1;


        TRACE("Jess is good")
        double EndTime = 4;
        double scale = 0.00006684491/1.29;
        double FSI_Iterations = 3000;

        double SamplingStep = 100;
        double dt = 0.0001;
        double RemeshingTime = 1000;
        double EdgeLength =0.00045;
        
        std::string Archieved = "DeformingPlexus/FlatForceFINAL9/";
        // std::string output_dir = "DeformingPlexus/Grow2Equi/";
        std::string output_dir = "DeformingPlexus_HemeLB/SecondCollapse/";
     
    
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
 

        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetPrintRemeshedIC(1);


        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
        p_simulator->RemoveAllForces();
        p_simulator->RemoveAllCellPopulationBoundaryConditions();

        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_simulator->AddForce(p_shear_force);

        /*
        -----------------------------
        Add the HemeLB Force
        ----------------------------
        */
       // Inlet1
        c_vector<double, 3> Point1 = Create_c_vector(0.02297977879048401, 0.037075481876763385,  0.0004635334404919405);
        c_vector<double, 3> PlaneNormal1 = Create_c_vector(0.999910929986149, 0.0020912833516396026,-0.013181753607835568  );
        // Inlet2
        c_vector<double, 3> Point2 = Create_c_vector(0.036837793060229,0.055384952301980456,-0.0007597519518940717)  ;
        c_vector<double, 3> PlaneNormal2 = Create_c_vector(0.6732805342808156,-0.7380737547778258,  -0.04405059212657047);

        // Outlet1
        c_vector<double, 3> Point3 = Create_c_vector(0.051243051697026636,0.05386889597979771, 0.00016345376729906492 ) ;
        c_vector<double, 3> PlaneNormal3 = Create_c_vector( -0.8734008505817445, -0.4862907639633924, 0.026310588875685135  );

        // Outlet2
        c_vector<double, 3> Point4 = Create_c_vector( 0.05849786867183286,0.04003892834739773, -8.346812953121241e-5) ;
        c_vector<double, 3> PlaneNormal4 = Create_c_vector( -0.7771752850914773, -0.6286228836915505, 0.028841746519592842);

        // Inlet3
        c_vector<double, 3> Point5 = Create_c_vector(0.056296622912376706, 0.020105116941221777, 0.0002243854912119816) ;
        c_vector<double, 3> PlaneNormal5 = Create_c_vector(-0.7208815333497773, 0.6929575165100672, -0.011819271867331964);

        // Outlet3
        c_vector<double, 3> Point6 = Create_c_vector(0.039492636709086544, 0.013468141696478432, -0.00030703284641268327) ;
        c_vector<double, 3> PlaneNormal6 = Create_c_vector(0.14577215937619503, 0.9879552859046741, -0.0519117578571096);

         // Inlet4
        c_vector<double, 3> Point7 = Create_c_vector(0.027954386139121126, 0.02122180981896879, 0.0008357837352478671) ;
        c_vector<double, 3> PlaneNormal7 = Create_c_vector(0.6878807670924608, 0.7247343474980099, -0.03975142539484216);



        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg

        double InletPressure = P_blood * 1.001; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later
        double OutletPressure = P_blood * (0.999);

        boost::shared_ptr<HemeLBForce<2, 3> > p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, InletPressure*1.1, "Inlet");
        p_ForceOut->Inlets(PlaneNormal3, Point3, OutletPressure, "Outlet");
        p_ForceOut->Inlets(PlaneNormal4, Point4, InletPressure*1.05, "Inlet");
        p_ForceOut->Inlets(PlaneNormal5, Point5, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal6, Point6, OutletPressure*1.01, "Outlet");
        p_ForceOut->Inlets(PlaneNormal7, Point7, OutletPressure*0.95, "Outlet");
        p_ForceOut->SetStartTime(EndTime);
        p_ForceOut->SetFluidSolidIterations(FSI_Iterations);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir+"HemeLBForce/", p_simulator->rGetCellPopulation());
        p_simulator->AddForce(p_ForceOut);

       

        /* 
        -----------------------------
        Edit  RemeshingTriggerOnStepHeteroModifier
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        // p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 
        p_Mesh_modifier->TurnOffRemeshing();   
        //AreaConstant           AreaDilationModulus        ShearModulus
        std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), 0) },
                                                                    {0,  Create_c_vector(pow(10, -5), pow(10, -4), pow(10, -4), 0)}    };
 

        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
        p_Mesh_modifier->SetmSetUpSolve(1);
        p_Mesh_modifier->SetStepSize(pow(10, -8));

        // Second collapse option 
        // Upstream 
        c_vector<double, 3> UpperPlanePoint =  Create_c_vector(0.030848605199273658, 0.034704038458335855, 0.001356825613262895);//Create_c_vector(0.03889835415936754,0.023284157207353464, 0.002357881097576096 );
        c_vector<double, 3> UpperPlaneNormal = Create_c_vector(0.18206343627104957,  -0.974167748852415, 0.13360427489042281); // Create_c_vector(-0.9111321813105736, 0.41098845246144305, -0.030440764175433815);
        // Down stream
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0.03209735497750507, 0.028022346707068654, 0.002273200213552209);// Create_c_vector(0.03473129104301778, 0.023894561960013404,0.0007601136405928878 );
        c_vector<double, 3> LowerPlaneNormal = -Create_c_vector(0.18206343627104957, -0.974167748852415, 0.13360427489042281);;//-Create_c_vector(-0.850657197084237, 0.5256842037978231, 0.006200881085641978);
        p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);
        p_Mesh_modifier->SetRadius(0.01);
        p_Mesh_modifier->SetUpdateFrequency(10/dt);

          /*
        -----------------------------
        Boundary conditions
        ----------------------------
        */

        std::vector<c_vector<double, 3> > boundary_plane_points;
        std::vector<c_vector<double, 3> > boundary_plane_normals;



        boundary_plane_points.push_back(Create_c_vector( 0.023667025594574218, 0.03580448541106122, 1.0883377925970793e-6));
        boundary_plane_normals.push_back(-Create_c_vector( 0.9992096159136262, -0.037477299160519,-0.01325124572923273 ));

        boundary_plane_points.push_back(Create_c_vector( 0.02842355040448236, 0.04796528990124343, -0.0012984792809731791 ));
        boundary_plane_normals.push_back(Create_c_vector( -0.6540763785122455,  0.7518206257291297, 0.08336568718943009 ));


        boundary_plane_points.push_back(Create_c_vector(0.06117298902054698,0.03654224598923696,0.0009160147066410285 ));
        boundary_plane_normals.push_back(Create_c_vector(0.7792901394890018, 0.6252088966438896,  -0.04266982601959587  ));

        boundary_plane_points.push_back(Create_c_vector(0.055625630272326206, 0.0200783962665125, -0.000722688871900608 ));
        boundary_plane_normals.push_back(Create_c_vector( 0.7212385318941745, -0.6926132859491577, 0.010090403254884016  ));


        boundary_plane_points.push_back(Create_c_vector( 0.026287008430771083, 0.023675441622844417, -0.0007261644860741899 ));
        boundary_plane_normals.push_back(Create_c_vector( -0.6727411476463672, -0.7373576645038901,  0.061016578573512954 ));


        boundary_plane_points.push_back(Create_c_vector( 0.03952667347394293,0.01380981593118016, 0.00035914153313716104  ));
        boundary_plane_normals.push_back(Create_c_vector( -0.14188750875173126, -0.9842863935116474,0.10511056278063648));


        boundary_plane_points.push_back(Create_c_vector( 0.051916367697719554, 0.05396743908086633, 0.0004919967151319857 ));
        boundary_plane_normals.push_back(Create_c_vector(0.8831316656538849, 0.4690231691065414, 0.009784066672481936 ));


        boundary_plane_points.push_back(Create_c_vector(0.036804371668431334, 0.053919748549890005, -0.0007634162035095087 ));
        boundary_plane_normals.push_back(Create_c_vector(-0.6416261436122902, 0.7612012690416773,  0.09427894697418274   ));

        boundary_plane_points.push_back(Create_c_vector(0.03751768498473936,0.050175026722357664, -0.002624281078237435 ));
        boundary_plane_normals.push_back(Create_c_vector(-0.6625540165601127,  0.7473601137814084, 0.04974972832792599 ));

        boundary_plane_points.push_back(Create_c_vector(0.026864281281697915, 0.037668141898196575, 0.00018974398398148207 ));
        boundary_plane_normals.push_back(Create_c_vector(-0.9998551753153383, -0.015809159893423953, -0.006300703024145976 ));

        boundary_plane_points.push_back(Create_c_vector(0.040547082072423954,0.018275901698485374, 0.00039994117540888903 ));
        boundary_plane_normals.push_back(Create_c_vector(-0.13760377878446087, -0.9886290666284141, 0.060644609667043325));
  
        boundary_plane_points.push_back(Create_c_vector(0.04660359474772499, 0.05104084087227261,0.0007679575729855695  ));
        boundary_plane_normals.push_back(Create_c_vector(0.8871018711566134,0.4562209461172769, -0.07009078765638171));


        boundary_plane_points.push_back(Create_c_vector( 0.029862726176558368, 0.022150802525023005,  0.0007074621356822784   ));
        boundary_plane_normals.push_back(Create_c_vector(-0.6983037235236903,  -0.7144025416789002 ,0.044731623664661664  ));


        boundary_plane_points.push_back(Create_c_vector(  0.053216101460717424,0.022875560649274096,  0.00035741534715402995   ));
        boundary_plane_normals.push_back(Create_c_vector(0.7104538656628817, -0.703732491551614,  0.003985611524662334   ));

        boundary_plane_points.push_back(Create_c_vector(0.05715768054001756,  0.03813699788011809,  0.0006839495043087638));
        boundary_plane_normals.push_back(Create_c_vector(0.7640527331012904, 0.6448510930472797, -0.01976079037329169 ));


     
        for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
      
             boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(&(p_simulator->rGetCellPopulation()) , boundary_plane_points[boundary_id], boundary_plane_normals[boundary_id], 0.008));
             p_simulator->AddCellPopulationBoundaryCondition(p_condition);
        }


        ///////////////////////////////////////////////////////////////////////////////////////

        std::vector<c_vector<double, 3> > boundary_plane_points1;
        std::vector<c_vector<double, 3> > boundary_plane_normals1;

        std::vector<c_vector<double, 3> > boundary_plane_points2;
        std::vector<c_vector<double, 3> > boundary_plane_normals2;



        boundary_plane_points1.push_back(Create_c_vector(0.0344839553141627, 0.039166859323769385,  -7.321768764989628e-5   ));
        boundary_plane_normals1.push_back(Create_c_vector( 0.7474390621015472,  0.662857583536715,  0.04421167710714476  ));
        /// ----------------------------------------
        boundary_plane_points2.push_back(Create_c_vector(0.038633725599634766 , 0.042834650295782095, 0.00014677123886970464 ));
        boundary_plane_normals2.push_back(Create_c_vector( 0.7713091680800281,0.6338501379897751,  -0.057586194578546857  ));
        /* */
        boundary_plane_points1.push_back(Create_c_vector(0.03408644744248119, 0.02948610276941591,1.1431840682180057e-5   ));
        boundary_plane_normals1.push_back(Create_c_vector( -0.13710180943275943, 0.9904392690600551,0 ));
        /// ----------------------------------------
        boundary_plane_points2.push_back(Create_c_vector(0.031214753375716997,0.03384381942477147, -0.00020100804806112238  ));
        boundary_plane_normals2.push_back(Create_c_vector( -0.1365186273602051,  0.9906313098812551, 0.00350317950355177 ));
        /* */
        boundary_plane_points1.push_back(Create_c_vector(0.05115507250002135,0.029210202823000585, 0.0013139290500819793  ));
        boundary_plane_normals1.push_back(Create_c_vector( 0.25675931195415885, 0.9657125663824311, 0.03839133829245353 ));
        /// ----------------------------------------
        boundary_plane_points2.push_back(Create_c_vector(0.052237407320433346, 0.033906865036683775, 0.0012929371808254831 ));
        boundary_plane_normals2.push_back(Create_c_vector( 0.23648919545312502, 0.9712857281725686, 0.026017199738188904 ));
        /* */
        boundary_plane_points1.push_back(Create_c_vector(0.04877608409834386,0.03986943093515008, 0.00376937928511186   ));
        boundary_plane_normals1.push_back(Create_c_vector(-0.7187448933270177, 0.6950578067249691 ,   0.017332732820003426));
        /// ----------------------------------------
        boundary_plane_points2.push_back(Create_c_vector(0.04524453117254053, 0.043228289032675365,0.003720997732043512  ));
        boundary_plane_normals2.push_back(Create_c_vector( -0.6624867475157061, 0.7438111690784874, 0.08863551274830604 ));
        /* */
        boundary_plane_points1.push_back(Create_c_vector(0.03663965409412324, 0.020155302012919925, -0.0009673119597267015     ));
        boundary_plane_normals1.push_back(Create_c_vector( 0.894371990942059,  -0.44368042276915565,  0.05697740139494638 ));
        /// ----------------------------------------
        boundary_plane_points2.push_back(Create_c_vector(0.037791247760704444, 0.01958401895618103, -0.0008939478391491307));
        boundary_plane_normals2.push_back(Create_c_vector( 0.894371990942059,  -0.44368042276915565,  0.05697740139494638 ));
        /* */
        boundary_plane_points1.push_back(Create_c_vector(0.04535333451092035,  0.01827253994971488, -0.000826775696576238   ));
        boundary_plane_normals1.push_back(Create_c_vector( 0.8800927838947658,  0.47304154554359323,  0.040845904642813985   ));
        /// ----------------------------------------
        boundary_plane_points2.push_back(Create_c_vector(0.04777378462420049,  0.02435646728406592, -0.0015296906390526728  ));
        boundary_plane_normals2.push_back(Create_c_vector(0.8999000530240061,0.43609399958642875, -0.001384952021894634    ));


        /* */
        boundary_plane_points1.push_back(Create_c_vector(0.03308530881950349, 0.02149761438499426, -0.00018009428504821387));
        boundary_plane_normals1.push_back(Create_c_vector(0.9261098741471802 ,  -0.376954441509806, -0.015028307726570638 ));
        /// ----------------------------------------
        boundary_plane_points2.push_back(Create_c_vector(0.034300824847584066,0.021164260640778476, -0.0002909367534368054 ));
        boundary_plane_normals2.push_back(Create_c_vector(0.9691730622969827, -0.2311539334429166, -0.08527270590145869 ));

       


        unsigned counter =0;
        for (unsigned boundary_id = 0; boundary_id < boundary_plane_points1.size()-1; boundary_id++)
        {
            counter +=1;
            boost::shared_ptr<EnclosedRegionBoundaryCondition<2, 3> > p_condition(new EnclosedRegionBoundaryCondition<2, 3>(&(p_simulator->rGetCellPopulation()) , boundary_plane_points1[boundary_id], boundary_plane_normals1[boundary_id], 0.01)); //0.01));

            p_condition->SetPointOnPlane2( boundary_plane_points2[boundary_id]);
            p_condition->SetNormalToPlane2(boundary_plane_normals2[boundary_id]);


             p_simulator->AddCellPopulationBoundaryCondition(p_condition);
        }

        boost::shared_ptr<EnclosedRegionBoundaryCondition<2, 3> > p_condition1(new EnclosedRegionBoundaryCondition<2, 3>(&(p_simulator->rGetCellPopulation()) , boundary_plane_points1[counter], boundary_plane_normals1[counter], 0.004));

        p_condition1->SetPointOnPlane2( boundary_plane_points2[counter]);
        p_condition1->SetNormalToPlane2(boundary_plane_normals2[counter]);

        p_simulator->AddCellPopulationBoundaryCondition(p_condition1);

        ///////////////////////////////////////////////////////////////////////////////////////



     for (int i =1; i<=10; i++)
        { 
            


            EndTime +=1;
            p_simulator->SetEndTime(EndTime);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);

        }

        
            // p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 
 
            // p_ForceOut->SetFluidSolidIterations(FSIIterations);
            p_simulator->SetSamplingTimestepMultiple(SamplingStep);
            p_simulator->SetDt(dt);
            p_Mesh_modifier->SetUpdateFrequency(2/dt);
    }


//////////////////////////////////////////////////////////////////////////////////////////////////////////
















    void offTestFSICylinder_HeteroContinued() throw(Exception)
    {
        std::string output_dir = "FSICylinder/Medium/Hetro8";
        std::string Archieve = "FSICylinder/Medium/Hetro8";
    
        double SamplingTimestepMultiple = 50000;
        double EndTime = 60;
        double scale = 1e3;
        double Length = 50e-6 * scale;
        double Radius = 0.5e-6 * scale;
        double dt = 0.00005;
        double FSIIterations = 200000;


        // Load and fix any settings in the simulator
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieve, EndTime);

        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetWriteVtkAsPoints(false);
   
        p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
        p_simulator->RemoveAllForces();

        /*
        -----------------------------
        Add the HemeLB Force
        ----------------------------
        */


        c_vector<double, 3> PlaneNormal1 = Create_c_vector(0, 0, 1);
        c_vector<double, 3> Point1 = Create_c_vector(0, 0, 0.0002e-6 * scale);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(0, 0, -1);
        c_vector<double, 3> Point2 = Create_c_vector(0, 0, Length - 0.0002e-6 * scale);

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = (0.002133152) * 1.001; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later
        double OutletPressure = (0.002133152) * (0.999);

        boost::shared_ptr<HemeLBForce<2, 3> > p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(EndTime);
        p_ForceOut->SetFluidSolidIterations(FSIIterations);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir+"/HemeLBForce/", p_simulator->rGetCellPopulation());
        p_simulator->AddForce(p_ForceOut);


        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_simulator->AddForce(p_shear_force);


        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_simulator->AddForce(p_membrane_force);
 
        /* 
        -----------------------------
        Edit  RemeshingTriggerOnStepHeteroModifier
        ----------------------------
         */

        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2,3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);

        
        std::map<double, c_vector<long double, 4> >  GrowthMaps =  { {1, Create_c_vector(pow(10, -8), pow(10, -9), pow(10, -8.05), 1e-9 ) },
                                                                     {0,  Create_c_vector(pow(10, -5), pow(10, -4), pow(10, -5), 1e-9 )}    };
 

        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);

        // p_Mesh_modifier->SetStepSize(pow(10, -8));

        // Upstream 
        // c_vector<double, 3> UpperPlanePoint = Create_c_vector(0,0,20e-6* scale);
        // c_vector<double, 3> UpperPlaneNormal = Create_c_vector(0,0,1);
        // // Down stream
        // c_vector<double, 3> LowerPlanePoint = Create_c_vector(0,0,30e-6 * scale);
        // c_vector<double, 3> LowerPlaneNormal = Create_c_vector(0,0,-1);
        // p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);
        p_Mesh_modifier->SetUpdateFrequency(1/dt);
        p_Mesh_modifier->SetmSetUpSolve(0);

       
        for (int j =1; j<=40; j++)
        {

            // for (int i =1; i<=1; i++)
            // { 
                // static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
                EndTime +=10;
                p_simulator->SetEndTime(EndTime);

                p_simulator->Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
            // }

            dt/=2 ;  SamplingTimestepMultiple*= 2; 
            FSIIterations*=2;
 
            p_ForceOut->SetFluidSolidIterations(FSIIterations);
            p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
            p_simulator->SetDt(dt);
            p_Mesh_modifier->SetUpdateFrequency(1/dt);
        }


       
    }

  
    
    
    void TestFSICylinder_HeteroContinued() throw(Exception)
    {
        std::string output_dir = "FSICylinder/Medium/Hetro10";
        std::string Archieve =  "FSICylinder/Medium/Hetro9";
    
        double SamplingTimestepMultiple = 500;
        double EndTime = 180;
        double scale = 1e3;
        double Length = 50e-6 * scale;
        double Radius = 0.5e-6 * scale;
        double dt = 0.0001;
        double FSIIterations = 2000;


        // Load and fix any settings in the simulator
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieve, EndTime);

        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

        /* Remove the constant pressure force   */
        p_simulator->RemoveAllForces();
        // // p_simulator->RemoveAllCellPopulationBoundaryConditions();
        
        p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);

        /*
        -----------------------------
        Add the HemeLB Force
        ----------------------------
        */


        c_vector<double, 3> PlaneNormal1 = Create_c_vector(0, 0, 1);
        c_vector<double, 3> Point1 = Create_c_vector(0, 0, 0.0002e-6 * scale);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(0, 0, -1);
        c_vector<double, 3> Point2 = Create_c_vector(0, 0, Length - 0.0002e-6 * scale);

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = (0.002133152) * 1.001; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later
        double OutletPressure = (0.002133152) * (0.999);

        boost::shared_ptr<HemeLBForce<2, 3> > p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(EndTime);
        p_ForceOut->SetFluidSolidIterations(FSIIterations);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir+"/HemeLBForce/", p_simulator->rGetCellPopulation());
        p_simulator->AddForce(p_ForceOut);


        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_simulator->AddForce(p_shear_force);


        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_simulator->AddForce(p_membrane_force);
 
        /* 
        -----------------------------
        Edit  RemeshingTriggerOnStepHeteroModifier
        ----------------------------
         */

        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2,3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);

        
        std::map<double, c_vector<long double, 4> >  GrowthMaps =  { {1, Create_c_vector(pow(10, -8), pow(10, -9), pow(10, -8.05), 1e-9 ) },
                                                                     {0,  Create_c_vector(pow(10, -5), pow(10, -4), pow(10, -5), 1e-9 )}    };
 

        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);

        // p_Mesh_modifier->SetStepSize(pow(10, -9));

        // Upstream 
        // c_vector<double, 3> UpperPlanePoint = Create_c_vector(0,0,20e-6* scale);
        // c_vector<double, 3> UpperPlaneNormal = Create_c_vector(0,0,1);
        // // Down stream
        // c_vector<double, 3> LowerPlanePoint = Create_c_vector(0,0,30e-6 * scale);
        // c_vector<double, 3> LowerPlaneNormal = Create_c_vector(0,0,-1);
        // p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);
        p_Mesh_modifier->SetUpdateFrequency(1/dt);
        p_Mesh_modifier->SetmSetUpSolve(0);
       
        for (int j =1; j<=40; j++)
        {

            for (int i =1; i<=5; i++)
            { 
                // static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
                EndTime +=10;
                p_simulator->SetEndTime(EndTime);

                p_simulator->Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
            }

            dt/=2 ;  SamplingTimestepMultiple*= 2; 
            FSIIterations*=2;
 
            // p_ForceOut->SetFluidSolidIterations(FSIIterations);
            p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
            p_simulator->SetDt(dt);
            p_Mesh_modifier->SetUpdateFrequency(2/dt);
        }


       
    }

  

    void offTestFSICylinder_Hetero() throw(Exception)
    {
        std::string output_dir = "FSICylinder/Medium/Hetro4";
        std::string Archieve = "FSICylinder/Medium";
    
        double SamplingTimestepMultiple = 100;
        double EndTime = 20;
        double scale = 1e3;
        double Length = 50e-6 * scale;
        double Radius = 0.5e-6 * scale;
        double dt = 0.005;
        double FSIIterations = 200;


        // Load and fix any settings in the simulator
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieve, EndTime);

        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(EndTime);

        /* Remove the constant pressure force   */
        p_simulator->RemoveAllForces();
        p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);

        /*
        -----------------------------
        Add the HemeLB Force
        ----------------------------
        */


        c_vector<double, 3> PlaneNormal1 = Create_c_vector(0, 0, 1);
        c_vector<double, 3> Point1 = Create_c_vector(0, 0, 0.0002e-6 * scale);

        c_vector<double, 3> PlaneNormal2 = Create_c_vector(0, 0, -1);
        c_vector<double, 3> Point2 = Create_c_vector(0, 0, Length - 0.0002e-6 * scale);

        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg

        double InletPressure = (0.002133152 - 0.001466542) * 1.001; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later
        double OutletPressure = (0.002133152 - 0.001466542) * (0.999);

        boost::shared_ptr<HemeLBForce<2, 3> > p_ForceOut(new HemeLBForce<2, 3>());
        p_ForceOut->Inlets(PlaneNormal1, Point1, InletPressure, "Inlet");
        p_ForceOut->Inlets(PlaneNormal2, Point2, OutletPressure, "Outlet");
        p_ForceOut->SetStartTime(EndTime);
        p_ForceOut->SetFluidSolidIterations(FSIIterations);
        p_ForceOut->SetUpHemeLBConfiguration(output_dir+"/HemeLBForce/", p_simulator->rGetCellPopulation());
        p_simulator->AddForce(p_ForceOut);

        // boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        // p_ForceOut->SetPressure(P_blood - P_tissue);
        // // p_ForceOut->SetRadiusThreshold(10 * Radius);
        // p_simulator->AddForce(p_ForceOut);


        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_simulator->AddForce(p_shear_force);


        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_simulator->AddForce(p_membrane_force);
 
        /* 
        -----------------------------
        Edit  RemeshingTriggerOnStepHeteroModifier
        ----------------------------
         */

        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2,3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);

        
        std::map<double, c_vector<long double, 4> >  GrowthMaps =  { {1, Create_c_vector(pow(10, -8), pow(10, -9), pow(10, -8.05), 1e-10 ) },
                                                                     {0,  Create_c_vector(pow(10, -5), pow(10, -4), pow(10, -5), 1e-10 )}    };
 

        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);

        p_Mesh_modifier->SetStepSize(pow(10, -8));

        // Upstream 
        c_vector<double, 3> UpperPlanePoint = Create_c_vector(0,0,20e-6* scale);
        c_vector<double, 3> UpperPlaneNormal = Create_c_vector(0,0,1);
        // Down stream
        c_vector<double, 3> LowerPlanePoint = Create_c_vector(0,0,30e-6 * scale);
        c_vector<double, 3> LowerPlaneNormal = Create_c_vector(0,0,-1);
        p_Mesh_modifier->Boundaries( UpperPlaneNormal,  UpperPlanePoint,  LowerPlaneNormal,  LowerPlanePoint);
        p_Mesh_modifier->SetUpdateFrequency(2/dt);
        p_Mesh_modifier->SetmSetUpSolve(1);

       
        for (int j =1; j<=40; j++)
        {

            
            for (int i =1; i<=10; i++)
            { 
                // static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetStartTime(NewEndTime);
                EndTime +=5;
                p_simulator->SetEndTime(EndTime);

                p_simulator->Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
            }

            dt/=2 ;  SamplingTimestepMultiple*= 2; 
            FSIIterations*=2;
 
            p_ForceOut->SetFluidSolidIterations(FSIIterations);
            p_simulator->SetSamplingTimestepMultiple(SamplingTimestepMultiple);
            p_simulator->SetDt(dt);
            p_Mesh_modifier->SetUpdateFrequency(2/dt);


        }


       
    }


  
    
 


};

#endif /*TESTRELAXATION_HPP_*/