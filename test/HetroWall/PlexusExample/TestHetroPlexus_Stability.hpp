#ifndef TESTRELAXATION_HPP_
#define TESTRELAXATION_HPP_

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "PetscSetupAndFinalize.hpp"

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"

#include "Debug.hpp"
#include "SmartPointers.hpp"
// #include "VtkMeshReader.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"

#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "HistoryDepMutableMesh.hpp"

#include "FixedRegionBoundaryCondition.hpp"
#include "MembraneBendingForce.hpp"
#include "MembraneDeformationForce.hpp"
#include "OutwardsPressure.hpp"

#include "RemeshingTriggerOnStepHeteroModifier.hpp"

#include "MembraneBendingForce0TargetAngle.hpp"
class TestRemeshing : public AbstractCellBasedTestSuite
{
public:



   void offTestContinuingHomoArchieve1() throw(Exception)
   {

        TRACE("Jess is good")
        double EndTime = 4;
        double scale = 0.00006684491 / 1.29;

        double SamplingStep = 50;
        double dt = 0.001;
        double RemeshingTime = 500;
        double EdgeLength =0.0004;
        

        std::string Archieved =  "DeformingPlexus/Grow2Equi/";
        std::string output_dir = "DeformingPlexus/Grow2Equi/";
     
    
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);

        // static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).ExecuteHistoryDependentRemeshing();

        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
    
 
        /* 
        -----------------------------
        Edit  RemeshingTriggerOnStepHeteroModifier
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        // p_Mesh_modifier->SetRemeshingInterval(RemeshingTime);                                                              
 
            
        for (int i =1; i<=1; i++)
        { 
            EndTime +=0.5;
            p_simulator->SetEndTime(EndTime);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }

            // dt/=2 ;  SamplingStep*= 2; 
            // // FSIIterations*=2;
            // p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 
 
            // // p_ForceOut->SetFluidSolidIterations(FSIIterations);
            // p_simulator->SetSamplingTimestepMultiple(SamplingStep);
            // p_simulator->SetDt(dt);
            // p_Mesh_modifier->SetUpdateFrequency(2/dt);
    }




  void TestContinuingHomoArchieve2() throw(Exception)
   {

        TRACE("Jess is good")
        double EndTime = 4;
        double scale = 0.00006684491 / 1.29;

        double SamplingStep = 50;
        double dt = 0.001;
        double RemeshingTime = 500;
        double EdgeLength =0.0004;
        

        std::string Archieved = "DeformingPlexus/Grow2Equi/";
        // std::string Archieved = "DeformingPlexus/FlatForceFINAL5/";
     
        // std::string output_dir = "DeformingPlexus/Grow2Equi/";
        std::string output_dir = "DeformingPlexus/TestParameters/";
     
    
        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);

        // static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).ExecuteHistoryDependentRemeshing();

        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);
        p_simulator->RemoveAllForces();
        p_simulator->RemoveAllCellPopulationBoundaryConditions();


        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_membrane_force->SetMembraneStiffness(pow(10, -9));
        p_simulator->AddForce(p_membrane_force);
        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_simulator->AddForce(p_shear_force);

        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure((P_blood - P_tissue));
        p_simulator->AddForce(p_ForceOut);


        /* 
        -----------------------------
        Edit  RemeshingTriggerOnStepHeteroModifier
        ----------------------------
        */

        double DilationParameter = -7;
        double AreaParameter = -7;
        double DeformationParamter = -8;


        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
        // p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 
        p_Mesh_modifier->TurnOffRemeshing();   
        //AreaConstant           AreaDilationModulus        ShearModulus
        std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), 1e-10) } };

        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
        p_Mesh_modifier->SetmSetUpSolve(1);


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

        boundary_plane_points.push_back(Create_c_vector(0.051817210775777105,0.05387570153154082, -0.00019851668457055458 ));
        boundary_plane_normals.push_back(Create_c_vector(0.8769643127867986, 0.4790312122440197,  0.038245153868211076 ));

        boundary_plane_points.push_back(Create_c_vector(0.06117298902054698,0.03654224598923696,0.0009160147066410285 ));
        boundary_plane_normals.push_back(Create_c_vector(0.7792901394890018, 0.6252088966438896,  -0.04266982601959587  ));

        boundary_plane_points.push_back(Create_c_vector(0.055625630272326206, 0.0200783962665125, -0.000722688871900608 ));
        boundary_plane_normals.push_back(Create_c_vector( 0.7212385318941745, -0.6926132859491577, 0.010090403254884016  ));

        boundary_plane_points.push_back(Create_c_vector(0.05266362792368705, 0.011971267579691799,-0.0014918860363199145   ));
        boundary_plane_normals.push_back(Create_c_vector( -0.1352823453802107, -0.9855392137170218, 0.10203501974549702  ));


        boundary_plane_points.push_back(Create_c_vector( 0.026287008430771083, 0.023675441622844417, -0.0007261644860741899 ));
        boundary_plane_normals.push_back(Create_c_vector( -0.6727411476463672, -0.7373576645038901,  0.061016578573512954 ));


        boundary_plane_points.push_back(Create_c_vector( 0.039700467012298214, 0.015344441098080423, -0.0004565231516265349 ));
        boundary_plane_normals.push_back(Create_c_vector( -0.1706692384045456, -0.9809531883324618,  0.09275156798022458   ));



        boundary_plane_points.push_back(Create_c_vector(0.036804371668431334, 0.053919748549890005, -0.0007634162035095087 ));
        boundary_plane_normals.push_back(Create_c_vector(-0.6416261436122902, 0.7612012690416773,  0.09427894697418274   ));

  
     
        for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
      
             boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(&(p_simulator->rGetCellPopulation()) , boundary_plane_points[boundary_id], boundary_plane_normals[boundary_id], 0.01));
             p_simulator->AddCellPopulationBoundaryCondition(p_condition);
        }

            
        for (int i =1; i<=10; i++)
        { 
            EndTime +=1;
            p_simulator->SetEndTime(EndTime);

            p_simulator->Solve();
            CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
        }

            // dt/=2 ;  SamplingStep*= 2; 
            // // FSIIterations*=2;
            // p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 
 
            // // p_ForceOut->SetFluidSolidIterations(FSIIterations);
            // p_simulator->SetSamplingTimestepMultiple(SamplingStep);
            // p_simulator->SetDt(dt);
            // p_Mesh_modifier->SetUpdateFrequency(2/dt);
    }



        



 void offTestMembraneParameters2() throw(Exception)
    {

        double DilationParameter = -7;
        double AreaParameter = -7;
        double DeformationParamter = -8;

        std::stringstream out;
        out << "DilationParameter_" << DilationParameter << "AreaParameter" << AreaParameter << "DeformationParamter" << DeformationParamter;
        std::string ParameterSet = out.str();
        std::string output_dir = "DeformingPlexus/TestingParameter/";

        TRACE("Jess is good")
        double EndTime = 5;
        double SamplingStep = 1;
        double dt = 0.001;
        double NewEndTime = EndTime+1;
        std::string Archieved = "DeformingPlexus/FlatForceFINAL5/";

        OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
        /* Update the ouput directory for the population  */
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
        // (p_simulator->rGetCellPopulation()).SetWriteVtkAsPoints(false);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetWriteVtkAsPoints(false);
        double EdgeLength =0.00045;
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);
        static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).ExecuteHistoryDependentRemeshing();
        // (p_simulator->rGetCellPopulation()).ExecuteHistoryDependentRemeshing();

        p_simulator->RemoveAllCellPopulationBoundaryConditions();
        p_simulator->RemoveAllForces();
        p_simulator->SetSamplingTimestepMultiple(SamplingStep);
        p_simulator->SetDt(dt);
        p_simulator->SetOutputDirectory(output_dir);

        boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
        p_membrane_force->SetMembraneStiffness(pow(10, -9));
        p_simulator->AddForce(p_membrane_force);
        /*
        -----------------------------
        Membrane forces
        ----------------------------
        */
        boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
        p_simulator->AddForce(p_shear_force);

        /*
        -----------------------------
        Constant Compressive tissue pressure
        ----------------------------
        */
        double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
        double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

        boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
        p_ForceOut->SetPressure((P_blood - P_tissue));
        p_simulator->AddForce(p_ForceOut);

        /*
        -----------------------------
        Update membrane properties
        ----------------------------
        */
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
        boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);
        p_Mesh_modifier->TurnOffRemeshing();
        // p_Mesh_modifier->SetRemeshingInterval(10);

        //AreaConstant           AreaDilationModulus        ShearModulus
        std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), 1e-10) } };

        p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
        p_Mesh_modifier->SetmSetUpSolve(1);


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

        boundary_plane_points.push_back(Create_c_vector(0.051817210775777105,0.05387570153154082, -0.00019851668457055458 ));
        boundary_plane_normals.push_back(Create_c_vector(0.8769643127867986, 0.4790312122440197,  0.038245153868211076 ));

        boundary_plane_points.push_back(Create_c_vector(0.06117298902054698,0.03654224598923696,0.0009160147066410285 ));
        boundary_plane_normals.push_back(Create_c_vector(0.7792901394890018, 0.6252088966438896,  -0.04266982601959587  ));

        boundary_plane_points.push_back(Create_c_vector(0.055625630272326206, 0.0200783962665125, -0.000722688871900608 ));
        boundary_plane_normals.push_back(Create_c_vector( 0.7212385318941745, -0.6926132859491577, 0.010090403254884016  ));

        boundary_plane_points.push_back(Create_c_vector(0.05266362792368705, 0.011971267579691799,-0.0014918860363199145   ));
        boundary_plane_normals.push_back(Create_c_vector( -0.1352823453802107, -0.9855392137170218, 0.10203501974549702  ));


        boundary_plane_points.push_back(Create_c_vector( 0.026287008430771083, 0.023675441622844417, -0.0007261644860741899 ));
        boundary_plane_normals.push_back(Create_c_vector( -0.6727411476463672, -0.7373576645038901,  0.061016578573512954 ));
        
  
     
        for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
        {
      
             boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(&(p_simulator->rGetCellPopulation()) , boundary_plane_points[boundary_id], boundary_plane_normals[boundary_id], 0.01));
             p_simulator->AddCellPopulationBoundaryCondition(p_condition);
        }




        p_simulator->SetEndTime(NewEndTime);
        p_simulator->Solve();
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
    }
   




  void offTestMembraneParameters() throw(Exception)
    {

        double N_D[2] = { -8, -7 };

        for (unsigned i = 1; i < 2; i++)
        {
            for (unsigned j = 0; j < 2; j++)
            {
                for (unsigned k = 0; k < 2; k++)
                {

                    double DilationParameter = N_D[i];
                    double AreaParameter = N_D[j];
                    double DeformationParamter = N_D[k];

                    std::stringstream out;
                    out << "DilationParameter_" << DilationParameter << "AreaParameter" << AreaParameter << "DeformationParamter" << DeformationParamter;
                    std::string ParameterSet = out.str();
                    std::string output_dir = "DeformingPlexus/ParameterSweep/" + ParameterSet;

                    TRACE("Jess is good")
                    double EndTime = 8;
                    double SamplingStep = 500;
                    double dt = 0.0001;
                    double NewEndTime = EndTime+0.05;
                    std::string Archieved = "DeformingPlexus/FlatForce4";

                    OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
                    /* Update the ouput directory for the population  */
                    static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);

                    p_simulator->RemoveAllForces();
                    p_simulator->SetSamplingTimestepMultiple(SamplingStep);
                    p_simulator->SetDt(dt);
                    p_simulator->SetOutputDirectory(output_dir);

                    boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
                    p_membrane_force->SetMembraneStiffness(pow(10, -9));
                    p_simulator->AddForce(p_membrane_force);
                    /*
                    -----------------------------
                    Membrane forces
                    ----------------------------
                    */
                    boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
                    p_simulator->AddForce(p_shear_force);

                    /*
                    -----------------------------
                    Constant Compressive tissue pressure
                    ----------------------------
                    */
                    double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg
                    double P_tissue = 0.001466542; // Pa == 1.5000e-05 mmHg , need to set up some collasping force for this -- this should be taken into consideration for the membrane properties :)

                    boost::shared_ptr<OutwardsPressure> p_ForceOut(new OutwardsPressure());
                    p_ForceOut->SetPressure((P_blood - P_tissue));
                    p_simulator->AddForce(p_ForceOut);

                    /*
                    -----------------------------
                    Update membrane properties
                    ----------------------------
                    */
                    std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
                    boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);
                    p_Mesh_modifier->TurnOffRemeshing();
                    // p_Mesh_modifier->SetRemeshingInterval(10);

                    //AreaConstant           AreaDilationModulus        ShearModulus
                    std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), 1e-10) } };

                    p_Mesh_modifier->SetMembranePropeties(GrowthMaps, 1);
                    p_Mesh_modifier->SetmSetUpSolve(1);

                    p_simulator->SetEndTime(NewEndTime);
                    p_simulator->Solve();
                    CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
                }
            }
        }
    }



};

#endif /*TESTRELAXATION_HPP_*/
