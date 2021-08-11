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

#include "EnclosedRegionBoundaryCondition.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"

#include "HistoryDepMeshBasedCellPopulation.hpp"
#include "HistoryDepMutableMesh.hpp"

#include "FixedRegionBoundaryCondition.hpp"
#include "MembraneBendingForce.hpp"
#include "MembraneBendingForce0TargetAngle.hpp"
#include "MembraneDeformationForce.hpp"
#include "OutwardsPressure.hpp"
#include "OutwardsPressureWithBreaks.hpp"

#include "HemeLBForce.hpp"

#include "RemeshingTriggerOnStepHeteroModifier.hpp"

class TestRemeshing : public AbstractCellBasedTestSuite
{
public:

  void TestMembraneParameters() throw(Exception)
    {

            double DilationParameter = -7.3;
            double AreaParameter = -7;
            double DeformationParamter = -8;
            double BendingParameter = -11;
            double FSI_Iterations = 3000;

            //AreaConstant           AreaDilationModulus        ShearModulus
            std::map<double, c_vector<long double, 4> > GrowthMaps = { { 1, Create_c_vector(pow(10, AreaParameter), pow(10, DilationParameter), pow(10, DeformationParamter), pow(10, BendingParameter)) },
                                                                        {0,  Create_c_vector(pow(10, -5), pow(10, -4), pow(10, -4), pow(10, BendingParameter))}    };
    



            std::stringstream out;
            out << "DilationParameter_" << DilationParameter << "AreaParameter" << AreaParameter << "DeformationParamter" << DeformationParamter;
            std::string ParameterSet = out.str();
            std::string output_dir = "HoneyComb_HemeLBTHird/Collapse" ;//+ ParameterSet;

            TRACE("Jess is good")
            double EndTime = 11;
            double SamplingStep = 50;
            double dt = 0.001;
            double NewEndTime = EndTime+10;

            std::string Archieved = "DeformingHoneyComb/RemeshingStep/";

            OffLatticeSimulation<2, 3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Load(Archieved, EndTime);
            /* Update the ouput directory for the population  */
            static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetChasteOutputDirectory(output_dir, EndTime);
            static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()). SetWriteVtkAsPoints(false);
            double EdgeLength = 0.0003;//(2e-6 * scale); 0.00045/2;//(2e-6 * scale);
            static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).SetTargetRemeshingEdgeLength(EdgeLength);
            // static_cast<HistoryDepMeshBasedCellPopulation<2, 3>&>(p_simulator->rGetCellPopulation()).ExecuteHistoryDependentRemeshing();

            p_simulator->RemoveAllForces();
            p_simulator->SetSamplingTimestepMultiple(SamplingStep);
            p_simulator->SetDt(dt);
            p_simulator->SetOutputDirectory(output_dir);

            boost::shared_ptr<MembraneBendingForce> p_membrane_force(new MembraneBendingForce());
            p_membrane_force->SetMembraneStiffness(pow(10, -8));
            p_simulator->AddForce(p_membrane_force);
            /*
            -----------------------------
            Membrane forces
            ----------------------------
            */
            boost::shared_ptr<MembraneDeformationForce> p_shear_force(new MembraneDeformationForce());
            p_simulator->AddForce(p_shear_force);


            ///////////////////////////////////////////////////////////////////////////////////////

            std::vector<c_vector<double, 3> > boundary_plane_points1;
            std::vector<c_vector<double, 3> > boundary_plane_normals1;

            std::vector<c_vector<double, 3> > boundary_plane_points2;
            std::vector<c_vector<double, 3> > boundary_plane_normals2;


            /* */
            boundary_plane_points1.push_back(Create_c_vector(0.017985684246494982, 0,0));
            boundary_plane_normals1.push_back(Create_c_vector(1,0,0));
            /// ----------------------------------------
            boundary_plane_points2.push_back(Create_c_vector(0.026139036127872058, 0,0 ));
            boundary_plane_normals2.push_back(Create_c_vector(1,0,0));

            /* */
            boundary_plane_points1.push_back(Create_c_vector(0.034963365591332625, 0,0));
            boundary_plane_normals1.push_back(Create_c_vector(1,0,0));
            /// ----------------------------------------
            boundary_plane_points2.push_back(Create_c_vector(0.04307533991933138, 0,0 ));
            boundary_plane_normals2.push_back(Create_c_vector(1,0,0));


            /* */
            boundary_plane_points1.push_back(Create_c_vector(0.05203772956749085, 0,0));
            boundary_plane_normals1.push_back(Create_c_vector(1,0,0));
            /// ----------------------------------------
            boundary_plane_points2.push_back(Create_c_vector(0.06008748970025353, 0,0 ));
            boundary_plane_normals2.push_back(Create_c_vector(1,0,0));

    
            unsigned counter =0;
            for (unsigned boundary_id = 0; boundary_id < boundary_plane_points1.size(); boundary_id++)
            {
                counter +=1;
                boost::shared_ptr<EnclosedRegionBoundaryCondition<2, 3> > p_condition(new EnclosedRegionBoundaryCondition<2, 3>(&(p_simulator->rGetCellPopulation()) , boundary_plane_points1[boundary_id], boundary_plane_normals1[boundary_id], 1)); //0.01));

                p_condition->SetPointOnPlane2( boundary_plane_points2[boundary_id]);
                p_condition->SetNormalToPlane2(boundary_plane_normals2[boundary_id]);


                p_simulator->AddCellPopulationBoundaryCondition(p_condition);
            }



            std::vector<c_vector<double, 3> > boundary_plane_points1D1;
            std::vector<c_vector<double, 3> > boundary_plane_normals1D1;

            std::vector<c_vector<double, 3> > boundary_plane_points2D1;
            std::vector<c_vector<double, 3> > boundary_plane_normals2D1;


            
            /* Diagonal 1 */
            boundary_plane_points1D1.push_back(Create_c_vector(0.012427020567408235, -0.014536134484497733 , -3.0128174549181702e-5 ));
            boundary_plane_normals1D1.push_back(Create_c_vector(0.7071067811865476,  -0.7071067811865476, 0  ));
            /// ----------------------------------------
            boundary_plane_points2D1.push_back(Create_c_vector( 0.017158353759139884,  -0.01926746767622938  , -3.0128174549181702e-5  ));
            boundary_plane_normals2D1.push_back(Create_c_vector( 0.7071067811865476,  -0.7071067811865476, 0   ));

            /* Diagonal 2 */
            boundary_plane_points1D1.push_back(Create_c_vector(0.011316836157625475, -0.0018692261147451787, 0.001297028002651899 ));
            boundary_plane_normals1D1.push_back(Create_c_vector(0.7071067811865476,  -0.7071067811865476, 0  ));
            /// ----------------------------------------
            boundary_plane_points2D1.push_back(Create_c_vector(0.01511884615100136, -0.005671236108121062,0.001297028002651899   ));
            boundary_plane_normals2D1.push_back(Create_c_vector( 0.7071067811865476,  -0.7071067811865476, 0   ));

            /* Diagonal 3 */ 
            boundary_plane_points1D1.push_back(Create_c_vector(0.027451585483880808,  -0.009569147737977221, 0.0024318538097381783   ));
            boundary_plane_normals1D1.push_back(Create_c_vector(0.7071067811865476,  -0.7071067811865476, 0  ));
            /// ----------------------------------------
            boundary_plane_points2D1.push_back(Create_c_vector(0.03145119971027805,-0.013568761964374502,  0.0024318538097381783  ));
            boundary_plane_normals2D1.push_back(Create_c_vector( 0.7071067811865476,  -0.7071067811865476, 0   ));


            /* Diagonal 4 */
            boundary_plane_points1D1.push_back(Create_c_vector( 0.02708341835052748,0.005187797143975152,  0.00235468249113942   ));
            boundary_plane_normals1D1.push_back(Create_c_vector(0.7071067811865476,  -0.7071067811865476, 0  ));
            /// ----------------------------------------
            boundary_plane_points2D1.push_back(Create_c_vector(0.031625346773351694, 0.000645868721150932, 0.00235468249113942  ));
            boundary_plane_normals2D1.push_back(Create_c_vector( 0.7071067811865476,  -0.7071067811865476, 0   ));


            /* Diagonal 5 */
            boundary_plane_points1D1.push_back(Create_c_vector( 0.04497342599254667, -0.0024329349981689874, 0.003183050459649167   ));
            boundary_plane_normals1D1.push_back(Create_c_vector(0.7071067811865476,  -0.7071067811865476, 0  ));
            /// ----------------------------------------
            boundary_plane_points2D1.push_back(Create_c_vector(0.04858572474811474,  -0.006045233753737065, 0.003183050459649167 ));
            boundary_plane_normals2D1.push_back(Create_c_vector( 0.7071067811865476,  -0.7071067811865476, 0   ));


             /* Diagonal 5 */
            boundary_plane_points1D1.push_back(Create_c_vector( 0.04645169728971372, -0.014822167982229684, 0.0035170578062890282  ));
            boundary_plane_normals1D1.push_back(Create_c_vector(0.7071067811865476,  -0.7071067811865476, 0  ));
            /// ----------------------------------------
            boundary_plane_points2D1.push_back(Create_c_vector( 0.050977788945146305,-0.019348259637662277,  0.0035170578062890282 ));
            boundary_plane_normals2D1.push_back(Create_c_vector( 0.7071067811865476,  -0.7071067811865476, 0   ));


             /* Diagonal 5 */
            // boundary_plane_points1D1.push_back(Create_c_vector( 0.04537342524530815, -0.0013781875058484624, 0.004500704378491893   ));
            // boundary_plane_normals1D1.push_back(Create_c_vector(0.7071067811865476,  -0.7071067811865476, 0  ));
            // /// ----------------------------------------
            // boundary_plane_points2D1.push_back(Create_c_vector(0.04945107627070632, -0.005455838531246647, 0.004500704378491893  ));
            // boundary_plane_normals2D1.push_back(Create_c_vector( 0.7071067811865476,  -0.7071067811865476, 0   ));


             /* Diagonal 6 */
            boundary_plane_points1D1.push_back(Create_c_vector(0.06296600825580362, -0.00829049002434537, 0.005023984536790793    ));
            boundary_plane_normals1D1.push_back(Create_c_vector(0.7071067811865476,  -0.7071067811865476, 0  ));
            /// ----------------------------------------
            boundary_plane_points2D1.push_back(Create_c_vector(0.06665564222165832, -0.011980123990200122, -0.011980123990200122 ));
            boundary_plane_normals2D1.push_back(Create_c_vector( 0.7071067811865476,  -0.7071067811865476, 0   ));


            /* Diagonal 7 */
            boundary_plane_points1D1.push_back(Create_c_vector(0.06282271772073537,0.006561715785187714,0  ));
            boundary_plane_normals1D1.push_back(Create_c_vector(0.7071067811865476,  -0.7071067811865476, 0  ));
            /// ----------------------------------------
            boundary_plane_points2D1.push_back(Create_c_vector(0.066891117006989,  0.0024945562814045054,0 ));
            boundary_plane_normals2D1.push_back(Create_c_vector( 0.7071067811865476,  -0.7071067811865476, 0   ));

       


            for (unsigned boundary_id = 0; boundary_id < boundary_plane_points2D1.size(); boundary_id++)
            {
                boost::shared_ptr<EnclosedRegionBoundaryCondition<2, 3> > p_condition(new EnclosedRegionBoundaryCondition<2, 3>(&(p_simulator->rGetCellPopulation()) , boundary_plane_points1D1[boundary_id], boundary_plane_normals1D1[boundary_id], 0.009)); //0.01));

                p_condition->SetPointOnPlane2( boundary_plane_points2D1[boundary_id]);
                p_condition->SetNormalToPlane2(boundary_plane_normals2D1[boundary_id]);


                p_simulator->AddCellPopulationBoundaryCondition(p_condition);
            }
            /////////////////////////////

            std::vector<c_vector<double, 3> > boundary_plane_points1D2;
            std::vector<c_vector<double, 3> > boundary_plane_normals1D2;

            std::vector<c_vector<double, 3> > boundary_plane_points2D2;
            std::vector<c_vector<double, 3> > boundary_plane_normals2D2;

            /* Diagonal 1 */
            boundary_plane_points1D2.push_back(Create_c_vector( 0.012230017099932914,  -0.012839892138548806, 0    ));
            boundary_plane_normals1D2.push_back(Create_c_vector( 0.7062437429245811,   0.7062437429245811,0     ));
            /// ----------------------------------------
            boundary_plane_points2D2.push_back(Create_c_vector( 0.015920150882942945,  -0.009149758355538776, 0       ));
            boundary_plane_normals2D2.push_back(Create_c_vector(0.7062437429245811,   0.7062437429245811,0    ));

            /* Diagonal 1 */
            boundary_plane_points1D2.push_back(Create_c_vector( 0.011729039598103213, 0.001567785414569594,0  ));
            boundary_plane_normals1D2.push_back(Create_c_vector( 0.7062437429245811,   0.7062437429245811,0     ));
            /// ----------------------------------------
            boundary_plane_points2D2.push_back(Create_c_vector( 0.015984337363639935, 0.005823083180106322,  0));
            boundary_plane_normals2D2.push_back(Create_c_vector(0.7062437429245811,   0.7062437429245811,0 ));

            /* Diagonal 1 */
            boundary_plane_points1D2.push_back(Create_c_vector(0.028682097853496537, -0.020756420678830635, 0 ));
            boundary_plane_normals1D2.push_back(Create_c_vector( 0.7062437429245811,   0.7062437429245811,0     ));
            /// ----------------------------------------
            boundary_plane_points2D2.push_back(Create_c_vector( 0.03300036879360806,  -0.01643845333505734,  0));
            boundary_plane_normals2D2.push_back(Create_c_vector(0.7062437429245811,   0.7062437429245811,0 ));

        
            /* Diagonal 1 */
            boundary_plane_points1D2.push_back(Create_c_vector(0.029580721775208946, -0.006450219329190718,  0 ));
            boundary_plane_normals1D2.push_back(Create_c_vector( 0.7062437429245811,   0.7062437429245811,0     ));
            /// ----------------------------------------
            boundary_plane_points2D2.push_back(Create_c_vector( 0.0333188092694944, -0.002712131834905262,   0));
            boundary_plane_normals2D2.push_back(Create_c_vector(0.7062437429245811,   0.7062437429245811,0 ));


            /* Diagonal 1 */
            boundary_plane_points1D2.push_back(Create_c_vector(0.045683197960433165, -0.012307360481108595,  0 ));
            boundary_plane_normals1D2.push_back(Create_c_vector( 0.7062437429245811,   0.7062437429245811,0     ));
            /// ----------------------------------------
            boundary_plane_points2D2.push_back(Create_c_vector( 0.049269597023970886,-0.008720961417570875,   0));
            boundary_plane_normals2D2.push_back(Create_c_vector(0.7062437429245811,   0.7062437429245811,0 ));

            /* Diagonal 1 */
            boundary_plane_points1D2.push_back(Create_c_vector(0.046871459363420374,0.0002045585330829235 ,  0 ));
            boundary_plane_normals1D2.push_back(Create_c_vector( 0.7062437429245811,   0.7062437429245811,0     ));
            /// ----------------------------------------
            boundary_plane_points2D2.push_back(Create_c_vector( 0.051270092497221566,0.004603191666884118 ,   0));
            boundary_plane_normals2D2.push_back(Create_c_vector(0.7062437429245811,   0.7062437429245811,0 ));


            /* Diagonal 1 */
            boundary_plane_points1D2.push_back(Create_c_vector(0.06266027233519558 ,-0.02018472214581361 ,  0 ));
            boundary_plane_normals1D2.push_back(Create_c_vector( 0.7062437429245811,   0.7062437429245811,0     ));
            /// ----------------------------------------
            boundary_plane_points2D2.push_back(Create_c_vector(0.06686927268562991,-0.015975721795379316,   0));
            boundary_plane_normals2D2.push_back(Create_c_vector(0.7062437429245811,   0.7062437429245811,0 ));

            /* Diagonal 1 */
            boundary_plane_points2D2.push_back(Create_c_vector(0.061500680003425276, -0.004213131041115538, 0 ));
            boundary_plane_normals2D2.push_back(-Create_c_vector( 0.7062437429245811,   0.7062437429245811,0     ));
            /// ----------------------------------------
            boundary_plane_points1D2.push_back(Create_c_vector( 0.06519758411107461, -0.0005173535084617061 ,   0));
            boundary_plane_normals1D2.push_back(-Create_c_vector(0.7062437429245811,   0.7062437429245811,0 ));



            /* Diagonal 1 */
            boundary_plane_points2D2.push_back(Create_c_vector(0.061500680003425276, -0.004213131041115538, 0 ));
            boundary_plane_normals2D2.push_back(-Create_c_vector( 0.7062437429245811,   0.7062437429245811,0     ));
            /// ----------------------------------------
            boundary_plane_points1D2.push_back(Create_c_vector( 0.06519758411107461, -0.0005173535084617061 ,   0));
            boundary_plane_normals1D2.push_back(-Create_c_vector(0.7062437429245811,   0.7062437429245811,0 ));




            for (unsigned boundary_id = 0; boundary_id < boundary_plane_points2D2.size(); boundary_id++)
            {
                boost::shared_ptr<EnclosedRegionBoundaryCondition<2, 3> > p_condition(new EnclosedRegionBoundaryCondition<2, 3>(&(p_simulator->rGetCellPopulation()) , boundary_plane_points1D2[boundary_id], boundary_plane_normals1D2[boundary_id], 0.008)); //0.01));

                p_condition->SetPointOnPlane2( boundary_plane_points2D2[boundary_id]);
                p_condition->SetNormalToPlane2(boundary_plane_normals2D2[boundary_id]);


                p_simulator->AddCellPopulationBoundaryCondition(p_condition);
            }


            /*
            -----------------------------
            Boundary conditions
            ----------------------------
            */

            std::vector<c_vector<double, 3> > boundary_plane_points;
            std::vector<c_vector<double, 3> > boundary_plane_normals;
            boundary_plane_points.push_back(Create_c_vector(0.008952984699510745,0,0 ));
            boundary_plane_normals.push_back(Create_c_vector(-1, 0, 0));

            boundary_plane_points.push_back(Create_c_vector(0.06911882536638361,0,0));
            boundary_plane_normals.push_back(Create_c_vector(1,0,0));
            for (unsigned boundary_id = 0; boundary_id < boundary_plane_points.size(); boundary_id++)
            {  
                boost::shared_ptr<FixedRegionBoundaryCondition<2, 3> > p_condition(new FixedRegionBoundaryCondition<2, 3>(&(p_simulator->rGetCellPopulation()), boundary_plane_points[boundary_id], boundary_plane_normals[boundary_id],2));
                p_simulator->AddCellPopulationBoundaryCondition(p_condition);
            }

            ///////////////////////////////////////////////////////////////////////////////////////


            double P_blood = 0.002133152; // Pa ==   1.6004e-05 mmHg

            double InletPressure = P_blood * 1.001; // Fluid - Tissue pressure, think about adding a negative tissue force in the HemeLB force. but do this later
            double OutletPressure = P_blood * (0.999);

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

        

            /* 
            -----------------------------
            Edit  RemeshingTriggerOnStepHeteroModifier
            ----------------------------
            */
            std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2, 3> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();
            boost::shared_ptr<RemeshingTriggerOnStepHeteroModifier<2, 3> > p_Mesh_modifier = boost::static_pointer_cast<RemeshingTriggerOnStepHeteroModifier<2, 3> >(*iter);     
            // p_Mesh_modifier->SetRemeshingInterval(RemeshingTime); 
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
            p_Mesh_modifier->SetUpdateFrequency(0.1/dt);
            p_Mesh_modifier->SetmSetUpSolve(1);


           ///////////////////////////////////////////////////////////////////////////////////////


            for (int i = 0; i < 35; i++)
            {
                PRINT_VARIABLE(EndTime)
                EndTime +=1;
                p_simulator->SetEndTime(EndTime);
                p_simulator->Solve();
                CellBasedSimulationArchiver<2, OffLatticeSimulation<2, 3>, 3>::Save(p_simulator);
            }
            
    }



};

#endif /*TESTRELAXATION_HPP_*/
