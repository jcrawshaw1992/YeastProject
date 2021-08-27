#include "HemeLBForce.hpp"
#include "EmptyBasementMatrix.hpp"
/*
    Chaste force that runs HemeLB and uses the forces from the HemeLB output.
    This is basically the linking script -- Go Jess
    Need to 
    1) SO far I think this will only work on mac, im going to need to check the bash scripts 
    will run on linx, 
    2) Need to put in a bunch of asserts that HemeLB exists on the machine 
    3) Templating has been removed, it dosnt make sense to have anything other than <2,3>
*/

HemeLBForce::HemeLBForce()
        : AbstractForce<2, 3>()
{
}


HemeLBForce::~HemeLBForce()
{
}


void HemeLBForce::AddForceContribution(AbstractCellPopulation<2, 3>& rCellPopulation)
{
    double P_tissue = 0.001466542;
    assert(2 ==2); assert(3 ==3);
    HistoryDepMeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<HistoryDepMeshBasedCellPopulation<2, 3>*>(&rCellPopulation);
    // PRINT_2_VARIABLES(mExecuteHemeLBCounter, mTriggerHemeLB)
    if (mExecuteHemeLBCounter == mTriggerHemeLB)
    {
        TRACE("mExecuteHemeLBCounter == mTriggerHemeLB");
        // /* Update mesh */
        MutableMesh<2, 3>& Mesh = p_cell_population->rGetMesh();
        mMesh = static_cast<HistoryDepMutableMesh<2, 3>*>(&Mesh); 

        // WriteOutVtuFile(mOutputDirectory);

        // /* Run HemeLB */
        ExecuteHemeLB();
        // /* Get the traction  */
        LoadTractionFromFile();
        UpdateCellData(rCellPopulation);
        mExecuteHemeLBCounter = 0;
    }
    else
    {
        mExecuteHemeLBCounter += 1;
    }

    // I have an edges issues in the hemelB force. I think it is hemelb, but not sure, so here is what I am doing. 
    // I also dont know if the cylinder should be going in or out, but we will find out soo 
    // THis bit takes care of the edges ...

    // TRACE("Adding HemeLB Force to the nodes :) ")
	
	for (typename AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
		 cell_iter != rCellPopulation.End();
		 ++cell_iter)
	{
		unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
		Node<3>* pNode = rCellPopulation.rGetMesh().GetNode(node_index);
  
         CellPtr p_cell = p_cell_population->GetCellUsingLocationIndex(node_index);
        // if (p_cell->GetMutationState()->IsType<EmptyBasementMatrix> ())
        if ( p_cell->GetCellData()->GetItem("FixedBoundary") !=2  || mCollapseType == 1)
        {
            c_vector<double, 3> NormalVector = Create_c_vector(0,0,0);
            std::set<unsigned>& containing_elements = pNode->rGetContainingElementIndices();
            assert(containing_elements.size() > 0);
            // Finding the normal to the node -- need to check this 
            for (std::set<unsigned>::iterator iter = containing_elements.begin();
                iter != containing_elements.end();
                ++iter)
            {
                
                Node<3>* pNode0 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(0));
                Node<3>* pNode1 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(1));
                Node<3>* pNode2 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(2));

                c_vector<double, 3> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
                c_vector<double, 3> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3

                NormalVector  += VectorProduct(vector_12, vector_13);

                CellPtr p_cell1 = p_cell_population->GetCellUsingLocationIndex(pNode0->GetIndex());
                CellPtr p_cell2 = p_cell_population->GetCellUsingLocationIndex(pNode1->GetIndex());
                CellPtr p_cell3 = p_cell_population->GetCellUsingLocationIndex(pNode2->GetIndex());

                // double Counter = p_cell1->GetMutationState()->IsType<EmptyBasementMatrix> () +  p_cell2->GetMutationState()->IsType<EmptyBasementMatrix> () + p_cell3->GetMutationState()->IsType<EmptyBasementMatrix> (); 

                // if (Counter == 3)
                // {
                //     HasBasementElement = 1;
                // }


            }
            NormalVector /=norm_2(NormalVector); // I think the normal is inwards facing 


            double Pressure;
            if (cell_iter->GetCellData()->GetItem("Boundary") == 1)
            {
            c_vector<double, 3> AverageForce = Create_c_vector(0,0,0);
            c_vector<unsigned, 2> NearestNodes =  p_cell_population->GetNearestInternalNodes(node_index);

                for ( int i = 0; i <2; i++)
                {  
                    AverageForce += mForceMap[NearestNodes[i]];
                }
                Pressure = norm_2(AverageForce)/2;

            }else
            {
                c_vector<long double,3> force = mForceMap[node_index];
                Pressure = norm_2(force);

            }

            c_vector<long double,3> HemeLBForce = Pressure * NormalVector; 
            c_vector<long double,3> TissueForce = P_tissue * NormalVector; 
            c_vector<long double,3> Force =  (Pressure - P_tissue)* NormalVector; //                 HemeLBForce -TissueForce;
            pNode->AddAppliedForceContribution(Force); 
        }

       
    }




      
            // mMinSS = norm_2(*MinShearStress);
            // mMaxSS = norm_2(*MaxShearStress);
     

    // double Nc =45;
	// for (typename AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
	// 	 cell_iter != rCellPopulation.End();
	// 	 ++cell_iter)
	// {
    //     unsigned ReferenceNode = 0;
    //     unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
    //     Node<3>* pNode = p_cell_population->rGetMesh().GetNode(node_index);

    //     if (cell_iter->GetCellData()->GetItem("Boundary") == 1)
    //     {     
    //         if (node_index < Nc + 1) // if on lower edge
    //         {
    //             ReferenceNode = node_index + (2 * Nc); // select node from two rows up
    //         }
    //         else if (node_index > Nc) // if on upper edge
    //         {
    //             ReferenceNode = node_index - (2 * Nc); // select node from two rows down
    //         }
    //         pNode->AddAppliedForceContribution(mForceMap[ReferenceNode]); // Add the new force
    //         cell_iter->GetCellData()->SetItem("HemeLBForce", norm_2(mForceMap[ReferenceNode]));
    //     }
    //     else
    //     {
    //         pNode->AddAppliedForceContribution(mForceMap[node_index]); 
    //     }


    // }


    // double P_tissue = 0.001466542;

}




void HemeLBForce::SetCollapseType(double CollapseType)
{
    mCollapseType = CollapseType;
}


// void HemeLBForce::SetUpHemeLBConfiguration(std::string outputDirectory,  AbstractCellPopulation<2, 3>& rCellPopulation)
// {
// //    TRACE("SetUpHemeLBConfiguration -- only 2 inputs") -- Hit here
//     MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);
//     MutableMesh<2, 3>& Mesh = p_cell_population->rGetMesh();
//     mMesh = static_cast<HistoryDepMutableMesh<2, 3>*>(&Mesh); 

//     SetUpFilePaths(outputDirectory, 1,0);
//     // WriteHemeLBBashScript();  
//     ExecuteHemeLB();
//     LoadTractionFromFile();
//     UpdateCellData(rCellPopulation);
// }


void HemeLBForce::SetUpHemeLBConfiguration(std::string outputDirectory, AbstractCellPopulation<2, 3>& rCellPopulation, bool RunInitalHemeLB)
{  
    MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);
    MutableMesh<2, 3>& Mesh = p_cell_population->rGetMesh();
    mMesh = static_cast<HistoryDepMutableMesh<2, 3>*>(&Mesh); 
    bool RenamePriorResults =0;
    SetUpFilePaths(outputDirectory,RunInitalHemeLB,RenamePriorResults);
    /*  Need to generate the HemeLB bash script first */

    TRACE("A WRITE HEMELB RUN FILE")
    // WriteHemeLBBashScript();   
    if (RunInitalHemeLB ==1)
    {
        TRACE("RunInitalHemeLB")
        ExecuteHemeLB();
    }
    // TRACE("ExecuteHemeLB")
    LoadTractionFromFile();
    // TRACE("Done LoadTractionFromFile")
    UpdateCellData(rCellPopulation);
}



void HemeLBForce::ExecuteHemeLB()
{       


    // if (mRemeshingCounter =2 )
    // {
    //      HistoryDepMeshBasedCellPopulation<2, 3>* pCellPopulation = static_cast<HistoryDepMeshBasedCellPopulation<2, 3>*>(&rCellPopulation);
    //     pCellPopulation->ExecuteHistoryDependentRemeshing();
    // }
    // else
    // {
    //     mRemeshingCounter+=1;
    // }


    int SystemOutput; 
    /*  Rename prior results directory :)    */
    std::string OldResultsDirectory = mHemeLBDirectory + "results_PriorTimeStep/";
    if (boost::filesystem::exists(mResultsDirectory))
    {      
        if (boost::filesystem::exists(OldResultsDirectory))
         {   
            std::string RemoveOldResults ="rm -r " + OldResultsDirectory +" >nul";
            SystemOutput = system(RemoveOldResults.c_str());
         }

        std::rename(mResultsDirectory.c_str(), OldResultsDirectory.c_str());
        SystemOutput = system(mRemoveResultsDirectory.c_str());
        
    }
    else
    {
        OldResultsDirectory = mResultsDirectory;
    }

    // Set up bash script to run HemeLB
    WriteHemeLBBashScript();  


    // Setup HemeLB 
    /*  Step -1: Generate the inital vtu and stl files  ... Need mesh here, call it config file and save somewhere paralle   */
    WriteOutVtuFile(mOutputDirectory);

    /*  Step 0: Create the HemeLB config.pr2 file */
    double HemeLBSimulationTime = 5000; //
    int Period = HemeLBSimulationTime*0.95;
    Writepr2File(mHemeLBDirectory,HemeLBSimulationTime);
      

    // Step 1: Run HemeLB setup
    std::string run_hemelb_setup = mhemelb_setup_exe + ' ' + mHemeLBDirectory + "config.pr2";
    SystemOutput = std::system(run_hemelb_setup.c_str());

    // Step 2: Update xml file
    PRINT_2_VARIABLES(double_to_string(mExpectedVelocity, 6), double_to_string(mEstimatedIC, 11)  )
    std::string update_xml_file = "python projects/VascularRemodelling/apps/update_xml_file.py -period "+std::to_string(Period) +" -directory " + mHemeLBDirectory + " -InitalConditions " + double_to_string(mEstimatedIC, 11) + " -ConvergenceTermination false -AverageVelocity " + double_to_string(mExpectedVelocity, 20); 
    SystemOutput = std::system(update_xml_file.c_str());

    /*  Step 3: run HemeLB simulation
        This command will open up a new terminal and run the bash script RunHemeLB (which is in apps/). 
        Running HemeLB will take some time, so In the mean time I will set up for some things I will need to 
        make the flow.vtus and then run the wait_file bash script, which will wait until HemeLB has finished 
        (as definded by the generation of myfile.txt)
    */

    // Run HemeLB
    TRACE(" Step 3: run HemeLB simulation")
    std::string HemeLBCommand =  "cd "+mChasteOutputDirectory + mOutputDirectory + ";./RunHemeLB";
    // std::string HemeLBCommand =  "cd  /data/vascrem/testoutput/FSICylinder/Medium/Hetro5/HemeLBForce;./RunHemeLB";

    PRINT_VARIABLE(HemeLBCommand)
    if(mMachine =="server")
    {
        SystemOutput = std::system(HemeLBCommand.c_str());  //////SystemOutput = std::system("./projects/VascularRemodelling/apps/RunHemeLB");
    }else
    {
        SystemOutput = std::system("open ./projects/VascularRemodelling/apps/RunHemeLB");
    }

    TRACE("GmyUnstructuredGridReader")
    
    std::string GmyUnstructuredGridReader = "python " +mHemeLBPath+ "Tools/hemeTools/converters/GmyUnstructuredGridReader.py " + mHemeLBDirectory + "config.xml >nul"; 
    SystemOutput = std::system(GmyUnstructuredGridReader.c_str());


    // PRINT_VARIABLE(boost::filesystem::exists(mHemeLBDirectory + "WaitFile.txt"))
    while(! boost::filesystem::exists(mHemeLBDirectory + "WaitFile.txt"))
    {
        TRACE("waiting within C")
        sleep(15); 
    }


    WriteOpenVtus(Period, mCenterlinesNumber);
    std::string GetVUtus =   "cd "+mChasteOutputDirectory + mOutputDirectory + ";./OpenVtus";
    SystemOutput = std::system(GetVUtus.c_str() );
    if ( boost::filesystem::exists(mHemeLB_output + "Centerlines_"+std::to_string(mCenterlinesNumber)+".vtp") )
    {
        CopyFile(mHemeLBDirectory + "centerlines.vtp", mHemeLB_output + "Centerlines_"+std::to_string(mCenterlinesNumber)+".vtp");
    }
    
    mCenterlinesNumber +=1;

    /*  Step 3a: 
        While the current HemeLB simulation is running, sort some things out 
        i) get the files from the last time step moved and duplication (only the last file is duplicated, this one is and easy & stable flow result to grab)
    */

    // if (mFlowVtus)
    // {
        // TRACE("GmyUnstructuredGridReader")
        // std::string GmyUnstructuredGridReader = "python " +mHemeLBPath+ "Tools/hemeTools/converters/GmyUnstructuredGridReader.py " + mHemeLBDirectory + "config.xml >nul"; 
        // SystemOutput = std::system(GmyUnstructuredGridReader.c_str());

        // For the not first ones here is what I will do, this one is the set up so I wont bother here, but in future reps have the vtu sorting when HemelB is going
        // std::string GenerateFlowVtus = "python " +mHemeLBPath+ "Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory + "results/Extracted/surface-pressure.xtr " + mHemeLBDirectory + "results/Extracted/wholegeometry-velocity.xtr " + mHemeLBDirectory + "results/Extracted/surface-traction.xtr >nul";
        // SystemOutput = std::system(GenerateFlowVtus.c_str());
    // }

 
    // Sort vtus files  --- still need to properly set up the file directories, but I think this will happen later
    // if (mCenterlinesNumber >1)
    // {
    //     TRACE("mCenterlinesNumber HemeLB")
    //     std::ostringstream strs1;
    //     strs1 << mStartTime;
    //     std::string StartTime = strs1.str();
    //     PRINT_VARIABLE(StartTime)
    //     std::string vtuFileSorting = "python projects/VascularRemodelling/apps/SortVtuFiles.py -Directory " + mChasteOutputDirectory+mOutputDirectory + " -CurrentNumberOfFiles " + std::to_string(mLatestFinialHemeLBVTU) + " -Time " + StartTime+ " >nul";
    //     SystemOutput =  std::system(vtuFileSorting.c_str());
    //     if (mFlowVtus)
    //     {
    //         UpdateCurrentyFlowVtuCount();  
    //     }
        
    // // }
    // CopyFile(mHemeLBDirectory + "results/Extracted/wholegeometry-velocity_"+std::to_string(Period)+".vtu", mHemeLB_output + "wholegeometry-velocity_"+std::to_string(mCenterlinesNumber)+".vtu");
    // CopyFile(mHemeLBDirectory + "results/Extracted/surface-pressure_"+std::to_string(Period)+".vtu", mHemeLB_output + "surface-pressure_"+std::to_string(mCenterlinesNumber)+".vtu");
    
    // mCenterlinesNumber +=1;

    // ---- I can other things Chaste needs running in the background here Maybe have some potts things going on
    /* Now wait*/
    // std::string WaitCommand = "./projects/VascularRemodelling/apps/wait_file " + mHemeLBDirectory + "WaitFile.txt >nul";
    // SystemOutput =  std::system(WaitCommand.c_str());
    // // TRACE("Have waited long enough :)") 

    
    
    // /* Generate a new stl file from the vtu while HemeLB is going*/
    // std::string ConvertVTUtoSTL = "python projects/VascularRemodelling/apps/vtuTostl.py -Directory " + mHemeLBDirectory + " >nul";
    // SystemOutput = std::system(ConvertVTUtoSTL.c_str());


    // if (CheckIfSteadyStateAchieved() ==0)
    // {    ReRunHemeLB();  }

    // Make your own HemeLB sim :) 
    // ~/hemelb-dev/Tools/setuptool/scripts/hemelb-setup-nogui
    // python ~/Chaste/projects/VascularRemodelling/apps/update_xml_file.py -period 1000 -directory /data/vascrem/testoutput/TestHemeLBOnNetwork/CollapsingMiddelBranch/HemeLBFluid/ -InitalConditions 0 -ConvergenceTermination false -AverageVelocity 0 
    // python ~/hemelb-dev/Tools/hemeTools/converters/GmyUnstructuredGridReader.py config.xml     
    // python ~/hemelb-dev/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py config.vtu  results2/Extracted/surface-pressure.xtr results2/Extracted/wholegeometry-velocity.xtr results2/Extracted/surface-traction.xtr 
        
    // Set up to generate the vtu files -> SHould have a vtu generated here, i think from the gmy file 
    // if (mFlowVtus)
    // {
        // TRACE("GmyUnstructuredGridReader")
        // std::string GmyUnstructuredGridReader = "python " +mHemeLBPath+ "Tools/hemeTools/converters/GmyUnstructuredGridReader.py " + mHemeLBDirectory + "config.xml >nul"; 
        // SystemOutput = std::system(GmyUnstructuredGridReader.c_str());

        // // For the not first ones here is what I will do, this one is the set up so I wont bother here, but in future reps have the vtu sorting when HemelB is going
        // std::string GenerateFlowVtus = "python " +mHemeLBPath+ "Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory + "results/Extracted/surface-pressure.xtr " + mHemeLBDirectory + "results/Extracted/wholegeometry-velocity.xtr " + mHemeLBDirectory + "results/Extracted/surface-traction.xtr >nul";
        // SystemOutput = std::system(GenerateFlowVtus.c_str());
    // }

}




void HemeLBForce::WriteOpenVtus(int Period, int mCenterlinesNumber)
{  
        ofstream bash_script;

        std::string BashFile =  mChasteOutputDirectory + mOutputDirectory + "OpenVtus";
        PRINT_VARIABLE(BashFile)
        bash_script.open(BashFile);
        bash_script << "#!/bin/bash\n# chmod 700 OpenVtus\n";
        // bash_script << "python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/GmyUnstructuredGridReader.py " + mHemeLBDirectory + "config.xml >" + mHemeLBDirectory + "OpenVTUS.txt\n";
        std::cout<<" python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory +"results/Extracted/wholegeometry-velocity.xtr " + mHemeLBDirectory +"results/Extracted/surface-pressure.xtr " + mHemeLBDirectory +"results/Extracted/surface-traction.xtr  > " + mHemeLBDirectory + "OpenVTUS.txt\n";
        
        bash_script << " python /home/vascrem/hemelb-dev/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory +"results/Extracted/wholegeometry-velocity.xtr " + mHemeLBDirectory +"results/Extracted/surface-pressure.xtr " + mHemeLBDirectory +"results/Extracted/surface-traction.xtr  > " + mHemeLBDirectory + "OpenVTUS.txt\n";
        
        
        
         // std::string GenerateFlowVtus = "python " +mHemeLBPath+ "Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory + "results/Extracted/surface-pressure.xtr " + mHemeLBDirectory + "results/Extracted/wholegeometry-velocity.xtr " + mHemeLBDirectory + "results/Extracted/surface-traction.xtr >nul";
  
        
        // bash_script << "echo 'HemeLB has finished' > " + mHemeLBDirectory + "WaitFile.txt \n";
        bash_script << "cp " + mHemeLBDirectory +"results/Extracted/wholegeometry-velocity_"+std::to_string(Period)+".vtu " +mHemeLB_output + "wholegeometry-velocity_"+std::to_string(mCenterlinesNumber)+".vtu \n";
        // bash_script << "cp " + mHemeLBDirectory +"results/Extracted/surface-pressure_"+std::to_string(Period)+".vtu " +mHemeLB_output + "surface-pressure_"+std::to_string(mCenterlinesNumber)+".vtu \n";
        bash_script << "cp " + mHemeLBDirectory +"results/Extracted/surface-traction_"+std::to_string(Period)+".vtu " +mHemeLB_output + "surface-traction_"+std::to_string(mCenterlinesNumber)+".vtu \n";
        bash_script << "echo 'HemeLB simulation complete' \n";
        // bash_script << "osascript -e 'tell application \"Terminal\" to close first window' & exit";  Need to think about with with application to linux 
        bash_script.close();
        std::string compileBashScript = "chmod 700 " + BashFile + " >nul";
        int SystemOutput = std::system(compileBashScript.c_str());

        


}








void HemeLBForce::ReRunHemeLB()
{
    int SystemOutput; 
    TRACE("HemeLB had initially failed, running again with longer simulation time ")
    bool SimulationSuccess =0;
    double SimulationDuration[3] = {1800,1800, 18000};
    for(unsigned i=0; i< 3;i++)
    {
        SystemOutput = system(mRemoveResultsDirectory.c_str());

        Writepr2File(mHemeLBDirectory,SimulationDuration[i]);
        
        /*  Step 1: Run HemeLB setup */
        std::string run_hemelb_setup = mhemelb_setup_exe + ' ' + mHemeLBDirectory + "config.pr2 >nul";
         SystemOutput = std::system(run_hemelb_setup.c_str());

        /*  Step 2: Update xml file */
        std::string update_xml_file = "python projects/VascularRemodelling/apps/update_xml_file.py -period 1111 -directory " + mHemeLBDirectory + " >nul"; 
        SystemOutput = std::system(update_xml_file.c_str());

        SystemOutput = std::system("open ./projects/VascularRemodelling/apps/RunHemeLB");

        /* Now wait */
        std::string WaitCommand = "./projects/VascularRemodelling/apps/wait_file " + mHemeLBDirectory + "WaitFile.txt";
        SystemOutput = std::system(WaitCommand.c_str());
        bool Check = CheckIfSteadyStateAchieved();
        if (Check)
        {
            SimulationSuccess = 1;
            break;
        }
    }
    // if we get to here the HemeLB simulation still hasnt run after four attempts  -- Need to do a bunch of other things here, like generate flow vtus, but I will worry about that later next time hi have a break in the code 
    assert(SimulationSuccess ==1);
  
}


void HemeLBForce::SetFluidSolidIterations(double Iterations)
{
    mTriggerHemeLB = Iterations;
    // TRACE("Have set fluid iterations") -- This was fine 
}


void HemeLBForce::Inlets(c_vector<double, 3> PlaneNormal, c_vector<double, 3> Point, double pressure, std::string FlowDirection)
{
    // Set the boundary planes for this hetro region, set an upper and a lower bound.
    std::vector<c_vector<double, 3> > CurrentBoundary;
    CurrentBoundary.push_back(PlaneNormal * mHemeLBScalling);
    CurrentBoundary.push_back(Point * mHemeLBScalling);

    mIolets.push_back(CurrentBoundary);
    mPressure.push_back(pressure); //*mHemeLBScalling?
    if (FlowDirection == "o" || FlowDirection == "O" || FlowDirection == "outlet")
    {
        FlowDirection = "Outlet";
    }
    else if (FlowDirection == "i" || FlowDirection == "I" || FlowDirection == "inlet")
    {
        FlowDirection = "Inlet";
    }
    assert(FlowDirection == "Inlet" || FlowDirection == "Outlet");
    mType.push_back(FlowDirection);
    // TRACE("I think the inlets should not be an issue ") -- Fine
}





void HemeLBForce::WriteOutVtuFile(std::string outputDirectory)
{

    VtkMeshWriter<2, 3> mesh_writer(outputDirectory + "HemeLBFluid/", "Chaste", false);
    mesh_writer.WriteFilesUsingMesh(*mMesh);
 
    // Not sure why this is called chaste.vtu, but it could be a bloody issue
    std::string VtuToStl = "meshio-convert " + mHemeLBDirectory + "Chaste.vtu " + mHemeLBDirectory + "config.stl  >nul";
    int SystemOutput = std::system(VtuToStl.c_str());
    // VtuToStl = "meshio-convert " + mHemeLBDirectory + "Chaste.vtu " + mHemeLBDirectory + "config.vtu  >nul";
    // int SystemOutput2 = std::system(VtuToStl.c_str());


}



void HemeLBForce::SetHemeLBPath(std::string HemeLBPath)
{

    std::string mhemelb_setup_exe = "env PYTHONPATH=" + HemeLBPath +"/Tools/setuptool:$PYTHONPATH " + HemeLBPath +"/Tools/setuptool/scripts/hemelb-setup-nogui";
    mHemeLBPath = HemeLBPath;

}


void HemeLBForce::Writepr2File(std::string outputDirectory, double SimulationDuration)
{

// scp linalg/src/UblasCustomFunctions.hpp vascrem@josborne.science.unimelb.edu.au:/home/vascrem/Chaste/linalg/src/UblasCustomFunctions.hpp
    std::string CenterlinesFile = outputDirectory + "centerlines.vtp";
    /*  Need radius of mesh so I can specify HemeLB discretisation, and iolet cap sizes, to do this I am going to generate and sort the centerlines */
    std::string GenerateCenterlinesFile = "vmtk vmtknetworkextraction -ifile " + outputDirectory + "config.stl -ofile " + CenterlinesFile +" >nul";
    int SystemOutput = std::system(GenerateCenterlinesFile.c_str());

    /* Read the centerlines file and get the max radius */
    vtkSmartPointer<vtkXMLPolyDataReader> Reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    Reader->SetFileName(CenterlinesFile.c_str());
    Reader->Update();

    vtkPointData* point_data = Reader->GetOutput()->GetPointData();
    double NumberOfDataPoints = Reader->GetOutput()->GetNumberOfPoints();
    // assert(point_data->HasArray(mRadiusdataName.c_str())); // No Radi point data
    vtkDataArray* p_scalars = point_data->GetArray(mRadiusdataName.c_str());
    assert(p_scalars->GetNumberOfComponents() == 1); // Radi are scalars, so should only have one component for each data point, otherwise there is a problem

    std::vector<double> RadiVector;
    double MinRadius = 1000000;
    double MaxRadius = 0;
    for (unsigned i = 0; i < NumberOfDataPoints; i++)
    {
        double* data = p_scalars->GetTuple(i); //RadiVector.push_back(*data);
        double Radius_t = *data; 

        if (*data < MinRadius & *data  >0)
        {
            MinRadius = *data;
        }
        else if(*data > MaxRadius & *data  >0)
        {
            MaxRadius = *data;
        }


    }
    mRadius = MaxRadius * mHemeLBScalling;
    // mRadius = MinRadius * mHemeLBScalling;
    PRINT_2_VARIABLES(mRadius, mHemeLBScalling)

    // Here I want to be able to say if there has been a descrese from the initials ---
    
    /* I have the max radius, this will be important for the discretisation and the cap sizes 
        https://royalsocietypublishing.org/doi/pdf/10.1098/rsif.2014.0543   &&&&    https://journals.aps.org/pre/pdf/10.1103/PhysRevE.89.023303

        double maxMA = 0.2; -- I think it would be <0.1, but lets see
        Blood viscosity = 0.004 Pa s
        Blood density = 1000 kg m-3
        Period 10000 -- oscillation = 60/70 s.

        double C = 1/3; // C^2_s -- Bernabeu et al 
        double tau = 0.8;// Bernabeu et al. demonstrated that the relaxation constant must be between t = 0.5-1 for stable HemeLB simulations, with minimal error at t = 0:8 Bernabeu et al 2014
        double nV = C*(tau -0.5) ;// Nondenominational kinematic viscosity
     */

    double V = 4; // Kinematic viscosity -- 4 mm^2/s  V = eta/rho
    // double deltaX = 2*mRadius/41;//15; // Diameter/15 This will need thinking about later -- Need to talk to someone about 
 
 
    double deltaX = 2*mRadius/21;//15; // Diameter/15 This will need thinking about later -- Need to talk to someone about 
    double deltaT = 0.1 * deltaX * deltaX / V;
   
    double MaxPressure = *std::min_element(mPressure.begin(), mPressure.end());
    double MinPressure  = *std::max_element(mPressure.begin(), mPressure.end());

    // This is for calculating the the velocity through the vessel 
    mExpectedVelocity = fabs((MaxPressure-MinPressure)/(2*V) * mRadius* mRadius);
    PRINT_VARIABLE(mExpectedVelocity)

    //
    int InletNumber = 1;
    int OutletNumber = 1;
    double HemeLBSimulationDuration = SimulationDuration * deltaT;
    std::string HemeLBRunTime = std::to_string(HemeLBSimulationDuration);
    if (HemeLBSimulationDuration < 1e-4)
    {
      std::stringstream ss;
      ss << setprecision(16) << HemeLBSimulationDuration;
      std::string str;
      ss >> str;
      HemeLBRunTime =  str;
    }
    ofstream config_pr2;
    std::string ConfigFile = outputDirectory + "config.pr2";
    config_pr2.open(ConfigFile);
    config_pr2 << "DurationSeconds: " + HemeLBRunTime + "\nIolets:\n";
    mEstimatedIC = 0;
    for (unsigned i = 0; i < mIolets.size(); i++)
    {
        std::vector<c_vector<double, 3> > Iolets = mIolets[i];
        c_vector<double, 3> Normal = Iolets[0];
        c_vector<double, 3> Point = Iolets[1];
        config_pr2 << "- Centre: {x: " + std::to_string(Point[0]) + ", y: " + std::to_string(Point[1]) + ", z: " + std::to_string(Point[2]) + "}\n";
        if (mType[i] == "Inlet")
        {
            config_pr2 << "  Name: Inlet" + std::to_string(InletNumber) + "\n";
            InletNumber += 1;
        }
        else if (mType[i] == "Outlet")
        {
            config_pr2 << "  Name: Outlet" + std::to_string(OutletNumber) + "\n";
            OutletNumber += 1;
        }
        config_pr2 << "  Normal: {x: " + std::to_string(Normal[0]) + ", y: " + std::to_string(Normal[1]) + ", z: " + std::to_string(Normal[2]) + "}\n";
        config_pr2 << "  Pressure: {x: " + std::to_string(mPressure[i]) + ", y: 0.0, z: 0.0}\n";
        mEstimatedIC+=mPressure[i];
        config_pr2 << "  Radius: " + std::to_string(mRadius * 2) + "\n";
        config_pr2 << "  Type: " + mType[i] + "\n";
    }


    config_pr2 << "OutputGeometryFile: config.gmy\n";
    config_pr2 << "OutputXmlFile: config.xml\n";


    mEstimatedIC/=mIolets.size();
    
    /* Get seed point */
    c_vector<double, 3> Seed = mMesh->GetNode(200)->rGetLocation();     
  
    config_pr2 << "SeedPoint: {x: " + std::to_string(Seed[0]) + ", y: " + std::to_string(Seed[1]) + ", z: " + std::to_string(Seed[2]) + "}\n";
    config_pr2 << "StlFile: config.stl\n";
    config_pr2 << "StlFileUnitId: 1\n";
    std::ostringstream strs;
    strs << deltaT;
    std::string dt = strs.str();
    config_pr2 << "TimeStepSeconds: " + dt + "\n";

    std::ostringstream strs1;
    strs1 << deltaX;
    std::string dx = strs1.str();

    config_pr2 << "VoxelSize: " + dx;
    config_pr2.close();
}


void HemeLBForce::WriteHemeLBBashScript()
{
    PRINT_VARIABLE(mMachine)
    if(mMachine =="server")
    {
            // Need to write bash scrip .... issue here 
            TRACE("Have set to 20 cores, will decrease later")
            int Cores =21;
            ofstream bash_script;

            std::string BashFile =  mChasteOutputDirectory + mOutputDirectory + "RunHemeLB";
            PRINT_VARIABLE(BashFile)
            bash_script.open(BashFile);
            bash_script << "#!/bin/bash\n# chmod 700 RunHemeLB\n";
            bash_script << "mpirun -np " + std::to_string(Cores) + " hemelb -in " + mHemeLBDirectory + "config.xml -out " + mHemeLBDirectory +"results/ >" + mHemeLBDirectory + "HemeLBTerminalOutput.txt\n";
            bash_script << "echo 'HemeLB has finished' > " + mHemeLBDirectory + "WaitFile.txt \n";
            bash_script << "echo 'HemeLB simulation complete' \n";
            // bash_script << "osascript -e 'tell application \"Terminal\" to close first window' & exit";  Need to think about with with application to linux 
            bash_script.close();
            std::string compileBashScript = "chmod 700 " + BashFile + " >nul";
            int SystemOutput = std::system(compileBashScript.c_str());

    }else
    {
        // Need to write bash scrip .... issue here 
        int Cores = 2;
        ofstream bash_script;

        std::string BashFile = "projects/VascularRemodelling/apps/RunHemeLB";
        bash_script.open(BashFile);
        bash_script << "#!/bin/bash\n# chmod 700 RunHemeLB\n";
        bash_script << "echo " + mHemeLBDirectory + "config.xml\n";
        bash_script << "mpirun -np " + std::to_string(Cores) + " hemelb -in " + mHemeLBDirectory + "config.xml > " + mHemeLBDirectory + "HemeLBTerminalOutput.txt\n";
        bash_script << "echo 'HemeLB has finished' > " + mHemeLBDirectory + "WaitFile.txt \n";
        bash_script << "echo 'HemeLB simulation complete' \n";
        bash_script << "osascript -e 'tell application \"Terminal\" to close first window' & exit";
        bash_script.close();
        std::string compileBashScript = "chmod 700 " + BashFile + " >nul";
        int SystemOutput = std::system(compileBashScript.c_str());
    }
    

}

// /data/vascrem/HemeLBSimulaitons/TestingHemeLB
// /home/vascrem/hemelb-dev



void HemeLBForce::SetStartTime(double StartTime)
{
    mStartTime = StartTime;
    // TRACE("Here the start time is set") -- This was fine
}


double HemeLBForce::GetStartTime()
{
    return mStartTime;
}




void HemeLBForce::SetUpFilePaths(std::string outputDirectory, bool CreateFiles, bool RenamePriorResults)
{
    int SystemOutput;
    assert(3 == 3);
    if (!boost::algorithm::ends_with(outputDirectory, "/"))
    {   
        outputDirectory = outputDirectory + "/";
    }
    char* C = getenv("CHASTE_TEST_OUTPUT");
    std::string directory;
    directory += C ;
    mChasteOutputDirectory = directory +"/" ;
    
    mOutputDirectory = outputDirectory;
    
    mHemeLBDirectory = mChasteOutputDirectory + mOutputDirectory + "HemeLBFluid/";
    // mHemeLB_output = mChasteOutputDirectory + mOutputDirectory + "HemeLB_results_from_time_190/";

    mHemeLB_output = mChasteOutputDirectory + mOutputDirectory + "HemeLB_results_from_time_"+std::to_string(int(mStartTime))+"/";
    // PINRT_VARIABLE(mHemeLB_output)
    // PINRT_VARIABLE(mHemeLBDirectory)
    if (CreateFiles ==1)
    {  
        if (boost::filesystem::exists(mHemeLBDirectory))
        {
            /* HemeLB directory already exists  */
            if (RenamePriorResults == 1)
            { /* Rename Old one, I want an option to go back and have a look at old stuff  */
                time_t now = time(0);
                std::string OldDirectory = mHemeLBDirectory + "_Circa_" + ctime(&now);
                std::rename(mHemeLBDirectory.c_str(), OldDirectory.c_str());
            }
            else //if (boost::filesystem::exists(mHemeLBDirectory))
            {
                std::string RemoveOldHemeLBDirectory = "rm -r " + mHemeLBDirectory;
                SystemOutput = system(RemoveOldHemeLBDirectory.c_str());
            }
        }
         if (boost::filesystem::exists(mHemeLB_output))
        {
             if (RenamePriorResults == 1)
            { /* Rename Old one, I want an option to go back and have a look at old stuff  */
                time_t now = time(0);
                std::string OldDirectory = mHemeLB_output + "_Circa_" + ctime(&now);
                std::rename(mHemeLB_output.c_str(), OldDirectory.c_str());
            } 
            else
            { 
                std::string RemoveOldHemeLBDirectory = "rm -r " + mHemeLB_output +" >nul";
                SystemOutput =  system(RemoveOldHemeLBDirectory.c_str());
            }
        }
        /*  Create HemeLB Paths, parrallel to the results_0 folder, going to need the file name here   */
        boost::filesystem::path Redir(mHemeLB_output.c_str());
        boost::filesystem::create_directories(Redir);

        boost::filesystem::path dir(mHemeLBDirectory.c_str());
        boost::filesystem::create_directories(dir);
    }

    // Set up some more path names 
    mResultsDirectory = mHemeLBDirectory + "results/";
    if (boost::filesystem::exists(mResultsDirectory))
    {
         mRemoveResultsDirectory = "rm -r " + mResultsDirectory + " >nul";
    }
   

}


void HemeLBForce::LoadTractionFromFile()
{

    TRACE("LoadTractionFromFile")
    std::string TractionFile = mHemeLBDirectory + "results/Extracted/surface-tractions.xtr";
 
 //    std::string TractionFile = "/data/vascrem/testoutput/TestxtrDirectory/surface-tractions.xtr";

	// PRINT_VARIABLE(mTractionFile);
	FILE* traction_file = fopen((char*)TractionFile.c_str(), "r");
	assert(traction_file != NULL);
	
	hemelb::io::writers::xdr::XdrFileReader reader(traction_file);

	// File format described in http://pauli.chem.ucl.ac.uk/trac/hemelb/wiki/ExtractionFiles

    unsigned hemelb_magic_number, extraction_magic_number, extraction_version_number;
	// ding the traction file");
    reader.readUnsignedInt(hemelb_magic_number);

    reader.readUnsignedInt(extraction_magic_number);

    reader.readUnsignedInt(extraction_version_number);
	

    assert(hemelb_magic_number == 0x686c6221);
    assert(extraction_magic_number == 0x78747204);
    assert(extraction_version_number == 4);


    double voxel_size;
    
    reader.readDouble(voxel_size);

    c_vector<double,3> origin;
    reader.readDouble(origin[0]);
    reader.readDouble(origin[1]);
    reader.readDouble(origin[2]);

    // unsigned long long number_fluid_sites;
    uint64_t number_fluid_sites; // Changed this so It would work on the server
    reader.readUnsignedLong(number_fluid_sites);
   
    unsigned field_count;
    reader.readUnsignedInt(field_count);
    
    assert(field_count == 2); // Traction and tangetial component of traction
    unsigned header_length;
    reader.readUnsignedInt(header_length);
    

    // Traction field header
    std::string field_name;
    unsigned number_floats;
    double traction_offset;
    

    reader.readString(field_name, header_length);
    reader.readUnsignedInt(number_floats);
    assert(number_floats == 3);
    reader.readDouble(traction_offset);

    // Tangential traction field header
    double tanget_traction_offset;

    reader.readString(field_name, header_length);
    reader.readUnsignedInt(number_floats);
    assert(number_floats == 3);
    reader.readDouble(tanget_traction_offset);

    // Data section (we are reading a single timestep)
    // unsigned long long timestep_num;
    uint64_t timestep_num;  // Jess changed so would work on server
    reader.readUnsignedLong(timestep_num);

    mAppliedPosition.clear();
    mAppliedTractions.clear();
    mAppliedTangentTractions.clear();

    double MinimumShearStress =10000;
    double MaximumShearStress = -1000;

    for (unsigned fluid_site_index = 0; fluid_site_index <  number_fluid_sites; fluid_site_index++)
    {
        
    	{
            
    		c_vector<unsigned,3> coords;
    		reader.readUnsignedInt(coords[0]);
    		reader.readUnsignedInt(coords[1]);
    		reader.readUnsignedInt(coords[2]);

    		mAppliedPosition.push_back(origin+voxel_size*coords);
    	}

    	{        
			c_vector<float,3> traction;
			reader.readFloat(traction[0]);
			traction[0] += traction_offset;
			reader.readFloat(traction[1]);
			traction[1] += traction_offset;
			reader.readFloat(traction[2]);
			traction[2] += traction_offset;
			// PRINT_VECTOR(traction);

			assert(fabs(traction[0])<1e10);
			assert(fabs(traction[1])<1e10);
			assert(fabs(traction[2])<1e10);

			mAppliedTractions.push_back(traction);
            // mAppliedTractionsMap[fluid_site_index] = traction;
    	}

    	{
        	c_vector<float,3> tangent_traction;
			reader.readFloat(tangent_traction[0]);
			tangent_traction[0] += tanget_traction_offset;
			reader.readFloat(tangent_traction[1]);
			tangent_traction[1] += tanget_traction_offset;
			reader.readFloat(tangent_traction[2]);
			tangent_traction[2] += tanget_traction_offset;

			assert(fabs(tangent_traction[0])<1e10);
			assert(fabs(tangent_traction[1])<1e10);
			assert(fabs(tangent_traction[2])<1e10);

			mAppliedTangentTractions.push_back(tangent_traction);
            // mAppliedTangentTractionsMap[fluid_site_index] = tangent_traction;
            //  if (mCenterlinesNumber <2)
            //  {
                if (MinimumShearStress > norm_2(tangent_traction))
                {
                    MinimumShearStress = norm_2(tangent_traction);
                }
                else if (MaximumShearStress < norm_2(tangent_traction))
                {
                    MaximumShearStress = norm_2(tangent_traction);
                }
            //  }
    	}
        // if (mCenterlinesNumber <2)
        //     {
        //     mMaxSS = MaximumShearStress;
        //     mMinSS = MinimumShearStress;
        //     }
    }

    assert(mAppliedPosition.size() == number_fluid_sites);
    assert(mAppliedTractions.size() == number_fluid_sites);
    assert(mAppliedTangentTractions.size() == number_fluid_sites);
}



void HemeLBForce::Network(std::string Network)
{
    mNetwork = Network;
    if (Network == "Honeycomb" || Network == "Honycomb"||Network == "Honey" ||Network == "H" )
    {
        mRegionOfForceCollection = 0.001;
        mMinSS =  0.000424725;
        mMaxSS = 0.00211783; 
    }   
    else if(Network == "Plexus" || Network == "plexus"||Network == "P" ||Network == "p" )
    {
        mRegionOfForceCollection = 0.0007;
        mMinSS =  5.69087e-06;
        mMaxSS =  0.0016764; 
    }
    else if(Network == "VascularNetwork" || Network == "vascularnetwork"||Network == "VN" ||Network == "vn" )
    {
        mRegionOfForceCollection =0.002;// 0.0007;
        mMinSS =  5.69087e-06;
        mMaxSS =  0.0016764; 
    }
    
}


void HemeLBForce::UpdateCellData(AbstractCellPopulation<2,3>& rCellPopulation)
{
	
    double MaximumShearStress = 0;
    double MinimumShearStress = 1000;
    
    TRACE("UpdateCellData")
	MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);
  
  
	for (typename AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
		 cell_iter != rCellPopulation.End();
		 ++cell_iter)
	{

		c_vector<double, 3> location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
		unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
		Node<3>* pNode = rCellPopulation.rGetMesh().GetNode(node_index);
		unsigned nearest_fluid_site = 100;
		double distance_to_fluid_site = 10e5;

        double counter =0;
       
        c_vector<float,3> shear_stress = Create_c_vector(0,0,0);
		for (unsigned fluid_site_index = 0; fluid_site_index <  mAppliedPosition.size(); fluid_site_index++)
		{
			// Find the closest fluid site 
            //  TRACE("F")
			double distance = norm_2(location - mAppliedPosition[fluid_site_index]*1e3);
			if (distance < distance_to_fluid_site)
			{
				distance_to_fluid_site = distance;	
				nearest_fluid_site = fluid_site_index;
			}
            if (distance <  mRegionOfForceCollection)
            {
                // shear_stress +=mAppliedTangentTractions[fluid_site_index];
                shear_stress[0] +=mAppliedTangentTractions[fluid_site_index][0];
                shear_stress[1] +=mAppliedTangentTractions[fluid_site_index][1];
                shear_stress[2] +=mAppliedTangentTractions[fluid_site_index][2];
            
                counter+=1;
			}
		}
        // shear_stress/=counter;
        // PRINT_3_VARIABLES(shear_stress, counter,WallShearStress)
        // PRINT_3_VARIABLES(shear_stress[0], shear_stress[1],shear_stress[2])
        double WallShearStress2 = sqrt(shear_stress[0]*shear_stress[0]+shear_stress[1]*shear_stress[1]+shear_stress[2]*shear_stress[2])/counter;


           if (mCenterlinesNumber <=1)
             {
                if (MinimumShearStress > norm_2(shear_stress))
                {
                    MinimumShearStress = norm_2(shear_stress);
                }
                else if (MaximumShearStress < norm_2(shear_stress))
                {
                    MaximumShearStress = norm_2(shear_stress);
                }
             }
    
        
		// assert(nearest_fluid_site != UNSIGNED_UNSET);
	
		// Get the HemeLB force at the closest lattice site 
		c_vector<double,3> force = mAppliedTractions[nearest_fluid_site]/133.3223874;//;  Convert to Pas
		double Pressure = norm_2(force); 
        


        mForceMap[node_index] = force;//mAppliedTractions[nearest_fluid_site]/133.3223874;//;  Convert to Pas
		// Store the force in CellData
		cell_iter->GetCellData()->SetItem("HemeLBForce", Pressure);
        cell_iter->GetCellData()->SetItem("shear_stress",1);




        if (mCenterlinesNumber >=2)
             {
                PRINT_3_VARIABLES(mMinSS,  0.3*(mMaxSS- mMinSS) ,mMaxSS )
                if ( norm_2(shear_stress)>= 1.1*mMaxSS )
                {
                    cell_iter->GetCellData()->SetItem("WallShearStressExtremes", 1);
                }else if ( norm_2(shear_stress) <= 0.9*mMinSS)
                {
                    cell_iter->GetCellData()->SetItem("WallShearStressExtremes", -1);
                }
                else if ( norm_2(shear_stress) <= 0.95*mMinSS)
                {
                    cell_iter->GetCellData()->SetItem("WallShearStressExtremes", -2);
                }
                else if ( norm_2(shear_stress) <= 0.98*mMinSS)
                {
                    cell_iter->GetCellData()->SetItem("WallShearStressExtremes", -3);
                }
                else if  ( norm_2(shear_stress)<1.1*mMaxSS  &&  norm_2(shear_stress) > 0.9*mMinSS  ) 
                {
                         cell_iter->GetCellData()->SetItem("WallShearStressExtremes", 0);  
                         assert(norm_2(shear_stress) > 0.9*mMinSS );

                }
                else 
                {
                        cell_iter->GetCellData()->SetItem("WallShearStressExtremes", 0);  
                }
             }
	}


        if (mCenterlinesNumber <2)
            {
                mMaxSS = MaximumShearStress;
                mMinSS = MinimumShearStress;
                // PRINT_2_VARIABLES(mMaxSS, MaximumShearStress);
                // PRINT_2_VARIABLES(mMinSS, MinimumShearStress);
                // PRINT_2_VARIABLES(mMinSS , mMaxSS )
                // PRINT_2_VARIABLES(0.9*mMinSS ,1.1*mMaxSS )
                PRINT_2_VARIABLES(MinimumShearStress , MaximumShearStress )
              } 
              


}



void  HemeLBForce::SetCenterlinesNumber(double CenterlinesNumber)
{
    mCenterlinesNumber = CenterlinesNumber;
}




void  HemeLBForce::CopyFile(std::string InputDirectory, std::string OutputDirectory)
{
  boost::filesystem::copy_file(InputDirectory.c_str(), OutputDirectory.c_str()); // This was changed so the code could work on linux and mac
}



void  HemeLBForce::UpdateCurrentyFlowVtuCount()
{
        ifstream inFile;
        inFile.open(mHemeLB_output + "CurrentLastFluidOutput.txt");
        std::string line;
        std::getline(inFile, line);
        mLatestFinialHemeLBVTU = std::stod(line);
        inFile.close(); 
}



bool HemeLBForce::CheckIfSteadyStateAchieved()
{

    /* Check if steady state was ever achieved by HemeLB*/
    std::ifstream inFile;
    std::string file = mHemeLBDirectory + "HemeLBTerminalOutput.txt";
    inFile.open(file);

    if (!inFile) {
    std::cout << "Unable to open file \n";
    exit(1); // terminate with error
    }
    bool SSAchieved =0;
    std::string line;
    std::string phrase = "steady flow simulation converged.";

    while(std::getline(inFile, line)) // read the stream line by line
    {
        // std::cout<< line<< std::endl;
        // check if the phrase appears somewhere in the line (pos)
        std::string::size_type pos = line.find(phrase);
        

        if(pos != std::string::npos) // phrase found pos = position of phrase beginning
        {
            // turn the part of the line after the phrase into an input-stream
            // std::istringstream iss(line.substr(pos + phrase.size()));
            inFile.close(); 
            SSAchieved =1;
            break;
    
        }

    }
    inFile.close(); 
    return SSAchieved;
}


std::pair<std::string, int> HemeLBForce::exec(const char* cmd)
{
    std::array<char, 128> buffer;
    std::string result;
    int return_code = -1;
    auto pclose_wrapper = [&return_code](FILE* cmd) { return_code = pclose(cmd); };
    { // scope is important, have to make sure the ptr goes out of scope first
        const std::unique_ptr<FILE, decltype(pclose_wrapper)> pipe(popen(cmd, "r"), pclose_wrapper);
        if (pipe)
        {
            while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
            {
                result += buffer.data();
            }
        }
    }
    return make_pair(result, return_code);
}



std::string HemeLBForce::double_to_string(double Number, long double precision)
{

    std::string NumberString = std::to_string(Number);
    if (Number < 1e-4)
    {
      std::stringstream ss;
      ss << setprecision(precision) << Number;
      std::string str;
      ss >> str;
      NumberString =  str;
    }

    return NumberString;
}



void HemeLBForce::OutputForceParameters(out_stream& rParamsFile)
{
    // No parameters to include

    // Call method on direct parent class
    AbstractForce<2, 3>::OutputForceParameters(rParamsFile);
}



void HemeLBForce::SetMachine(std::string Machine)
{
    assert( Machine == "mac" || Machine == "server");
    mMachine = Machine;
}



void HemeLBForce::SetGenerateFlowVtus(bool FlowVtus)
{

    bool mFlowVtus = FlowVtus;
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////
// template class HemeLBForce;

// // Serialization for Boost >= 1.36
// #include "SerializationExportWrapperForCpp.hpp"
// EXPORT_TEMPLATE_CLASS_ALL_DIMS(HemeLBForce)
// // CHASTE_CLASS_EXPORT(HemeLBForce)

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(HemeLBForce)





// ----------------------------------------------------------------------------------------------------------

// if(mCollapsedRegion.size() !=0)
//         {

//             c_vector<double, 3> Node_location = location;
//             c_vector<double, 3> UpperPlane = -mCollapsedRegion[0] ;
//             c_vector<double, 3> UpperPoint = mCollapsedRegion[1];
//             c_vector<double, 3> LowerPlane = -mCollapsedRegion[2];
//             c_vector<double, 3> LowerPoint = mCollapsedRegion[3];

//              // Vector connecting the node to upper plane
//             c_vector<double, 3> NodeToUpperPlane = Node_location - UpperPoint;
//             c_vector<double, 3> NodeToLowerPlane = Node_location - LowerPoint;

//             double DotToUpperPlane = inner_prod(NodeToUpperPlane,UpperPlane );
//             double DotToLowerPlane = inner_prod(NodeToLowerPlane,LowerPlane );

//             double radius = 0.002; // XXX TODO Radius threshold needs fixing
//             if (DotToLowerPlane >= 0 && DotToUpperPlane >= 0)
//             {
//                 if (norm_2(NodeToUpperPlane) <radius ||  norm_2(NodeToLowerPlane)<radius  )
//                 {
//                     cell_iter->GetCellData()->SetItem("Pressure", 0);
//                    	cell_iter->GetCellData()->SetItem("applied_force_x", 0);
//                     cell_iter->GetCellData()->SetItem("applied_force_y", 0);
//                     cell_iter->GetCellData()->SetItem("applied_force_z", 0);
//                     cell_iter->GetCellData()->SetItem("applied_force_mag", 0);
//                     cell_iter->GetCellData()->SetItem("voronoi_cell_area", voronoi_cell_area);
//                     cell_iter->GetCellData()->SetItem("applied_shear_stress_x", 0);
//                     cell_iter->GetCellData()->SetItem("applied_shear_stress_y", 0);
//                     cell_iter->GetCellData()->SetItem("applied_shear_stress_z", 0);
//                     cell_iter->GetCellData()->SetItem("applied_shear_stress_mag", 0);
//                     mForceMap[node_index] = Create_c_vector(0,0,0);
//                 }
                
//             }
//         }



// double Y =0;
//     // determine if I need aditional 0 inlets ect 
//     if (mNewInlets ==1)
//     {   
//         if ( XValue.size()>0)
//         {// Additional inlets to be added 
//             if ( XValue.size()==1)
//             { // Just one inlet 
//                 double MaxX = *XValue.begin();
//                 c_vector<double, 3> UpperPlaneNormal = Create_c_vector(1,0,0);
//                 c_vector<double, 3> LowerPlaneNormal = Create_c_vector(-1,0,0);
//                 c_vector<double, 3> UpperPoint = Create_c_vector(MaxX,0,0);
//                 c_vector<double, 3> LowerPoint = Create_c_vector(MaxX,0,0);
//                 Inlets(UpperPlaneNormal , UpperPoint, 0, "Outlet");
//                 Inlets(LowerPlaneNormal , LowerPoint, 0, "Outlet");
//                 // CollapsedRegions(UpperPlaneNormal , UpperPoint, LowerPlaneNormal , LowerPoint);
//             }
//             if ( XValue.size()>1)
//             {   //Going to have a start and end to the inlets -- if there are more than one region needing cutting off, im not sure what to do, but will deal with later
//                 double MaxX = *std::min_element(XValue.begin(), XValue.end());
//                 double MinX  = *std::max_element(XValue.begin(), XValue.end());
//                 // Now I want to have the Max pointing towards pos and 
//                 // min pointing towards neg
//                 c_vector<double, 3> UpperPlaneNormal = Create_c_vector(1,0,0);
//                 c_vector<double, 3> LowerPlaneNormal = Create_c_vector(-1,0,0);
//                 c_vector<double, 3> UpperPoint = Create_c_vector(MaxX,0,0);
//                 c_vector<double, 3> LowerPoint = Create_c_vector(MinX,0,0);
//                 Inlets(UpperPlaneNormal , UpperPoint, 0, "Outlet");
//                 Inlets(LowerPlaneNormal , LowerPoint, 0, "Outlet");
//                 // CollapsedRegions(UpperPlaneNormal , UpperPoint, LowerPlaneNormal , LowerPoint);
//             }  
//            mNewInlets = 0; 
//         }
//     }



// 
// void HemeLBForce::UpdateCellData(AbstractCellPopulation<2,3>& rCellPopulation)
// {
	
// 	assert(3==3); // Currently assumes that 3 = 3
	
// 	// std::map<unsigned, c_vector<unsigned, 2>  > LatticeToNodeMap;
	

// 	MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);

// 	for (typename AbstractCellPopulation<2, 3>::Iterator cell_iter = rCellPopulation.Begin();
// 		 cell_iter != rCellPopulation.End();
// 		 ++cell_iter)
// 	{

// 		c_vector<double, 3> location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

// 		unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
// 		Node<3>* pNode = rCellPopulation.rGetMesh().GetNode(node_index);
// 		unsigned nearest_fluid_site = UNSIGNED_UNSET;
// 		double distance_to_fluid_site = DBL_MAX;
	
	
// 		for (unsigned fluid_site_index = 0; fluid_site_index <  mAppliedPosition.size(); fluid_site_index++)
// 		{
// 			// Find the closest fluid site 
// 			double distance = norm_2(location - mAppliedPosition[fluid_site_index]*1e3);
// 			if (distance < distance_to_fluid_site)
// 			{
// 				distance_to_fluid_site = distance;	
// 				nearest_fluid_site = fluid_site_index;
// 			}
// 		}
//         double counter =0;
//         c_vector<double,3> shear_stress ;
//         for (unsigned fluid_site_index = 0; fluid_site_index <  mAppliedPosition.size(); fluid_site_index++)
// 		{
// 			// Find the closest fluid site 
// 			double distance = norm_2(location - mAppliedPosition[fluid_site_index]*1e3);
// 			if (distance <  mRegionOfForceCollection)
//             {
//                 shear_stress +=mAppliedTangentTractions[fluid_site_index];
//                 mAppliedTractions[nearest_fluid_site]/133.3223874;//*voronoi_cell_area;  Convert to Pas
//                 counter+=1;
// 			}
// 		}
//         shear_stress/=counter;

// 		// PRINT_VARIABLE(distance_to_fluid_site);
// 		// LatticeToNodeMap[node_index] = Create_c_vector(nearest_fluid_site, 0 );
// 		assert(nearest_fluid_site != UNSIGNED_UNSET);
	
// 		// c_vector<double, 3> NormalVector = Create_c_vector(0,0,0);
// 		// std::set<unsigned>& containing_elements = pNode->rGetContainingElementIndices();
//         // assert(containing_elements.size() > 0);
// 		// Finding the normal to the node -- need to check this 
//         // for (std::set<unsigned>::iterator iter = containing_elements.begin();
//         //     iter != containing_elements.end();
//         //     ++iter)
//         // {
//         //     Node<3>* pNode0 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(0));
//         //     Node<3>* pNode1 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(1));
//         //     Node<3>* pNode2 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(2));

//         //     c_vector<double, 3> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
//         //     c_vector<double, 3> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3

//         //     NormalVector  += VectorProduct(vector_12, vector_13);
//         // }
// 		// NormalVector /=norm_2(NormalVector); // I think the normal is inwards facing 

// 		// Get the HemeLB force at the closest lattice site 
// 		c_vector<double,3> force = mAppliedTractions[nearest_fluid_site]/133.3223874;//*voronoi_cell_area;  Convert to Pas
// 		double Pressure = norm_2(force); 
        

// 		// Check I have maintained outward pointing norm by getting the direction of the force 
// 		// Vector and making sure the normal and the force vector are goin in the same direcvtion 
// 		// DO this by checking the angle betweem these two vectors is below a certain value -- basic dot proudct thing

// 		// c_vector<long double,3> forceDirection = force / Pressure;
//         // c_vector<long double,3> Force = Pressure * NormalVector; 
//         mForceMap[node_index] = force;
// 		// Store the force in CellData
// 		cell_iter->GetCellData()->SetItem("HemeLBForce", Pressure);
//         cell_iter->GetCellData()->SetItem("shear_stress", norm_2(shear_stress));
//         // pNode->AddAppliedForceContribution(Force); 


//         // unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
//         //     Node<3>* pNode = p_cell_population->rGetMesh().GetNode(node_index);
//         //     c_vector<double, 3> ForceOnNode = mForceMap[node_index]/1;
//         //     rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(ForceOnNode); 


// 		// double Angle = abs(acos(inner_prod(forceDirection, NormalVector) )) ;

// 		// if (Angle > M_PI/2)
// 		// {
// 		// 	// TRACE(" Normal was the wrong way");
// 		// 	NormalVector = -NormalVector;

// 		// }

// 		// Calculate the approximate area of the voronoi region around the cell by including a third of the area
// 		// of all surrounding triangles. Useful for turning stresses into forces.
// 		// double voronoi_cell_area = rCellPopulation.GetVolumeOfCell(*cell_iter);
		

// 		// PRINT_VARIABLE(Pressure);
// 		// // XXXX  -Replaced the direction of the force with the normal -- did this because the force acts normal to the lattice site it was selected from, which is not the identical position to the node, but is slightly off
// 		// location = location - centroid;
// 		// location /= norm_2(location);
		
		
       

// 		// c_vector<double,3> shear_stress = mAppliedTangentTractions[nearest_fluid_site];

// 		// assert(fabs(Force[0])<1e10);
// 		// assert(fabs(Force[1])<1e10);
// 		// assert(fabs(Force[2])<1e10);
// 		// assert(fabs(shear_stress[0])<1e10);
// 		// assert(fabs(shear_stress[1])<1e10);
// 		// assert(fabs(shear_stress[2])<1e10);

       
// 		// cell_iter->GetCellData()->SetItem("applied_force_x", Force[0]);
// 		// PRINT_VECTOR(Force);
// 		// cell_iter->GetCellData()->SetItem("applied_force_y", Force[1]);
// 		// cell_iter->GetCellData()->SetItem("applied_force_z", Force[2]);
// 		// cell_iter->GetCellData()->SetItem("applied_force_mag", norm_2(Force));
// 		// cell_iter->GetCellData()->SetItem("voronoi_cell_area", voronoi_cell_area);
// 		// cell_iter->GetCellData()->SetItem("applied_shear_stress_x", shear_stress[0]);
// 		// cell_iter->GetCellData()->SetItem("applied_shear_stress_y", shear_stress[1]);
// 		// cell_iter->GetCellData()->SetItem("applied_shear_stress_z", shear_stress[2]);
		
// 	}
// }
