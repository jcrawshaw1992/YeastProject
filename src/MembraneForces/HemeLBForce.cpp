#include "HemeLBForce.hpp"
/*
    Chaste force that runs HemeLB and uses the forces from the HemeLB output.
    This is basically the linking script -- Go Jess
    Need to 
    1) SO far I think this will only work on mac, im going to need to check the bash scripts 
    will run on linx, 
    2) Need to put in a bunch of asserts that HemeLB exists on the machine 


         TO DO 
    ----------------
    1) For the inital time step, have some measure of complexity of the vessel for how long we let the first time step run for, maybe this will
    be the spread of the radii 
    2) See if I can get the code to run for long time in more complicated geometry i.e check for reliability
*/
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HemeLBForce<ELEMENT_DIM, SPACE_DIM>::HemeLBForce()
        : AbstractForce<ELEMENT_DIM, SPACE_DIM>()
{
    TRACE("The constructor is called ")
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
HemeLBForce<ELEMENT_DIM, SPACE_DIM>::~HemeLBForce()
{
    TRACE("The distructor is called ")
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HemeLBForce<ELEMENT_DIM, SPACE_DIM>::AddForceContribution(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);
    if (mExecuteHemeLBCounter == mTriggerHemeLB)
    {
        // /* Update mesh */
        MutableMesh<ELEMENT_DIM, SPACE_DIM>& Mesh = p_cell_population->rGetMesh();
        mMesh = static_cast<HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>*>(&Mesh); 

        WriteOutVtuFile(mOutputDirectory);

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

   MAKE_PTR(EmptyBasementMatrix, p_Basement);
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        Node<SPACE_DIM>* pNode = p_cell_population->rGetMesh().GetNode(node_index);
        c_vector<double, 3> ForceOnNode = mForceMap[node_index]/1;
        rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(ForceOnNode); 
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HemeLBForce<ELEMENT_DIM, SPACE_DIM>::SetUpHemeLBConfiguration(std::string outputDirectory,  AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
   TRACE("SetUpHemeLBConfiguration -- only 2 inputs")
    MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);
    MutableMesh<ELEMENT_DIM, SPACE_DIM>& Mesh = p_cell_population->rGetMesh();
    mMesh = static_cast<HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>*>(&Mesh); 
    

    SetUpFilePaths(outputDirectory, 1,0);
    TRACE("B WRITE HEMELB RUN FILE")
    WriteHemeLBBashScript();  
    ExecuteHemeLB();
    LoadTractionFromFile();
    UpdateCellData(rCellPopulation);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HemeLBForce<ELEMENT_DIM, SPACE_DIM>::SetUpHemeLBConfiguration(std::string outputDirectory, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation, bool RunInitalHemeLB)
{  
    MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);
    MutableMesh<ELEMENT_DIM, SPACE_DIM>& Mesh = p_cell_population->rGetMesh();
    mMesh = static_cast<HistoryDepMutableMesh<ELEMENT_DIM, SPACE_DIM>*>(&Mesh); 
    bool RenamePriorResults =0;
    SetUpFilePaths(outputDirectory,RunInitalHemeLB,RenamePriorResults);
    /*  Need to generate the HemeLB bash script first */

    TRACE("A WRITE HEMELB RUN FILE")
    WriteHemeLBBashScript();   
    if (RunInitalHemeLB ==1)
    {
        ExecuteHemeLB();
    }
    LoadTractionFromFile();
    UpdateCellData(rCellPopulation);
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HemeLBForce<ELEMENT_DIM, SPACE_DIM>::ExecuteHemeLB()
{       
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

    // Setup HemeLB 
    /*  Step -1: Generate the inital vtu and stl files  ... Need mesh here, call it config file and save somewhere paralle   */
    WriteOutVtuFile(mOutputDirectory);

    /*  Step 0: Create the HemeLB config.pr2 file */
    double HemeLBSimulationTime = 500; // Too short right now, but who cares
    int Period = HemeLBSimulationTime/2;
    Writepr2File(mHemeLBDirectory,HemeLBSimulationTime);

    // Step 1: Run HemeLB setup
    std::string run_hemelb_setup = mhemelb_setup_exe + ' ' + mHemeLBDirectory + "config.pr2 >nul";
    SystemOutput = std::system(run_hemelb_setup.c_str());

    // Step 2: Update xml file
    std::string update_xml_file = "python /Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/apps/update_xml_file.py -period "+std::to_string(Period) +" -directory " + mHemeLBDirectory + " -InitalConditions " + std::to_string(mEstimatedIC) +" >nul"; 
    SystemOutput = std::system(update_xml_file.c_str());


    /*  Step 3: run HemeLB simulation
        This command will open up a new terminal and run the bash script RunHemeLB (which is in apps/). 
        Running HemeLB will take some time, so In the mean time I will set up for some things I will need to 
        make the flow.vtus and then run the wait_file bash script, which will wait until HemeLB has finished 
        (as definded by the generation of myfile.txt)
    */

    // Run HemeLB
    SystemOutput = std::system("open ./projects/VascularRemodelling/apps/RunHemeLB");

    /*  Step 3a: 
        While the current HemeLB simulation is running, sort some things out 
        i) get the files from the last time step moved and duplication (only the last file is duplicated, this one is and easy & stable flow result to grab)
    */


    // Sort vtus files  --- still need to properly set up the file directories, but I think this will happen later
    if (mCenterlinesNumber >1)
    {

        std::ostringstream strs1;
        strs1 << mStartTime;
        std::string StartTime = strs1.str();
        PRINT_VARIABLE(StartTime)
        std::string vtuFileSorting = "python /Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/apps/SortVtuFiles.py -Directory " + mOutputDirectory + " -CurrentNumberOfFiles " + std::to_string(mLatestFinialHemeLBVTU) + " -Time " + StartTime+ " >nul";
        SystemOutput =  std::system(vtuFileSorting.c_str());
        UpdateCurrentyFlowVtuCount();
    }
    
    CopyFile(mHemeLBDirectory + "centerlines.vtp", mHemeLB_output + "Centerlines_"+std::to_string(mCenterlinesNumber)+".vtp");
    mCenterlinesNumber +=1;


    // ---- I can other things Chaste needs running in the background here Maybe have some potts things going on
    /* Now wait*/
    std::string WaitCommand = "./projects/VascularRemodelling/apps/wait_file " + mHemeLBDirectory + "WaitFile.txt >nul";
    SystemOutput =  std::system(WaitCommand.c_str());

        /* Generate a new stl file from the vtu while HemeLB is going*/
    // std::string ConvertVTUtoSTL = "python /Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/apps/vtuTostl.py -Directory " + mHemeLBDirectory + " >nul";
    // SystemOutput = std::system(ConvertVTUtoSTL.c_str());

    // Set up to generate the vtu files 
    std::string GmyUnstructuredGridReader = "python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/GmyUnstructuredGridReader.py " + mHemeLBDirectory + "config.xml >nul"; 
    SystemOutput = std::system(GmyUnstructuredGridReader.c_str());

    if (CheckIfSteadyStateAchieved() ==0)
    { 
        ReRunHemeLB();
    }

    std::cout <<" Continue Chaste "<< std::endl;

    // For the not first ones here is what I will do, this one is the set up so I wont bother here, but in future reps have the vtu sorting when HemelB is going
    std::string GenerateFlowVtus = "python /Users/jcrawshaw/Documents/HemeLB/hemelb/Tools/hemeTools/converters/ExtractedPropertyUnstructuredGridReader.py " + mHemeLBDirectory + "config.vtu " + mHemeLBDirectory + "results/Extracted/surface-pressure.xtr " + mHemeLBDirectory + "results/Extracted/wholegeometry-velocity.xtr " + mHemeLBDirectory + "results/Extracted/surface-traction.xtr >nul";
    SystemOutput = std::system(GenerateFlowVtus.c_str());

}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HemeLBForce<ELEMENT_DIM, SPACE_DIM>::ReRunHemeLB()
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
        std::string update_xml_file = "python /Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/apps/update_xml_file.py -period 1111 -directory " + mHemeLBDirectory + " >nul"; 
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

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HemeLBForce<ELEMENT_DIM, SPACE_DIM>::SetFluidSolidIterations(double Iterations)
{
    mTriggerHemeLB = Iterations;
    TRACE("Have set fluid iterations")
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HemeLBForce<ELEMENT_DIM, SPACE_DIM>::Inlets(c_vector<double, 3> PlaneNormal, c_vector<double, 3> Point, double pressure, std::string FlowDirection)
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
    TRACE("I think the inlets should not be an issue ")
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HemeLBForce<ELEMENT_DIM, SPACE_DIM>::WriteOutVtuFile(std::string outputDirectory)
{

    VtkMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(outputDirectory + "HemeLBFluid/", "Chaste", false);
    mesh_writer.WriteFilesUsingMesh(*mMesh);

    std::string VtuToStl = "meshio-convert " + mHemeLBDirectory + "Chaste.vtu " + mHemeLBDirectory + "config.stl  >nul";
    int SystemOutput = std::system(VtuToStl.c_str());


}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HemeLBForce<ELEMENT_DIM, SPACE_DIM>::Writepr2File(std::string outputDirectory, double SimulationDuration)
{

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
    assert(point_data->HasArray(mRadiusdataName.c_str())); // No Radi point data
    vtkDataArray* p_scalars = point_data->GetArray(mRadiusdataName.c_str());
    assert(p_scalars->GetNumberOfComponents() == 1); // Radi are scalars, so should only have one component for each data point, otherwise there is a problem

    std::vector<double> RadiVector;
    double MinRadius = 1000000;
    for (unsigned i = 0; i < NumberOfDataPoints; i++)
    {
        double* data = p_scalars->GetTuple(i); //RadiVector.push_back(*data);
        if (*data < MinRadius & *data  >0)
        {
            MinRadius = *data;
        }
    }
    mRadius = MinRadius * mHemeLBScalling;

    /* I have the max radius, this will be important for the discretisation and the cap sizes 
        https://royalsocietypublishing.org/doi/pdf/10.1098/rsif.2014.0543   &&&&    https://journals.aps.org/pre/pdf/10.1103/PhysRevE.89.023303

        double maxMA = 0.2; -- I think it would be <0.1, but lets see
        Blood viscosity = 0.004 Pa s
        Blood density = 1000 kg m-3
        Period 10000oscillation = 60/70 s.

        double C = 1/3; // C^2_s -- Bernabeu et al 
        double tau = 0.8;// Bernabeu et al. demonstrated that the relaxation constant must be between t = 0.5-1 for stable HemeLB simulations, with minimal error at t = 0:8 Bernabeu et al 2014
        double nV = C*(tau -0.5) ;// Nondenominational kinematic viscosity
     */

    double V = 4; // Kinematic viscosity -- 4 mm^2/s  V = eta/rho
    double deltaX = mRadius / 15; // Diameter/15
    double deltaT = 0.1 * deltaX * deltaX / V;

    //
    int InletNumber = 1;
    int OutletNumber = 1;
    double HemeLBSimulationDuration = SimulationDuration * deltaT;
    ofstream config_pr2;
    std::string ConfigFile = outputDirectory + "/config.pr2";
    config_pr2.open(ConfigFile);
    config_pr2 << "DurationSeconds: " + std::to_string(HemeLBSimulationDuration) + "\nIolets:\n";
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
        config_pr2 << "  Radius: " + std::to_string(mRadius * 3) + "\n";
        config_pr2 << "  Type: " + mType[i] + "\n";
    }
    config_pr2 << "OutputGeometryFile: config.gmy\n";
    config_pr2 << "OutputXmlFile: config.xml\n";


    mEstimatedIC/=mIolets.size();
    // /* Get seed point */
    // typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = mMesh->GetNodeIteratorBegin();
    // std::advance(node_iter, 100);

    c_vector<double, SPACE_DIM> Seed = mMesh->GetNode(200)->rGetLocation();     
  
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

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HemeLBForce<ELEMENT_DIM, SPACE_DIM>::WriteHemeLBBashScript()
{
    int Cores = 2;
    ofstream bash_script;

    std::string BashFile = "/Users/jcrawshaw/Documents/Chaste/projects/VascularRemodelling/apps/RunHemeLB";
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

// /data/vascrem/HemeLBSimulaitons/TestingHemeLB
// /home/vascrem/hemelb-dev


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HemeLBForce<ELEMENT_DIM, SPACE_DIM>::SetStartTime(double StartTime)
{
    mStartTime = StartTime;
    TRACE("Here the start time is set")
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double HemeLBForce<ELEMENT_DIM, SPACE_DIM>::GetStartTime()
{
    return mStartTime;
}



template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HemeLBForce<ELEMENT_DIM, SPACE_DIM>::SetUpFilePaths(std::string outputDirectory, bool CreateFiles, bool RenamePriorResults)
{
    int SystemOutput;
    assert(SPACE_DIM == 3);
    if (!boost::algorithm::ends_with(outputDirectory, "/"))
    {   
        outputDirectory = outputDirectory + "/";
    }
    char* C = getenv("CHASTE_TEST_OUTPUT");
    std::string directory;
    directory += C ;
     std::string mChasteOutputDirectory = directory +"/" ;
    
    mOutputDirectory = outputDirectory;
    mHemeLBDirectory = mChasteOutputDirectory + mOutputDirectory + "HemeLBFluid/";
    mHemeLB_output = mChasteOutputDirectory + mOutputDirectory + "HemeLB_results_from_time_0/";

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
            else
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
            else{ 
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
    mRemoveResultsDirectory = "rm -r " + mResultsDirectory + " >nul";

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HemeLBForce<ELEMENT_DIM, SPACE_DIM>::LoadTractionFromFile()
{
    TRACE("This needs to be uncommented and fixed later")
    //  scons b=GccOpt projects/VascularRemodelling/test/SetUpHemeLBChasteLinkage/TestHemeLBForce.hpp 
    // alias vascser="ssh vascrem@josborne.science.unimelb.edu.au"
//     std::string TractionFile = mHemeLBDirectory + "results/Extracted/surface-tractions.xtr";
//     // std::string TractionFile = "/Users/jcrawshaw/Documents/ChasteWorkingDirectory/ShrunkPlexus/results/Extracted/surface-tractions.xtr";
// 	TRACE("Load tracrtion file");
// 	// PRINT_VARIABLE(TractionFile); 
// 	FILE* traction_file = fopen((char*)TractionFile.c_str(), "r");
// 	assert(traction_file != NULL);
	

//     // I need to figure out how to effecrtivly point to the HemeLB reader 
// 	hemelb::io::writers::xdr::XdrFileReader reader(traction_file);
	
// 	// File format described in http://pauli.chem.ucl.ac.uk/trac/hemelb/wiki/ExtractionFiles

//     unsigned hemelb_magic_number, extraction_magic_number, extraction_version_number;
// 	// TRACE("Reading the traction file");
//     reader.readUnsignedInt(hemelb_magic_number);
//     reader.readUnsignedInt(extraction_magic_number);
//     reader.readUnsignedInt(extraction_version_number);
	

//     assert(hemelb_magic_number == 0x686c6221);
//     assert(extraction_magic_number == 0x78747204);
//     assert(extraction_version_number == 4);

//     double voxel_size;
//     reader.readDouble(voxel_size);
//     // PRINT_VARIABLE(voxel_size)

//     c_vector<double,3> origin;
//     reader.readDouble(origin[0]);
//     reader.readDouble(origin[1]);
//     reader.readDouble(origin[2]);


//     unsigned long long number_fluid_sites;
//     // uint64_t& number_fluid_sites; // Changed this so It would work on the server
//     reader.readUnsignedLong(number_fluid_sites);
//     unsigned field_count;
//     reader.readUnsignedInt(field_count);
//     assert(field_count == 2); // Traction and tangetial component of traction
//     unsigned header_length;
//     reader.readUnsignedInt(header_length);

//     // Traction field header
//     std::string field_name;
//     unsigned number_floats;
//     double traction_offset;

//     reader.readString(field_name, header_length);
//     reader.readUnsignedInt(number_floats);
//     assert(number_floats == 3);
//     reader.readDouble(traction_offset);

//     // Tangential traction field header
//     double tanget_traction_offset;

//     reader.readString(field_name, header_length);
//     reader.readUnsignedInt(number_floats);
//     assert(number_floats == 3);
//     reader.readDouble(tanget_traction_offset);

//    // Data section (we are reading a single timestep)
//     unsigned long long timestep_num;
//     // uint64_t&  timestep_num; // Changed this so It would work on the server
//     reader.readUnsignedLong(timestep_num);
//     // PRINT_VARIABLE(timestep_num)

//     mAppliedPosition.clear();
//     mAppliedTractions.clear();
//     mAppliedTangentTractions.clear();

//     for (unsigned fluid_site_index = 0; fluid_site_index <  number_fluid_sites; fluid_site_index++)
//     {
//     	{
//     		c_vector<unsigned,3> coords;
//     		reader.readUnsignedInt(coords[0]);
//     		reader.readUnsignedInt(coords[1]);
//     		reader.readUnsignedInt(coords[2]);
    	
//             mAppliedPosition.push_back((origin+coords*voxel_size)*1e3); // mAppliedPosition.push_back(origin+voxel_size*coords); <---- this was the original line, its not good, give completluy wrong results, it looks like a time based problem, may need to revist 

//     	}

//     	{
// 			c_vector<float,3> traction;
// 			reader.readFloat(traction[0]);
// 			traction[0] += traction_offset;
// 			reader.readFloat(traction[1]);
// 			traction[1] += traction_offset;
// 			reader.readFloat(traction[2]);
// 			traction[2] += traction_offset;
// 			// PRINT_VECTOR(traction);

// 			assert(fabs(traction[0])<1e10);
// 			assert(fabs(traction[1])<1e10);
// 			assert(fabs(traction[2])<1e10);

// 			mAppliedTractions.push_back(traction);
//     	}

//     	{
// 			c_vector<float,3> tangent_traction;
// 			reader.readFloat(tangent_traction[0]);
// 			tangent_traction[0] += tanget_traction_offset;
// 			reader.readFloat(tangent_traction[1]);
// 			tangent_traction[1] += tanget_traction_offset;
// 			reader.readFloat(tangent_traction[2]);
// 			tangent_traction[2] += tanget_traction_offset;

// 			assert(fabs(tangent_traction[0])<1e10);
// 			assert(fabs(tangent_traction[1])<1e10);
// 			assert(fabs(tangent_traction[2])<1e10);

// 			mAppliedTangentTractions.push_back(tangent_traction);
//     	}
//     }

//     assert(mAppliedPosition.size() == number_fluid_sites);
//     assert(mAppliedTractions.size() == number_fluid_sites);
//     assert(mAppliedTangentTractions.size() == number_fluid_sites);




    // 
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HemeLBForce<ELEMENT_DIM, SPACE_DIM>::UpdateCellData(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
	
	assert(SPACE_DIM==3); // Currently assumes that SPACE_DIM = 3
	
	std::map<unsigned, c_vector<unsigned, 2>  > LatticeToNodeMap;
	

	MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation);

	for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
		 cell_iter != rCellPopulation.End();
		 ++cell_iter)
	{
		c_vector<double, SPACE_DIM> location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
		unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
		Node<SPACE_DIM>* pNode = rCellPopulation.rGetMesh().GetNode(node_index);
		unsigned nearest_fluid_site = UNSIGNED_UNSET;
		double distance_to_fluid_site = DBL_MAX;
	
	
		for (unsigned fluid_site_index = 0; fluid_site_index <  mAppliedPosition.size(); fluid_site_index++)
		{
			// Find the closest fluid site 
			double distance = norm_2(location - mAppliedPosition[fluid_site_index]);
			if (distance < distance_to_fluid_site)
			{
				distance_to_fluid_site = distance;	
				nearest_fluid_site = fluid_site_index;
			}
		}
		// PRINT_VARIABLE(distance_to_fluid_site);
		LatticeToNodeMap[node_index] = Create_c_vector(nearest_fluid_site, 0 );
		assert(nearest_fluid_site != UNSIGNED_UNSET);

	
		c_vector<double, 3> NormalVector = Create_c_vector(0,0,0);
		std::set<unsigned>& containing_elements = pNode->rGetContainingElementIndices();
        assert(containing_elements.size() > 0);
		// Finding the normal to the node -- need to check this 
        for (std::set<unsigned>::iterator iter = containing_elements.begin();
            iter != containing_elements.end();
            ++iter)
        {
            Node<SPACE_DIM>* pNode0 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(0));
            Node<SPACE_DIM>* pNode1 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(1));
            Node<SPACE_DIM>* pNode2 = p_cell_population->rGetMesh().GetNode(p_cell_population->rGetMesh().GetElement(*iter)->GetNodeGlobalIndex(2));

            c_vector<double, 3> vector_12 = pNode1->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 2
            c_vector<double, 3> vector_13 = pNode2->rGetLocation() - pNode0->rGetLocation(); // Vector 1 to 3

            NormalVector  += VectorProduct(vector_12, vector_13);
        }
		NormalVector /=norm_2(NormalVector);

		// Get the HemeLB force at the closest lattice site 
		c_vector<double,3> force = mAppliedTractions[nearest_fluid_site]/133.3223874;//*voronoi_cell_area;  Convert to Pas
		double Pressure = norm_2(force); 
        
        // PRINT_3_VARIABLES(distance_to_fluid_site, nearest_fluid_site, Pressure)    

		// Check I have maintained outward pointing norm by getting the direction of the force 
		// Vector and making sure the normal and the force vector are goin in the same direcvtion 
		// DO this by checking the angle betweem these two vectors is below a certain value -- basic dot proudct thing

		c_vector<long double,3> forceDirection = force / Pressure;

		double Angle = abs(acos(inner_prod(forceDirection, NormalVector) )) ;

		// if (Angle > M_PI/2)
		// {
		// 	// TRACE(" Normal was the wrong way");
		// 	NormalVector = -NormalVector;

		// }

		// Calculate the approximate area of the voronoi region around the cell by including a third of the area
		// of all surrounding triangles. Useful for turning stresses into forces.
		double voronoi_cell_area = rCellPopulation.GetVolumeOfCell(*cell_iter);
		

		// PRINT_VARIABLE(Pressure);
		// // XXXX  -Replaced the direction of the force with the normal -- did this because the force acts normal to the lattice site it was selected from, which is not the identical position to the node, but is slightly off
		// location = location - centroid;
		// location /= norm_2(location);
		
		c_vector<long double,3> Force = Pressure * NormalVector; 
       

		c_vector<double,3> shear_stress = mAppliedTangentTractions[nearest_fluid_site];

		assert(fabs(Force[0])<1e10);
		assert(fabs(Force[1])<1e10);
		assert(fabs(Force[2])<1e10);
		assert(fabs(shear_stress[0])<1e10);
		assert(fabs(shear_stress[1])<1e10);
		assert(fabs(shear_stress[2])<1e10);

        mForceMap[node_index] = Force;


		// Store the force in CellData
		cell_iter->GetCellData()->SetItem("Pressure", Pressure);
		// PRINT_VARIABLE(Pressure);


		// TRACE("Setting Force");
		cell_iter->GetCellData()->SetItem("applied_force_x", Force[0]);
		// PRINT_VECTOR(Force);
		cell_iter->GetCellData()->SetItem("applied_force_y", Force[1]);
		cell_iter->GetCellData()->SetItem("applied_force_z", Force[2]);
		cell_iter->GetCellData()->SetItem("applied_force_mag", norm_2(Force));
		cell_iter->GetCellData()->SetItem("voronoi_cell_area", voronoi_cell_area);
		cell_iter->GetCellData()->SetItem("applied_shear_stress_x", shear_stress[0]);
		cell_iter->GetCellData()->SetItem("applied_shear_stress_y", shear_stress[1]);
		cell_iter->GetCellData()->SetItem("applied_shear_stress_z", shear_stress[2]);
		cell_iter->GetCellData()->SetItem("applied_shear_stress_mag", norm_2(shear_stress));
	}
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void  HemeLBForce<ELEMENT_DIM, SPACE_DIM>::CopyFile(std::string InputDirectory, std::string OutputDirectory)
{
//   copyfile(InputDirectory.c_str(), OutputDirectory.c_str(), NULL, COPYFILE_DATA | COPYFILE_XATTR); 
//   std::filesystem::copy_file(InputDirectory.c_str(), OutputDirectory.c_str()); // This was changed so the code could work on linux and mac
  boost::filesystem::copy_file(InputDirectory.c_str(), OutputDirectory.c_str()); // This was changed so the code could work on linux and mac
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void  HemeLBForce<ELEMENT_DIM, SPACE_DIM>::UpdateCurrentyFlowVtuCount()
{
        ifstream inFile;
        inFile.open(mHemeLB_output + "CurrentLastFluidOutput.txt");
        std::string line;
        std::getline(inFile, line);
        mLatestFinialHemeLBVTU = std::stod(line);
        inFile.close(); 
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool HemeLBForce<ELEMENT_DIM, SPACE_DIM>::CheckIfSteadyStateAchieved()
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

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<std::string, int> HemeLBForce<ELEMENT_DIM, SPACE_DIM>::exec(const char* cmd)
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

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void HemeLBForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // No parameters to include

    // Call method on direct parent class
    AbstractForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////
template class HemeLBForce<2, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(HemeLBForce)
// CHASTE_CLASS_EXPORT(HemeLBForce)
