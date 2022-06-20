#include "PressureForce.hpp"
#include "UblasCustomFunctions.hpp"
#include "Debug.hpp"
#include "VtkMeshWriter.hpp" 

// Area constrant on a membrane 


PressureForce::PressureForce()
   : AbstractForce<2,3>()
{
}



void PressureForce::AddForceContribution(AbstractCellPopulation<2,3>& rCellPopulation)
{
/*
 * Itterate over all the cells, get the force, then add the force .
*/
        // Helper variables
        MeshBasedCellPopulation<2,3>* p_cell_population = static_cast<MeshBasedCellPopulation<2,3>*>(&rCellPopulation);

        // TRACE("Writing the current chaste mesh");
        // /* Update mesh */
        // MutableMesh<2, 3>& Mesh = p_cell_population->rGetMesh();
        // MutableMesh<2, 3> *mMesh;
        // mMesh = static_cast<MutableMesh<2, 3>*>(&Mesh); 


        // std::string OutputFile = " /data/vascrem/testoutput/YeastDeformation/FirstIterationTest/";

        // VtkMeshWriter<2, 3> mesh_writer(OutputFile, "ChasteMesh", false);
        // mesh_writer.WriteFilesUsingMesh(*mMesh);

        // std::string InitialChasteMesh = " /data/vascrem/testoutput/YeastDeformation/FirstIterationTest/ChasteMesh.vtu";
        // std::string Initial_U = " reaction_diffusion_deforming_membrane/Output/Updated_u.xml";
        // std::string Initial_V = " reaction_diffusion_deforming_membrane/Output/Updated_v.xml";


        // mCurrentTime+= mCounter*FEMTimeStep;
        // mCurrentEndTime += mCounter*FEMTimeStep;
        // mCounter +=1;

        

        

        // TRACE(" Step 1: Attain initial concentration profile ")
        // std::string HemeLBCommand =  "cd projects/YeastProject/;./YeastBash";
        // HemeLBCommand +=InitialChasteMesh+Initial_U+Initial_V+OutputFile +" "+to_string(mCurrentTime)+" "+to_string(mCurrentEndTime)+" "+to_string(FEMTimeStep)+" 1 ";

        // int SystemOutput = std::system(HemeLBCommand.c_str()); 

        // mMesh_file=  OutputFile+"ReactionDiffusionForChaste_u.vtu"
    
        // VtkMeshReader<2,3> mesh_reader(mMesh_file);
        // MutableMesh<2,3> YeastMesh;
        // YeastMesh.ConstructFromMeshReader(mesh_reader);
        // // std::vector<double> mConcentrationVector;
        // mesh_reader.GetPointData("Concentration_u", mConcentrationVector);

        

        if (mFRMIterationCounter ==50)
        {
            GetConcentrationProfile(rCellPopulation);
            mFRMIterationCounter=0;
        }
        else
        {
            mFRMIterationCounter+=1;   
        }
        // 

        // Calculate midpoint
        c_vector<double,3> centroid = zero_vector<double>(3);
        for (AbstractCellPopulation<2,3>::Iterator cell_iter = rCellPopulation.Begin();
                     cell_iter != rCellPopulation.End();
                     ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            centroid += rCellPopulation.GetNode(node_index)->rGetLocation();
        }
        centroid /= rCellPopulation.GetNumRealCells();

        for (AbstractCellPopulation<2,3>::Iterator cell_iter = rCellPopulation.Begin();
                     cell_iter != rCellPopulation.End();
                     ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<3>* p_node = rCellPopulation.GetNode(node_index);
            c_vector<double,3> cell_location = p_node->rGetLocation()-centroid;
            cell_location(2) = 0.0;
            
            c_vector<double, 3> Force = zero_vector<double>(3);
            c_vector<double,3> normal = zero_vector<double>(3);

            std::set<unsigned>&  containing_elements = p_node->rGetContainingElementIndices();
            assert(containing_elements.size()>0);
            for (std::set<unsigned>::iterator iter = containing_elements.begin();
                    iter != containing_elements.end();
                    ++iter)
            {
                // Negative as normals point inwards for these surface meshes
                normal += - p_cell_population->rGetMesh().GetElement(*iter)->CalculateNormal();
            }
            double Pressure = 1e-3;
            normal /= norm_2(normal);
            double CellArea = rCellPopulation.GetVolumeOfCell(*cell_iter);
            Force =  mConcentrationVector[node_index]*Pressure* normal;//CellArea; //mStrength * CellArea * normal;


            cell_iter->GetCellData()->SetItem("area", CellArea);
            
            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(Force);

            cell_iter->GetCellData()->SetItem("applied_force_x", Force[0]);
            cell_iter->GetCellData()->SetItem("applied_force_y", Force[1]);
            cell_iter->GetCellData()->SetItem("applied_force_z", Force[2]);
            cell_iter->GetCellData()->SetItem("applied_force_mag", norm_2(Force));

            cell_iter->GetCellData()->SetItem("norm_x", normal[0]);
            cell_iter->GetCellData()->SetItem("norm_y", normal[1]);
            cell_iter->GetCellData()->SetItem("norm_z", normal[2]);

            cell_iter->GetCellData()->SetItem("Radius", norm_2(cell_location));
        }
}

// /home/vascrem/anaconda3/bin/activate
    //  conda config --set auto_activate_base false


 void PressureForce::SetReactionDiffusionResultsFile(std::string mesh_file)
{

    mMesh_file = mesh_file;
}

 void PressureForce::GetConcentrationProfile(AbstractCellPopulation<2,3>& rCellPopulation)
{

        /* 1) Write out current mesh */
        TRACE("GetConcentrationProfile")

        MeshBasedCellPopulation<2, 3>* p_cell_population = static_cast<MeshBasedCellPopulation<2, 3>*>(&rCellPopulation);
        MutableMesh<2, 3>& Mesh = p_cell_population->rGetMesh();
        mMesh = static_cast<MutableMesh<2, 3>*>(&Mesh); 

        mOutputDirectory = "/data/vascrem/testoutput/YeastDeformation/FirstIterationTest/";
        WriteOutVtuFile("YeastDeformation/FirstIterationTest/");


        /* 1) Set up conditions for FEM */
        std::string InitialChasteMesh = mOutputDirectory+"ChasteMesh.vtu";
        std::string Initial_U = mOutputDirectory+ "Updated_u.xml";
        std::string Initial_V = mOutputDirectory+ "Updated_V.xml"; 


        mCurrentTime+= mCounter*FEMTimeStep;
        mCurrentEndTime += mCounter*FEMTimeStep;
        mCounter +=1;

        /* 3) Orgonise meshes for fenix input */
        std::string HemeLBCommand =  "cd projects/YeastProject/;./YeastBash ";

        std::string waitFile = mOutputDirectory + "WaitFile.txt";
        HemeLBCommand +=InitialChasteMesh+" "+Initial_U+" "+Initial_V+" "+mOutputDirectory +" "+std::to_string(mCurrentTime)+" "+std::to_string(mCurrentEndTime)+" "+std::to_string(FEMTimeStep)+" "+std::to_string(mLoadPreviousResult)+" "+waitFile;
        mLoadPreviousResult = 1;
        int SystemOutput = std::system(HemeLBCommand.c_str()); 

        while(! boost::filesystem::exists(waitFile))
            {
                TRACE("waiting within C")
                sleep(2); 
            }
        // remove(waitFile.c_str())
        boost::filesystem::remove(waitFile.c_str());




        mMesh_file=  mOutputDirectory+"ReactionDiffusionForChaste_u000000.vtu";
    
        VtkMeshReader<2,3> mesh_reader(mMesh_file);
        MutableMesh<2,3> YeastMesh;
        YeastMesh.ConstructFromMeshReader(mesh_reader);
        // std::vector<double> mConcentrationVector;
        mesh_reader.GetPointData("Concentration_u", mConcentrationVector);

        

        /* 4) Orgonise meshes for fenix input */

        /* 5) Run fenics code */

        /* 6) Read in fenics output  */


        // std::string mesh_file = "projects/YeastProject/reaction_diffusion_deforming_membrane/Output/ReactionDiffusionResults.vtu";

        // Read in the data from fenix output 
        // VtkMeshReader<2,3> mesh_reader(mMesh_file);
        // MutableMesh<2,3> YeastMesh;
        // YeastMesh.ConstructFromMeshReader(mesh_reader);
        // std::vector<double> mConcentrationVector;
        // mesh_reader.GetPointData("Concentration", mConcentrationVector);

        // I now have everything in ConcentrationVector, I  need to associate this with the appropirate cells in my mesh. assume everything lines up nicely 


        for (AbstractCellPopulation<2,3>::Iterator cell_iter = rCellPopulation.Begin();
                        cell_iter != rCellPopulation.End();
                        ++cell_iter)
            {
                unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                Node<3>* p_node = rCellPopulation.GetNode(node_index);
                c_vector<double,3> cell_location = p_node->rGetLocation();
                // Check mesh and yeast mesh nodes match locations

                // This should only work on the first one
                // I dont know if this is working right... probs isnt .... are the cell counters right?????
                assert(std::abs(YeastMesh.GetNode(node_index)->GetPoint()[0]-cell_location[0])<0.01);
                // PRINT_3_VARIABLES(YeastMesh.GetNode(node_index)->GetPoint()[0],YeastMesh.GetNode(node_index)->GetPoint()[1],YeastMesh.GetNode(node_index)->GetPoint()[2]) ;
                // PRINT_3_VARIABLES(cell_location[0],cell_location[1],cell_location[2]);

                // PRINT_VECTOR(YeastMesh.GetNode(node_index)->GetPoint());
                // PRINT_VECTOR(cell_location);


                double ConcentrationForThisNode = mConcentrationVector[node_index];

                cell_iter->GetCellData()->SetItem("Concentration", ConcentrationForThisNode);
                mConcentrationAcrossYeast[node_index]=ConcentrationForThisNode;
            }
    }


void PressureForce::SetOutputDirectory(std::string outputDirectory)
{
     mOutputDirectory=outputDirectory;
}

void PressureForce::WriteOutVtuFile(std::string OutputDirectory)
{

    VtkMeshWriter<2, 3> mesh_writer(OutputDirectory , "ChasteMesh", false);
    mesh_writer.WriteFilesUsingMesh(*mMesh);
 
}



void PressureForce::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<2,3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PressureForce)


