#include "PressureForce.hpp"
#include "UblasCustomFunctions.hpp"
#include "Debug.hpp"
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
   void AddForceContribution(AbstractCellPopulation<2,3>& rCellPopulation)
    {


        VtkMeshReader<2,3> mesh_reader(mMesh_file);
        MutableMesh<2,3> YeastMesh;
        YeastMesh.ConstructFromMeshReader(mesh_reader);
        // std::vector<double> mConcentrationVector;
        mesh_reader.GetPointData("Concentration", mConcentrationVector);



        // Helper variables
        MeshBasedCellPopulation<2,3>* p_cell_population = static_cast<MeshBasedCellPopulation<2,3>*>(&rCellPopulation);

        GetConcentrationProfile(rCellPopulation)

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
            // cell_location(2) = 0.0;
            
            c_vector<double, 3> force = zero_vector<double>(3);
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
            normal /= norm_2(normal);
            double CellArea = rCellPopulation.GetVolumeOfCell(*cell_iter);
            Force =  mConcentrationVector[node_index]*0.001 * normal/CellArea; //mStrength * CellArea * normal;


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
} 

 void PressureForce::SetReactionDiffusionResultsFile(std::string mesh_file)
{

    mMesh_file = std::string mesh_file;
}

 void PressureForce::GetConcentrationProfile(AbstractCellPopulation<2,3>& rCellPopulation)
{

    // std::string mesh_file = "projects/YeastProject/reaction_diffusion_deforming_membrane/Output/ReactionDiffusionResults.vtu";

    // Read in the data from fenix output 
    VtkMeshReader<2,3> mesh_reader(mMesh_file);
    MutableMesh<2,3> YeastMesh;
    YeastMesh.ConstructFromMeshReader(mesh_reader);
    // std::vector<double> mConcentrationVector;
    mesh_reader.GetPointData("Concentration", mConcentrationVector);

    // I now have everything in ConcentrationVector, I  need to associate this with the appropirate cells in my mesh. assume everything lines up nicely 


    for (AbstractCellPopulation<2,3>::Iterator cell_iter = rCellPopulation.Begin();
                     cell_iter != rCellPopulation.End();
                     ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<3>* p_node = rCellPopulation.GetNode(node_index);
            c_vector<double,3> cell_location = p_node->rGetLocation();
            // Check mesh and yeast mesh nodes match locations
            TS_ASSERT_DELTA(YeastMesh.GetNode(node_index)->GetPoint()[0]==cell_location[0]);

            ConcentrationForThisNode = ConcentrationVector[node_index];

            cell_iter->GetCellData()->SetItem("Concentration", ConcentrationForThisNode);
            mConcentrationAcrossYeast[node_index]=ConcentrationForThisNode;
        }
}



void PressureForce::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<2,3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PressureForce)


