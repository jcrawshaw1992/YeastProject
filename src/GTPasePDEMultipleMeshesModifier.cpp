/*

Copyright (c) 2005-2013, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "GTPasePDEMultipleMeshesModifier.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "TetrahedralMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "GTPasePDESystemSolver.hpp"
#include <set>
#include "PottsArbitrarySurfaceIn3DMesh.hpp"
#include <numeric>
#include "petscvec.h"
#include "FlowBasedCDC42StimulusProtocol.hpp"
#include "CPMGTPaseEventHandler.hpp"

template<unsigned DIM>
GTPasePDEMultipleMeshesModifier<DIM>::GTPasePDEMultipleMeshesModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mSolution(0),
      mpFeMesh(0),
      mPreviousFeMesh(0),
      mDependentVariableName("GTPasePDESolution"),
      mOutputDirectory("")
{
    assert(DIM==2);
}

template<unsigned DIM>
GTPasePDEMultipleMeshesModifier<DIM>::~GTPasePDEMultipleMeshesModifier()
{
    OutputFileHandler output_file_handler(mOutputDirectory, false);

    for (std::map<std::string, std::vector<std::string> >::const_iterator iter_cell = mFilesCollection.begin();
         iter_cell != mFilesCollection.end();
         ++iter_cell)
    {
        out_stream mpVtkMetaFile = output_file_handler.OpenOutputFile(iter_cell->first + ".pvd");
        *mpVtkMetaFile << "<?xml version=\"1.0\"?>\n";
        *mpVtkMetaFile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
        *mpVtkMetaFile << "    <Collection>\n";

        unsigned timestep = 0;
        for (std::vector<std::string>::const_iterator iter_file = iter_cell->second.begin();
             iter_file != iter_cell->second.end();
             ++iter_file, ++timestep)
        {
            *mpVtkMetaFile << "<DataSet timestep=\"" << timestep << "\" group=\"\" part=\"0\"\n";
            *mpVtkMetaFile << "         file=\"" << *iter_file << "\"/>\n";
        }

        *mpVtkMetaFile << "    </Collection>\n";
        *mpVtkMetaFile << "</VTKFile>\n";
    }
}


template<unsigned DIM>
void GTPasePDEMultipleMeshesModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Get the FeMesh from the cell population
    SetupFeMesh(rCellPopulation);

    // Construct the solution vector from cell data (takes care of cells dividing);
    UpdateSolutionVector(rCellPopulation);

    CPMGTPaseEventHandler::BeginEvent(CPMGTPaseEventHandler::SOLVE_PDE);
    for (unsigned cell_number=0; cell_number<rCellPopulation.GetNumRealCells(); ++cell_number)
    {
        FlowBasedCDC42StimulusProtocol<DIM,DIM> stimulus(*mpFeMesh.at(cell_number), 10.0, Create_c_vector(-1, 0));

        // Use previous solution as initial condition. In the first timestep, UpdateSolutionVector() already created a suitable initialisation.
        // The noise centres are computed new if the constructor gets a NULL value in the third argument.
    	GTPasePDESystemSolver<DIM, DIM> solver(mpFeMesh.at(cell_number), &stimulus, mSolution.at(cell_number));

        // Current time and time step of the Potts simulation
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        double current_time = p_simulation_time->GetTime();
        double dt = p_simulation_time->GetTimeStep();

        // Solve the PDE for a whole Potts timestep
        solver.SetTimes(current_time,current_time + dt);

        // PDE timestep needs to be kept small for stability and accuracy
        assert(dt > 0.025);
        solver.SetTimeStep(0.025);

        // Solve the current time step and store the solution for the next one
        Vec current_solution = solver.Solve();
        VecCopy(current_solution, mSolution[cell_number]);
    }
    CPMGTPaseEventHandler::EndEvent(CPMGTPaseEventHandler::SOLVE_PDE);

    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void GTPasePDEMultipleMeshesModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    CPMGTPaseEventHandler::BeginEvent(CPMGTPaseEventHandler::IO);

    if (DIM>1)
    {
        for (unsigned cell_number=0; cell_number<rCellPopulation.GetNumRealCells(); ++cell_number)
        {
            std::ostringstream cell_id_string;
            cell_id_string << "cell" << cell_number;

            std::ostringstream time_string;
            time_string << "time" << SimulationTime::Instance()->GetTimeStepsElapsed();

            std::string results_file = "pde_results_" + cell_id_string.str() + time_string.str();
            VtkMeshWriter<DIM,DIM>* p_vtk_mesh_writer = new VtkMeshWriter<DIM,DIM>(mOutputDirectory, results_file, false);

            ReplicatableVector solution_repl(mSolution[cell_number]);
            std::vector<std::vector<double> > each_variable_for_all_nodes(GTPASE_PROBLEM_DIM);
            for (unsigned i=0; i<GTPASE_PROBLEM_DIM*mpFeMesh[cell_number]->GetNumNodes(); i+=GTPASE_PROBLEM_DIM)
            {
                for (unsigned variable_id=0; variable_id<GTPASE_PROBLEM_DIM; ++variable_id)
                {
                    each_variable_for_all_nodes[variable_id].push_back(solution_repl[i+variable_id]);
                }
            }

            std::string variable_names[GTPASE_PROBLEM_DIM] = {"CDC42 Active", "RAC Active", "RHO Active",
                                                              "CDC42 Inactive", "RAC Inactive", "RHO Inactive",
                                                              "PIP1", "PIP2", "PIP3", "ARP",
                                                              "F0", "F1", "F2", "F3", "F4",
                                                              "B0", "B1", "B2", "B3", "B4",
                                                              "P0", "P1", "P2", "P3", "P4"};

            for (unsigned variable_id=0; variable_id<GTPASE_PROBLEM_DIM; ++variable_id)
            {
                p_vtk_mesh_writer->AddPointData(variable_names[variable_id], each_variable_for_all_nodes[variable_id]);
            }

            p_vtk_mesh_writer->WriteFilesUsingMesh(*mpFeMesh[cell_number]);
            delete p_vtk_mesh_writer;

            mFilesCollection["pde_results_" + cell_id_string.str()].push_back(results_file + ".vtu");
        }
    }

    CPMGTPaseEventHandler::EndEvent(CPMGTPaseEventHandler::IO);
}

template<unsigned DIM>
void GTPasePDEMultipleMeshesModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    // Cache the Output Directory
    mOutputDirectory = outputDirectory;

    //Setup a FeMesh to save the ics to
    SetupFeMesh(rCellPopulation);

    // Copy the cell data to mSolution (this is the initial conditions)
    UpdateSolutionVector(rCellPopulation);

    // Output the initial conditions on FeMesh
    UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template<unsigned DIM>
void GTPasePDEMultipleMeshesModifier<DIM>::SetupFeMesh(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    CPMGTPaseEventHandler::BeginEvent(CPMGTPaseEventHandler::SETUP_MESH);

    // Cast rCellPopulation to the appropriate type
    PottsBasedCellPopulation<DIM>* potts_cell_population = dynamic_cast<PottsBasedCellPopulation<DIM>*>(&rCellPopulation);
    assert(potts_cell_population);

    // Keep a copy of the previous meshes. This will be needed for solution interporlation when the cell changes shape.
    // Nothing will happen here in the first time step when mpFeMesh.size() is 0
    std::for_each(mPreviousFeMesh.begin(), mPreviousFeMesh.end(), std::default_delete<TetrahedralMesh<DIM,DIM> >());
    mPreviousFeMesh.resize(0);
    mPreviousFeMesh.insert(mPreviousFeMesh.begin(), mpFeMesh.begin(), mpFeMesh.end());

    // Pointer to the original mesh. It will be used as a template for submesh creation
    TetrahedralMesh<DIM, DIM>* original_triangular_mesh = static_cast<PottsArbitrarySurfaceIn3DMesh<DIM>*>(&potts_cell_population->rGetMesh())->GetDelaunayMesh();

    // Create the new set of submeshes to be used for FEM in the current time step.
    assert(potts_cell_population->GetNumRealCells() == potts_cell_population->GetNumElements());
    mpFeMesh.resize(rCellPopulation.GetNumRealCells());
    for (unsigned cell_number=0; cell_number<potts_cell_population->GetNumRealCells(); ++cell_number)
    {
        PottsElement<DIM>* potts_element = potts_cell_population->GetElement(cell_number);

        // Loop over the nodes in the Potts element, they correspond to triangles/tets in the original mesh
        std::vector<unsigned> element_subset;
        for (unsigned potts_node_number=0; potts_node_number<potts_element->GetNumNodes(); ++potts_node_number)
        {
            element_subset.push_back(potts_element->GetNode(potts_node_number)->GetIndex());
        }

        mpFeMesh[cell_number] = new TetrahedralSubsetMesh<DIM,DIM>(*original_triangular_mesh, element_subset);
    }

    CPMGTPaseEventHandler::EndEvent(CPMGTPaseEventHandler::SETUP_MESH);
}

template<unsigned DIM>
void GTPasePDEMultipleMeshesModifier<DIM>::UpdateSolutionVector(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    CPMGTPaseEventHandler::BeginEvent(CPMGTPaseEventHandler::UPDATE_SOLUTION);

    bool first_timestep = (mSolution.size() == 0);

    // Holds true in the first time step
    if (first_timestep)
    {
        // No previous meshes should be available.
        assert(mPreviousFeMesh.size() == 0);

        for (unsigned cell_number=0; cell_number<rCellPopulation.GetNumRealCells(); ++cell_number)
        {
            mSolution.push_back(GTPasePDESystemSolver<DIM, DIM>::CreateInitialCondition(mpFeMesh[cell_number]));
        }
    }
    else
    {
        // Previous meshes should be available.
        assert(mPreviousFeMesh.size() != 0);

        for (unsigned cell_number=0; cell_number<rCellPopulation.GetNumRealCells(); ++cell_number)
        {
            if (HasCellChangedShape(cell_number))
            {
                Vec previous_solution = mSolution[cell_number];
                mSolution[cell_number] = InterpolateSolutionBetweenMeshes(mpFeMesh[cell_number], mPreviousFeMesh[cell_number], previous_solution);
                PetscTools::Destroy(previous_solution);
            }
        }
    }

    CPMGTPaseEventHandler::EndEvent(CPMGTPaseEventHandler::UPDATE_SOLUTION);
}

template<unsigned DIM>
bool GTPasePDEMultipleMeshesModifier<DIM>::HasCellChangedShape(unsigned cellNumber)
{
    return mpFeMesh[cellNumber]->GetOriginalNodeIndices() != mPreviousFeMesh[cellNumber]->GetOriginalNodeIndices();
}

template<unsigned DIM>
Vec GTPasePDEMultipleMeshesModifier<DIM>::InterpolateSolutionBetweenMeshes(TetrahedralSubsetMesh<DIM,DIM>* newMesh, TetrahedralSubsetMesh<DIM,DIM>* oldMesh, Vec previousSolution)
{
    // Initialise the vector to DBL_MAX. This is required for the interporlation algorithm below to work properly.
    Vec interpolated_solution = PetscTools::CreateAndSetVec(GTPASE_PROBLEM_DIM*newMesh->GetNumNodes(), DBL_MAX);

    const std::map<unsigned, unsigned>& new_global_to_local_map = ComputeGlobalToLocalMap(newMesh);
    const std::map<unsigned, unsigned>& old_global_to_local_map = ComputeGlobalToLocalMap(oldMesh);

    std::vector<unsigned> common_nodes;
    std::vector<unsigned> new_nodes;
    GetNodeSetIntersectionAndDifference(newMesh, oldMesh, common_nodes, new_nodes);

    // For all nodes that are shared between the old and new mesh, copy values over
    for (std::vector<unsigned>::const_iterator node_index_iter = common_nodes.begin(); node_index_iter != common_nodes.end(); ++node_index_iter)
    {
        unsigned local_index_new = new_global_to_local_map.at(*node_index_iter);
        unsigned local_index_old = old_global_to_local_map.at(*node_index_iter);

        bool same_boundary_status = (newMesh->GetNode(local_index_new)->IsBoundaryNode() == oldMesh->GetNode(local_index_old)->IsBoundaryNode());

        unsigned first_entry_new = local_index_new * GTPASE_PROBLEM_DIM;
        unsigned first_entry_old = local_index_old * GTPASE_PROBLEM_DIM;

        for (unsigned var_index=0; var_index<GTPASE_PROBLEM_DIM; ++var_index)
        {
            if ((var_index < P0) || (var_index > P4) || same_boundary_status)
            {
                // It isn't a Pi variable or if it is, its node has the same boundary status in both the old and new meshes
                PetscVecTools::SetElement(interpolated_solution, first_entry_new+var_index, PetscVecTools::GetElement(previousSolution, first_entry_old+var_index));
            }
            else
            {
                // It's a Pi variable from a node with different boundary status in both the old and new meshes

                // Node stopped being boundary, set Pi to zero
                if (!newMesh->GetNode(local_index_new)->IsBoundaryNode())
                {
                    PetscVecTools::SetElement(interpolated_solution, first_entry_new+var_index, 0.0);
                }

                // Node became boundary, interpolate a suitable value from neighbouring boundary nodes
                if (!oldMesh->GetNode(local_index_old)->IsBoundaryNode())
                {
                    const std::set<unsigned>& neighs = oldMesh->GetNeighbouringNodes(local_index_old);

                    double new_value = 0.0;
                    for (std::set<unsigned>::const_iterator neigh = neighs.begin(); neigh != neighs.end(); ++neigh)
                    {
                        if (oldMesh->GetNode(*neigh)->IsBoundaryNode())
                        {
                            new_value = PetscVecTools::GetElement(previousSolution, *neigh*GTPASE_PROBLEM_DIM + var_index);
                            break;
                        }
                    }
                    //assert(new_value != 0.0);
                    if (new_value == 0.0)
                    {
                        new_value = BasalEdgeBarbedEndDensity;
                    }

                    PetscVecTools::SetElement(interpolated_solution, first_entry_new+var_index, new_value);
                }
            }
        }
    }

    // For all the nodes that are new in the new mesh, perform interpolation while taking care of conserving mass

    // new_nodes uses node numbering according to the template mesh, we need it to obtain the equivalent local indices in the new mesh.
    // We need to visit the local nodes in ascending order to avoid having to deal with a node that has no information in the neighbourhood
    // to interpolate from. Otherwise we would trip the assertion in line 319. Creating a set with the local indices will allow us to
    // easily iterate in ascending order.
    std::set<unsigned> new_nodes_local_indices;
    for (std::vector<unsigned>::const_iterator node_index_iter = new_nodes.begin(); node_index_iter != new_nodes.end(); ++node_index_iter)
    {
        new_nodes_local_indices.insert(new_global_to_local_map.at(*node_index_iter));
    }

    for (std::set<unsigned>::const_iterator node_index_iter = new_nodes_local_indices.begin(); node_index_iter != new_nodes_local_indices.end(); ++node_index_iter)
    {
        unsigned local_index = *node_index_iter;
        unsigned first_entry_new = local_index * GTPASE_PROBLEM_DIM;

        const std::set<unsigned>& neighbour_nodes = newMesh->GetNeighbouringNodes(local_index);

        for (unsigned var_index=0; var_index<GTPASE_PROBLEM_DIM; ++var_index)
        {
            // For barbed ends and pushing barbed ends the new nodes are set to 0
            if ((var_index >= B0 && var_index <= B4)
                    || (var_index >= P0 && var_index <= P4))
            {
                PetscVecTools::SetElement(interpolated_solution, first_entry_new+var_index, 0.0);
                continue;
            }

            std::vector<double> non_zero_neigh_values;

            for (std::set<unsigned>::const_iterator neighbour = neighbour_nodes.begin(); neighbour != neighbour_nodes.end(); ++neighbour)
            {
                unsigned neigh_first_entry = (*neighbour) * GTPASE_PROBLEM_DIM;

                double vec_entry = PetscVecTools::GetElement(interpolated_solution, neigh_first_entry+var_index);

                // A node may have as a neighbour other nodes without available values (because they have been interpolated yet).
                if ((vec_entry != DBL_MAX) && (vec_entry != 0.0))
                {
                    non_zero_neigh_values.push_back(vec_entry);
                }
            }

            if (non_zero_neigh_values.size() != 0)
            {
                double mean = std::accumulate(non_zero_neigh_values.begin(), non_zero_neigh_values.end(), 0.0) / non_zero_neigh_values.size();
                PetscVecTools::SetElement(interpolated_solution, first_entry_new+var_index, mean);
            }
            else
            {
                std::cout << "WARNING: node " << local_index  << " (of " << new_nodes_local_indices.size() << " new nodes) has no valid neighboring nodes to interpolate from." << std::endl;
                PetscVecTools::SetElement(interpolated_solution, first_entry_new+var_index, 0.0);
            }
        }

    }

    EnforceMassConservation(oldMesh, newMesh, previousSolution, interpolated_solution);

    PetscVecTools::Finalise(interpolated_solution);

    return interpolated_solution;
}

template<unsigned DIM>
void GTPasePDEMultipleMeshesModifier<DIM>::EnforceMassConservation(TetrahedralSubsetMesh<DIM,DIM>* previousMesh, TetrahedralSubsetMesh<DIM,DIM>* newMesh, const Vec previousSolution, Vec newSolution) const
{
    DistributedVectorFactory* p_previous_factory = previousMesh->GetDistributedVectorFactory();
    DistributedVector distributed_previous_solution = p_previous_factory->CreateDistributedVector(previousSolution, true);

    DistributedVectorFactory* p_new_factory = newMesh->GetDistributedVectorFactory();
    DistributedVector distributed_new_solution = p_new_factory->CreateDistributedVector(newSolution);

    for (unsigned variable = 0; variable < GTPASE_PROBLEM_DIM; ++variable)
    {
        // The stripe objects will allow access to each variable independently using an iterator
        DistributedVector::Stripe previous_solution_variable(distributed_previous_solution, variable);
        DistributedVector::Stripe new_solution_variable(distributed_new_solution, variable);

        double previous_total_mass = 0.0;
        double new_total_mass = 0.0;

        // Each iteration corresponds to one mesh node with the GTPASE_PROBLEM_DIM dof accessible via the stripe objects
        for (DistributedVector::Iterator index = distributed_previous_solution.Begin();
             index!= distributed_previous_solution.End();
             ++index)
        {
            previous_total_mass += previous_solution_variable[index];
        }

        for (DistributedVector::Iterator index = distributed_new_solution.Begin();
             index!= distributed_new_solution.End();
             ++index)
        {
            new_total_mass += new_solution_variable[index];
        }

        double normalisation_factor = previous_total_mass / new_total_mass;

        for (DistributedVector::Iterator index = distributed_new_solution.Begin();
             index!= distributed_new_solution.End();
             ++index)
        {
            new_solution_variable[index] *= normalisation_factor;
        }
    }

    distributed_new_solution.Restore();
}

template<unsigned DIM>
const std::map<unsigned, unsigned> GTPasePDEMultipleMeshesModifier<DIM>::ComputeGlobalToLocalMap(TetrahedralSubsetMesh<DIM,DIM>* mesh)
{
    const std::vector<unsigned>& original_node_indices = mesh->GetOriginalNodeIndices();

    // Construct a map that is the inverse of original_node_indices
    std::map<unsigned, unsigned> global_to_local_map;
    for (unsigned node_index = 0; node_index < original_node_indices.size(); ++node_index)
    {
        global_to_local_map[original_node_indices[node_index]] = node_index;
    }

    return global_to_local_map;
}

template<unsigned DIM>
void GTPasePDEMultipleMeshesModifier<DIM>::GetNodeSetIntersectionAndDifference(TetrahedralSubsetMesh<DIM,DIM>* newMesh, TetrahedralSubsetMesh<DIM,DIM>* oldMesh, std::vector<unsigned>& commonNodes, std::vector<unsigned>& newNodes)
{
    // We have to make copies of these vectors since we need to sort them
    std::vector<unsigned> new_mesh_set = newMesh->GetOriginalNodeIndices();
    std::vector<unsigned> old_mesh_set = oldMesh->GetOriginalNodeIndices();

    // Set algebra below requires sorted STL containers
    std::sort(new_mesh_set.begin(), new_mesh_set.end());
    std::sort(old_mesh_set.begin(), old_mesh_set.end());

    // Upper bound for intersection and difference cardinality
    unsigned max_set_size = std::max(new_mesh_set.size(), old_mesh_set.size());

    commonNodes.resize(max_set_size);
    std::vector<unsigned>::iterator iter = std::set_intersection(new_mesh_set.begin(), new_mesh_set.end(), old_mesh_set.begin(), old_mesh_set.end(), commonNodes.begin());
    commonNodes.resize(iter - commonNodes.begin());

    newNodes.resize(max_set_size);
    iter = std::set_difference(new_mesh_set.begin(), new_mesh_set.end(), old_mesh_set.begin(), old_mesh_set.end(), newNodes.begin());
    newNodes.resize(iter - newNodes.begin());
}


template<unsigned DIM>
void GTPasePDEMultipleMeshesModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    CPMGTPaseEventHandler::BeginEvent(CPMGTPaseEventHandler::UPDATE_CELL_DATA);

    PottsBasedCellPopulation<DIM>* potts_cell_population = dynamic_cast<PottsBasedCellPopulation<DIM>*>(&rCellPopulation);
    assert(potts_cell_population);

    // local cell index used by the CA simulation
    unsigned cell_index = 0;

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                 cell_iter != rCellPopulation.End();
                 ++cell_iter)
    {
        if (!cell_iter->HasCellVecData())
        {
            boost::shared_ptr<AbstractCellProperty> p_vec_data(new CellVecData());
            cell_iter->AddCellProperty(p_vec_data);
        }

        Vec solution_at_each_element = AverageNodeWiseSolutionAtElements(mSolution[cell_index], mpFeMesh[cell_index]);
        cell_iter->GetCellVecData()->SetItem(mDependentVariableName, solution_at_each_element);

        const std::vector<unsigned>& element_indices = mpFeMesh[cell_index]->GetOriginalElementIndices();
        cell_iter->GetCellVecData()->SetItem("GTPasePDEElementIndices", PetscTools::CreateVec(std::vector<double>(element_indices.begin(), element_indices.end())));

        cell_index++;
    }

    CPMGTPaseEventHandler::EndEvent(CPMGTPaseEventHandler::UPDATE_CELL_DATA);
}

template<unsigned DIM>
void GTPasePDEMultipleMeshesModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}


template<unsigned DIM>
Vec GTPasePDEMultipleMeshesModifier<DIM>::AverageNodeWiseSolutionAtElements(const Vec nodeWiseSolution, TetrahedralSubsetMesh<DIM, DIM>* pMesh) const
{
    assert(PetscVecTools::GetSize(nodeWiseSolution)/GTPASE_PROBLEM_DIM == pMesh->GetNumNodes());

    Vec element_wise_solution = PetscTools::CreateVec(pMesh->GetNumElements()*GTPASE_PROBLEM_DIM);

    for(typename TetrahedralSubsetMesh<DIM, DIM>::ElementIterator elem_iter = pMesh->GetElementIteratorBegin();
        elem_iter != pMesh->GetElementIteratorEnd(); ++elem_iter)
    {
        Element<DIM, DIM>& element = *elem_iter;

        for (unsigned var_index=0; var_index<GTPASE_PROBLEM_DIM; ++var_index)
        {
            std::vector<double> node_values;
            for (unsigned node_local_index = 0; node_local_index<element.GetNumNodes(); ++node_local_index)
            {
                unsigned node_global_index = element.GetNode(node_local_index)->GetIndex();
                node_values.push_back(PetscVecTools::GetElement(nodeWiseSolution, node_global_index*GTPASE_PROBLEM_DIM + var_index));
            }

            double mean = std::accumulate(node_values.begin(), node_values.end(), 0.0) / node_values.size();
            PetscVecTools::SetElement(element_wise_solution, element.GetIndex()*GTPASE_PROBLEM_DIM + var_index, mean);
        }
    }

    return element_wise_solution;
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

//template class GTPasePDEMultipleMeshesModifier<1>;
template class GTPasePDEMultipleMeshesModifier<2>;
//template class GTPasePDEMultipleMeshesModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
//EXPORT_TEMPLATE_CLASS_SAME_DIMS(GTPasePDEMultipleMeshesModifier)
EXPORT_TEMPLATE_CLASS1(GTPasePDEMultipleMeshesModifier, 2)
