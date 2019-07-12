#include "InitialConditionsGenerator.hpp"
#include "PetscVecTools.hpp"
#include "PetscMatTools.hpp"
#include "OutputFileHandler.hpp"
#include "Debug.hpp"
#include "PetscException.hpp"

InitialConditionsGenerator::InitialConditionsGenerator()
{

}

void InitialConditionsGenerator::SetInitialRestLengths(MeshBasedCellPopulation<2,3>& rCellPopulation, double springConstant)
{
    // Store the edges and current rest lengths
    std::vector< std::pair<unsigned, unsigned > > edges;
    std::vector< double > rest_lengths;

    for (MeshBasedCellPopulation<2,3>::SpringIterator spring_iterator = rCellPopulation.SpringsBegin();
                 spring_iterator != rCellPopulation.SpringsEnd();
                 ++spring_iterator)
    {
        // Note that nodeA_global_index is always less than nodeB_global_index
        Node<3>* p_nodeA = spring_iterator.GetNodeA();
        Node<3>* p_nodeB = spring_iterator.GetNodeB();

        unsigned nodeA_global_index = p_nodeA->GetIndex();
        unsigned nodeB_global_index = p_nodeB->GetIndex();

        double current_rest_length = rCellPopulation.GetRestLength(nodeA_global_index,nodeB_global_index);
        rest_lengths.push_back(current_rest_length);

        std::pair<unsigned,unsigned> node_pair = rCellPopulation.CreateOrderedPair(nodeA_global_index,nodeB_global_index);

        edges.push_back(node_pair);
    }
    unsigned num_edges = edges.size();


    //Set up the linear system to calculate the least squares fit
    unsigned num_nodes = rCellPopulation.GetNumNodes();
    unsigned num_equations = 3 * num_nodes;
    unsigned num_variables = num_edges;

    Mat lsq_matrix;
    // \todo: currently preallocating enough memory for a dense matrix
    unsigned max_num_nonzeros_per_row = num_variables;
    PetscTools::SetupMat(lsq_matrix, num_equations, num_variables, max_num_nonzeros_per_row);

    Vec lsq_rhs = PetscTools::CreateAndSetVec(num_equations, 0.0);

    for (unsigned edge_index = 0; edge_index<num_edges; edge_index++)
    {
        std::pair<unsigned,unsigned> edge =  edges[edge_index];

        Node<3>* p_nodeA = rCellPopulation.GetNode(edge.first);
        Node<3>* p_nodeB = rCellPopulation.GetNode(edge.second);

        double voronoi_cell_area_A = rCellPopulation.GetVolumeOfCell(rCellPopulation.GetCellUsingLocationIndex(edge.first));
        double voronoi_cell_area_B = rCellPopulation.GetVolumeOfCell(rCellPopulation.GetCellUsingLocationIndex(edge.second));

        c_vector<double, 3> r_ij = p_nodeB->rGetLocation() -  p_nodeA->rGetLocation();
        c_vector<double, 3> r_ij_hat = r_ij/norm_2(r_ij);

        {
            PetscInt rows[3] = {edge.first, num_nodes + edge.first, 2 * num_nodes + edge.first};
            c_vector<double, 3> matrix_coefficients =  springConstant * r_ij_hat;
            MatSetValues(lsq_matrix, 3, rows, 1, (PetscInt*) &edge_index, matrix_coefficients.data(), INSERT_VALUES);
        }

        {
            PetscInt rows[3] = {edge.second, num_nodes + edge.second, 2 * num_nodes + edge.second};
            c_vector<double, 3> matrix_coefficients = - springConstant * r_ij_hat;
            MatSetValues(lsq_matrix, 3, rows, 1, (PetscInt*) &edge_index, matrix_coefficients.data(), INSERT_VALUES);
        }
    }

    //Loop over nodes (i.e cells) to construct the rhs vector
    for (AbstractCellPopulation<2,3>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        unsigned rows[3] = {node_index, num_nodes + node_index, 2 * num_nodes + node_index};
        c_vector<double,3> applied_force;
        applied_force(0) = -cell_iter->GetCellData()->GetItem("applied_force_x");
        applied_force(1) = -cell_iter->GetCellData()->GetItem("applied_force_y");
        applied_force(2) = -cell_iter->GetCellData()->GetItem("applied_force_z");
        PetscVecTools::AddMultipleValues(lsq_rhs, rows, applied_force);
    }

    PetscMatTools::Finalise(lsq_matrix);
    PetscVecTools::Finalise(lsq_rhs);

    // Solve the least squares problem
    PetscOptionsSetValue("-ksp_monitor","");
    // PetscOptionsSetValue("-pc_hypre_boomeramg_max_iter", "1");
    // PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold", "0.0");
    KSP mKspSolver;
    KSPCreate(PETSC_COMM_WORLD, &mKspSolver);
    KSPSetType(mKspSolver, KSPLSQR);
    KSPSetTolerances(mKspSolver, 1e-4, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetFromOptions(mKspSolver);

    // With the PETSc built-in preconditioners, one should call KSPSetOperators(ksp,A,A'*A)) since the preconditioner needs to work for the normal equations A'*A.
    Mat lsq_precond_matrix;
    MatTransposeMatMult(lsq_matrix, lsq_matrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &lsq_precond_matrix);
    KSPSetOperators(mKspSolver, lsq_matrix, lsq_precond_matrix);

    // Use point Jacobi since the lhs matrix of the normal equations may be singular and factorisation-based preconditioners will fail.
    PC prec;
    KSPGetPC(mKspSolver, &prec);
    PCSetType(prec, PCHYPRE);

    KSPSetUp(mKspSolver);
    Vec lsq_solution = PetscTools::CreateVec(num_variables);
    KSPSolve(mKspSolver, lsq_rhs, lsq_solution);

    // KSPConvergedReason is an enum where positive values correspond to different convergence criteria being met
    KSPConvergedReason reason;
    KSPGetConvergedReason(mKspSolver, &reason);
    PRINT_VARIABLE(reason);
    assert(reason>0);

    // Save the rest lengths
    for (unsigned edge_index = 0; edge_index<num_edges; edge_index++)
    {
        double new_rest_length = rest_lengths[edge_index] - PetscVecTools::GetElement(lsq_solution, edge_index);
        assert(new_rest_length>0);

        rCellPopulation.SetRestLength(edges[edge_index].first, edges[edge_index].second, new_rest_length);
    }

    PetscTools::Destroy(lsq_matrix);
    PetscTools::Destroy(lsq_precond_matrix);
    PetscTools::Destroy(lsq_solution);
    PetscTools::Destroy(lsq_rhs);
}

void InitialConditionsGenerator::SetApproximateInitialRestLengths(MeshBasedCellPopulation<2,3>& rCellPopulation, double springConstant)
{
    // Store the edges and current rest lengths
    std::vector< std::pair<unsigned, unsigned > > edges;
    std::vector< double > rest_lengths;

    for (MeshBasedCellPopulation<2,3>::SpringIterator spring_iterator = rCellPopulation.SpringsBegin();
                 spring_iterator != rCellPopulation.SpringsEnd();
                 ++spring_iterator)
    {
        // Note that nodeA_global_index is always less than nodeB_global_index
        Node<3>* p_nodeA = spring_iterator.GetNodeA();
        Node<3>* p_nodeB = spring_iterator.GetNodeB();

        unsigned nodeA_global_index = p_nodeA->GetIndex();
        unsigned nodeB_global_index = p_nodeB->GetIndex();

        double current_rest_length = rCellPopulation.GetRestLength(nodeA_global_index,nodeB_global_index);
        rest_lengths.push_back(current_rest_length);

        std::pair<unsigned,unsigned> node_pair = rCellPopulation.CreateOrderedPair(nodeA_global_index,nodeB_global_index);

        edges.push_back(node_pair);
    }
    unsigned num_edges = edges.size();

    //Set up the linear system to calculate the extensions
    unsigned num_nodes = rCellPopulation.GetNumNodes();

    Mat matrix;
    // \todo: currently preallocating enough memory for a dense matrix
    unsigned max_num_nonzeros_per_row = num_nodes;
    PetscTools::SetupMat(matrix, num_nodes, num_nodes, max_num_nonzeros_per_row);

    Vec rhs = PetscTools::CreateAndSetVec(num_nodes, 0.0);

    for (unsigned edge_index = 0; edge_index<num_edges; edge_index++)
    {
        std::pair<unsigned,unsigned> edge =  edges[edge_index];

        Node<3>* p_nodeA = rCellPopulation.GetNode(edge.first);
        Node<3>* p_nodeB = rCellPopulation.GetNode(edge.second);

        double voronoi_cell_area_A = rCellPopulation.GetVolumeOfCell(rCellPopulation.GetCellUsingLocationIndex(edge.first));
        double voronoi_cell_area_B = rCellPopulation.GetVolumeOfCell(rCellPopulation.GetCellUsingLocationIndex(edge.second));

        c_vector<double, 3> r_ij = p_nodeB->rGetLocation() -  p_nodeA->rGetLocation();
        c_vector<double, 3> r_ij_hat = r_ij/norm_2(r_ij);


        c_vector<double,3> normal_A;
        normal_A(0) = rCellPopulation.GetCellUsingLocationIndex(edge.first)->GetCellData()->GetItem("applied_force_x");
        normal_A(1) = rCellPopulation.GetCellUsingLocationIndex(edge.first)->GetCellData()->GetItem("applied_force_y");
        normal_A(2) = rCellPopulation.GetCellUsingLocationIndex(edge.first)->GetCellData()->GetItem("applied_force_z");
        normal_A /= norm_2(normal_A);

        c_vector<double,3> normal_B;
        normal_B(0) = rCellPopulation.GetCellUsingLocationIndex(edge.second)->GetCellData()->GetItem("applied_force_x");
        normal_B(1) = rCellPopulation.GetCellUsingLocationIndex(edge.second)->GetCellData()->GetItem("applied_force_y");
        normal_B(2) = rCellPopulation.GetCellUsingLocationIndex(edge.second)->GetCellData()->GetItem("applied_force_z");
        normal_B /= norm_2(normal_B);

        {
            PetscInt rows[2] = {edge.first,edge.second};
            c_vector<double, 2> matrix_coefficients;
            matrix_coefficients(0) =  -springConstant * inner_prod(r_ij_hat,normal_A);
            matrix_coefficients(1) =  -springConstant * inner_prod(r_ij_hat,normal_A);
            MatSetValues(matrix, 2, rows, 1, (PetscInt*) &(edge.first), matrix_coefficients.data(), ADD_VALUES);
        }
        {
            PetscInt rows[2] = {edge.second,edge.first};
            c_vector<double, 2> matrix_coefficients;
            matrix_coefficients(0) =  springConstant * inner_prod(r_ij_hat,normal_B);
            matrix_coefficients(1) =  springConstant * inner_prod(r_ij_hat,normal_B);
            MatSetValues(matrix, 2, rows, 1, (PetscInt*) &(edge.second), matrix_coefficients.data(), ADD_VALUES);
        }
    }

    //Loop over nodes (i.e cells) to construct the rhs vector
    for (AbstractCellPopulation<2,3>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        unsigned rows[1] = {node_index};
        c_vector<double,3> applied_force;
        applied_force(0) = cell_iter->GetCellData()->GetItem("applied_force_x");
        applied_force(1) = cell_iter->GetCellData()->GetItem("applied_force_y");
        applied_force(2) = cell_iter->GetCellData()->GetItem("applied_force_z");
        c_vector<double,1> vector_coeficients;
        vector_coeficients(0) = norm_2(applied_force);

        PetscVecTools::AddMultipleValues(rhs, rows, vector_coeficients);
    }

    PetscMatTools::Finalise(matrix);
    PetscVecTools::Finalise(rhs);

    // Solve the linear system
    PetscOptionsSetValue("-ksp_monitor_true_residual","");
    KSP mKspSolver;
    KSPCreate(PETSC_COMM_WORLD, &mKspSolver);
    KSPSetTolerances(mKspSolver, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetFromOptions(mKspSolver);

    KSPSetOperators(mKspSolver, matrix, matrix);

    PC prec;
    KSPGetPC(mKspSolver, &prec);
    PCSetType(prec, PCHYPRE);

    KSPSetUp(mKspSolver);
    Vec solution = PetscTools::CreateVec(num_nodes);
    KSPSolve(mKspSolver, rhs, solution);

    // KSPConvergedReason is an enum where positive values correspond to different convergence criteria being met
    KSPConvergedReason reason;
    KSPGetConvergedReason(mKspSolver, &reason);
    PRINT_VARIABLE(reason);
    assert(reason>0);

    // Save the rest lengths
    for (unsigned edge_index = 0; edge_index<num_edges; edge_index++)
    {
        std::pair<unsigned,unsigned> edge =  edges[edge_index];

        double new_rest_length = rest_lengths[edge_index] - PetscVecTools::GetElement(solution, edge.first)
                                                          - PetscVecTools::GetElement(solution, edge.second);

        PRINT_4_VARIABLES(rest_lengths[edge_index],PetscVecTools::GetElement(solution, edge.first),PetscVecTools::GetElement(solution, edge.second), new_rest_length);
        PRINT_4_VARIABLES(edge.first,edge.second, rest_lengths[edge_index],new_rest_length);

        //assert(new_rest_length>0);

        rCellPopulation.SetRestLength(edge.first, edge.second, new_rest_length);
    }

    PetscViewer mviewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD,"file.mat",&mviewer);
    PetscViewerSetFormat(mviewer, PETSC_VIEWER_ASCII_MATLAB);
    VecView(rhs, mviewer);
    MatView(matrix, mviewer);
    //VecView(solution, mviewer);


    PetscTools::Destroy(matrix);
    PetscTools::Destroy(solution);
    PetscTools::Destroy(rhs);
}

void InitialConditionsGenerator::SetSimpleInitialRestLengths(MeshBasedCellPopulation<2,3>& rCellPopulation, double springConstant)
{
    // Store the edges and current rest lengths
    std::vector< std::pair<unsigned, unsigned > > edges;
    std::vector< double > rest_lengths;

    for (MeshBasedCellPopulation<2,3>::SpringIterator spring_iterator = rCellPopulation.SpringsBegin();
                 spring_iterator != rCellPopulation.SpringsEnd();
                 ++spring_iterator)
    {
        // Note that nodeA_global_index is always less than nodeB_global_index
        Node<3>* p_nodeA = spring_iterator.GetNodeA();
        Node<3>* p_nodeB = spring_iterator.GetNodeB();

        unsigned nodeA_global_index = p_nodeA->GetIndex();
        unsigned nodeB_global_index = p_nodeB->GetIndex();

        if (nodeA_global_index == 4)
        {
            PRINT_2_VARIABLES(nodeA_global_index, nodeB_global_index);
            PRINT_3_VARIABLES(p_nodeA->rGetLocation()[0],p_nodeA->rGetLocation()[1],p_nodeA->rGetLocation()[2]);
            PRINT_3_VARIABLES(p_nodeB->rGetLocation()[0],p_nodeB->rGetLocation()[1],p_nodeB->rGetLocation()[2]);
        }
        if (nodeB_global_index == 4)
        {
            PRINT_2_VARIABLES(nodeA_global_index, nodeB_global_index);
            PRINT_3_VARIABLES(p_nodeA->rGetLocation()[0],p_nodeA->rGetLocation()[1],p_nodeA->rGetLocation()[2]);
            PRINT_3_VARIABLES(p_nodeB->rGetLocation()[0],p_nodeB->rGetLocation()[1],p_nodeB->rGetLocation()[2]);
        }


        double current_rest_length = rCellPopulation.GetRestLength(nodeA_global_index,nodeB_global_index);
        rest_lengths.push_back(current_rest_length);

        std::pair<unsigned,unsigned> node_pair = rCellPopulation.CreateOrderedPair(nodeA_global_index,nodeB_global_index);

        edges.push_back(node_pair);
    }
    unsigned num_edges = edges.size();

    //Set up the linear system to calculate the extensions
    unsigned num_nodes = rCellPopulation.GetNumNodes();

    // store the calculated extension
    double extensions[num_nodes];
    double extension_coficient[num_nodes];
    for (unsigned i=0; i<num_nodes; i++)
    {
        extension_coficient[i] = 0.0;
    }

    for (unsigned edge_index = 0; edge_index<num_edges; edge_index++)
    {
        std::pair<unsigned,unsigned> edge =  edges[edge_index];

        Node<3>* p_nodeA = rCellPopulation.GetNode(edge.first);
        Node<3>* p_nodeB = rCellPopulation.GetNode(edge.second);

        c_vector<double, 3> r_ij = p_nodeB->rGetLocation() -  p_nodeA->rGetLocation();
        c_vector<double, 3> r_ij_hat = r_ij/norm_2(r_ij);


        c_vector<double,3> normal_A;
        normal_A(0) = rCellPopulation.GetCellUsingLocationIndex(edge.first)->GetCellData()->GetItem("applied_force_x");
        normal_A(1) = rCellPopulation.GetCellUsingLocationIndex(edge.first)->GetCellData()->GetItem("applied_force_y");
        normal_A(2) = rCellPopulation.GetCellUsingLocationIndex(edge.first)->GetCellData()->GetItem("applied_force_z");
        normal_A /= norm_2(normal_A);

        c_vector<double,3> normal_B;
        normal_B(0) = rCellPopulation.GetCellUsingLocationIndex(edge.second)->GetCellData()->GetItem("applied_force_x");
        normal_B(1) = rCellPopulation.GetCellUsingLocationIndex(edge.second)->GetCellData()->GetItem("applied_force_y");
        normal_B(2) = rCellPopulation.GetCellUsingLocationIndex(edge.second)->GetCellData()->GetItem("applied_force_z");
        normal_B /= norm_2(normal_B);

        extension_coficient[edge.first] += -springConstant * inner_prod(r_ij_hat,normal_A);
        extension_coficient[edge.second] += springConstant * inner_prod(r_ij_hat,normal_B);
    }

    //Loop over nodes (i.e cells) to construct the rhs vector
    for (AbstractCellPopulation<2,3>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        c_vector<double,3> applied_force;
        applied_force(0) = cell_iter->GetCellData()->GetItem("applied_force_x");
        applied_force(1) = cell_iter->GetCellData()->GetItem("applied_force_y");
        applied_force(2) = cell_iter->GetCellData()->GetItem("applied_force_z");

        extensions[node_index] = norm_2(applied_force)/extension_coficient[node_index];

        if (node_index == 4)
        {
            PRINT_3_VARIABLES(applied_force(0),applied_force(1),applied_force(2));
            PRINT_VARIABLE(extensions[node_index]);
        }


    }

    // Save the rest lengths
    for (unsigned edge_index = 0; edge_index<num_edges; edge_index++)
    {
        std::pair<unsigned,unsigned> edge =  edges[edge_index];

        double common_extension =  0.5*extensions[edge.first] + 0.5*extensions[edge.second];

        double new_rest_length = rest_lengths[edge_index]- common_extension;

        PRINT_4_VARIABLES(rest_lengths[edge_index],extensions[edge.first],extensions[edge.second], new_rest_length);
        PRINT_4_VARIABLES(edge.first,edge.second, rest_lengths[edge_index],new_rest_length);

        assert(new_rest_length>0);

        rCellPopulation.SetRestLength(edge.first, edge.second, new_rest_length);
    }
}
