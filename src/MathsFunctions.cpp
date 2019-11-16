#include "MathsFunctions.hpp"
#include <math.h>
#include "Debug.hpp"


// Need to include this to use the functions in this class
// MathsFunctions maths_functions = MathsFunctions();
// c_vector<c_vector<long double, 3>, 3> NewVec =  maths_functions.TheNewMatrixMultiplication(Matrix2, Matrix1);
 

c_vector<c_vector<long double, 3>, 3> MathsFunctions::MatrixMultiplication(c_vector<c_vector<long double, 3>, 3> Matrix1, c_vector<c_vector<long double, 3>, 3> Matrix2)
{

    c_vector<c_vector<long double, 3>, 3> MatrixTranspose;
    c_vector<c_vector<long double, 3>, 3> Answer;

    // This will give us the determinat
    for (int i = 0; i < 3; i++)
    {
        MatrixTranspose[i] = Create_c_vector(Matrix2[0][i], Matrix2[1][i], Matrix2[2][i]);
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Answer[i][j] = inner_prod(Matrix1[i], MatrixTranspose[j]);
        }
    }

    return Answer;
}

c_vector<long double, 3> MathsFunctions::MatrixMultiplication(c_vector<c_vector<long double, 3>, 3> Matrix, c_vector<long double, 3> Vector)
{

    c_vector<long double, 3> Answer;


        for (int i = 0; i < 3; i++)
        {
            Answer[i] = inner_prod(Matrix[i], Vector);
        }

    return Answer;
}

std::pair<double, double> MathsFunctions::Create_pair(double x, double y)
{
    std::pair<double, double> NodePair = std::make_pair(std::min(x, y), std::max(x, y));
    return NodePair;
}

std::pair<unsigned, unsigned> MathsFunctions::Create_pair(unsigned x, unsigned y)
{
    std::pair<unsigned, unsigned> NodePair = std::make_pair(std::min(x, y), std::max(x, y));
    return NodePair;
}


double MathsFunctions::MaintainOutwardsPointingNormal(c_vector<double, 3> Normal, c_vector<double, 3> x1)
{

    double direction = 1;
    c_vector<double, 2> Normal2D = Create_c_vector(Normal[0], Normal[1]);

    c_vector<double, 2> x12D = Create_c_vector(x1[0], x1[1]);

    c_vector<double, 2> Extension = Normal2D + x12D;

    double absExtensin = norm_2(Extension);
    double absx1 = norm_2(x12D);
    if (absExtensin < absx1) // meaning the normal points towards the center
    {
        direction = -1;
        //  Normal = -Normal;// reverse the normal so it points out
    }
    return direction;
}


void MathsFunctions::PRINT_PAIR(std::pair<unsigned, unsigned> Pair)
{
    std::cout << "DEBUG: Pair = [" << Pair.first << ", " << Pair.second << "] " <<  std::endl;
}




double MathsFunctions::AreVectorsSame(std::vector<unsigned> Vector1, std::vector<unsigned> Vector2)
{
    sort(Vector1.begin(), Vector1.end());
    sort(Vector2.begin(), Vector2.end());

    std::vector<unsigned> DifferenceVector;
    std::set_difference(Vector1.begin(), Vector1.end(), Vector2.begin(), Vector2.end(),
                        std::inserter(DifferenceVector, DifferenceVector.begin()));

    double difference = DifferenceVector.size();
    return difference;
}


double MathsFunctions::AreVectorsSame(std::vector<std::pair<unsigned, unsigned> > Vector1, std::vector<std::pair<unsigned, unsigned> > Vector2)
{
    sort(Vector1.begin(), Vector1.end());
    sort(Vector2.begin(), Vector2.end());

    std::vector<std::pair<unsigned, unsigned> > DifferenceVector;
    std::set_difference(Vector1.begin(), Vector1.end(), Vector2.begin(), Vector2.end(),
                        std::inserter(DifferenceVector, DifferenceVector.begin()));

    double difference = DifferenceVector.size();
    return difference;
}


std::vector<unsigned> MathsFunctions::Intersection(std::vector<unsigned> Vector1, std::vector<unsigned> Vector2)
{
    std::vector<unsigned> IntersectionVector;

    for (std::vector<unsigned>::iterator it = Vector1.begin(); it != Vector1.end(); ++it)
    {
        for (std::vector<unsigned>::iterator it2 = Vector2.begin(); it2 != Vector2.end(); ++it2)
        {
            if (*it == *it2)
            {
                IntersectionVector.push_back(*it);
            }
        }
    }

    return IntersectionVector;
}


std::vector<unsigned> MathsFunctions::RemoveElement(std::vector<unsigned> Vector1, unsigned number)
{

    unsigned ElementToRemove;
    for (std::vector<unsigned>::iterator it = Vector1.begin(); it != Vector1.end(); ++it)
    {
        if (*it == number)
        {
            ElementToRemove = *it;
        }
    }
    Vector1.erase(Vector1.begin() + ElementToRemove);

    return Vector1;
}


std::vector<double> MathsFunctions::RemoveElement(std::vector<double> Vector1, double number)
{
    unsigned ElementToRemove;
    for (std::vector<double>::iterator it = Vector1.begin(); it != Vector1.end(); ++it)
    {
        if (*it == number)
        {
            ElementToRemove = *it;
        }
    }
    Vector1.erase(Vector1.begin() + ElementToRemove);

    return Vector1;
}


bool MathsFunctions::IsPairInVector(std::vector<std::pair<unsigned, unsigned> > Vector, std::pair<unsigned, unsigned> Pair)
{
    bool IsInVector = 0;

    for (std::vector<std::pair<unsigned, unsigned> >::iterator it = Vector.begin(); it != Vector.end(); ++it)
    {
        if (*it == Pair)
        {
            IsInVector = 1;
            return IsInVector;
        }
    
    }


    return IsInVector;
}

bool MathsFunctions::IsNumberInVector(std::vector<double> Vector, double number)
{
    bool IsInVector = 0;
    for (std::vector<double>::iterator it = Vector.begin(); it != Vector.end(); ++it)
    {
        if (*it == number)
        {
            IsInVector = 1;
        } 
    }
    return IsInVector;
}

bool MathsFunctions::IsNumberInVector(std::vector<unsigned> Vector, unsigned number)
{
  
    bool IsInVector = 0;
    for (std::vector<unsigned>::iterator it = Vector.begin(); it != Vector.end(); ++it)
    { 
        if (*it == number)
        {
            IsInVector = 1;
            return IsInVector;
        } 
    }
    return IsInVector;
}


double MathsFunctions::min_value(std::vector<double> Vector)
{
    double min = Vector[0];
    for (std::vector<double>::iterator it = Vector.begin(); it != Vector.end(); ++it)
    {
        /* std::cout << *it; ... */
        if (*it < min)
        {
            min = *it;
        } //else if (min + *it < -M_PI)//AddAngles(*it,min)
        // {
        //  // Have crossed the periodic boundary
        //     min = *it;
        // }
    }
    return min;
}


double MathsFunctions::max_value(std::vector<double> Vector)
{
    double max = Vector[0];
    for (std::vector<double>::iterator it = Vector.begin(); it != Vector.end(); ++it)
    {
        /* std::cout << *it; ... */
        if (*it > max)
        {
            max = *it;
        }
    }
    return max;
}



double MathsFunctions::AddAngles(double alpha, double beta)
{
    double eta = alpha + beta;

    if (eta > M_PI)
    { // This anlge is now negative
        eta = eta - 2 * M_PI;
    }
    else if (eta < -M_PI)
    { // This anlge is now positive
        eta = eta + 2 * M_PI;
    }
    return eta;
}


double MathsFunctions::PeriodicAngle(double eta)
{
    if (eta > M_PI)
    { // This anlge is now negative
        eta = eta - 2 * M_PI;
    }
    else if (eta < -M_PI)
    { // This anlge is now positive
        eta = eta + 2 * M_PI;
    }
    return eta;
}


bool MathsFunctions::IsVectorInVector(std::vector<c_vector<double, 3> > Vector, c_vector<double, 3> Location)
{
  
    bool IsInVector = 0;
    for (typename std::vector<c_vector<double, 3> >::iterator it = Vector.begin(); it != Vector.end(); ++it)
    {
        c_vector<double, 3> AVector = *it;
        if ((AVector[0] == Location[0]) && (AVector[1] == Location[1]))
        {
            //  IsInVector = 1;
            return 1;
        }
    }
    return IsInVector;
}












c_vector<c_vector<long double, 3>, 3> MathsFunctions::MatrixTranspose(c_vector<c_vector<long double, 3>, 3> Matrix)
{
    c_vector<c_vector<long double, 3>, 3> MatrixTranspose;
    for (int i = 0; i < 3; i++)
    { // Transporse of the matrix
        MatrixTranspose[i] = Create_c_vector(Matrix[0][i], Matrix[1][i], Matrix[2][i]);
    }
    return MatrixTranspose;
}

c_vector<c_vector<long double, 3>, 3> MathsFunctions::Inverse(c_vector<c_vector<long double, 3>, 3> Matrix)
{

    c_vector<c_vector<long double, 3>, 3> MTranspose = MatrixTranspose(Matrix);

    long double det = Matrix[0][0] * (Matrix[1][1] * Matrix[2][2] - Matrix[1][2] * Matrix[2][1]) - Matrix[0][1] * (Matrix[1][0] * Matrix[2][2] - Matrix[1][2] * Matrix[2][0]) + Matrix[0][2] * (Matrix[1][0] * Matrix[2][1] - Matrix[1][1] * Matrix[2][0]);

    double l = 0, m = 1;
    c_vector<c_vector<long double, 3>, 3> InverseMatrix;
    for (int k = 0; k < 3; k++)
    {
        l += 1, m += 1;
        if (l == 3)
        {
            l = 0;
        }
        if (m == 3)
        {
            m = 0;
        }
        InverseMatrix[k] = (1 / det) * VectorProduct(MTranspose[l], MTranspose[m]);
    }
    return InverseMatrix;
}

c_vector<c_vector<long double, 3>, 3> MathsFunctions::Inverse(c_vector<c_vector<long double, 3>, 3> Matrix, double elem)
{

    long double det = Matrix[0][0] * (Matrix[1][1] * Matrix[2][2] - Matrix[1][2] * Matrix[2][1]) - Matrix[0][1] * (Matrix[1][0] * Matrix[2][2] - Matrix[1][2] * Matrix[2][0]) + Matrix[0][2] * (Matrix[1][0] * Matrix[2][1] - Matrix[1][1] * Matrix[2][0]);

    c_vector<c_vector<long double, 3>, 3> MTranspose = MatrixTranspose(Matrix);

    double l = 0, m = 1;
    c_vector<c_vector<long double, 3>, 3> InverseMatrix;
    for (int k = 0; k < 3; k++)
    {
        l += 1, m += 1;
        if (l == 3)
        {
            l = 0;
        }
        if (m == 3)
        {
            m = 0;
        }
        InverseMatrix[k] = (1 / det) * VectorProduct(MTranspose[l], MTranspose[m]);
    }
    return InverseMatrix;
}

c_vector<c_vector<long double, 3>, 3> MathsFunctions::RowReduction(c_vector<c_vector<long double, 3>, 3> Matrix)
{
    // row reduction to reduced row eshioln form
    c_vector<c_vector<long double, 3>, 3> Identity;
    Identity[0] = Create_c_vector(1, 0, 0);
    Identity[1] = Create_c_vector(0, 1, 0);
    Identity[2] = Create_c_vector(0, 0, 1);
    for (int i = 0; i < 3; i++)
    {

        Identity[i] /= Matrix[i][i];
        Matrix[i] /= Matrix[i][i]; // sets all the diagonal elements to 1

        for (int j = 0; j < 3; j++)
        {
            if (j != i) //prevents taking out the t row that we just make 1
            {
                Identity[j] -= (Identity[i] * Matrix[j][i]);
            }
            if (j != i) //prevents taking out the t row that we just make 1
            {
                Matrix[j] -= (Matrix[i] * Matrix[j][i]);
            }
        }
    }

    // by this point the oridinal identity matrix is now the inverse
    return Identity;
}

// c_vector<c_vector<long double, 3>, 3> MathsFunctions::MappingMatrix(MeshBasedCellPopulation<2, 3>* p_cell_population, typename AbstractTetrahedralMesh<2, 3>::ElementIterator elem_iter, double a, double b, double theta)
// {
//     // Inverse Mapping Matrix
//     c_vector<c_vector<long double, 3>, 3> NewPositionVector;
//     c_vector<c_vector<long double, 3>, 3> PositionVector;

//     Node<3>* pNode0 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(0));
//     Node<3>* pNode1 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(1));
//     Node<3>* pNode2 = p_cell_population->rGetMesh().GetNode(elem_iter->GetNodeGlobalIndex(2));

//     int SpConst = 450; // prevents matrix begin singular
//     PositionVector[0] = Create_c_vector(0, 0, SpConst);
//     PositionVector[1] = pNode1->rGetLocation() - pNode0->rGetLocation() + Create_c_vector(0, 0, SpConst);
//     PositionVector[2] = pNode2->rGetLocation() - pNode0->rGetLocation() + Create_c_vector(0, 0, SpConst);
//     PositionVector = MatrixTranspose(PositionVector);

//     NewPositionVector[0] = Create_c_vector(0, 0, SpConst);
//     NewPositionVector[1] = Create_c_vector(a, 0, SpConst);
//     NewPositionVector[2] = Create_c_vector(b * cos(theta), b * sin(theta), SpConst);
//     NewPositionVector = MatrixTranspose(NewPositionVector);

//     c_vector<c_vector<long double, 3>, 3> InverseE = MatrixMultiplication(PositionVector, Inverse(NewPositionVector));

//     return InverseE;
// }

// c_vector<c_vector<long double, 3>, 3> MathsFunctions::Mapping(c_vector<c_vector<long double, 3>, 3> PositionVector,   c_vector<c_vector<long double, 3>, 3> NewPositionVector)
// {

//         int SpConst =450; // prevents matrix begin singular
//         for (int i=0; i<3; i++)
//         {
//            PositionVector[i][2] += SpConst;
//            NewPositionVector[i][2] = SpConst;
//         }
//           PositionVector = MatrixTranspose(PositionVector);
//           NewPositionVector = MatrixTranspose(NewPositionVector);

//         c_vector<c_vector<long double, 3>, 3>  InverseE  = MatrixMultiplication(PositionVector, Inverse(NewPositionVector));

// return InverseE;
// }

long double MathsFunctions::det(c_vector<c_vector<long double, 2>, 2> Matrix)
{

    long double Determinate = Matrix[0][0] * Matrix[1][1] - Matrix[0][1] * Matrix[1][0];
    return Determinate;
}

long double MathsFunctions::tr(c_vector<c_vector<long double, 2>, 2> Matrix)
{

    long double Trace = Matrix[0][0] + Matrix[1][1];
    return Trace;
}
