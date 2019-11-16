

#ifndef MathsFunctions_HPP_
#define MathsFunctions_HPP_

#include <math.h>
#include <cstdio>
#include "Debug.hpp"
#include <iostream>
#include <cassert>
#include <sstream>
#include <string>
#include "PetscTools.hpp"
#include "UblasCustomFunctions.hpp"




/**
 * Membrane Surface Force
 * Force tyring to minimising surface area to the relaxed state 
 */
class MathsFunctions 
{
public:


    // std::map<unsigned, c_vector< c_vector< double, 3>, 3> > mMapping;
  
    // std::map<unsigned, c_vector<double, 3> > mForceMap;
  
    c_vector<c_vector<long double, 3>, 3> Inverse(c_vector<c_vector<long double, 3>, 3> Matrix);
  
    c_vector<c_vector<long double, 3>, 3> Inverse(c_vector<c_vector<long double, 3>, 3> Matrix, double elem);
  
    c_vector<c_vector<long double, 3>, 3> MatrixMultiplication(c_vector<c_vector<long double, 3>, 3> Matrix1, c_vector<c_vector<long double, 3>, 3> Matrix2);
  
    c_vector<long double, 3> MatrixMultiplication(c_vector<c_vector<long double, 3>, 3> Matrix, c_vector<long double, 3> Vector);
  
    c_vector<c_vector<long double, 3>, 3> MatrixTranspose(c_vector<c_vector<long double, 3>, 3> Matrix1);

    c_vector<c_vector<long double, 3>, 3> RowReduction(c_vector<c_vector<long double, 3>, 3> Matrix);

    // c_vector<c_vector<long double, 3>, 3> MappingMatrix(MeshBasedCellPopulation<2, 3>* p_cell_population, typename AbstractTetrahedralMesh<2, 3>::ElementIterator elem_iter, double a, double b, double theta);
    
    long double det(c_vector<c_vector<long double, 2>, 2> Matrix);

    long double tr(c_vector<c_vector<long double, 2>, 2> Matrix);




    std::pair<double, double> Create_pair(double x, double y);

    std::pair<unsigned, unsigned> Create_pair(unsigned x, unsigned y);

    double MaintainOutwardsPointingNormal(c_vector<double, 3> Normal, c_vector<double, 3> x1);

    void PRINT_PAIR(std::pair<unsigned, unsigned> Pair);

    double AreVectorsSame(std::vector<unsigned> Vector1, std::vector<unsigned> Vector2);

    double AreVectorsSame(std::vector<std::pair<unsigned, unsigned> > Vector1, std::vector<std::pair<unsigned, unsigned> > Vector2);

    std::vector<unsigned> Intersection(std::vector<unsigned> Vector1, std::vector<unsigned> Vector2);

    std::vector<unsigned> RemoveElement(std::vector<unsigned> Vector1, unsigned number);

    std::vector<double> RemoveElement(std::vector<double> Vector1, double number);

    bool IsPairInVector(std::vector<std::pair<unsigned, unsigned> > Vector, std::pair<unsigned, unsigned> Pair);

    bool IsNumberInVector(std::vector<unsigned> Vector, unsigned number);

    bool IsNumberInVector(std::vector<double> Vector, double number);

    double min_value(std::vector<double> Vector);

    double max_value(std::vector<double> Vector);
    
    double AddAngles(double alpha, double beta);

    double PeriodicAngle(double eta);

    bool IsVectorInVector(std::vector< c_vector<double,3>  > Vector, c_vector<double,3> Locaiton );

};

// // Declare identifier for the serializer
// #include "SerializationExportWrapper.hpp"
// CHASTE_CLASS_EXPORT(MathsFunctions)

#endif /*MathsFunctions_HPP_*/
