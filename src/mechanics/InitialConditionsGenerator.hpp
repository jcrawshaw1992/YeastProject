#ifndef INITIALCONDITIONSGENERATOR_HPP_
#define INITIALCONDITIONSGENERATOR_HPP_

#include <cmath>
#include <vector>

#include "MutableMesh.hpp"
#include "MeshBasedCellPopulation.hpp"

class InitialConditionsGenerator
{
public:

    InitialConditionsGenerator();

    /*
     * Calculate the rest lengths to give zero force on each node.
     *
     * Note this calculates the least squares approximation to a 3*NumNodes by NumEdges Linear system.
     */
    void SetInitialRestLengths(MeshBasedCellPopulation<2,3>& rCellPopulation, double springConstant);

    /*
     * Calculate the rest lengths to give zero force in the direction of the normal to the surface
     *
     * Note this is solving a NumNodes by NumNodes linear system.
     */
    void SetApproximateInitialRestLengths(MeshBasedCellPopulation<2,3>& rCellPopulation, double springConstant);

    /*
     * Calculate the rest lengths to give zero force in the direction of the normal to the surface
     *
     * Note this done by assuming only local interactions so gives a diagonal system.
     */
    void SetSimpleInitialRestLengths(MeshBasedCellPopulation<2,3>& rCellPopulation, double springConstant);


};

#endif /*INITIALCONDITIONSGENERATOR_HPP_*/
