

#ifndef LINEARSPRINGFORCE_HPP_
#define LINEARSPRINGFORCE_HPP_

#include "AbstractForce.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "AbstractCellPopulation.hpp"

#include "UblasCustomFunctions.hpp"
/**
 * An abstract class for two-body force laws.
 */
template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class LinearSpringForce : public AbstractForce<ELEMENT_DIM, SPACE_DIM>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<ELEMENT_DIM,SPACE_DIM> >(*this);
        archive & mNearestNodesMap;
        archive & mMeinekeSpringStiffness;

    }

protected:

    /** Whether to have zero force if the cells are far enough apart. */
    // bool mUseCutOffLength;

    // /** Mechanics cut off length. */
    // double mMechanicsCutOffLength;
    double mMeinekeSpringStiffness;

public:

    /**
     * Constructor.
     */
    LinearSpringForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);


    // When the linear spring force is corrected for drag, the edges break everything, so I need to be able to remove the spring force on these regions and replace it with the force in the neighbourhood 
    void SetNearestNeighboursMap( std::map<unsigned, c_vector<unsigned, 2> > NearestNodesMap);
    std::map<unsigned, c_vector<unsigned, 2> > mNearestNodesMap;


    /**
     * Overridden CalculateForceBetweenNodes() method.
     *
     * Calculates the force between two nodes.
     *
     * Note that this assumes they are connected and is called by AddForceContribution()
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     * @return The force exerted on Node A by Node B.
     */

    c_vector<double, SPACE_DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                     unsigned nodeBGlobalIndex,
                                                     AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * Set mMeinekeSpringStiffness.
     *
     * @param springStiffness the new value of mMeinekeSpringStiffness
     */
    void SetMeinekeSpringStiffness(double springStiffness);



    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);

    /**
     * Overridden WriteDataToVisualizerSetupFile() method.
     * Write any data necessary to a visualization setup file.
     * Used by AbstractCellBasedSimulation::WriteVisualizerSetupFile().
     *
     * @param pVizSetupFile a visualization setup file
     */
    virtual void WriteDataToVisualizerSetupFile(out_stream& pVizSetupFile);
};

#endif /*LinearSpringForce_HPP_*/
