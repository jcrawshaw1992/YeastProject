
#ifndef MembraneStiffnessForceForTrianglePair_HPP_
#define MembraneStiffnessForceForTrianglePair_HPP_

#include "AbstractForce.hpp"
#include "MeshBasedCellPopulation.hpp"


#include "projects/VascularRemodelling/src/MutationStates/EmptyBasementMatrix.hpp"
#include "projects/VascularRemodelling/src/MutationStates/LostEndothelialCell.hpp"
#include "projects/VascularRemodelling/src/MutationStates/HasEndothelialCell.hpp"

/**
 * Membrane Stiffness Force!
 */
class MembraneStiffnessForceForTrianglePair : public AbstractForce<2, 3>
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
        archive & boost::serialization::base_object<AbstractForce<2,3> >(*this);
        archive & mMembraneStiffness;
        archive & mOriginalAngles;
        archive & mMembraneStiffnessMap;
    }

protected:

    double mMembraneStiffness;
    std::map<std::pair<unsigned, unsigned>, double> mOriginalAngles;
    std::map<std::pair<Node<3>*, Node<3>*>, double> mMembraneStiffnessMap;


    bool CalculateElementNormals(MutableMesh<2, 3>& rMesh, std::pair<Node<3>*, Node<3>*> edge,
                                 std::pair<c_vector<double, 3>, c_vector<double, 3> >& nonUnitNormals,
                                 std::pair<Node<3>*,  Node<3>*>& otherNodes);

public:

    /**
     * Constructor.
     */
    MembraneStiffnessForceForTrianglePair();

    void SetMembraneStiffness(double membraneStiffnes);

    double GetMembraneStiffness() const;

    double GetOriginalAngle(std::pair<Node<3>*, Node<3>*> edge);

    double SetOriginalAngle(std::pair<Node<3>*, Node<3>*> edge, double angle);

    void SetupInitialMembrane(MutableMesh<2,3>& rMesh);

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<2,3>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};


// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MembraneStiffnessForceForTrianglePair)

#endif /*MembraneStiffnessForceForTrianglePair_HPP_*/
