

#ifndef ElementQualityWriter_HPP_
#define ElementQualityWriter_HPP_

#include "AbstractCellPopulationWriter.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"


/**
 * A class written using the visitor pattern for writing Voronoi data from a cell population to file.
 *
 * The output file is called voronoi.dat by default.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ElementQualityWriter : public AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    ElementQualityWriter();

    /**
     * Visit the MeshBasedCellPopulation and write the index and location of each node, as well
     * as the volume and surface area (area and perimeter in 2 dimensions) of the corresponding
     * element in the dual Voronoi tessellation.
     *
     * Outputs a line of space-separated values of the form:
     * ...[node index] [node x-pos] [node y-pos] [node z-pos] [elem volume] [elem surface area] ...
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * @param pCellPopulation a pointer to the MeshBasedCellPopulation to visit.
     */
    virtual void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the data.
     *
     * This is an empty dummy function, since this class is defined for use with a MeshBasedCellPopulation only.
     *
     * @param pCellPopulation a pointer to the CaBasedCellPopulation to visit.
     */
    virtual void Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the data.
     *
     * This is an empty dummy function, since this class is defined for use with a MeshBasedCellPopulation only.
     *
     * @param pCellPopulation a pointer to the NodeBasedCellPopulation to visit.
     */
    virtual void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the data.
     *
     * This is an empty dummy function, since this class is defined for use with a MeshBasedCellPopulation only.
     *
     * @param pCellPopulation a pointer to the PottsBasedCellPopulation to visit.
     */
    virtual void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the data.
     *
     * This is an empty dummy function, since this class is defined for use with a MeshBasedCellPopulation only.
     *
     * @param pCellPopulation a pointer to the VertexBasedCellPopulation to visit.
     */
    virtual void Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ElementQualityWriter)

#endif /*VORONOIDATAWRITER_HPP_*/
