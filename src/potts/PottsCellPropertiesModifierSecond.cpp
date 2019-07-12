

#include "PottsCellPropertiesModifier.hpp"
#include "MeshBasedCellPopulation.hpp"

template<unsigned DIM>
PottsCellPropertiesModifier<DIM>::PottsCellPropertiesModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
PottsCellPropertiesModifier<DIM>::~PottsCellPropertiesModifier()
{
}

template<unsigned DIM>
void PottsCellPropertiesModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PottsCellPropertiesModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    // UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PottsCellPropertiesModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    // rCellPopulation.Update();

    // /**
    //  * This hack is needed because in the case of a MeshBasedCellPopulation in which
    //  * multiple cell divisions have occurred over one time step, the Voronoi tessellation
    //  * (while existing) is out-of-date. Thus, if we did not regenerate the Voronoi
    //  * tessellation here, an assertion may trip as we try to access a Voronoi element
    //  * whose index exceeds the number of elements in the out-of-date tessellation.
    //  *
    //  * \todo work out how to properly fix this (#1986)
    //  */
    // if (bool(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation)))
    // {
    //     static_cast<MeshBasedCellPopulation<DIM>*>(&(rCellPopulation))->CreateVoronoiTessellation();
    // }

    // // Iterate over cell population
    // for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
    //      cell_iter != rCellPopulation.End();
    //      ++cell_iter)
    // {
    //     // Get the volume of this cell
    //     double cell_volume = rCellPopulation.GetVolumeOfCell(*cell_iter);

    //     // Store the cell's volume in CellData
    //     cell_iter->GetCellData()->SetItem("volume", cell_volume);
    // }
}

template<unsigned DIM>
void PottsCellPropertiesModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class PottsCellPropertiesModifier<1>;
template class PottsCellPropertiesModifier<2>;
template class PottsCellPropertiesModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PottsCellPropertiesModifier)

