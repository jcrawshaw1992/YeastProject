
#include "CellMajorAxisAngleWriter.hpp"
#include "AbstractCellPopulation.hpp"
using namespace std;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellMajorAxisAngleWriter<ELEMENT_DIM, SPACE_DIM>::CellMajorAxisAngleWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("CelllMajorAxisAngle.dat")
{
    
    this->mVtkCellDataName ="/Major Axis Angle";
    // ofstream myfile;
    // myfile.open ("CelllMajorAxisAngle.dat", ios::app );
    // myfile << "Writing this to a file.\n";
    // myfile.close();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellMajorAxisAngleWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
     double ApsectRatio  =  pCell->GetCellData()->GetItem("MajorAxisAngle");
    return ApsectRatio;
}


/* I have several parameters writen into the CellParameteres.dat file 
At each time step I have 
time cell i.d AspectRatio MajorAxisAngle id AspectRatio MajorAxisAngle
*/

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellMajorAxisAngleWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
                                                                                        // AbstractCellPopulation<DIM,DIM>& rCellPopulation
{
    unsigned cell_id = pCell->GetCellId();
    

    double MajorAxisAngle  =  pCell->GetCellData()->GetItem("MajorAxisAngle");
    *this->mpOutStream << " " << cell_id << " " << MajorAxisAngle;
}

// Explicit instantiation
template class CellMajorAxisAngleWriter<1,1>;
template class CellMajorAxisAngleWriter<1,2>;
template class CellMajorAxisAngleWriter<2,2>;
template class CellMajorAxisAngleWriter<1,3>;
template class CellMajorAxisAngleWriter<2,3>;
template class CellMajorAxisAngleWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellMajorAxisAngleWriter)
