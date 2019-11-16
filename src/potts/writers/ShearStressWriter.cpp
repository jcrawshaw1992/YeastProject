

#include "ShearStressWriter.hpp"
#include "AbstractCellPopulation.hpp"
using namespace std;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ShearStressWriter<ELEMENT_DIM, SPACE_DIM>::ShearStressWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("WallShearStress.dat")
{
    
    this->mVtkCellDataName ="/Wall shear stress";
    // ofstream myfile;
    // myfile.open ("CelllMajorAxisAngle.dat", ios::app );
    // myfile << "Writing this to a file.\n";
    // myfile.close();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ShearStressWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
     double WallShearStress  =  pCell->GetCellData()->GetItem("WallShearStress");

    return WallShearStress;
}


/* I have several parameters writen into the CellParameteres.dat file 
At each time step I have 
time cell i.d AspectRatio MajorAxisAngle id AspectRatio MajorAxisAngle
*/

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ShearStressWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
                                                                                        // AbstractCellPopulation<DIM,DIM>& rCellPopulation
{
    unsigned cell_id = pCell->GetCellId();
    

    double WallShearStress  =  pCell->GetCellData()->GetItem("WallShearStress");
    *this->mpOutStream << " " << cell_id << " " << WallShearStress;
}

// Explicit instantiation
template class ShearStressWriter<1,1>;
template class ShearStressWriter<1,2>;
template class ShearStressWriter<2,2>;
template class ShearStressWriter<1,3>;
template class ShearStressWriter<2,3>;
template class ShearStressWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ShearStressWriter)
