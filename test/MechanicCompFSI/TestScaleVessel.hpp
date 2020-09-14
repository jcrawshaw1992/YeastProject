
#ifndef TestBendingForceParameterSweep_HPP_
#define TestBendingForceParameterSweep_HPP_

#include <cmath>
#include <cstdio>
#include <ctime>
#include <cxxtest/TestSuite.h>

#include "Debug.hpp"
#include "SmartPointers.hpp"

#include "MutableMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "VtkMeshWriter.hpp"
#include "AbstractCellBasedTestSuite.hpp"

class SurfaceAreaForceOnCylinder : public AbstractCellBasedTestSuite
{


public:
    void TestScaleMesh() throw(Exception)
    {
       std::string mesh_file = "projects/VascularRemodelling/test/data/bifurcation_cut/config.vtu";
       std::string output_dir = "ScalledMesh/";

   
        VtkMeshReader<2, 3> mesh_reader(mesh_file);
        MutableMesh<2, 3> mutable_mesh;
        mutable_mesh.ConstructFromMeshReader(mesh_reader);
        double scale = 1e-3; // so distances are in mm

        mutable_mesh.Scale(4.8*scale , 4.8 *scale , 4.8*scale ); // so distances are back in original scale
        VtkMeshWriter<2,3> mesh_writer(output_dir, "config", false);
		mesh_writer.WriteFilesUsingMesh( mutable_mesh);


// Still need to covert to stl :( 



    }
};

#endif /*TESTCYLINDRICALGROWTHDEFORMABLEMEMBRANE_HPP_*/
