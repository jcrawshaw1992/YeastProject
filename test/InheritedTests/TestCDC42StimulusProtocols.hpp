#ifndef TESTCDC42STIMULUSPROTOCOLS_HPP
#define TESTCDC42STIMULUSPROTOCOLS_HPP

#include <cxxtest/TestSuite.h>
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "GradientCDC42StimulusProtocol.hpp"
#include "VShapedCDC42StimulusProtocol.hpp"
#include "RandomActivationCDC42StimulusProtocol.hpp"
#include "FlowBasedCDC42StimulusProtocol.hpp"
#include "UblasCustomFunctions.hpp"

class TestCDC42StimulusProtocols : public CxxTest::TestSuite
{
public:

    // Define a 2D simulation
    static const unsigned DIM = 2;

    void TestGradientProtocol()
    {
        // Circular mesh centred at (0,0) with radius 1.
        TrianglesMeshReader<DIM,DIM> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<DIM,DIM> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        GradientCDC42StimulusProtocol<DIM,DIM> stimulus_protocol(mesh, 10.0);

        // Node 51 is (-1, 0), left-hand-side edge of the circle
        TS_ASSERT_EQUALS(stimulus_protocol.GetCDC42ActivationRate(mesh.GetNode(51), 1), 0.0);
        // Node 1 is (1, 0), right-hand-side edge of the circle
        TS_ASSERT_EQUALS(stimulus_protocol.GetCDC42ActivationRate(mesh.GetNode(1), 1), 2.0 * stimulus_protocol.GetCDC42StimulusGradient());

        // The stimulus is turned off after 10s.
        TS_ASSERT_EQUALS(stimulus_protocol.GetCDC42ActivationRate(mesh.GetNode(1), 11), 0.0);
    }

    void TestVShapedProtocol()
    {
        // Circular mesh centred at (0,0) with radius 1.
        TrianglesMeshReader<DIM,DIM> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<DIM,DIM> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        VShapedCDC42StimulusProtocol<DIM,DIM> stimulus_protocol(mesh, 10.0);

        // Node 51 is (-1, 0), left-hand-side edge of the circle
        TS_ASSERT_EQUALS(stimulus_protocol.GetCDC42ActivationRate(mesh.GetNode(51), 1), 0.8*stimulus_protocol.GetCDC42StimulusGradient());
        // Node 1 is (1, 0), right-hand-side edge of the circle
        TS_ASSERT_EQUALS(stimulus_protocol.GetCDC42ActivationRate(mesh.GetNode(1), 1), stimulus_protocol.GetCDC42StimulusGradient());

        // The stimulus is turned off after 10s.
        TS_ASSERT_EQUALS(stimulus_protocol.GetCDC42ActivationRate(mesh.GetNode(1), 11), 0.0);
    }

    void TestRandomProtocol()
    {
        // Circular mesh centred at (0,0) with radius 1.
        TrianglesMeshReader<DIM,DIM> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<DIM,DIM> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        RandomActivationCDC42StimulusProtocol<DIM,DIM> stimulus_protocol(mesh, 20.0);

        // Create a mock node located at one of the noise centres
        Node<DIM> node(0, stimulus_protocol.GetNoiseCentre()[0]);
        TS_ASSERT_EQUALS(stimulus_protocol.GetCDC42ActivationRate(&node, 1), 0.6195);

        // The stimulus is turned off after 20s.
        TS_ASSERT_EQUALS(stimulus_protocol.GetCDC42ActivationRate(&node, 21), 0.0);
    }

    void TestFlowBasedProtocol()
    {
        // Circular mesh centred at (0,0) with radius 1.
        TrianglesMeshReader<DIM,DIM> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<DIM,DIM> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Flow in the negative direction of x
        FlowBasedCDC42StimulusProtocol<DIM,DIM> stimulus_protocol(mesh, 10.0, Create_c_vector(-1, 0));

        // Node 51 is (-1, 0), left-hand-side edge of the circle
        TS_ASSERT_EQUALS(stimulus_protocol.GetCDC42ActivationRate(mesh.GetNode(51), 1), 0.0);
        // Node 1 is (1, 0), right-hand-side edge of the circle
        TS_ASSERT_EQUALS(stimulus_protocol.GetCDC42ActivationRate(mesh.GetNode(1), 1), 1.0);
        // Node 131 is not at edge of the mesh
        TS_ASSERT_EQUALS(stimulus_protocol.GetCDC42ActivationRate(mesh.GetNode(131), 1), 0.0);

        // The stimulus is turned off after 10s.
        TS_ASSERT_EQUALS(stimulus_protocol.GetCDC42ActivationRate(mesh.GetNode(1), 11), 0.0);
    }
};

#endif
