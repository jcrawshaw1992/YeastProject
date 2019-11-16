

#include "Honeycomb3DMeshGenerator.hpp"

#include <boost/foreach.hpp>
#include "TrianglesMeshReader.hpp"
#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
#include "MathsCustomFunctions.hpp"
#include "ChasteSyscalls.hpp"
#include "Debug.hpp"

Honeycomb3DMeshGenerator::Honeycomb3DMeshGenerator(unsigned numNodesLongWidth, unsigned numNodesAlongLength, double width, double length)
  : mpMesh(NULL),
    mMeshFilename("mesh")
{
	TRACE("DING");
    // The code below won't work in parallel
    assert(PetscTools::IsSequential());

    // Get a unique temporary foldername
    std::stringstream pid;
    pid << getpid();
    OutputFileHandler output_file_handler("cylinder_temporary_honeycomb_mesh_" + pid.str());
    std::string output_dir = output_file_handler.GetOutputDirectoryFullPath();

    // double angle_spacing = 2.0*M_PI / (double)num_nodes_along_width;
    // double length_spacing = length / ((double)numNodesAlongLength -1);

    unsigned num_nodes_along_width  = numNodesLongWidth;
    unsigned num_nodes_along_length = numNodesAlongLength;

    double horizontal_spacing = width/ (double)num_nodes_along_width;
	 
    double vertical_spacing = length/ (double)num_nodes_along_length;
	 
    // This line is needed to define ghost nodes later
    // mDomainDepth = (double)(num_nodes_along_length)*vertical_spacing;

    unsigned num_nodes = num_nodes_along_width * num_nodes_along_length;
    unsigned num_elem_along_width = num_nodes_along_width - 1;
    unsigned num_elem_along_length = num_nodes_along_length - 1;
    unsigned num_elem = 2 * num_elem_along_width * num_elem_along_length;
    unsigned num_edges = 3 * num_elem_along_width * num_elem_along_length + num_elem_along_width + num_elem_along_length;

    double x0 = -horizontal_spacing * 0;
    double y0 = -vertical_spacing * 0;
    double z = 0;

	// Write node file
	out_stream p_node_file = output_file_handler.OpenOutputFile(mMeshFilename+".node");
	(*p_node_file) << std::scientific;
	//(*p_node_file) << std::setprecision(20);
	
	(*p_node_file) << num_nodes << "\t3\t0\t1" << std::endl;

	unsigned node = 0;
	for (unsigned i=0; i<num_nodes_along_length; i++)
	{
	 
		for (unsigned j=0; j<num_nodes_along_width; j++)
		{
 
			unsigned boundary = 0;
			if ((i==0) || (i==numNodesAlongLength-1))
			{
				boundary = 1;
				 
			}
 
            double x = x0 + horizontal_spacing * ((double)j + 0.25 * (1.0 + SmallPow(-1.0, i + 1)));
            double y = y0 + vertical_spacing * (double)i;
	 

			(*p_node_file) << node++ << "\t" << x << "\t" << y << "\t" << z << "\t" << boundary << std::endl;
		}
	}
	 
	p_node_file->close();
	

	// Write faceent file and edge file
	out_stream p_face_file = output_file_handler.OpenOutputFile(mMeshFilename+".face");
	(*p_face_file) << std::scientific;

	out_stream p_edge_file = output_file_handler.OpenOutputFile(mMeshFilename+".edge");
	(*p_node_file) << std::scientific;

	(*p_face_file) << num_elem << "\t0" << std::endl;
	(*p_edge_file) << num_edges << "\t1" << std::endl;

	unsigned face = 0;
	unsigned edge = 0;
	for (unsigned i=0; i<num_elem_along_length; i++)
	{
		for (unsigned j=0; j < num_elem_along_width; j++)
		{
			unsigned node0 =     i*num_nodes_along_width + j;
			unsigned node1 =     i*num_nodes_along_width + j+1;
			unsigned node2 = (i+1)*num_nodes_along_width + j;

			if (j == num_nodes_along_width -1)
			{
				node1 = node1 - num_nodes_along_width;

			}
			if (i%2 != 0)
			{
				node2 = node2 + 1;
				if (j == num_nodes_along_width -1)
				{
					node2 = node2 - num_nodes_along_width;

				}
			}

			(*p_face_file) << face++ << "\t" << node0 << "\t" << node1 << "\t" << node2 << std::endl;

			unsigned horizontal_edge_is_boundary_edge = 0;
			unsigned vertical_edge_is_boundary_edge = 0;
			if (i==0)
			{
				horizontal_edge_is_boundary_edge = 1;
			}

			(*p_edge_file) << edge++ << "\t" << node0 << "\t" << node1 << "\t" << horizontal_edge_is_boundary_edge << std::endl;
			(*p_edge_file) << edge++ << "\t" << node1 << "\t" << node2 << "\t" << 0 << std::endl;
			(*p_edge_file) << edge++ << "\t" << node2 << "\t" << node0 << "\t" << vertical_edge_is_boundary_edge << std::endl;

			node0 = i*num_nodes_along_width + j + 1;
			node1 = (i+1)*num_nodes_along_width + j+1;
			node2 = (i+1)*num_nodes_along_width + j;

			if (i%2 != 0)
			{
				node0 = node0 - 1;
				if (j == num_nodes_along_width -1)
				{
					node1 = node1 - num_nodes_along_width;
				}
			}
			else
			{
				if (j == num_nodes_along_width -1)
				{
					node0 = node0 - num_nodes_along_width;
					node1 = node1 - num_nodes_along_width;
				}
			}

			(*p_face_file) << face++ << "\t" << node0 << "\t" << node1 << "\t" << node2 << std::endl;
		}
	}
	for (unsigned i=0; i<num_elem_along_length; i++)
	{
		
		unsigned node0, node1;

		if (i%2==0)
		{
			 node0 = (i+1)*num_nodes_along_width - 1;
			 node1 = (i+2)*num_nodes_along_width - 1;
		}
		else
		{
			node0 = (i+1)*num_nodes_along_width;
			node1 = (i)*num_nodes_along_width;
		}
		
		(*p_edge_file) << edge++ << "\t" << node0 << "\t" << node1 << "\t" << 1 << std::endl;
	}

	for (unsigned j=0; j<num_nodes_along_width; j++)
	{
		
		unsigned node0 = num_nodes_along_width*(numNodesAlongLength-1) + j;
		unsigned node1 = num_nodes_along_width*(numNodesAlongLength-1) + j+1;

		if( j== num_nodes_along_width -1 )
		{
			node1 = node1 -num_nodes_along_width;
		}
		(*p_edge_file) << edge++ << "\t" << node1 << "\t" << node0 << "\t" << 1 << std::endl;
		
	}

	p_face_file->close();
	p_edge_file->close();
	

	// Having written the mesh to file, now construct it using TrianglesMeshReader.
	// Nested scope so the reader closes files before we delete them below.
	{
		
		TrianglesMeshReader<2,3> mesh_reader(output_file_handler.GetOutputDirectoryFullPath() + mMeshFilename);
		
		mpMesh = new MutableMesh<2,3>();

		mpMesh->ConstructFromMeshReader(mesh_reader);

	}


//	// Make the mesh cylindrical (we use Triangle library mode inside this ReMesh call)
//	mpMesh->ReMesh();
//
//	// Delete the temporary files
//	output_file_handler.FindFile("").Remove();

	// The original files have been deleted, it is better if the mesh object forgets about them
	TRACE("The original files have been delete")
	mpMesh->SetMeshHasChangedSinceLoading();
}

Honeycomb3DMeshGenerator::~Honeycomb3DMeshGenerator()
{
    delete mpMesh;
}

MutableMesh<2,3>* Honeycomb3DMeshGenerator::GetMesh()
{

    return mpMesh;
}
