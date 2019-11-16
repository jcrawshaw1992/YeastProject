
/*

PeriodicRectangleMeshGenerator.cpp

*/

#include "PeriodicRectangleMeshGenerator.hpp"

#include <boost/foreach.hpp>
#include "TrianglesMeshReader.hpp" 
#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
#include "MathsCustomFunctions.hpp"
#include "ChasteSyscalls.hpp"
#include "Debug.hpp"

PeriodicRectangleMeshGenerator::PeriodicRectangleMeshGenerator(unsigned numNodesLongWidth, unsigned numNodesAlongLength, double width, double length)
  : mpMesh(NULL),
    mMeshFilename("mesh")
{
    // The code below won't work in parallel
    assert(PetscTools::IsSequential());

    // Get a unique temporary foldername
    std::stringstream pid;
    pid << getpid();
    OutputFileHandler output_file_handler("Periodic_rectangle_mesh" + pid.str());
    std::string output_dir = output_file_handler.GetOutputDirectoryFullPath();

   
	unsigned num_nodes =  numNodesLongWidth * numNodesAlongLength;
    unsigned num_face_around       = numNodesLongWidth;
	unsigned num_face_along_length = numNodesAlongLength-1;
	unsigned num_face              = 2*num_face_around*num_face_along_length;
	unsigned num_edges             = 3*num_face_around*num_face_along_length + num_face_around + num_face_along_length;

    double horizontal_spacing = width/ (double)numNodesLongWidth;
    double vertical_spacing = length/ (double)numNodesAlongLength;
	 
    // This line is needed to define ghost nodes later
    // mDomainDepth = (double)(num_nodes_along_length)*vertical_spacing;

    double x0 = -horizontal_spacing * 0;
    double y0 = -vertical_spacing * 0;
    double z = 0;
/////

	// Write node file
	out_stream p_node_file = output_file_handler.OpenOutputFile(mMeshFilename+".node");
	(*p_node_file) << std::scientific;
	//(*p_node_file) << std::setprecision(20);
	(*p_node_file) << num_nodes << "\t3\t0\t1" << std::endl;

	unsigned node = 0;
	for (unsigned i=0; i<numNodesAlongLength; i++)
	{
		for (unsigned j=0; j<numNodesLongWidth; j++)
		{
			unsigned boundary = 0;
			if ((i==0) || (i==numNodesAlongLength-1))
			{
			  boundary = 2;
			}
			if ((j==0) ||(j==numNodesLongWidth-1) )
			{
			  boundary = 1;
		
			}	
		    if (( j==0 && i==0) || (j==0 && i==numNodesAlongLength-1) || ( j==numNodesLongWidth-1 && i==0) || (j==numNodesLongWidth-1 && i==numNodesAlongLength-1)  ) 
			{
			  boundary = 3;
			}	
			double x = x0 + horizontal_spacing * ((double)j + 0.25 * (1.0 + SmallPow(-1.0, i + 1)));
            double y = y0 + vertical_spacing * (double)i;

			(*p_node_file) << node++ << "\t" << x << "\t" << y << "\t" << z << "\t" << boundary << std::endl;
			mBoundaryVector.push_back(boundary);
		}
	}
	p_node_file->close();



	// Write faceent file and edge file
	out_stream p_face_file = output_file_handler.OpenOutputFile(mMeshFilename+".face");
	(*p_face_file) << std::scientific;

	out_stream p_edge_file = output_file_handler.OpenOutputFile(mMeshFilename+".edge");
	(*p_node_file) << std::scientific;

	(*p_face_file) << num_face << "\t0" << std::endl;
	(*p_edge_file) << num_edges << "\t1" << std::endl;

	unsigned face = 0;
	unsigned edge = 0;
	for (unsigned i=0; i<num_face_along_length; i++)
	{
		for (unsigned j=0; j < num_face_around; j++)
		{
			unsigned node0 =     i*numNodesLongWidth + j;
			unsigned node1 =     i*numNodesLongWidth + j+1;
			unsigned node2 = (i+1)*numNodesLongWidth + j;

			if (j == num_face_around -1)
			{
				node1 = node1 - numNodesLongWidth;

			}
			if (i%2 != 0)
			{
				node2 = node2 + 1;
				if (j == num_face_around -1)
				{
					node2 = node2 - numNodesLongWidth;

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

			node0 = i*numNodesLongWidth + j + 1;
			node1 = (i+1)*numNodesLongWidth + j+1;
			node2 = (i+1)*numNodesLongWidth + j;

			if (i%2 != 0)
			{
				node0 = node0 - 1;
				if (j == num_face_around -1)
				{
					node1 = node1 - numNodesLongWidth;
				}
			}
			else
			{
				if (j == num_face_around -1)
				{
					node0 = node0 - numNodesLongWidth;
					node1 = node1 - numNodesLongWidth;
				}
			}

			(*p_face_file) << face++ << "\t" << node0 << "\t" << node1 << "\t" << node2 << std::endl;
		}
	}

	for (unsigned i=0; i<num_face_along_length; i++)
	{
		unsigned node0, node1;

		if (i%2==0)
		{
			 node0 = (i+1)*numNodesLongWidth - 1;
			 node1 = (i+2)*numNodesLongWidth - 1;
		}
		else
		{
			node0 = (i+1)*numNodesLongWidth;
			node1 = (i)*numNodesLongWidth;
		}
		(*p_edge_file) << edge++ << "\t" << node0 << "\t" << node1 << "\t" << 1 << std::endl;
	}

	for (unsigned j=0; j<num_face_around; j++)
	{
		unsigned node0 = numNodesLongWidth*(numNodesAlongLength-1) + j;
		unsigned node1 = numNodesLongWidth*(numNodesAlongLength-1) + j+1;

		if( j== num_face_around -1 )
		{
			node1 = node1 -numNodesLongWidth;
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
	mpMesh->SetMeshHasChangedSinceLoading();
}



std::vector<unsigned > PeriodicRectangleMeshGenerator::GetBoundaryVector()
{
		return mBoundaryVector;
}

PeriodicRectangleMeshGenerator::~PeriodicRectangleMeshGenerator()
{
    delete mpMesh;
}

MutableMesh<2,3>* PeriodicRectangleMeshGenerator::GetMesh()
{

	// bool bound;
	//  for (typename MutableMesh<2,3>::NodeIterator node_iter = mpMesh->GetNodeIteratorBegin();
	//          node_iter != mpMesh->GetNodeIteratorEnd();
	//          ++node_iter)
	//     {
			
	//          bound = node_iter->IsBoundaryNode();
	//          PRINT_VARIABLE(bound)
	//     }



    return mpMesh;



}
