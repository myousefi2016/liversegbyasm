#include <itkTriangleMeshToSimplexMeshFilter.h>
#include <itkSimplexMesh.h>

#include "kmUtility.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

using namespace std;
using namespace km;

const unsigned Dimensions = 3;
//typedef itk::Mesh<float, Dimensions> MeshType;
typedef itk::SimplexMesh<float, Dimensions> MeshType;

int main(int argc, char* argv[]) {

  if (argc < 2) {
    std::cout << "usage " << argv[0] << " inputMesh [outputMesh]" << std::endl;
    exit(-1);
  }

  itk::TimeProbesCollectorBase chronometer;
  itk::MemoryProbesCollectorBase memorymeter;
  memorymeter.Start( "converting" );
  chronometer.Start( "converting" );
  
  MeshType::Pointer inputMesh = readMesh<MeshType>( argv[1] );
  std::cout<<"Number of points: "<<inputMesh->GetNumberOfPoints()<<std::endl;
  std::cout<<"Number of cells: "<<inputMesh->GetNumberOfCells()<<std::endl;
  
	if(argc>2)
	{
		//MeshType::Pointer outputMesh = triangleMeshToSimplexMesh<MeshType, MeshType>( inputMesh );
		writeMesh<MeshType>(argv[2], inputMesh);
	}

  chronometer.Stop( "converting" );
  memorymeter.Stop( "converting" );
  chronometer.Report( std::cout );
  memorymeter.Report( std::cout );

  return EXIT_SUCCESS;
}

