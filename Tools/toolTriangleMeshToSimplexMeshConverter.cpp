#include <itkTriangleMeshToSimplexMeshFilter.h>
#include <itkSimplexMesh.h>
#include "itkSimplexMeshAdaptTopologyFilter.h"
#include "vtkDelaunay3D.h"
#include "vtkDataSetSurfaceFilter.h"
#include <vtkImageWriter.h>

#include "kmUtility.h"
#include "kmVtkItkUtility.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

using namespace std;
using namespace km;

const unsigned Dimensions = 3;
typedef itk::Image<unsigned char, Dimensions> UCharImageType;
typedef itk::Mesh<float, Dimensions> TriangleMeshType;
typedef itk::SimplexMesh<float, Dimensions> SimplexMeshType;

int main(int argc, char* argv[]) {

  if (argc < 3) {
    std::cout << "usage " << argv[0] << " inputMesh outputMesh geoImage " << std::endl;
    exit(-1);
  }

  itk::TimeProbesCollectorBase chronometer;
  itk::MemoryProbesCollectorBase memorymeter;
  memorymeter.Start( "converting" );
  chronometer.Start( "converting" );
  
  SimplexMeshType::Pointer inputMesh = readMesh<SimplexMeshType>( argv[1] );

	KM_DEBUG_INFO( "Converting from triangle mesh to simplex mesh!" );
	SimplexMeshType::Pointer outputMesh = triangleMeshToSimplexMesh<SimplexMeshType, SimplexMeshType>( inputMesh );
	
	writeMesh<SimplexMeshType>(argv[2], outputMesh);

	//km::writeSimplexMeshGeometryData<SimplexMeshType>( argv[3], outputMesh );

  chronometer.Stop( "converting" );
  memorymeter.Stop( "converting" );
  chronometer.Report( std::cout );
  memorymeter.Report( std::cout );

  return EXIT_SUCCESS;
}

