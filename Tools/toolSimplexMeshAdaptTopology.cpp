#include "itkRegularSphereMeshSource.h"
#include "itkTriangleMeshToSimplexMeshFilter.h"
#include "itkSimplexMeshAdaptTopologyFilter.h"
#include "itkDefaultDynamicMeshTraits.h"

#include "kmUtility.h"
#include "kmVtkItkUtility.h"

int main( int argc, char * argv[] )
{
	if (argc < 4) {
		std::cout << "usage " << argv[0] << " inputMesh geometryImage outputMesh" << std::endl;
		exit(-1);
	}

	itk::TimeProbesCollectorBase chronometer;
	itk::MemoryProbesCollectorBase memorymeter;
	memorymeter.Start( "converting" );
	chronometer.Start( "converting" );

	// Declare the type of the input and output mesh
	typedef itk::DefaultDynamicMeshTraits<double, 3, 3,double,double> TriangleMeshTraits;
	typedef itk::DefaultDynamicMeshTraits<double, 3, 3, double,double> SimplexMeshTraits;
	typedef itk::Mesh<double,3, TriangleMeshTraits> TriangleMeshType;
	typedef itk::SimplexMesh<double,3, SimplexMeshTraits> SimplexMeshType;

	SimplexMeshType::Pointer simplexMesh = km::readMesh<SimplexMeshType>( argv[1] );
	km::readSimplexMeshGeometryData<SimplexMeshType>( argv[2], simplexMesh );
	km::ComputeGeometry<SimplexMeshType>( simplexMesh, true );

	km::writeMesh<SimplexMeshType>("computedMesh.vtk", simplexMesh);

	km::adaptMesh<SimplexMeshType>( simplexMesh, 0.05, 50 );

	km::writeMesh<SimplexMeshType>( argv[3], simplexMesh );

	chronometer.Stop( "converting" );
	memorymeter.Stop( "converting" );
	chronometer.Report( std::cout );
	memorymeter.Report( std::cout );

	return EXIT_SUCCESS;
}

