#include <iostream>
#include "itkSimplexMesh.h"

#include "kmUtility.h"
#include "kmVtkItkUtility.h"

int main(int argc, char* argv[])
{
	if (argc<2)
	{
		std::cerr<<"Input:"<<std::endl;
		std::cerr<<"       InputMesh GeometryFile OutputMesh"<<std::endl;
		return -1;
	}

	typedef itk::SimplexMesh<double, 3> MeshType;
	MeshType::Pointer inputMesh = km::readMesh<MeshType>( argv[1] );

	km::readSimplexMeshGeometryData<MeshType>(argv[2], inputMesh);
	km::ComputeGeometry<MeshType>(inputMesh, true);

	MeshType::PointDataContainerPointer allpointdata = inputMesh->GetPointData();
	MeshType::GeometryMapPointer        allgeodata   = inputMesh->GetGeometryData();

	allpointdata->Reserve( inputMesh->GetNumberOfPoints() );

	MeshType::GeometryMapIterator        geodataIt = allgeodata->Begin();
	MeshType::PointDataContainerIterator pointdataIt = allpointdata->Begin();

	SimplexMeshGeometry *geodata;

	while ( geodataIt != allgeodata->End() && pointdataIt != allpointdata->End() )
	{
		geodata = geodataIt.Value();

		allpointdata->InsertElement( geodataIt.Index(), geodata->phi );

		geodataIt++;
		pointdataIt++;
	}

	km::writeMesh<MeshType>( argv[3], inputMesh );

	return 0;
}