#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"
#include "itkSimplexMesh.h"
#include "itkListSample.h"
#include "itkKdTreeGenerator.h"
#include "itkWeightedCentroidKdTreeGenerator.h"

#include "kmUtility.h"
#include "kmVtkItkUtility.h"

using namespace km;

int main(int argc, char* argv[])
{
	if(argc<5)
	{
		std::cout<<"Usage:"<<std::endl;
		std::cout<<"inputMesh geometryData outputMesh iterations";
		return EXIT_FAILURE;
	}

	const char* inputMeshFile = argv[1];
	const char* geometryDataFile = argv[2];
	const char* outputMeshFile = argv[3];
	const int iterations = atoi(argv[4]);
	
	typedef itk::SimplexMesh<double, 3> MeshType;
	typedef MeshType::PointType PointType;
	typedef MeshType::VectorType VectorType;
	typedef MeshType::PointsContainer PointsContainer;
	typedef MeshType::PointsContainerConstPointer PointsContainerConstPointer;
	typedef MeshType::PointsContainerPointer PointsContainerPointer;
	typedef MeshType::PointsContainerIterator PointsContainerIterator;
	typedef MeshType::PointDataContainer PointDataContainer;
	typedef MeshType::GeometryMapType GeometryMapType;
	typedef GeometryMapType::Pointer GeometryMapPointer;
	typedef GeometryMapType::Iterator GeometryMapIterator;
	
	typedef Statistics::ListSample<PointType>                       SampleType;
	typedef Statistics::WeightedCentroidKdTreeGenerator<SampleType> TreeGeneratorType;
	typedef TreeGeneratorType::KdTreeType                           KdTreeType;
	typedef KdTreeType::InstanceIdentifierVectorType                NeighborhoodIdentifierType;
	
	SampleType::Pointer samples = SampleType::New();
	samples->SetMeasurementVectorSize( PointType::PointDimension );
	
	TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
	treeGenerator = TreeGeneratorType::New();
	treeGenerator->SetBucketSize( 4 );
	
	GeometryImageType::Pointer geoImage = km::readImage<GeometryImageType>( geometryDataFile );
	MeshType::Pointer inputMesh = km::readMesh<MeshType>(inputMeshFile);
	loadSimplexMeshGeometryData<MeshType>(geoImage, inputMesh);

	GeometryMapIterator geoIt = inputMesh->GetGeometryData()->Begin();
	GeometryMapIterator geoItEnd = inputMesh->GetGeometryData()->End();
	while(geoIt!=geoItEnd)
	{
		samples->PushBack(inputMesh->GetPoint(geoIt.Index()));
		geoIt++;
	}

	itk::TimeProbesCollectorBase chronometer;
	itk::MemoryProbesCollectorBase memorymeter;
	memorymeter.Start( "starting" );
	chronometer.Start( "starting" );
	
	for(int iter=0;iter<iterations;iter++)
	{
		//std::cout<<"Iteration: "<<iter<<std::endl;
		treeGenerator->SetSample( samples );
		treeGenerator->Update();
		km::ComputeGeometry<MeshType>( inputMesh );

		geoIt = inputMesh->GetGeometryData()->Begin();
		
		itk::SimplexMeshGeometry *geodata;
		MeshType::PointIdentifier idx;
		VectorType normal;
		PointType pos_cur;
		PointType pos_new;
		PointType pos_closest;
		while(geoIt!=geoItEnd)
		{
			idx = geoIt.Index();
			geodata = geoIt.Value();
			normal.Set_vnl_vector(geodata->normal.Get_vnl_vector());
			pos_cur = geodata->pos;
			
			TreeGeneratorType::KdTreeType::InstanceIdentifierVectorType neighbors;
			treeGenerator->GetOutput()->Search( pos_cur, 1u, neighbors );
			
			if(neighbors[0]!=idx && neighbors[0]!=geodata->neighborIndices[0] && neighbors[0]!=geodata->neighborIndices[1] && neighbors[0]!=geodata->neighborIndices[2])
			{
				pos_new = pos_cur + normal*0;
				std::cout<<idx<<", "<<neighbors.size()<<", "<<neighbors[0]<<std::endl;
			}
			else
			{
				pos_new = pos_cur + normal*0;
			}

			//pos_new = pos_cur;
			
			geodata->pos = pos_new;
			samples->SetMeasurementVector(idx, pos_new);
			inputMesh->SetPoint(idx, pos_new);
			
			geoIt++;
		}
	}

	chronometer.Stop( "starting" );
	memorymeter.Stop( "starting" );
	chronometer.Report( std::cout );
	memorymeter.Report( std::cout );

	MeshType::Pointer outputMesh = inputMesh;
	km::writeMesh<MeshType>(outputMeshFile, outputMesh);
	
	return EXIT_SUCCESS;
}