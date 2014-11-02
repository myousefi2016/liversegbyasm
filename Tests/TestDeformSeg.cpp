
#include <iostream>
#include <fstream>

#include "itkMesh.h"
#include "itkSimplexMesh.h"
#include "itkRegularSphereMeshSource.h"
#include "itkDefaultDynamicMeshTraits.h"
//#include "itkSimplexMeshAdaptTopologyFilter2.h"
#include "itkSmoothingQuadEdgeMeshFilter.h"
#include "itkDeformableSimplexMesh3DGradientConstraintForceFilter.h"

#include "kmUtility.h"
#include "kmVtkItkConverter.h"
#include "kmSegmentation.h"

const unsigned Dimension = 3;
typedef itk::DefaultDynamicMeshTraits<double, Dimension, 3,double,double> MeshTraitsType;
typedef itk::Mesh<double,Dimension, MeshTraitsType> TriangleMeshType;
typedef itk::SimplexMesh<double,Dimension, MeshTraitsType> SimplexMeshType;
typedef itk::Image<unsigned char, Dimension> ImageType;
typedef itk::Image<float, Dimension> FloatImageType;

void segSimplexMesh( int argc, char* argv[] )
{
	//argv:
	//    1-input binary image
	//    2-initial mesh
	//    3-geometry image
	//    4-deformed mesh
	//    5-iterations

	for (int i=1;i<argc;i++)
	{
		KM_DEBUG_INFO(argv[i]);
	}

	ImageType::Pointer inputImage = km::readImage<ImageType>(argv[1]);

	SimplexMeshType::Pointer initialMesh = SimplexMeshType::New();

	if (false)
	{
		KM_DEBUG_INFO("Initial mesh from loading..");
		//SimplexMeshType::Pointer  initialMesh = km::readMesh<SimplexMeshType>( argv[2] );
		vtkSmartPointer<vtkPolyData> initPolydata = km::readPolyData(argv[2]);

		//KM_DEBUG_INFO( initPolydata->GetNumberOfPolys() );

		//initialMesh = km::polyData2Mesh<SimplexMeshType>( initPolydata );
		km::readSimplexMeshGeometryData<SimplexMeshType>( argv[3], initialMesh );
	}
	else
	{
		ImageType::PointType centroid = km::getCentroid<ImageType>( inputImage );
		SimplexMeshType::Pointer triangleMesh 
			= km::createRegularSphereMesh<SimplexMeshType>( centroid, 130, 4 );
		initialMesh 
			= km::triangleMeshToSimplexMesh<SimplexMeshType, SimplexMeshType>( triangleMesh );
		//km::writeMesh<SimplexMeshType>( argv[2], initialMesh );
	}

	KM_DEBUG_PRINT( "number of points: ", initialMesh->GetNumberOfPoints() );
	KM_DEBUG_PRINT( "number of cells : ", initialMesh->GetNumberOfCells() );

	KM_DEBUG_INFO("Generate initial mesh done");

	unsigned int iterations = 1000;
	if (argc>3)
	{
		iterations = atoi(argv[3]);
	}

	unsigned int adaptioniterations = 0;
	if (argc>4)
	{
		adaptioniterations = atoi(argv[4]);
	}


	typedef itk::DeformableSimplexMesh3DBalloonForceFilter<SimplexMeshType,SimplexMeshType> DeformFilterType;
	typedef DeformFilterType::FeatureImageType                                              FeatureImageType;
	typedef DeformFilterType::GradientImageType                                             GradientImageType;
	typedef itk::GradientRecursiveGaussianImageFilter<FeatureImageType,GradientImageType>   GradientFilterType;

	FeatureImageType::Pointer featureImage = km::gaussSmooth<FeatureImageType>( km::castImage<ImageType, FeatureImageType>( inputImage ), 2 );
	
	FeatureImageType::SpacingType spa;
	spa.Fill( 3 );
	featureImage = km::resampleImage<FeatureImageType>( featureImage, spa, 0 );

	//typedef itk::SobelEdgeDetectionImageFilter<FeatureImageType,FeatureImageType>   EdgeFilterType;

	//EdgeFilterType::Pointer edgeFilter = EdgeFilterType::New();
	//edgeFilter->SetInput( featureImage );
	//edgeFilter->Update();
	//FeatureImageType::Pointer edgeMap = edgeFilter->GetOutput();

	//km::writeImage<FeatureImageType>( "edgeMap.nii.gz", edgeMap );
	
	FeatureImageType::Pointer edgeMap
		= km::calculateGradientMagnitudeImage<FeatureImageType, FeatureImageType>( featureImage, 1.0 );
	edgeMap 
		= km::rescaleIntensity<FeatureImageType, FeatureImageType>( edgeMap, 0, 30 );

	KM_DEBUG_INFO( "Generate edge map done!" );
	GradientImageType::Pointer gradientMapS2 
		= km::calculateRecursiveGradientImage<FeatureImageType, GradientImageType>( edgeMap, 1.0 );

	typedef itk::SimplexMeshAdaptTopologyFilter<SimplexMeshType, SimplexMeshType> AdaptFilterType;

	std::vector<SimplexMeshType::Pointer> meshes;
	meshes.push_back( initialMesh );

	km::writeMesh<SimplexMeshType>( ".", "initialMesh", 0, ".vtk", initialMesh );

	unsigned k = 0;
	while(true)
	{
		SimplexMeshType::Pointer deformedMesh = meshes[k];
		km::deformSegSimplexMesh2<DeformFilterType, FeatureImageType, GradientImageType, SimplexMeshType>(
			featureImage,//potential image
			gradientMapS2,   //gradient image
			deformedMesh,  //initial mesh
			0.5,           //alpha
			0.1,           //beta
			-0.05,         //kamma
			iterations,           //iterations
			1,      //rigidity
			0.001,
			1.0);
		km::writeMesh<SimplexMeshType>( ".", "deformedMesh", k, ".vtk", deformedMesh );
		KM_DEBUG_PRINT("Deform doen ", k);

		if (k >= adaptioniterations)
		{
			break;
		}

		iterations *= 0.9;
		if (iterations<5)
		{
			break;
		}

		AdaptFilterType::Pointer adaptor = AdaptFilterType::New();
		adaptor->SetInput( deformedMesh );
		adaptor->SetThreshold( 0.8 );
		adaptor->Update();
		SimplexMeshType::Pointer adaptedMesh = adaptor->GetOutput();

		km::writeMesh<SimplexMeshType>( ".", "adaptedMesh", k, ".vtk", adaptedMesh );
		KM_DEBUG_PRINT("Adapt done ", k);

		meshes.push_back( adaptedMesh );

		k++;
	}

	//km::writeSimplexMeshGeometryData<SimplexMeshType>( argv[3], deformedMesh );
	//km::writeMesh<SimplexMeshType>( argv[2], deformedMesh );
}

int main(int argc, char * argv[] )
{
	if(argc<3)
	{
		std::cout<<"Usage:"<<std::endl;
		std::cout<<"InputImage InitialMesh OutputMesh";
		return EXIT_FAILURE;
	}

	itk::TimeProbesCollectorBase chronometer;
	itk::MemoryProbesCollectorBase memorymeter;
	memorymeter.Start( "starting" );
	chronometer.Start( "starting" );

	ImageType::Pointer inputImage = km::readImage<ImageType>(argv[1]);

	SimplexMeshType::Pointer initialMesh = km::readMesh<SimplexMeshType>(argv[2]);
	km::readSimplexMeshGeometryData<SimplexMeshType>( "geoImage.mha", initialMesh );

	typedef itk::DeformableSimplexMesh3DBalloonForceFilter<SimplexMeshType, SimplexMeshType> DeformableMeshFilterType;
	typedef DeformableMeshFilterType::FeatureImageType OriginalImageType;
	OriginalImageType::Pointer originalImage = km::castImage<ImageType, OriginalImageType>(inputImage);
	typedef DeformableMeshFilterType::GradientImageType GradientImageType;
	FloatImageType::Pointer edgeImage = km::sobelEdge<ImageType,FloatImageType>(inputImage);
	km::writeImage<FloatImageType>("edgeImage.nii.gz", edgeImage);
	GradientImageType::Pointer gradientImage = km::calculateRecursiveGradientImage<FloatImageType, GradientImageType>(edgeImage, 3);
	//km::deformSegSimplexMeshWithGradientConstaint<DeformableMeshFilterType, OriginalImageType, GradientImageType, SimplexMeshType>(
	//	originalImage,
	//	gradientImage,
	//	initialMesh
	//	);
	km::deformSegSimplexMesh2<DeformableMeshFilterType, OriginalImageType, GradientImageType, SimplexMeshType>(
		originalImage,
		gradientImage,
		initialMesh,
		0.3,
		0.05,
		0.05,
		300,
		0,
		0.01,
		1.0
		);

	km::writeMesh<SimplexMeshType>( argv[3], initialMesh );

	//segSimplexMesh( argc, argv );

	chronometer.Stop( "starting" );
	memorymeter.Stop( "starting" );
	chronometer.Report( std::cout );
	memorymeter.Report( std::cout );
}