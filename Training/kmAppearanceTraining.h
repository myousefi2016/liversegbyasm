#ifndef __kmAppearanceTraining_h
#define __kmAppearanceTraining_h

#include <iostream>
#include <fstream>
#include  <direct.h> 
#include <cmath>
#include <ctime>

#include "itkDefaultDynamicMeshTraits.h"
#include "itkSimplexMesh.h"
#include "itkCovariantVector.h"
#include "itkDifferenceOfGaussiansGradientImageFilter.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include "kmKNNProfileClassifier-FLANN.h"
#include "kmUtility.h"
#include "kmVtkItkUtility.h"
#include "kmProcessing.h"

#define KMeansFlag 0

using namespace std;
using std::ostream;
using std::istream;

namespace km
{
	typedef std::vector<std::string>    StringVectorType;

	const unsigned int Dimension = 3;
	typedef itk::Image<short, Dimension>           OriginalImageType;
	typedef itk::Image<unsigned char, Dimension>   SegImageType;
	typedef itk::CovariantVector<double, Dimension>   GradientVectorType;
	typedef itk::Image<GradientVectorType, Dimension> GradientImageType;

	typedef OriginalImageType::SpacingType SpacingType;

	typedef double MeshPixelType;
	typedef itk::DefaultDynamicMeshTraits<MeshPixelType, Dimension, Dimension, MeshPixelType, MeshPixelType, MeshPixelType> MeshTraitsType;
	typedef itk::SimplexMesh<MeshPixelType,Dimension, MeshTraitsType> SimplexMeshType;

	typedef itk::LinearInterpolateImageFunction<OriginalImageType> OriginalInterpolatorType;
	typedef itk::LinearInterpolateImageFunction<GradientImageType> GradientInterpolatorType;

	typedef SimplexMeshType::PointsContainerConstPointer PointsContainerConstPointer;
	typedef SimplexMeshType::PointsContainerConstIterator PointsContainerConstIterator;

	typedef SimplexMeshType::GeometryMapType GeometryMapType;
	typedef GeometryMapType::Pointer         GeometryMapPointer;
	typedef GeometryMapType::Iterator        GeometryMapIterator;

	typedef SimplexMeshType::PointType PointType;
	typedef PointType::VectorType      VectorType;
	//typedef itk::Vector< FLOATTYPE, Dimension > NormalVectorType;

	char sampleFileName[1024];

	std::ofstream sampleFile;

	void extractProfile( 
		GradientInterpolatorType * gradientInterpolator,
		OriginalInterpolatorType * intensityInterpolator,
		SimplexMeshType* livermesh,
		PROFILE_CATEGORY category)
	{
		PointsContainerConstPointer points  = livermesh->GetPoints();
		GeometryMapIterator geoIt = livermesh->GetGeometryData()->Begin();
		GeometryMapIterator geoItEnd = livermesh->GetGeometryData()->End();

		typedef GradientInterpolatorType::InputImageType GradientImageType;
		GradientImageType::ConstPointer gradImage = gradientInterpolator->GetInputImage();

		g_liverCentroid.CastFrom(km::getMeshCentroid<SimplexMeshType>( livermesh ));

		SimplexMeshGeometry *geodata;
		while (geoIt!=geoItEnd)
		{
			SimplexMeshType::PointIdentifier idx = geoIt.Index();
			geodata = geoIt.Value();

			SimplexMeshType::PointType mpoint = points->GetElement( idx );
			OriginalImageType::PointType ipoint;
			ipoint.CastFrom( mpoint );

			// compute normal
			VectorType normal;
			normal.Set_vnl_vector(geodata->normal.Get_vnl_vector());

			//提取向法线内部方向shift后的位置的Profile
			OriginalImageType::PointType ipoint_inside = ipoint;
			for (int ttt=0;ttt<NUMBER_OF_INSIDE_PER_POINT;ttt++)
			{
				ipoint_inside -= normal*SHIFT_INSIDE;

				std::vector<double> feature;
				km::extractFeature<GradientInterpolatorType, OriginalInterpolatorType>(
					gradientInterpolator, intensityInterpolator, geodata, ipoint_inside, feature, category, true);

				sampleFile << IPClass;
				for (int i=0;i<feature.size();i++)
				{
					sampleFile << " " << feature[i];
				}
				sampleFile << " " << std::endl;
			}

			//提取当前位置的Profile
			OriginalImageType::PointType ipoint_boundary = ipoint;
			for (int ttt=0;ttt<NUMBER_OF_BOUNDARY_PER_POINT;ttt++)
			{
				ipoint_boundary += normal*((float)rand()/RAND_MAX)*0.25;

				std::vector<double> feature;
				km::extractFeature<GradientInterpolatorType, OriginalInterpolatorType>(
					gradientInterpolator, intensityInterpolator, geodata, ipoint_boundary, feature, category, true);

				sampleFile << BPClass;
				for (int i=0;i<feature.size();i++)
				{
					sampleFile << " " << feature[i];
				}
				sampleFile << " " << std::endl;
			}

			//提取向法线外部方向shift后的位置的Profile
			OriginalImageType::PointType ipoint_outside = ipoint;
			for (int ttt=0;ttt<NUMBER_OF_OUTSIDE_PER_POINT;ttt++)
			{
				//提取向法线内部方向shift后的位置的Profile
				ipoint_outside += normal*SHIFT_OUTSIDE;

				std::vector<double> feature;
				km::extractFeature<GradientInterpolatorType, OriginalInterpolatorType>(
					gradientInterpolator, intensityInterpolator, geodata, ipoint_outside, feature, category, true);

				sampleFile << OPClass;
				for (int i=0;i<feature.size();i++)
				{
					sampleFile << " " << feature[i];
				}
				sampleFile << " " << std::endl;
			}

			geoIt++;
		}
	}

	int trainAppearance( 
		const char* origlistfile,
		const char* meshlistfile,
		const char* referencegeometryfile,
		const char* outputdir,
		PROFILE_CATEGORY profile_category)
	{
		KM_DEBUG_INFO( "======================================" );
		KM_DEBUG_PRINT( "original image list", origlistfile );
		KM_DEBUG_PRINT( "shape mesh list", meshlistfile );
		KM_DEBUG_PRINT( "reference geometry", referencegeometryfile );
		KM_DEBUG_PRINT( "output directory", outputdir );
		KM_DEBUG_PRINT( "profile category", category2string(profile_category) );

		srand(time(0));

		//////////////////////////////////////////////////////////////////////////
		//导入训练数据文件
		KM_DEBUG_INFO("Import Datalist");

		StringVectorType origDataFileNamesList;
		StringVectorType shapeMeshFileNamesList;

		int numberOfGrayData = km::getDataList(origlistfile, origDataFileNamesList);
		int numberOfShapeMesh= km::getDataList(meshlistfile, shapeMeshFileNamesList);

		int numberOfData = std::min( numberOfGrayData, numberOfShapeMesh );

		if(numberOfGrayData!=numberOfShapeMesh)
		{
			std::cout<<"number of gray data "<<numberOfGrayData<<" does not match the number of shape mesh "<<numberOfShapeMesh<< "!"<<std::endl;
		}
		KM_DEBUG_PRINT( "Number Of Data:" , numberOfData);

		if (numberOfData<1)
		{
			return EXIT_FAILURE;
		}

		int refidx = 7;

		refidx = std::min( numberOfData-1, refidx );

		SimplexMeshType::Pointer reflivermesh = km::readMesh<SimplexMeshType>( shapeMeshFileNamesList[refidx] );
		OriginalImageType::Pointer reforigimage = km::readImage<OriginalImageType>( origDataFileNamesList[refidx] );

		km::GeometryImageType::Pointer geometryImage = km::readImage<km::GeometryImageType>(referencegeometryfile);
		km::loadSimplexMeshGeometryData<SimplexMeshType>(geometryImage, reflivermesh);
		km::ComputeGeometry<SimplexMeshType>( reflivermesh, true );

		const unsigned int numberOfPoints = reflivermesh->GetNumberOfPoints();
		KM_DEBUG_PRINT( "Number of Points ", numberOfPoints );

		std::stringstream ss0;
		ss0 << outputdir << "\\samples_" << category2string(profile_category) << ".txt";

		sprintf(sampleFileName, "%s", ss0.str().c_str());
		
		std::cout<<"Output sample file: "<< sampleFileName << std::endl;

		sampleFile.open (sampleFileName, std::ofstream::out | std::ofstream::app);

		for (unsigned int i=0;i<numberOfData;i++)
		{
			KM_DEBUG_PRINT("Start to extract appearance: ", i+1);

			SimplexMeshType::Pointer livermesh = km::readMesh<SimplexMeshType>( shapeMeshFileNamesList[i] );
			km::loadSimplexMeshGeometryData<SimplexMeshType>(geometryImage, livermesh);

			km::ComputeGeometry<SimplexMeshType>( livermesh, true );

			//读入图像
			KM_DEBUG_INFO("Reading images..");
			std::string origimagefile = origDataFileNamesList[i];
			OriginalImageType::Pointer origimage = km::readImage<OriginalImageType>(origimagefile);

			KM_DEBUG_INFO("Down-resampling to speed up..");
			SpacingType downspac;
			downspac.Fill( RESAMPLE_SPACING );
			origimage = km::resampleImage<OriginalImageType>( origimage, downspac );

			//KM_DEBUG_INFO("Histogram matching..");
			//origimage = km::histogramMatch<OriginalImageType>( origimage, reforigimage );
			KM_DEBUG_INFO("Smoothing..");
			origimage = km::minMaxSmooth<OriginalImageType>( origimage, 3, 1.0, 1.0 );
			
			//KM_DEBUG_INFO("Thresholding..");
			//origimage = km::thresholdImage<OriginalImageType>( origimage, 0, 400, 0 );

			KM_DEBUG_INFO( "Estimate liver threshold value" );

			g_phase = TRAINING;

			KM_DEBUG_INFO("Detect liver intensity range including tumor..");
			g_liverThresholds.clear();
			//km::detectLiverIntensityRangeWithoutTumor<OriginalImageType, SimplexMeshType>( origimage, livermesh, g_liverThresholds );
			km::detectLiverIntensityRangeIncludingTumor<OriginalImageType, SimplexMeshType>( origimage, livermesh, g_liverThresholds );

			for(int k=0;k<g_liverThresholds.size();k++)
			{
				std::cout<<g_liverThresholds[k].first<<","<<g_liverThresholds[k].second<<std::endl;
			}

			KM_DEBUG_INFO( "Estimate liver threshold value DONE!" );

			GradientImageType::Pointer gradimage = km::calculateRecursiveGradientImage<OriginalImageType, GradientImageType>( origimage, SIGMA );
			GradientInterpolatorType::Pointer gradInterpolator = GradientInterpolatorType::New();
			gradInterpolator->SetInputImage( gradimage );

			OriginalInterpolatorType::Pointer intensityInterpolator = OriginalInterpolatorType::New();
			intensityInterpolator->SetInputImage( origimage );

			KM_DEBUG_INFO( "Extract profile.." );
			extractProfile( gradInterpolator, intensityInterpolator, livermesh, profile_category );

			KM_DEBUG_PRINT( "Appearance generating done: ", i+1 );
		}

		sampleFile.close();

		//////////////////////////////////////////////////////////////////////////
		if(true)
		{
			//Write profiles.
			std::stringstream ss;
			ss << outputdir << "\\profile_" << category2string(profile_category) << ".h5";
			KM_DEBUG_PRINT( "Output intensity profile file: ", ss.str().c_str() );
			km::KNNProfileClassifier knnClassifier(profile_category);
			knnClassifier.setShapeNumber(numberOfData);
			knnClassifier.setShapePointsNumber(numberOfPoints);
			knnClassifier.loadSamples( sampleFileName );
			knnClassifier.save( ss.str().c_str() );

			ss.str( "" );
			ss.clear();

			ss << outputdir << "\\AdaboostClassifier_" << category2string(profile_category) << ".h5";
			KM_DEBUG_PRINT( "Output non-clustered adaboost classifier file: ", ss.str().c_str() );
			km::AdaboostProfileClassifier adaboostClassifier(knnClassifier.profile_category);
			adaboostClassifier.train( knnClassifier, outputdir );
			adaboostClassifier.save(ss.str().c_str());

			//Generate classifier map mesh.
			KM_DEBUG_INFO("Generate classifier map mesh.");
			km::AdaboostProfileClassifier::FloatVectorType error_map;
			adaboostClassifier.generateErrorMap(knnClassifier, error_map);
			km::assigneMesh<SimplexMeshType>( reflivermesh, 0.0 );
			for (int pid=0;pid<reflivermesh->GetNumberOfPoints();pid++)
			{
				reflivermesh->SetPointData( pid, error_map[pid] );
			}

			km::smoothMeshData<SimplexMeshType>(reflivermesh, 2);

			ss.str( "" );
			ss.clear();

			ss << outputdir << "\\ErrorMap_" << category2string(profile_category) << ".vtk";
			km::writeMesh<SimplexMeshType>( ss.str().c_str(), reflivermesh );
		}

		return EXIT_SUCCESS;
	}
}




#endif