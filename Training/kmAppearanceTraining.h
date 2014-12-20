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

#include "kmCommon.h"

#define EXTRACT_FLAG true

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
	
	struct ProfileUnit
	{
		PROFILE_CATEGORY category;
		char sampleTxtFilename[1024];
		char sampleHd5Filename[1024];
		char classifierHd5FileName[1024];
		std::ofstream sampleFile;

		ProfileUnit(PROFILE_CATEGORY category_, const char* outputdir)
		{
			category = category_;
			std::stringstream ss0;
			ss0 << outputdir << "\\samples_" << ProfileCategoryUtils::category2string(category_) << ".txt";
			sprintf(sampleTxtFilename, "%s", ss0.str().c_str());

			ss0.str( "" );
			ss0.clear();

			ss0 << outputdir << "\\samples_" << ProfileCategoryUtils::category2string(category_) << ".h5";
			sprintf(sampleHd5Filename, "%s", ss0.str().c_str());

			ss0.str( "" );
			ss0.clear();

			ss0 << outputdir << "\\AdaboostClassifier_" << ProfileCategoryUtils::category2string(category_) << ".h5";
			sprintf(classifierHd5FileName, "%s", ss0.str().c_str());
		}

		~ProfileUnit()
		{
		}

		void openTxtStream()
		{
			if (!sampleFile.is_open()){
				sampleFile.open (sampleTxtFilename, std::ofstream::out | std::ofstream::app);
			}
		}

		void closeTxtSteam()
		{
			if (sampleFile.is_open()){
				sampleFile.close();
			}
		}

		void writeLine(int pointID, int classLabel, const std::vector<FLOATTYPE>& features)
		{
			sampleFile << pointID << " " << classLabel;
			for (int k=0;k<features.size();k++)
			{
				sampleFile << " " << features[k];
			}
			sampleFile << " " << std::endl;
		}
	};
	
	int trainAppearance( 
		const char* origlistfile,
		const char* meshlistfile,
		const char* referencegeometryfile,
		const char* outputdir)
	{
		KM_DEBUG_INFO( "======================================" );
		KM_DEBUG_PRINT( "original image list", origlistfile );
		KM_DEBUG_PRINT( "shape mesh list", meshlistfile );
		KM_DEBUG_PRINT( "reference geometry", referencegeometryfile );
		KM_DEBUG_PRINT( "output directory", outputdir );

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

		int refidx = 1;
		refidx = std::min( numberOfData-1, refidx );

		SimplexMeshType::Pointer reflivermesh = km::readMesh<SimplexMeshType>( shapeMeshFileNamesList[refidx] );
		OriginalImageType::Pointer reforigimage = km::readImage<OriginalImageType>( origDataFileNamesList[refidx] );

		km::GeometryImageType::Pointer geometryImage = km::readImage<km::GeometryImageType>(referencegeometryfile);
		km::loadSimplexMeshGeometryData<SimplexMeshType>(geometryImage, reflivermesh);
		km::ComputeGeometry<SimplexMeshType>( reflivermesh, true );

		const unsigned int numberOfPoints = reflivermesh->GetNumberOfPoints();
		KM_DEBUG_PRINT( "Number of Points ", numberOfPoints );

		ProfileUnit profileUnitPlain(PLAIN, outputdir);
		ProfileUnit profileUnitLiver(LIVER, outputdir);
		ProfileUnit profileUnitBoundary(BOUNDARY, outputdir);
		ProfileUnit profileUnitCoordinate(COORDINATE, outputdir);
		profileUnitPlain.openTxtStream();
		profileUnitLiver.openTxtStream();
		profileUnitBoundary.openTxtStream();
		profileUnitCoordinate.openTxtStream();

		for (unsigned int i=0;i<numberOfData;i++)
		{
			KM_DEBUG_PRINT("Start to extract appearance: ", i+1);

			SimplexMeshType::Pointer livermesh = km::readMesh<SimplexMeshType>( shapeMeshFileNamesList[i] );
			km::loadSimplexMeshGeometryData<SimplexMeshType>(geometryImage, livermesh);
			km::ComputeGeometry<SimplexMeshType>( livermesh, true );

			g_liverCentroid.CastFrom(km::getMeshCentroid<SimplexMeshType>( livermesh ));

			//读入图像
			KM_DEBUG_INFO("Reading images..");
			std::string origimagefile = origDataFileNamesList[i];
			OriginalImageType::Pointer origimage = km::readImage<OriginalImageType>(origimagefile);

			KM_DEBUG_INFO("Down-resampling to speed up..");
			SpacingType downspac;
			downspac.Fill( RESAMPLE_SPACING );
			origimage = km::resampleImage<OriginalImageType>( origimage, downspac );

			KM_DEBUG_INFO("Smoothing..");
			origimage = km::minMaxSmooth<OriginalImageType>( origimage, 3, 1.0, 1.0 );

			KM_DEBUG_INFO( "Estimate liver threshold value" );

			g_phase = TRAINING;

			KM_DEBUG_INFO("Detect liver intensity range..");
			g_liverThresholds.clear();
			km::detectLiverIntensityRangeWithoutTumor<OriginalImageType, SimplexMeshType>( origimage, livermesh, g_liverThresholds );

			for(int k=0;k<g_liverThresholds.size();k++)
			{
				std::cout<<g_liverThresholds[k].first<<","<<g_liverThresholds[k].second<<std::endl;
			}

			//KM_DEBUG_INFO( "Estimate liver threshold value DONE!" );
			km::ProfileExtractor<OriginalImageType> profileExtractor;
			profileExtractor.setImage(origimage);
			profileExtractor.enableCache(false);

			KM_DEBUG_INFO( "Extract profile.." );
			GeometryMapIterator geoIt = livermesh->GetGeometryData()->Begin();
			GeometryMapIterator geoItEnd = livermesh->GetGeometryData()->End();
			SimplexMeshGeometry *geodata;
			while (geoIt!=geoItEnd)
			{
				SimplexMeshType::PointIdentifier idx = geoIt.Index();
				geodata = geoIt.Value();

				SimplexMeshType::PointType mpoint = livermesh->GetPoints()->GetElement( idx );
				OriginalImageType::PointType ipoint;
				ipoint.CastFrom( mpoint );

				// compute normal
				VectorType normal;
				normal.Set_vnl_vector(geodata->normal.Get_vnl_vector());

				std::vector<PointType> points_inside;
				std::vector<PointType> points_boundary;
				std::vector<PointType> points_outside;
				for (int t=1;t<=NUMBER_OF_INSIDE_PER_POINT;t++){
					points_inside.push_back(ipoint-normal*SHIFT_INSIDE*t);
				}
				for (int t=0;t<NUMBER_OF_BOUNDARY_PER_POINT;t++){
					points_boundary.push_back(ipoint+normal*((rand()%10-5.0)/10.0)*SHIFT_BOUNDARY);
				}
				for (int t=1;t<=NUMBER_OF_OUTSIDE_PER_POINT;t++){
					points_outside.push_back(ipoint+normal*SHIFT_OUTSIDE*t);
				}

				std::vector<FLOATTYPE> features;
				//提取向法线内部方向shift后的位置的Profile
				for (int t=0;t<points_inside.size();t++)
				{
					PointType ipoint_inside = points_inside[t];

					//Liver
					profileExtractor.extractFeatureSet(features, profileUnitLiver.category, geodata, ipoint_inside);
					profileUnitLiver.writeLine(idx, IPClass, features);

					//Boundary
					profileExtractor.extractFeatureSet(features, profileUnitBoundary.category, geodata, ipoint_inside);
					profileUnitBoundary.writeLine(idx, IPClass, features);

					//Plain
					if (t<1)
					{
						profileExtractor.extractFeatureSet(features, profileUnitPlain.category, geodata, ipoint_inside);
						profileUnitPlain.writeLine(idx, IPClass, features);
					}
					
				}

				//提取当前位置的Profile
				for (int t=0;t<points_boundary.size();t++)
				{
					PointType ipoint_boundary = points_boundary[t];

					//Boundary
					profileExtractor.extractFeatureSet(features, profileUnitBoundary.category, geodata, ipoint_boundary);
					profileUnitBoundary.writeLine(idx, BPClass, features);

					//Plain
					if (t<1)
					{
						profileExtractor.extractFeatureSet(features, profileUnitPlain.category, geodata, ipoint_boundary);
						profileUnitPlain.writeLine(idx, BPClass, features);
					}

					//Coordinate
					profileExtractor.extractFeatureSet(features, profileUnitCoordinate.category, geodata, ipoint_boundary);
					profileUnitCoordinate.writeLine(idx, BPClass, features);
				}

				//提取向法线外部方向shift后的位置的Profile
				for (int t=0;t<points_outside.size();t++)
				{
					//提取向法线内部方向shift后的位置的Profile
					PointType ipoint_outside = points_outside[t];

					//Liver
					profileExtractor.extractFeatureSet(features, profileUnitLiver.category, geodata, ipoint_outside);
					profileUnitLiver.writeLine(idx, OPClass, features);

					//Plain
					if (t<1)
					{
						profileExtractor.extractFeatureSet(features, profileUnitPlain.category, geodata, ipoint_outside);
						profileUnitPlain.writeLine(idx, OPClass, features);
					}
				}

				geoIt++;
			}

			KM_DEBUG_PRINT( "Appearance generating done: ", i+1 );
		}

		profileUnitPlain.closeTxtSteam();
		profileUnitLiver.closeTxtSteam();
		profileUnitBoundary.closeTxtSteam();
		profileUnitCoordinate.closeTxtSteam();

		//Output plain sample.
		//Write profiles.
		km::ProfileContainer profilesPlain;
		profilesPlain.setShapeNumber(numberOfData);
		profilesPlain.setShapePointsNumber(numberOfPoints);
		profilesPlain.loadSamples( profileUnitPlain.sampleTxtFilename );
		profilesPlain.cluster(9);
		profilesPlain.save( profileUnitPlain.sampleHd5Filename );
		{
			km::ProfileContainerUtils<SimplexMeshType>::assignClusterLabels(profilesPlain, reflivermesh);
			std::stringstream ss;
			ss << outputdir << "\\ClusteredMeshPlain.vtk";
			km::writeMesh<SimplexMeshType>( ss.str().c_str(), reflivermesh );
		}

		////Output coordinate sample.
		////Write profiles.
		//km::ProfileContainer profilesCoordinate;
		//profilesCoordinate.setShapeNumber(numberOfData);
		//profilesCoordinate.setShapePointsNumber(numberOfPoints);
		//profilesCoordinate.loadSamples( profileUnitCoordinate.sampleTxtFilename );
		//profilesCoordinate.cluster(9);
		//profilesCoordinate.save( profileUnitCoordinate.sampleHd5Filename );
		//{
		//	km::ProfileContainerUtils<SimplexMeshType>::assignClusterLabels(profilesCoordinate, reflivermesh);
		//	std::stringstream ss;
		//	ss << outputdir << "\\ClusteredMeshCoordinate.vtk";
		//	km::writeMesh<SimplexMeshType>( ss.str().c_str(), reflivermesh );
		//}

		//Output liver sample & classifier.
		km::ProfileContainer profilesLiver;
		profilesLiver.setShapeNumber(numberOfData);
		profilesLiver.setShapePointsNumber(numberOfPoints);
		profilesLiver.loadSamples( profileUnitLiver.sampleTxtFilename );
		profilesLiver.copyCluster(profilesPlain);
		profilesLiver.save( profileUnitLiver.sampleHd5Filename );
		{
			km::ProfileClassifier adaboostClassifierLiver;
			adaboostClassifierLiver.train( profilesLiver, profileUnitLiver.category, outputdir );
			adaboostClassifierLiver.save( profileUnitLiver.classifierHd5FileName );
			adaboostClassifierLiver.test(profilesLiver);
		}

		//Output Boundary sample & classifier.
		km::ProfileContainer profilesBoundary;
		profilesBoundary.setShapeNumber(numberOfData);
		profilesBoundary.setShapePointsNumber(numberOfPoints);
		profilesBoundary.loadSamples( profileUnitBoundary.sampleTxtFilename );
		profilesBoundary.copyCluster(profilesPlain);
		profilesBoundary.save( profileUnitBoundary.sampleHd5Filename );
		{
			km::ProfileClassifier adaboostClassifierBoundary;
			adaboostClassifierBoundary.train( profilesBoundary, profileUnitBoundary.category, outputdir );
			adaboostClassifierBoundary.save( profileUnitBoundary.classifierHd5FileName );
			adaboostClassifierBoundary.test(profilesBoundary);
		}

		return EXIT_SUCCESS;
	}
}




#endif