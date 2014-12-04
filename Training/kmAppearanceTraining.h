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

#include "kmProfileClassifier.h"
#include "kmProfileExtractor.h"
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
			ss0 << outputdir << "\\samples_" << category2string(category_) << ".txt";
			sprintf(sampleTxtFilename, "%s", ss0.str().c_str());

			ss0.str( "" );
			ss0.clear();

			ss0 << outputdir << "\\samples_" << category2string(category_) << ".h5";
			sprintf(sampleHd5Filename, "%s", ss0.str().c_str());

			ss0.str( "" );
			ss0.clear();

			ss0 << outputdir << "\\AdaboostClassifier_" << category2string(category_) << ".h5";
			sprintf(classifierHd5FileName, "%s", ss0.str().c_str());
		}

		~ProfileUnit()
		{
			std::cout<<"De-constructor of ProfileUnit.."<<std::endl;
			//delete[] sampleTxtFilename;
			//delete[] sampleHd5Filename;
			//delete[] classifierHd5FileName;
			std::cout<<"De-constructor end.."<<std::endl;
		}

		void openTxtStream()
		{
			if (!sampleFile.is_open())
			{
				sampleFile.open (sampleTxtFilename, std::ofstream::out | std::ofstream::app);
			}
		}

		void closeTxtSteam()
		{
			if (sampleFile.is_open())
			{
				sampleFile.close();
			}
		}
	};
	
	int trainAppearance( 
		const char* origlistfile,
		const char* meshlistfile,
		const char* referencegeometryfile,
		const char* outputdir/*,
		PROFILE_CATEGORY profile_category*/)
	{
		KM_DEBUG_INFO( "======================================" );
		KM_DEBUG_PRINT( "original image list", origlistfile );
		KM_DEBUG_PRINT( "shape mesh list", meshlistfile );
		KM_DEBUG_PRINT( "reference geometry", referencegeometryfile );
		KM_DEBUG_PRINT( "output directory", outputdir );
		//KM_DEBUG_PRINT( "profile category", category2string(profile_category) );

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

		ProfileUnit profileUnitPlain(PLAIN, outputdir);
		ProfileUnit profileUnitLiver(LIVER, outputdir);
		ProfileUnit profileUnitBoundary(BOUNDARY, outputdir);
		ProfileUnit profileUnitCoordinate(COORDINATE, outputdir);

		profileUnitPlain.openTxtStream();
		profileUnitLiver.openTxtStream();
		profileUnitBoundary.openTxtStream();
		profileUnitCoordinate.openTxtStream();

		int refidx = 1;

		refidx = std::min( numberOfData-1, refidx );

		SimplexMeshType::Pointer reflivermesh = km::readMesh<SimplexMeshType>( shapeMeshFileNamesList[refidx] );
		OriginalImageType::Pointer reforigimage = km::readImage<OriginalImageType>( origDataFileNamesList[refidx] );

		km::GeometryImageType::Pointer geometryImage = km::readImage<km::GeometryImageType>(referencegeometryfile);
		km::loadSimplexMeshGeometryData<SimplexMeshType>(geometryImage, reflivermesh);
		km::ComputeGeometry<SimplexMeshType>( reflivermesh, true );

		const unsigned int numberOfPoints = reflivermesh->GetNumberOfPoints();
		KM_DEBUG_PRINT( "Number of Points ", numberOfPoints );

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

				std::vector<double> features;

				//提取向法线内部方向shift后的位置的Profile
				OriginalImageType::PointType ipoint_inside = ipoint;
				for (int ttt=0;ttt<NUMBER_OF_INSIDE_PER_POINT;ttt++)
				{
					ipoint_inside -= normal*SHIFT_INSIDE;

					//Plain
					profileExtractor.extractFeatureSet(features, profileUnitPlain.category, geodata, ipoint_inside);

					profileUnitPlain.sampleFile << IPClass;
					for (int i=0;i<features.size();i++)
					{
						profileUnitPlain.sampleFile << " " << features[i];
					}
					profileUnitPlain.sampleFile << " " << std::endl;

					//Liver
					profileExtractor.extractFeatureSet(features, profileUnitLiver.category, geodata, ipoint_inside);

					profileUnitLiver.sampleFile << IPClass;
					for (int i=0;i<features.size();i++)
					{
						profileUnitLiver.sampleFile << " " << features[i];
					}
					profileUnitLiver.sampleFile << " " << std::endl;

					//Boundary
					profileExtractor.extractFeatureSet(features, profileUnitBoundary.category, geodata, ipoint_inside);

					profileUnitBoundary.sampleFile << IPClass;
					for (int i=0;i<features.size();i++)
					{
						profileUnitBoundary.sampleFile << " " << features[i];
					}
					profileUnitBoundary.sampleFile << " " << std::endl;
				}

				//提取当前位置的Profile
				OriginalImageType::PointType ipoint_boundary = ipoint - normal*SHIFT_BOUNDARY;
				for (int ttt=0;ttt<NUMBER_OF_BOUNDARY_PER_POINT;ttt++)
				{
					ipoint_boundary = ipoint; //+= normal*(SHIFT_BOUNDARY*2/(NUMBER_OF_BOUNDARY_PER_POINT-1));

					//Plain
					profileExtractor.extractFeatureSet(features, profileUnitPlain.category, geodata, ipoint_boundary);

					profileUnitPlain.sampleFile << BPClass;
					for (int i=0;i<features.size();i++)
					{
						profileUnitPlain.sampleFile << " " << features[i];
					}
					profileUnitPlain.sampleFile << " " << std::endl;

					//Boundary
					profileExtractor.extractFeatureSet(features, profileUnitBoundary.category, geodata, ipoint_boundary);

					profileUnitBoundary.sampleFile << BPClass;
					for (int i=0;i<features.size();i++)
					{
						profileUnitBoundary.sampleFile << " " << features[i];
					}
					profileUnitBoundary.sampleFile << " " << std::endl;

					//Coordinate
					profileExtractor.extractFeatureSet(features, profileUnitCoordinate.category, geodata, ipoint_boundary);

					profileUnitCoordinate.sampleFile << BPClass;
					for (int i=0;i<features.size();i++)
					{
						profileUnitCoordinate.sampleFile << " " << features[i];
					}
					profileUnitCoordinate.sampleFile << " " << std::endl;
				}

				//提取向法线外部方向shift后的位置的Profile
				OriginalImageType::PointType ipoint_outside = ipoint;
				for (int ttt=0;ttt<NUMBER_OF_OUTSIDE_PER_POINT;ttt++)
				{
					//提取向法线内部方向shift后的位置的Profile
					ipoint_outside += normal*SHIFT_OUTSIDE;

					//Plain
					profileExtractor.extractFeatureSet(features, profileUnitPlain.category, geodata, ipoint_outside);

					profileUnitPlain.sampleFile << OPClass;
					for (int i=0;i<features.size();i++)
					{
						profileUnitPlain.sampleFile << " " << features[i];
					}
					profileUnitPlain.sampleFile << " " << std::endl;

					//Liver
					profileExtractor.extractFeatureSet(features, profileUnitLiver.category, geodata, ipoint_outside);

					profileUnitLiver.sampleFile << OPClass;
					for (int i=0;i<features.size();i++)
					{
						profileUnitLiver.sampleFile << " " << features[i];
					}
					profileUnitLiver.sampleFile << " " << std::endl;

					//Boundary
					profileExtractor.extractFeatureSet(features, profileUnitBoundary.category, geodata, ipoint_outside);

					profileUnitBoundary.sampleFile << OPClass;
					for (int i=0;i<features.size();i++)
					{
						profileUnitBoundary.sampleFile << " " << features[i];
					}
					profileUnitBoundary.sampleFile << " " << std::endl;
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
		km::KNNProfileClassifier knnClassifierPlain;
		knnClassifierPlain.setShapeNumber(numberOfData);
		knnClassifierPlain.setShapePointsNumber(numberOfPoints);
		knnClassifierPlain.loadSamples( profileUnitPlain.sampleTxtFilename );
		knnClassifierPlain.cluster(9);
		knnClassifierPlain.save( profileUnitPlain.sampleHd5Filename );
		{
			km::assigneMesh<SimplexMeshType>( reflivermesh, 0.0 );
			for (int pid=0;pid<reflivermesh->GetNumberOfPoints();pid++)
			{
				reflivermesh->SetPointData( pid, knnClassifierPlain.cluster_labels[pid][0] );
			}

			std::stringstream ss;
			ss << outputdir << "\\ClusteredMeshPlain.vtk";
			km::writeMesh<SimplexMeshType>( ss.str().c_str(), reflivermesh );
		}

		//Output liver sample & classifier.
		km::KNNProfileClassifier knnClassifierLiver;
		knnClassifierLiver.setShapeNumber(numberOfData);
		knnClassifierLiver.setShapePointsNumber(numberOfPoints);
		knnClassifierLiver.loadSamples( profileUnitLiver.sampleTxtFilename );
		knnClassifierLiver.copyCluster(knnClassifierPlain);
		knnClassifierLiver.save( profileUnitLiver.sampleHd5Filename );
		{
			km::ProfileClassifier adaboostClassifierLiver( profileUnitLiver.category );
			adaboostClassifierLiver.train( knnClassifierLiver, outputdir );
			adaboostClassifierLiver.save( profileUnitLiver.classifierHd5FileName );

			km::ProfileClassifier::IntFloatMapType errorMap;
			adaboostClassifierLiver.test(knnClassifierLiver, errorMap);
			km::assigneMesh<SimplexMeshType>( reflivermesh, 0.0 );
			for (int pid=0;pid<reflivermesh->GetNumberOfPoints();pid++)
			{
				reflivermesh->SetPointData( pid, errorMap[knnClassifierLiver.cluster_labels[pid][0]] );
			}

			std::stringstream ss;
			ss << outputdir << "\\ErrorMap_Liver.vtk";
			km::writeMesh<SimplexMeshType>( ss.str().c_str(), reflivermesh );
		}

		//Output coordinate sample.
		//Write profiles.
		km::KNNProfileClassifier knnClassifierCoordinate;
		knnClassifierCoordinate.setShapeNumber(numberOfData);
		knnClassifierCoordinate.setShapePointsNumber(numberOfPoints);
		knnClassifierCoordinate.loadSamples( profileUnitCoordinate.sampleTxtFilename );
		knnClassifierCoordinate.cluster(9);
		knnClassifierCoordinate.save( profileUnitCoordinate.sampleHd5Filename );
		{
			km::assigneMesh<SimplexMeshType>( reflivermesh, 0.0 );
			for (int pid=0;pid<reflivermesh->GetNumberOfPoints();pid++)
			{
				reflivermesh->SetPointData( pid, knnClassifierCoordinate.cluster_labels[pid][0] );
			}

			std::stringstream ss;
			ss << outputdir << "\\ClusteredMeshCoordinate.vtk";
			km::writeMesh<SimplexMeshType>( ss.str().c_str(), reflivermesh );
		}

		//Output Boundary sample & classifier.
		km::KNNProfileClassifier knnClassifierBoundary;
		knnClassifierBoundary.setShapeNumber(numberOfData);
		knnClassifierBoundary.setShapePointsNumber(numberOfPoints);
		knnClassifierBoundary.loadSamples( profileUnitBoundary.sampleTxtFilename );
		knnClassifierBoundary.copyCluster(knnClassifierPlain);
		knnClassifierBoundary.save( profileUnitBoundary.sampleHd5Filename );
		{
			km::ProfileClassifier adaboostClassifierBoundary( profileUnitBoundary.category );
			adaboostClassifierBoundary.train( knnClassifierBoundary, outputdir );
			adaboostClassifierBoundary.save( profileUnitBoundary.classifierHd5FileName );

			km::ProfileClassifier::IntFloatMapType errorMap;
			adaboostClassifierBoundary.test(knnClassifierBoundary, errorMap);
			km::assigneMesh<SimplexMeshType>( reflivermesh, 0.0 );
			for (int pid=0;pid<reflivermesh->GetNumberOfPoints();pid++)
			{
				reflivermesh->SetPointData( pid, errorMap[knnClassifierBoundary.cluster_labels[pid][0]] );
			}

			std::stringstream ss;
			ss << outputdir << "\\ErrorMap_Boundary.vtk";
			km::writeMesh<SimplexMeshType>( ss.str().c_str(), reflivermesh );
		}

		return EXIT_SUCCESS;
	}
}




#endif