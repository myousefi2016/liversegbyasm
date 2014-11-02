#include "LiverSegmentationAPI.h"

#include "kmKNNProfileClassifier-FLANN.h"
#include "AdaSegmentAPI.h"
#include "kmUtility.h"
#include "kmVtkItkUtility.h"
#include "kmModelFitting.h"
#include "kmRegistration.h"
#include "kmSegmentation.h"
#include "kmGlobal.h"
#include "kmProcessing.h"

void Notifier::notifyMesh(MeshType *mesh)
{
	//do nothing
}

void Notifier::notifyImage(UCharImageType *image)
{
	//do nothing
}


bool WRITE_MIDDLE_RESULT = true;

char* outputdir = NULL;

template<class FloatImageType, class MeshType, class GradientInterpolatorType, class IntensityInterpolatorType>
void
generateBestPointSet(
						 AdaboostProfileClassifier * classifier,
						 AdaboostProfileClassifier * classifier_liver,
						 const typename GradientInterpolatorType* gradientInterpolator, 
						 const typename IntensityInterpolatorType* intensityInterpolator,
						 typename MeshType* outputMesh,
						 const typename MeshType* liverMesh,
						 const typename MeshType* shapeMesh,
						 const typename MeshType* varianceMap,
						 const typename MeshType* errorMap,
						 const double spacingDistForTest = 1.5,
						 const int numberOfPointsForTest = 5)
{
	typedef MeshType::PointsContainer PointsContainer;
	typedef MeshType::PointsContainerConstPointer PointsContainerConstPointer;
	typedef MeshType::PointsContainerPointer PointsContainerPointer;
	typedef MeshType::PointsContainerIterator PointsContainerIterator;
	typedef MeshType::PointDataContainer PointDataContainer;
	typedef MeshType::GeometryMapType GeometryMapType;
	typedef GeometryMapType::Pointer GeometryMapPointer;
	typedef GeometryMapType::Iterator GeometryMapIterator;
	typedef itk::SimplexMeshGeometry::VectorType VectorType;
	typedef itk::SimplexMeshGeometry::PointType PointType;

	km::assigneMesh<MeshType>(outputMesh, 0.0);

	GeometryMapIterator geoIt = liverMesh->GetGeometryData()->Begin();
	GeometryMapIterator geoItEnd = liverMesh->GetGeometryData()->End();

	itk::SimplexMeshGeometry *geodata;
	while (geoIt!=geoItEnd)
	{
		MeshType::PointIdentifier idx = geoIt.Index();
		geodata = geoIt.Value();

		PointType mpoint = liverMesh->GetPoint(idx);
		VectorType normal;
		normal.Set_vnl_vector(geodata->normal.Get_vnl_vector());

		double best_offset = 0.0;
		double current_confidence = 0.0;
		//Check the direction we should search for the next best position.
		double direction = 1.0;
		{
			std::vector<FLOATTYPE> feature;
			km::extractFeature<GradientInterpolatorType, IntensityInterpolatorType>(
				gradientInterpolator,
				intensityInterpolator,
				geodata,
				mpoint,
				feature,
				classifier_liver->profile_category);

			double current_confidence = 2*classifier->classify( feature, idx ) - 1.0;
			//if ( current_confidence > 0 )
			//{
			//	direction = 1.0;
			//}
			//else
			//{
			//	direction = -1.0; //By default search towards inside.
			//}
			direction = current_confidence;
		}
		if (classifier->profile_category == LIVER)
		{
			double pointVariance = 1.0;
			varianceMap->GetPointData(idx, &pointVariance);

			double pointError = 0.0;
			//errorMap->GetPointData(idx, &pointError);
			//pointError = std::log( 3*pointError+1 );

			int tmpNumberOfPointsForTest = numberOfPointsForTest;//static_cast<int>(numberOfPointsForTest*pointVariance+3.0);

			//if (direction<0)
			//{
			//	tmpNumberOfPointsForTest = 0;
			//}

			std::vector<double> confidences_trace;
			confidences_trace.push_back(current_confidence);

			for(int i=0;i<tmpNumberOfPointsForTest;i++)
			{
				double offset = spacingDistForTest*(i+1)*direction;
				PointType pttest = mpoint + normal*offset;

				double confidence = 0.0;
				if (intensityInterpolator->IsInsideBuffer(pttest))
				{
					std::vector<FLOATTYPE> feature;
					km::extractFeature<GradientInterpolatorType, IntensityInterpolatorType>(
						gradientInterpolator,
						intensityInterpolator,
						geodata,
						pttest,
						feature,
						classifier->profile_category);

					confidence = classifier->classify(feature,idx);
					//let's bring in the error map.
					confidence = (1.0-pointError)*confidence+pointError*(1.0-confidence);

					confidence = 2*confidence - 1.0;
				}

				confidences_trace.push_back(confidence);
				double accumulated_confidences = 0.0;
				for (int k=confidences_trace.size()-2;k<confidences_trace.size();k++)
				{
					accumulated_confidences += confidences_trace[k];
				}
				if (accumulated_confidences*direction<0)
				{
					best_offset = offset - direction*spacingDistForTest;
					break;
				}
				else
				{
					best_offset = offset;		
				}

				if (confidence * direction > 0)
				{
					direction = confidence;
				}
			}
		}
		else if (classifier->profile_category == BOUNDARY) //Find best boundary point.
		{
			int tmpNumberOfPointsForTest = numberOfPointsForTest;

			std::vector<double> confidences_trace;
			std::vector<double> offset_trace;

			direction = 1.0;
			double offset = -1.0*spacingDistForTest*(tmpNumberOfPointsForTest/5.0);
			for(int i=0;i<tmpNumberOfPointsForTest;i++)
			{
				PointType pttest = mpoint + normal*offset;

				double confidence = 0.0;
				if (intensityInterpolator->IsInsideBuffer(pttest))
				{
					std::vector<FLOATTYPE> feature;
					km::extractFeature<GradientInterpolatorType, IntensityInterpolatorType>(
						gradientInterpolator,
						intensityInterpolator,
						geodata,
						pttest,
						feature,
						classifier->profile_category);

					confidence = classifier->classify( feature, idx );
				}

				confidences_trace.push_back(confidence);
				offset_trace.push_back(offset);

				offset = offset + spacingDistForTest;
			}

			double bestconfidence = 0.15;
			best_offset = 0.0;
			int N = 3; //accumulate confidence for local N confidences.
			for (int k=0;k<confidences_trace.size();k++)
			{
				double accumulated_confidences = 0.0;
				int local_idx_start = std::max(0, k-N/2);
				int local_idx_end = std::min(static_cast<int>(confidences_trace.size()-1), k+N/2);

				for (int t=local_idx_start;t<local_idx_end;t++)
				{
					accumulated_confidences += confidences_trace[t];
				}
				if (accumulated_confidences>bestconfidence)
				{
					bestconfidence = accumulated_confidences;
					best_offset = offset_trace[k];
				}
				//if (confidences_trace[k]>bestconfidence)
				//{
				//	bestconfidence = confidences_trace[k];
				//	best_offset = offset_trace[k];
				//}
			}
		}
		else
		{
			KM_DEBUG_ERROR( "Unknown profile category." );
			break;
		}

		outputMesh->SetPointData(idx, best_offset);

		geoIt++;
	}

	if (classifier->profile_category = LIVER)
	{
		//Remove noise point.
		km::smoothMeshData<MeshType>(outputMesh, 1);
	}

	geoIt = liverMesh->GetGeometryData()->Begin();
	geoItEnd = liverMesh->GetGeometryData()->End();

	while (geoIt!=geoItEnd)
	{
		unsigned int idx = geoIt.Index();
		geodata = geoIt.Value();

		VectorType normal;
		normal.Set_vnl_vector(geodata->normal.Get_vnl_vector());

		double offsetVal = 0.0;
		outputMesh->GetPointData(idx, &offsetVal);

		PointType oldPos = liverMesh->GetPoint(idx);
		outputMesh->SetPoint(idx, oldPos + normal*offsetVal);

		geoIt++;
	}
}

void LiverSeg( km::Notifier* notifier,
			  const char* inputImageFile,
			  const char* SSMFile,
			  const char* boundaryClassifierFile,
			  const char* liverClassifierFile,
			  const char* adaboostSegmentFile,
			  const char* geoFile,
			  const char* atlasImageFile,
			  const char* configFile,
			  const KM_STRATEGY strategy)
{
	if ( notifier == NULL )
	{
		std::cout<<"Notifier is NULL!"<<std::endl;
		return;
	}
	if ( inputImageFile == NULL )
	{
		std::cout<<"Notifier is NULL!"<<std::endl;
		return;
	}
	if ( SSMFile == NULL )
	{
		std::cout<<"SSMFile is NULL!"<<std::endl;
		return;
	}
	if ( boundaryClassifierFile == NULL )
	{
		std::cout<<"boundaryClassifierFile is NULL!"<<std::endl;
		return;
	}
	if ( liverClassifierFile == NULL )
	{
		std::cout<<"liverClassifierFile is NULL!"<<std::endl;
		return;
	}
	if ( adaboostSegmentFile == NULL )
	{
		std::cout<<"adaboostSegmentFile is NULL!"<<std::endl;
		return;
	}
	if ( geoFile == NULL )
	{
		std::cout<<"geoFile is NULL!"<<std::endl;
		return;
	}
	if ( atlasImageFile == NULL )
	{
		std::cout<<"atlasImageFile is NULL!"<<std::endl;
		return;
	}
	if ( configFile == NULL )
	{
		std::cout<<"configFile is NULL!"<<std::endl;
		return;
	}

	KM_DEBUG_INFO("Load config file...");
	km::loadConfig(configFile);

	KM_DEBUG_PRINT("Segmentation strategy", strategy);

	itk::TimeProbesCollectorBase chronometer;
	itk::MemoryProbesCollectorBase memorymeter;
	memorymeter.Start( "segmentation" );
	chronometer.Start( "segmentation" );

	g_phase = INITIALIZATION;

	memorymeter.Start( "initialization" );
	chronometer.Start( "initialization" );

	//���������

	typedef km::AdaboostProfileClassifier AdaboostProfileClassifierType;
	AdaboostProfileClassifierType adaboostProfileClassifier_Boundary(BOUNDARY);
	adaboostProfileClassifier_Boundary.load(boundaryClassifierFile);

	AdaboostProfileClassifierType adaboostProfileClassifier_Liver(LIVER);
	adaboostProfileClassifier_Liver.load(liverClassifierFile);

	//������ָ�ͼ��
	KM_DEBUG_INFO("read input image..");
	KM_DEBUG_INFO(inputImageFile);
	ShortImageType::Pointer inputImage = km::readImage<ShortImageType>( inputImageFile );
	km::shiftMinimum<ShortImageType>( inputImage, -1024 );

	KM_DEBUG_INFO("read atlas image..");
	ShortImageType::Pointer atlasImage = km::readImage<ShortImageType>( atlasImageFile );

	//Store the original image information for final segmentation mask generation.
	ShortImageType::Pointer inputCached = ShortImageType::New();
	inputCached->CopyInformation( inputImage );

	double minval = -1024;
	double maxval = 1024;

	SpacingType downspac;
	downspac.Fill( RESAMPLE_SPACING );

	inputImage = km::resampleImage<ShortImageType>( inputImage, downspac, minval );

	double zdist = RESAMPLE_SPACING * inputImage->GetLargestPossibleRegion().GetSize()[2];

	if ( zdist > 400 )
	{
		KM_DEBUG_INFO("roi locating by template matching..");
		inputImage = km::extractRoiByTemplateMatching<ShortImageType, ShortImageType>( inputImage, atlasImage );
	}

	if(WRITE_MIDDLE_RESULT)
	{
		km::writeImage<ShortImageType>( outputdir, "inputImage.nii.gz", inputImage );
	}

	//return;

	KM_DEBUG_INFO("histogram matching..");
	inputImage = km::histogramMatch<ShortImageType>( inputImage, atlasImage );

	inputImage = km::minMaxSmooth<ShortImageType>( inputImage, 3, 1.0, 1 );

	if(WRITE_MIDDLE_RESULT)
	{
		km::writeImage<ShortImageType>( outputdir, "inputImageSmoothed.nii.gz", inputImage );
	}

	//Calculate gradient image.
	GradientImageType::Pointer gradImage = km::calculateRecursiveGradientImage<ShortImageType, GradientImageType>( inputImage, SIGMA );

	km::calculateMinAndMax<ShortImageType>( inputImage, minval, maxval );
	KM_DEBUG_PRINT( "Min value", minval );
	KM_DEBUG_PRINT( "Max value", maxval );

	memorymeter.Stop( "initialization" );
	chronometer.Stop( "initialization" );

	FloatImageType::PointType liverCentroid;
	
	FloatImageType::Pointer probablityMap = FloatImageType::New();
	UCharImageType::Pointer liverMask  = UCharImageType::New();
	UCharImageType::SizeType radius;

	{
		g_phase = ADABOOST_REGION_SEGMENTATION;

		memorymeter.Start( "adboost segmentation" );
		chronometer.Start( "adboost segmentation" );

		KM_DEBUG_INFO("Start to generate liver probabilistic map by Adboost...");

		downspac.Fill( 4.0 );
		ShortImageType::Pointer downsampleInput = km::resampleImage<ShortImageType>( inputImage,downspac, minval );

		if(WRITE_MIDDLE_RESULT)
		{
			km::writeImage<ShortImageType>( outputdir, "downsampleInput.nii.gz", downsampleInput );
		}

		//Get abdominal mask
		ShortImageType::Pointer bodyImage = ShortImageType::New();
		UCharImageType::Pointer bodyMask = UCharImageType::New();

		KM_DEBUG_INFO("Generate body mask image...");
		km::removeAir<ShortImageType, UCharImageType>( downsampleInput, bodyMask, minval, true );

		if(WRITE_MIDDLE_RESULT)
		{
			km::writeImage<UCharImageType>( outputdir, "bodyMask.nii.gz", bodyMask );
		}
		
		//Locate liver
		KM_DEBUG_INFO( "Generate liver probilistic map.." );
		probablityMap = AdaSegment::adaSegment<FloatImageType, ShortImageType, UCharImageType>( downsampleInput, bodyMask, adaboostSegmentFile, 1, 300 );

		//probablityMap = km::medianSmooth<FloatImageType>( probablityMap, 1 );

		if(WRITE_MIDDLE_RESULT)
		{
			km::writeImage<FloatImageType>( outputdir, "probabilityMap.nii.gz", probablityMap );
		}

		KM_DEBUG_INFO( "Start to locate liver centroid" );

		radius.Fill(2);

		liverMask = km::binaryThresholdImage<FloatImageType, UCharImageType>(probablityMap, 0.3, 1.0, 1, 0);
		liverMask = km::binaryOpen<UCharImageType>( liverMask, radius );
		//liverMask = km::extractMaxConnectedComponent<UCharImageType>( liverMask );

		//liverMask = km::fillSliceHole<UCharImageType>(liverMask, 0, 0);
		
		FloatImageType::Pointer maskedProbablityMap = km::maskImage<FloatImageType, UCharImageType>(probablityMap, liverMask);
		km::locateLiverFromProbablityMap<FloatImageType>( maskedProbablityMap, liverCentroid);

		KM_DEBUG_PRINT( "Liver centroid ", liverCentroid );

		UCharImageType::SizeType padding;
		padding.Fill( 20 );
		liverMask = km::padImage<UCharImageType>(liverMask,padding, padding, 0);

		if(WRITE_MIDDLE_RESULT)
		{
			km::writeImage<UCharImageType>( outputdir, "liverMask.nii.gz", liverMask );
		}

		memorymeter.Stop( "adboost segmentation" );
		chronometer.Stop( "adboost segmentation" );
	}
	
	//����SSM
	KM_DEBUG_PRINT( "Loading statistical shape model..", SSMFile );
	StatisticalModelType::Pointer model = StatisticalModelType::New();
	model->Load( SSMFile );
	
	StatisticalModelType::VectorType paramweights = model->GetPCAVarianceVector();
	StatisticalModelType::VectorType variances = model->GetPCAVarianceVector();
	double totalVariances = 0.0;
	for (int i=0;i<model->GetNumberOfPrincipalComponents();i++)
	{
		totalVariances += std::abs(variances[i]);
	}
	
	//Since the first variance is alway the biggest one, so the weight will be based on this value.
	//Some times we might get a crazy shape model which generate zero variance(should never happen), so for preventing from crashing, add a small value.
	paramweights /= (variances[0]+0.000001);

	unsigned int numberOfMainShapeComponents = 1;
	unsigned int numberOfMainProfileComponents = 1;
	double totalMainVariances = 0.0;
	double variancesThresholdLow = 0.90;
	double variancesThresholdHigh = 0.98;
	for (int i=0;i<model->GetNumberOfPrincipalComponents();i++)
	{
		totalMainVariances += std::abs(variances[i]);
		double variancePercentage = totalMainVariances/totalVariances;
		if ( variancePercentage < variancesThresholdLow )
		{
			numberOfMainShapeComponents ++;
			numberOfMainProfileComponents ++;
		}
		else if( variancePercentage >= variancesThresholdLow && variancePercentage < variancesThresholdHigh )
		{
			numberOfMainProfileComponents ++;
		}
		else
		{
			break;
		}
	}
	KM_DEBUG_PRINT( "numberOfMainComponents for shape", numberOfMainShapeComponents );
	KM_DEBUG_PRINT( "numberOfMainComponents for profile", numberOfMainProfileComponents );

	//FIXEME
	//const char* geofile = "geoImage.mha";
	GeometryImageType::Pointer geoImage = km::readImage<GeometryImageType>( geoFile );
	MeshType::Pointer deformedMesh = model->DrawMean();

	unsigned int numberOfPointsOnMesh = deformedMesh->GetNumberOfPoints();
	deformedMesh->GetPointData()->Reserve( numberOfPointsOnMesh );
	MeshType::Pointer referenceShapeMesh = model->GetRepresenter()->GetReference();
	MeshType::Pointer meanShapeMesh = model->DrawMean();

	loadSimplexMeshGeometryData<MeshType>(geoImage, deformedMesh);

	//FIXED
	//MeshType::Pointer errorMap = km::readMesh<MeshType>("D:\\Workspace\\ASM\\projects\\LiverSegbyASM\\experiments\\training\\appearance\\output_20140815\\ErrorMap_LIVER.vtk");
	//km::copyMeshData<MeshType, MeshType>( errorMap, referenceShapeMesh );
	loadSimplexMeshGeometryData<MeshType>(geoImage, referenceShapeMesh);

	typedef itk::StatisticalShapeModelTransform<RepresenterType, double, Dimension> ShapeTransformType;
	ShapeTransformType::Pointer shapeTransform = ShapeTransformType::New();
	shapeTransform->SetStatisticalModel( model );
	shapeTransform->SetUsedNumberOfCoefficients( numberOfMainShapeComponents );
	shapeTransform->SetIdentity();

	/************************************************************************/
	/* Locate model                                                         */
	/************************************************************************/
	//////////////////////////////////////////////////////////////////////////
	MeshType::PointType meshCentroid = km::getMeshCentroid<MeshType>( deformedMesh );

	//typedef itk::Euler3DTransform<double> RigidTransformType;
	//typedef itk::VersorRigid3DTransform<double> RigidTransformType;
	//typedef itk::ScaleSkewVersor3DTransform<double> RigidTransformType;
	//typedef itk::ScaleVersor3DTransform<double> RigidTransformType;
	//typedef itk::AffineTransform<double, Dimension> RigidTransformType;
	//typedef itk::CenteredAffineTransform<double, Dimension> RigidTransformType;
	//typedef itk::Rigid3DTransform< double > RigidTransformType;
	typedef itk::Similarity3DTransform<double> RigidTransformType;

	RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
	rigidTransform->SetIdentity();
	rigidTransform->SetCenter( meshCentroid );

	typedef itk::CompositeTransform<double, Dimension> CompositeTransformType;
	CompositeTransformType::Pointer compositeTransform = CompositeTransformType::New();
	compositeTransform->AddTransform( rigidTransform );
	compositeTransform->AddTransform( shapeTransform );

	rigidTransform->Translate( liverCentroid - meshCentroid );

	km::transformMesh<MeshType, RigidTransformType>( deformedMesh, deformedMesh, rigidTransform );
	
	if(WRITE_MIDDLE_RESULT)
	{
		km::writeMesh<MeshType>( outputdir, "locatedMesh.vtk", deformedMesh );
	}
	
	//Nofify locate result
	notifier->notifyMesh( deformedMesh );

	memorymeter.Start( "fitting" );
	chronometer.Start( "fitting" );

	bool flag_fittingByDistMap = true;
	bool flag_fittingByLiverProfile = false;
	bool flag_deformingByLiverProfile = true;
	bool flag_deformingBoundaryProfile = true;

	std::vector<double> opt_scales;
	opt_scales.resize( compositeTransform->GetNumberOfParameters() );

	int p=0;
	for (int t=0;t<shapeTransform->GetNumberOfParameters();t++)
	{
		opt_scales[p++] = 1.0;
	}
	for (int t=0;t<3;t++)
	{
		opt_scales[p++] = 1.0;
	}
	for (int t=0;t<3;t++)
	{
		opt_scales[p++] = 1e-3;
	}
	for (int t=0;t<1;t++)
	{
		opt_scales[p++] = 1.0;
	}

	if(flag_fittingByDistMap)
	{
		std::cout<<"*******************************************************************************"<<std::endl;
		std::cout<<"**                                                                           **"<<std::endl;
		std::cout<<"**                                                                           **"<<std::endl;
		std::cout<<"**            Shape      Fitting         By           Distance Map           **"<<std::endl;
		std::cout<<"**                                                                           **"<<std::endl;
		std::cout<<"**                                                                           **"<<std::endl;
		std::cout<<"*******************************************************************************"<<std::endl;

		g_phase = MODEL_FITTING_BY_DISTMAP;

		FloatImageType::Pointer liverDistMap = km::calculateDistanceMap<UCharImageType, FloatImageType>(liverMask);
		//liverDistMap = km::powImage<FloatImageType>(liverDistMap, 2.0);
		liverDistMap = km::absImage<FloatImageType>(liverDistMap);
		liverDistMap = km::shiftScale<FloatImageType>(liverDistMap, 0, 100);

		if(WRITE_MIDDLE_RESULT)
		{
			km::writeImage<FloatImageType>( outputdir, "liverDistMap.nii.gz", liverDistMap );
		}

		int tmp_numberOfMainShapeComponents = numberOfMainShapeComponents;
		shapeTransform->SetUsedNumberOfCoefficients( tmp_numberOfMainShapeComponents );

		double ssmParamDiff = km::transformFittingToDistanceMap<FloatImageType, MeshType, CompositeTransformType>(
			liverDistMap,
			referenceShapeMesh,
			compositeTransform,
			tmp_numberOfMainShapeComponents,
			opt_scales,
			1500);

		std::cout<<compositeTransform->GetParameters()<<std::endl;

		km::transformMesh<MeshType, CompositeTransformType>( referenceShapeMesh, deformedMesh, compositeTransform );

		MeshType::Pointer rigidTransformedMeanMesh = km::transformMesh<MeshType, RigidTransformType>( meanShapeMesh, rigidTransform );

		if(WRITE_MIDDLE_RESULT)
		{
			km::writeMesh<MeshType>( outputdir, "rigidTransformedMeanMesh.vtk", rigidTransformedMeanMesh );
			km::writeMesh<MeshType>( outputdir, "shapeFittingMeshByDistMap", tmp_numberOfMainShapeComponents, ".vtk", deformedMesh );
		}

		//Nofify rough fitting result
		notifier->notifyMesh( deformedMesh );
	}
	
	ComputeGeometry<MeshType>( deformedMesh );

	if (true)
	{
		std::cout<<"*******************************************************************************"<<std::endl;
		std::cout<<"**                                                                           **"<<std::endl;
		std::cout<<"**                                                                           **"<<std::endl;
		std::cout<<"**            Tumor detection                                                **"<<std::endl;
		std::cout<<"**                                                                           **"<<std::endl;
		std::cout<<"**                                                                           **"<<std::endl;
		std::cout<<"*******************************************************************************"<<std::endl;

		g_phase = TUMOR_DETECTION;
		
		km::detectLiverIntensityRangeIncludingTumor<ShortImageType, MeshType>( inputImage, deformedMesh, g_liverThresholds );
		//km::detectLiverIntensityRangeWithoutTumor<ShortImageType, MeshType>( inputImage, deformedMesh, g_liverThresholds );
		KM_DEBUG_INFO( " Liver gray value range: " );
		for(int i=0;i<g_liverThresholds.size();i++)
		{
			std::cout<<g_liverThresholds[i].first<<","<<g_liverThresholds[i].second<<std::endl;
		}
		std::cout<<std::endl;
	}

	MeshType::Pointer errorMapLiver = km::readMesh<MeshType>("D:\\Workspace\\ASM\\projects\\LiverSegbyASM\\experiments\\training\\appearance\\output_20140815\\ErrorMap_LIVER.vtk");
	MeshType::Pointer varianceMap = km::readMesh<MeshType>("D:\\Workspace\\ASM\\projects\\LiverSegbyASM\\experiments\\buidModel\\shape\\meanMeshWithVariance.vtk");

	//TEST
	typedef itk::LinearInterpolateImageFunction<ShortImageType> IntensityInterpolatorType;
	typedef itk::LinearInterpolateImageFunction<GradientImageType> GradientInterpolatorType;
	IntensityInterpolatorType::Pointer intensityInterpolator = IntensityInterpolatorType::New();
	intensityInterpolator->SetInputImage( inputImage );
	GradientInterpolatorType::Pointer gradientInterpolator = GradientInterpolatorType::New();
	gradientInterpolator->SetInputImage( gradImage );

	if (flag_deformingByLiverProfile)
	{
		std::cout<<"*******************************************************************************"<<std::endl;
		std::cout<<"**                                                                           **"<<std::endl;
		std::cout<<"**                                                                           **"<<std::endl;
		std::cout<<"**            Deforming By Liver Profile                                     **"<<std::endl;
		std::cout<<"**                                                                           **"<<std::endl;
		std::cout<<"**                                                                           **"<<std::endl;
		std::cout<<"*******************************************************************************"<<std::endl;

		g_phase = DEFORMATION_BY_LIVER_PROFILE;

		shapeTransform->SetUsedNumberOfCoefficients( numberOfMainProfileComponents );

		typedef itk::DeformableSimplexMesh3DAdaboostClassifierForceFilter<MeshType, MeshType> DeformableMeshFilterType;
		DeformableMeshFilterType::Pointer deformableFilter = DeformableMeshFilterType::New();
		deformableFilter->SetGradient( gradImage );
		deformableFilter->SetIntensity( inputImage );
		deformableFilter->SetConfidenceThreshold( 1.3 );
		deformableFilter->SetAlpha( 0.3 );
		deformableFilter->SetBeta( 0.1 );
		deformableFilter->SetKappa( 0.3 );
		deformableFilter->SetRigidity( 2 );
		deformableFilter->SetGamma( 0.15 );
		deformableFilter->SetP2p(false);
		deformableFilter->SetInput( deformedMesh );
		deformableFilter->SetVarianceMap(varianceMap);
		deformableFilter->SetIterations( 50 );

		unsigned int maxIterations = 5;
		double rigidParaDiffTollerance = 0.002;
		double shapeParaDiffTollerance = 0.02;
		unsigned int iter_ellapsed = 0;

		bool fittingByComputing = false;
		bool fittingByOptimizing = false;
		bool flagLooping = true;

		MeshType::Pointer shapeMesh = km::transformMesh<MeshType, CompositeTransformType>( referenceShapeMesh, compositeTransform );
		km::loadSimplexMeshGeometryData<MeshType>(geoImage, shapeMesh);
		km::ComputeGeometry<MeshType>(shapeMesh, false);

		typedef itk::IdentityTransform<double, Dimension> IdentityTransformType;
		IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();
		MeshType::Pointer bestMesh = km::transformMesh<MeshType, IdentityTransformType>( deformedMesh, identityTransform );
		km::loadSimplexMeshGeometryData<MeshType>(geoImage, bestMesh);
		km::assigneMesh<MeshType>(bestMesh, 0.0);
		
		while ( flagLooping )
		{
			std::cout<<"********************************iteration: "<<iter_ellapsed<<"************************************"<<std::endl;

			generateBestPointSet<FloatImageType, MeshType, GradientInterpolatorType, IntensityInterpolatorType>
				( &adaboostProfileClassifier_Liver, &adaboostProfileClassifier_Liver, gradientInterpolator, intensityInterpolator, bestMesh, deformedMesh, shapeMesh, varianceMap, errorMapLiver, 1.5, 8);
			
			deformableFilter->SetBestMesh( bestMesh );
			deformableFilter->SetShapeMesh( shapeMesh );
			deformableFilter->Update();

			if (fittingByOptimizing && !fittingByComputing)
			{
				KM_DEBUG_INFO("Fitting model by optimizating.");
				//Store rigid parameters before fitting.
				RigidTransformType::ParametersType rigidParamPre = rigidTransform->GetParameters();
				//Store shape parameters before fitting.
				ShapeTransformType::ParametersType shapeParamPre = shapeTransform->GetParameters();

				km::transformFitting<MeshType, CompositeTransformType>(
					deformedMesh,
					referenceShapeMesh,
					compositeTransform,
					shapeTransform->GetUsedNumberOfCoefficients(),
					KdTree,
					opt_scales,
					200
					);

				//Store rigid parameters after fitting.
				RigidTransformType::ParametersType rigidParamPost = rigidTransform->GetParameters();
				//Store shape parameters after fitting.
				ShapeTransformType::ParametersType shapeParamPost = shapeTransform->GetParameters();

				for (int ttt=0;ttt<shapeTransform->GetUsedNumberOfCoefficients();ttt++)
				{
					if (shapeParamPost[ttt]<-3.0){
						shapeParamPost[ttt] = -3.0;
					}else if (shapeParamPost[ttt]>3.0){
						shapeParamPost[ttt] = 3.0;
					}
				}

				//Calculate shape parameters variance.
				double rigidParamDiff = 0;
				for (int p=0;p<rigidTransform->GetNumberOfParameters();p++)
				{
					rigidParamDiff += std::abs(rigidParamPost[p]-rigidParamPre[p])*(opt_scales[p+shapeTransform->GetNumberOfParameters()]);
				}
				rigidParamDiff /= rigidTransform->GetNumberOfParameters();
				KM_DEBUG_PRINT("Rigid parameters difference", rigidParamDiff);
				if (rigidParamDiff < rigidParaDiffTollerance)
				{
					KM_DEBUG_INFO( "Rigid parameters difference is smaller than tollerance. Switch to fitting by computing now.." );
					fittingByComputing = true;
				}
			}
			else
			{
				RigidTransformType::Pointer inversedRigidTransform = RigidTransformType::New();
				rigidTransform->GetInverse( inversedRigidTransform );
				
				MeshType::Pointer displacedShapeMesh = shapeMesh;
				km::transformMesh<MeshType, RigidTransformType>( deformedMesh, displacedShapeMesh, inversedRigidTransform );

				KM_DEBUG_INFO("Fitting model by computing.");
				//Store shape parameters before fitting.
				ShapeTransformType::ParametersType shapeParamPre = shapeTransform->GetParameters();
				//Store shape parameters after fitting.
				ShapeTransformType::ParametersType shapeParamPost = shapeTransform->GetParameters();

				StatisticalModelType::VectorType coefficient = model->ComputeCoefficientsForDataset( displacedShapeMesh );
				for (int ttt=0;ttt<shapeTransform->GetUsedNumberOfCoefficients();ttt++)
				{
					if (coefficient[ttt]<-3.0){
						shapeParamPost[ttt] = -3.0;
					}else if (coefficient[ttt]>3.0){
						shapeParamPost[ttt] = 3.0;
					}else{
						shapeParamPost[ttt] = coefficient[ttt];
					}
				}
				shapeTransform->SetParameters( shapeParamPost );

				//Calculate shape parameters variance.
				double shapeParamDiff = 0;
				for (int p=0;p<shapeTransform->GetUsedNumberOfCoefficients();p++)
				{
					shapeParamDiff += std::abs(shapeParamPost[p]-shapeParamPre[p]);
				}
				shapeParamDiff /= shapeTransform->GetUsedNumberOfCoefficients();
				KM_DEBUG_PRINT("SSM parameters difference", shapeParamDiff);
				if (shapeParamDiff < shapeParaDiffTollerance && iter_ellapsed > 0)
				{
					KM_DEBUG_INFO( "SSM parameters difference is smaller than tollerance. Stop fitting now.." );
					flagLooping = false;
				}
				else if (iter_ellapsed >= maxIterations)
				{
					KM_DEBUG_INFO( "Exceed maximum iteration number. Stop fitting now.." );
					flagLooping = false;
				}
			}

			std::cout<<"Rigid transform parameters: "<<rigidTransform->GetParameters()<<std::endl;
			std::cout<<"Shape transform parameters: "<<shapeTransform->GetParameters()<<std::endl;

			km::transformMesh<MeshType, CompositeTransformType>( referenceShapeMesh, shapeMesh, compositeTransform );
			km::ComputeGeometry<MeshType>( shapeMesh, false );

			if(WRITE_MIDDLE_RESULT)
			{
				km::writeMesh<MeshType>( outputdir, "LIVER-shapeMesh-iter", iter_ellapsed, ".vtk", shapeMesh );
				km::writeMesh<MeshType>( outputdir, "LIVER-bestMesh-iter", iter_ellapsed, ".vtk", bestMesh );
				km::writeMesh<MeshType>( outputdir, "LIVER-deformedMesh-iter", iter_ellapsed, ".vtk", deformedMesh );
			}

			iter_ellapsed++;
		}
	}

	if(flag_deformingBoundaryProfile)
	{
		std::cout<<"*******************************************************************************"<<std::endl;
		std::cout<<"**                                                                           **"<<std::endl;
		std::cout<<"**                                                                           **"<<std::endl;
		std::cout<<"**            Deforming By Boundary Profile                                  **"<<std::endl;
		std::cout<<"**                                                                           **"<<std::endl;
		std::cout<<"**                                                                           **"<<std::endl;
		std::cout<<"*******************************************************************************"<<std::endl;

		g_phase = DEFORMATION_BY_BOUNDARY_PROFILE;

		shapeTransform->SetUsedNumberOfCoefficients( numberOfMainProfileComponents );

		typedef itk::DeformableSimplexMesh3DAdaboostClassifierForceFilter<MeshType, MeshType> DeformableMeshFilterType;
		DeformableMeshFilterType::Pointer deformableFilter = DeformableMeshFilterType::New();
		deformableFilter->SetGradient( gradImage );
		deformableFilter->SetIntensity( inputImage );
		deformableFilter->SetConfidenceThreshold( 1.3 );
		deformableFilter->SetAlpha( 0.3 );
		deformableFilter->SetBeta( 0.1 );
		deformableFilter->SetKappa( 0.3 );
		deformableFilter->SetRigidity( 1 );
		deformableFilter->SetGamma( 0.15 );
		deformableFilter->SetP2p(false);
		deformableFilter->SetInput( deformedMesh );
		deformableFilter->SetVarianceMap(varianceMap);
		deformableFilter->SetIterations( 50 );

		unsigned int maxIterations = 0;
		double rigidParaDiffTollerance = 0.002;
		double shapeParaDiffTollerance = 0.02;
		unsigned int iter_ellapsed = 0;

		bool fittingByComputing = false;
		bool fittingByOptimizing = false;
		bool flagLooping = true;

		MeshType::Pointer shapeMesh = km::transformMesh<MeshType, CompositeTransformType>( referenceShapeMesh, compositeTransform );
		km::loadSimplexMeshGeometryData<MeshType>(geoImage, shapeMesh);
		km::ComputeGeometry<MeshType>(shapeMesh, false);

		typedef itk::IdentityTransform<double, Dimension> IdentityTransformType;
		IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();
		MeshType::Pointer bestMesh = km::transformMesh<MeshType, IdentityTransformType>( deformedMesh, identityTransform );
		km::loadSimplexMeshGeometryData<MeshType>(geoImage, bestMesh);
		km::assigneMesh<MeshType>(bestMesh, 0.0);

		while ( flagLooping )
		{
			std::cout<<"********************************iteration: "<<iter_ellapsed<<"************************************"<<std::endl;

			generateBestPointSet<FloatImageType, MeshType, GradientInterpolatorType, IntensityInterpolatorType>
				( &adaboostProfileClassifier_Boundary, &adaboostProfileClassifier_Liver, gradientInterpolator, intensityInterpolator, bestMesh, deformedMesh, shapeMesh, varianceMap, errorMapLiver, 1.5, 8);

			deformableFilter->SetBestMesh( bestMesh );
			deformableFilter->SetShapeMesh( shapeMesh );
			deformableFilter->Update();

			if (fittingByOptimizing && !fittingByComputing)
			{
				KM_DEBUG_INFO("Fitting model by optimizating.");
				//Store rigid parameters before fitting.
				RigidTransformType::ParametersType rigidParamPre = rigidTransform->GetParameters();
				//Store shape parameters before fitting.
				ShapeTransformType::ParametersType shapeParamPre = shapeTransform->GetParameters();

				km::transformFitting<MeshType, CompositeTransformType>(
					deformedMesh,
					referenceShapeMesh,
					compositeTransform,
					shapeTransform->GetUsedNumberOfCoefficients(),
					KdTree,
					opt_scales,
					200
					);

				//Store rigid parameters after fitting.
				RigidTransformType::ParametersType rigidParamPost = rigidTransform->GetParameters();
				//Store shape parameters after fitting.
				ShapeTransformType::ParametersType shapeParamPost = shapeTransform->GetParameters();

				for (int ttt=0;ttt<shapeTransform->GetUsedNumberOfCoefficients();ttt++)
				{
					if (shapeParamPost[ttt]<-3.0){
						shapeParamPost[ttt] = -3.0;
					}else if (shapeParamPost[ttt]>3.0){
						shapeParamPost[ttt] = 3.0;
					}
				}

				//Calculate shape parameters variance.
				double rigidParamDiff = 0;
				for (int p=0;p<rigidTransform->GetNumberOfParameters();p++)
				{
					rigidParamDiff += std::abs(rigidParamPost[p]-rigidParamPre[p])*(opt_scales[p+shapeTransform->GetNumberOfParameters()]);
				}
				rigidParamDiff /= rigidTransform->GetNumberOfParameters();
				KM_DEBUG_PRINT("Rigid parameters difference", rigidParamDiff);
				if (rigidParamDiff < rigidParaDiffTollerance)
				{
					KM_DEBUG_INFO( "Rigid parameters difference is smaller than tollerance. Switch to fitting by computing now.." );
					fittingByComputing = true;
				}
			}
			else
			{
				RigidTransformType::Pointer inversedRigidTransform = RigidTransformType::New();
				rigidTransform->GetInverse( inversedRigidTransform );
				MeshType::Pointer displacedShapeMesh = km::transformMesh<MeshType, RigidTransformType>( deformedMesh, inversedRigidTransform );

				KM_DEBUG_INFO("Fitting model by computing.");
				//Store shape parameters before fitting.
				ShapeTransformType::ParametersType shapeParamPre = shapeTransform->GetParameters();
				//Store shape parameters after fitting.
				ShapeTransformType::ParametersType shapeParamPost = shapeTransform->GetParameters();

				StatisticalModelType::VectorType coefficient = model->ComputeCoefficientsForDataset( displacedShapeMesh );
				for (int ttt=0;ttt<shapeTransform->GetUsedNumberOfCoefficients();ttt++)
				{
					if (coefficient[ttt]<-3.0){
						shapeParamPost[ttt] = -3.0;
					}else if (coefficient[ttt]>3.0){
						shapeParamPost[ttt] = 3.0;
					}else{
						shapeParamPost[ttt] = coefficient[ttt];
					}
				}
				shapeTransform->SetParameters( shapeParamPost );

				//Calculate shape parameters variance.
				double shapeParamDiff = 0;
				for (int p=0;p<shapeTransform->GetUsedNumberOfCoefficients();p++)
				{
					shapeParamDiff += std::abs(shapeParamPost[p]-shapeParamPre[p]);
				}
				shapeParamDiff /= shapeTransform->GetUsedNumberOfCoefficients();
				KM_DEBUG_PRINT("SSM parameters difference", shapeParamDiff);
				if (shapeParamDiff < shapeParaDiffTollerance && iter_ellapsed > 0)
				{
					KM_DEBUG_INFO( "SSM parameters difference is smaller than tollerance. Stop fitting now.." );
					flagLooping = false;
				}
				else if (iter_ellapsed >= maxIterations)
				{
					KM_DEBUG_INFO( "Exceed maximum iteration number. Stop fitting now.." );
					flagLooping = false;
				}
			}

			std::cout<<"Rigid transform parameters: "<<rigidTransform->GetParameters()<<std::endl;
			std::cout<<"Shape transform parameters: "<<shapeTransform->GetParameters()<<std::endl;

			km::transformMesh<MeshType, CompositeTransformType>( referenceShapeMesh, shapeMesh, compositeTransform );
			km::ComputeGeometry<MeshType>( shapeMesh, false );

			if(WRITE_MIDDLE_RESULT)
			{
				km::writeMesh<MeshType>( outputdir, "BOUNDARY-shapeMesh-iter", iter_ellapsed, ".vtk", shapeMesh );
				km::writeMesh<MeshType>( outputdir, "BOUNDARY-bestMesh-iter", iter_ellapsed, ".vtk", bestMesh );
				km::writeMesh<MeshType>( outputdir, "BOUNDARY-deformedMesh-iter", iter_ellapsed, ".vtk", deformedMesh );
			}

			iter_ellapsed++;
		}

		//Nofify profile based fitting result
		notifier->notifyMesh( deformedMesh );
	}

	memorymeter.Stop( "fitting" );
	chronometer.Stop( "fitting" );

	KM_DEBUG_INFO( "Start to genereate final segmentation image..." );

	km::loadSimplexMeshGeometryData<MeshType>(geoImage, deformedMesh);
	MeshType::Pointer triangleLiverMesh = km::simplexMeshToTriangleMesh<MeshType, MeshType>(deformedMesh);

	if(WRITE_MIDDLE_RESULT)
	{
		km::writeMesh<MeshType>( outputdir, "triangleLiverMesh.vtk", triangleLiverMesh );
	}

	UCharImageType::Pointer liversegmask = km::generateBinaryFromMesh<MeshType, UCharImageType, ShortImageType>(triangleLiverMesh, inputCached);

	KM_DEBUG_INFO("Fill slice hole...");
	km::fillSliceHole<UCharImageType>(liversegmask);

	notifier->notifyImage( liversegmask );

	//liversegmask->Print(std::cout);

	if (WRITE_MIDDLE_RESULT)
	{
		km::writeImage<UCharImageType>( outputdir, "finalResult.nii.gz", liversegmask );
	}

	memorymeter.Stop( "segmentation" );
	chronometer.Stop( "segmentation" );
	chronometer.Report( std::cout );
	memorymeter.Report( std::cout );
}