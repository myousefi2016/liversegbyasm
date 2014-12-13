#include "LiverSegmentationAPI.h"

#include "kmProfileClassifier.h"
#include "kmClassifierUtils.h"
#include "kmSSMUtils.h"
#include "kmModelFitting.h"
#include "kmProfileExtractor.h"
#include "kmGlobal.h"
#include "AdaSegmentAPI.h"
#include "kmUtility.h"
#include "kmVtkItkUtility.h"
#include "kmProcessing.h"

namespace km
{
	bool WRITE_MIDDLE_RESULT = true;

	void LiverSeg( 
		UCharImageType * segmentationResult,
		NotifierType * notifier,
		const char* outputdir,
		const char* inputImageFile,
		const char* SSMFile,
		const char* boundaryClassifierFile,
		const char* liverClassifierFile,
		const char* adaboostSegmentFile,
		const char* geoFile,
		const char* atlasImageFile,
		const char* configFile)
	{
		if ( notifier == NULL ){
			std::cout<<"Notifier is NULL!"<<std::endl;
			return;
		}
		if ( outputdir == NULL ){
			std::cout<<"outputdir is NULL!"<<std::endl;
			return;
		}
		if ( inputImageFile == NULL ){
			std::cout<<"inputImageFile is NULL!"<<std::endl;
			return;
		}
		if ( SSMFile == NULL ){
			std::cout<<"SSMFile is NULL!"<<std::endl;
			return;
		}
		if ( boundaryClassifierFile == NULL ){
			std::cout<<"boundaryClassifierFile is NULL!"<<std::endl;
			return;
		}
		if ( liverClassifierFile == NULL ){
			std::cout<<"liverClassifierFile is NULL!"<<std::endl;
			return;
		}
		if ( adaboostSegmentFile == NULL ){
			std::cout<<"adaboostSegmentFile is NULL!"<<std::endl;
			return;
		}
		if ( geoFile == NULL ){
			std::cout<<"geoFile is NULL!"<<std::endl;
			return;
		}
		if ( atlasImageFile == NULL ){
			std::cout<<"atlasImageFile is NULL!"<<std::endl;
			return;
		}
		if ( configFile == NULL ){
			std::cout<<"configFile is NULL!"<<std::endl;
			return;
		}

		KM_DEBUG_INFO("Load config file...");
		km::loadConfig(configFile);

		itk::TimeProbesCollectorBase chronometer;
		itk::MemoryProbesCollectorBase memorymeter;
		memorymeter.Start( "segmentation" );
		chronometer.Start( "segmentation" );

		g_phase = INITIALIZATION;

		memorymeter.Start( "initialization" );
		chronometer.Start( "initialization" );

		//读入分类器

		typedef km::ProfileClassifier ProfileClassifierType;
		ProfileClassifierType ProfileClassifier_Boundary;
		ProfileClassifier_Boundary.load(boundaryClassifierFile);

		ProfileClassifierType ProfileClassifier_Liver;
		ProfileClassifier_Liver.load(liverClassifierFile);

		//读入待分割图像
		KM_DEBUG_INFO("Read input image..");
		KM_DEBUG_INFO(inputImageFile);
		ShortImageType::Pointer inputImage = km::readImage<ShortImageType>( inputImageFile );
		km::shiftMinimum<ShortImageType>( inputImage, -1024 );

		KM_DEBUG_INFO("Read atlas image..");
		ShortImageType::Pointer atlasImage = km::readImage<ShortImageType>( atlasImageFile );

		KM_DEBUG_INFO("Store the original image information for final segmentation mask generation.");
		ShortImageType::Pointer inputCached = ShortImageType::New();
		inputCached->CopyInformation( inputImage );

		double minval = -1024;
		double maxval = 1024;

		KM_DEBUG_INFO("Down-sample input image..");
		SpacingType downspac;
		downspac.Fill( RESAMPLE_SPACING );
		inputImage = km::resampleImage<ShortImageType>( inputImage, downspac, minval );

		if(WRITE_MIDDLE_RESULT){
			km::writeImage<ShortImageType>( outputdir, "inputImage-downsampled.nii.gz", inputImage );
		}

		double zdist = RESAMPLE_SPACING * inputImage->GetLargestPossibleRegion().GetSize()[2];
		if ( zdist > 400 ){
			KM_DEBUG_INFO("ROI locating by template matching..");
			inputImage = km::extractRoiByTemplateMatching<ShortImageType, ShortImageType>( inputImage, atlasImage );
		}

		if(WRITE_MIDDLE_RESULT){
			km::writeImage<ShortImageType>( outputdir, "inputImage.nii.gz", inputImage );
		}

		KM_DEBUG_INFO("Histogram matching..");
		inputImage = km::histogramMatch<ShortImageType>( inputImage, atlasImage );

		KM_DEBUG_INFO("Smooth input image..");
		inputImage = km::minMaxSmooth<ShortImageType>( inputImage, 3, 1.0, 1 );

		if(WRITE_MIDDLE_RESULT){
			km::writeImage<ShortImageType>( outputdir, "inputImageSmoothed.nii.gz", inputImage );
		}

		KM_DEBUG_INFO("Calculate gradient image..");
		GradientImageType::Pointer gradImage = km::calculateRecursiveGradientImage<ShortImageType, GradientImageType>( inputImage, SIGMA );

		KM_DEBUG_INFO("Calculate minimum & maximum grey value..");
		km::calculateMinAndMax<ShortImageType>( inputImage, minval, maxval );
		KM_DEBUG_PRINT( "Min value", minval );
		KM_DEBUG_PRINT( "Max value", maxval );

		memorymeter.Stop( "initialization" );
		chronometer.Stop( "initialization" );

		FloatImageType::PointType liverCentroid;

		FloatImageType::Pointer probablityMap = FloatImageType::New();
		UCharImageType::Pointer liverMask  = UCharImageType::New();
		UCharImageType::SizeType radius;

		KM_DEBUG_INFO("Liver region initialization by region segmentation..");
		{
			g_phase = ADABOOST_REGION_SEGMENTATION;
			memorymeter.Start( "adboost segmentation" );
			chronometer.Start( "adboost segmentation" );

			KM_DEBUG_INFO("Start to generate liver probabilistic map by Adboost...");
			downspac.Fill( 4.0 );
			ShortImageType::Pointer downsampleInput = km::resampleImage<ShortImageType>( inputImage,downspac, minval );

			if(WRITE_MIDDLE_RESULT){
				km::writeImage<ShortImageType>( outputdir, "downsampleInput.nii.gz", downsampleInput );
			}

			//Get abdominal mask
			ShortImageType::Pointer bodyImage = ShortImageType::New();
			UCharImageType::Pointer bodyMask = UCharImageType::New();

			KM_DEBUG_INFO("Generate body mask image...");
			km::removeAir<ShortImageType, UCharImageType>( downsampleInput, bodyMask, minval, true );
			if(WRITE_MIDDLE_RESULT){
				km::writeImage<UCharImageType>( outputdir, "bodyMask.nii.gz", bodyMask );
			}

			//Locate liver
			KM_DEBUG_INFO( "Generate liver probilistic map.." );
			probablityMap = AdaSegment::adaSegment<FloatImageType, ShortImageType, UCharImageType>( downsampleInput, bodyMask, adaboostSegmentFile, 1, 300 );
			if(WRITE_MIDDLE_RESULT){
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

			if(WRITE_MIDDLE_RESULT){
				km::writeImage<UCharImageType>( outputdir, "liverMask.nii.gz", liverMask );
			}
			memorymeter.Stop( "adboost segmentation" );
			chronometer.Stop( "adboost segmentation" );
		}

		//载入SSM
		KM_DEBUG_PRINT( "Loading statistical shape model..", SSMFile );
		StatisticalModelType::Pointer model = StatisticalModelType::New();
		model->Load( SSMFile );

		StatisticalModelType::VectorType paramweights = model->GetPCAVarianceVector();
		StatisticalModelType::VectorType variances = model->GetPCAVarianceVector();
		double totalVariances = 0.0;
		for (int i=0;i<model->GetNumberOfPrincipalComponents();i++){
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
			if ( variancePercentage < variancesThresholdLow ){
				numberOfMainShapeComponents ++;
				numberOfMainProfileComponents ++;
			}else if( variancePercentage >= variancesThresholdLow && variancePercentage < variancesThresholdHigh ){
				numberOfMainProfileComponents ++;
			}else{
				break;
			}
		}
		KM_DEBUG_PRINT( "numberOfMainComponents for shape", numberOfMainShapeComponents );
		KM_DEBUG_PRINT( "numberOfMainComponents for profile", numberOfMainProfileComponents );

		//FIXEME
		//const char* geofile = "geoImage.mha";
		GeometryImageType::Pointer geoImage = km::readImage<GeometryImageType>( geoFile );

		MeshType::Pointer deformedMesh = model->DrawMean();
		MeshType::Pointer meanShapeMesh = model->DrawMean();
		MeshType::Pointer referenceShapeMesh = model->GetRepresenter()->GetReference();

		loadSimplexMeshGeometryData<MeshType>(geoImage, deformedMesh);
		loadSimplexMeshGeometryData<MeshType>(geoImage, referenceShapeMesh);
		km::ComputeGeometry<MeshType>( deformedMesh );

		km::initModelVariance<StatisticalModelType, MeshType>(model, deformedMesh, numberOfMainShapeComponents);
		if (WRITE_MIDDLE_RESULT){
			MeshType::Pointer varianceMesh = model->DrawMean();
			varianceMesh->GetPointData()->Reserve( varianceMesh->GetNumberOfPoints() );
			for (int i=0;i<varianceMesh->GetNumberOfPoints();i++){
				varianceMesh->SetPointData(i, g_varianceMap[i]);
			}
			km::writeMesh<MeshType>(outputdir,"varianceMesh.vtk", varianceMesh);
		}

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

		typedef itk::AffineTransform<double, Dimension> RigidTransformType;
		//typedef itk::Similarity3DTransform<double> RigidTransformType;
		RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
		rigidTransform->SetIdentity();
		rigidTransform->SetCenter( meshCentroid );
		rigidTransform->Translate( liverCentroid - meshCentroid );

		typedef itk::CompositeTransform<double, Dimension> CompositeTransformType;
		CompositeTransformType::Pointer compositeTransform = CompositeTransformType::New();
		compositeTransform->AddTransform( rigidTransform );
		compositeTransform->AddTransform( shapeTransform );

		km::transformMesh<MeshType, RigidTransformType>( deformedMesh, deformedMesh, rigidTransform );

		//Nofify locate result
		notifier->notify( deformedMesh );
		if (WRITE_MIDDLE_RESULT){
			writeMesh<MeshType>(outputdir, "locatedMesh.vtk", deformedMesh);
		}

		memorymeter.Start( "fitting" );
		chronometer.Start( "fitting" );

		bool flag_fittingByDistMap = true;
		bool flag_fittingByLiverProfile = false;
		bool flag_deformingByLiverProfile = true;
		bool flag_deformingBoundaryProfile = true;

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

			if(WRITE_MIDDLE_RESULT){
				km::writeImage<FloatImageType>( outputdir, "liverDistMap.nii.gz", liverDistMap );
			}

			std::vector<double> opt_scales;
			opt_scales.resize( compositeTransform->GetNumberOfParameters() );

			int p=0;
			for (int t=0;t<shapeTransform->GetNumberOfParameters();t++){
				opt_scales[p++] = 1.0;
			}
			for (int t=0;t<9;t++){
				opt_scales[p++] = 1.0;
			}
			for (int t=0;t<3;t++){
				opt_scales[p++] = 1e-3;
			}

			int tmp_numberOfMainShapeComponents = numberOfMainShapeComponents;
			shapeTransform->SetUsedNumberOfCoefficients( tmp_numberOfMainShapeComponents );

			typedef itk::Similarity3DTransform<double> Similarity3DTransformType;
			Similarity3DTransformType::Pointer similarityTransform = Similarity3DTransformType::New();
			similarityTransform->SetIdentity();
			similarityTransform->SetCenter( meshCentroid );
			similarityTransform->Translate( liverCentroid - meshCentroid );

			CompositeTransformType::Pointer compositeTransformTmp = CompositeTransformType::New();
			compositeTransformTmp->AddTransform( similarityTransform );
			compositeTransformTmp->AddTransform( shapeTransform );

			double ssmParamDiff = km::transformFittingToDistanceMap<FloatImageType, MeshType, CompositeTransformType>(
				liverDistMap,
				referenceShapeMesh,
				compositeTransformTmp,
				tmp_numberOfMainShapeComponents,
				opt_scales,
				300);

			ShapeTransformType::ParametersType shapeParamPost = shapeTransform->GetParameters();
			for (int ttt=0;ttt<shapeTransform->GetUsedNumberOfCoefficients();ttt++)
			{
				if (shapeParamPost[ttt]<-3.0){
					shapeParamPost[ttt] = 0;
				}else if (shapeParamPost[ttt]>3.0){
					shapeParamPost[ttt] = 0;
				}
			}
			shapeTransform->SetParameters( shapeParamPost );

			rigidTransform->SetMatrix(similarityTransform->GetMatrix());
			rigidTransform->SetOffset(similarityTransform->GetOffset());
			std::cout<<compositeTransform->GetParameters()<<std::endl;

			km::transformMesh<MeshType, CompositeTransformType>( referenceShapeMesh, deformedMesh, compositeTransform );

			//Nofify rough fitting result
			notifier->notify( deformedMesh );
			if (WRITE_MIDDLE_RESULT)
			{
				writeMesh<MeshType>(outputdir, "initializedMesh.vtk", deformedMesh);
			}
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

			//km::detectLiverIntensityRangeIncludingTumor<ShortImageType, MeshType>( inputImage, deformedMesh, g_liverThresholds );
			km::detectLiverIntensityRangeWithoutTumor<ShortImageType, MeshType>( inputImage, deformedMesh, g_liverThresholds );
			KM_DEBUG_INFO( " Liver gray value range: " );
			for(int i=0;i<g_liverThresholds.size();i++)
			{
				std::cout<<g_liverThresholds[i].first<<","<<g_liverThresholds[i].second<<std::endl;
			}
			std::cout<<std::endl;
		}

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

			unsigned int maxIterations = 20;
			double rigidParaDiffTollerance = 0.002;
			double shapeParaDiffTollerance = 0.03;
			unsigned int iter_ellapsed = 0;

			bool flagLooping = true;

			MeshType::Pointer compositeFittedMesh = km::transformMesh<MeshType, CompositeTransformType>( referenceShapeMesh, compositeTransform );
			km::loadSimplexMeshGeometryData<MeshType>(geoImage, compositeFittedMesh);
			km::ComputeGeometry<MeshType>(compositeFittedMesh, true);

			typedef itk::IdentityTransform<double, Dimension> IdentityTransformType;
			IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();
			MeshType::Pointer bestMesh = km::transformMesh<MeshType, IdentityTransformType>( deformedMesh, identityTransform );
			km::loadSimplexMeshGeometryData<MeshType>(geoImage, bestMesh);
			km::assigneMesh<MeshType>(bestMesh, 0.0);

			typedef km::ProfileExtractor<ShortImageType> ProfileExtractorType;
			ProfileExtractorType profileExtractor;
			profileExtractor.setImage(inputImage);
			profileExtractor.enableCache(false);

			typedef km::ClassifierUtils<MeshType, ProfileExtractorType> ClassifierUtilsType;
			ClassifierUtilsType classifierUtils;
			classifierUtils.SetBoundaryClassifier(&ProfileClassifier_Boundary);
			classifierUtils.SetRegionClassifier(&ProfileClassifier_Liver);
			classifierUtils.SetProfileExtractor(&profileExtractor);

			typedef km::SSMUtils<MeshType, StatisticalModelType, RigidTransformType, ShapeTransformType> SSMUtilsType;
			SSMUtilsType ssmUtils;
			ssmUtils.SetSSM(model);
			ssmUtils.SetRigidTransform(rigidTransform);
			ssmUtils.SetShapeTransform(shapeTransform);

			MeshType::Pointer tmpMesh = MeshType::New();

			while ( flagLooping )
			{
				std::cout<<"********************************iteration: "<<iter_ellapsed<<"************************************"<<std::endl;

				g_liverCentroid.CastFrom(km::getMeshCentroid<MeshType>(compositeFittedMesh));
				rigidTransform->SetCenter( g_liverCentroid );

				classifierUtils.deformByLiverClassification(bestMesh, compositeFittedMesh, 1.5, 20);

				//Store rigid parameters before fitting.
				RigidTransformType::ParametersType rigidParamPre = rigidTransform->GetParameters();
				//Store shape parameters before fitting.
				ShapeTransformType::ParametersType shapeParamPre = shapeTransform->GetParameters();

				km::transformMesh<MeshType, CompositeTransformType>( referenceShapeMesh, tmpMesh, compositeTransform );

				KM_DEBUG_PRINT("Composite transform fitting...", iter_ellapsed);
				ssmUtils.compositeTransformFitting(bestMesh);

				km::transformMesh<MeshType, CompositeTransformType>( referenceShapeMesh, compositeFittedMesh, compositeTransform );
				km::ComputeGeometry<MeshType>( compositeFittedMesh, true );

				classifierUtils.updateShapeNormals(tmpMesh, compositeFittedMesh);

				//Store rigid parameters before fitting.
				RigidTransformType::ParametersType rigidParamPost = rigidTransform->GetParameters();
				//Store shape parameters before fitting.
				ShapeTransformType::ParametersType shapeParamPost = shapeTransform->GetParameters();

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

				std::cout<<"Rigid transform parameters: "<<rigidTransform->GetParameters()<<std::endl;
				std::cout<<"Shape transform parameters: "<<shapeTransform->GetParameters()<<std::endl;

				notifier->notify( compositeFittedMesh );
				notifier->notify( bestMesh );

				if (WRITE_MIDDLE_RESULT)
				{
					km::writeMesh<MeshType>(outputdir, "bestMesh", iter_ellapsed, ".vtk", bestMesh);
					km::writeMesh<MeshType>(outputdir, "compositeFittedMesh", iter_ellapsed, ".vtk", compositeFittedMesh);
				}

				iter_ellapsed++;
			}
		}

		//if(flag_deformingBoundaryProfile)
		//{
		//	std::cout<<"*******************************************************************************"<<std::endl;
		//	std::cout<<"**                                                                           **"<<std::endl;
		//	std::cout<<"**                                                                           **"<<std::endl;
		//	std::cout<<"**            Deforming By Boundary Profile                                  **"<<std::endl;
		//	std::cout<<"**                                                                           **"<<std::endl;
		//	std::cout<<"**                                                                           **"<<std::endl;
		//	std::cout<<"*******************************************************************************"<<std::endl;

		//	g_phase = DEFORMATION_BY_BOUNDARY_PROFILE;

		//	shapeTransform->SetUsedNumberOfCoefficients( numberOfMainProfileComponents );

		//	typedef itk::DeformableSimplexMesh3DWithShapePriorFilter<MeshType, MeshType> DeformableMeshFilterType;
		//	DeformableMeshFilterType::Pointer deformableFilter = DeformableMeshFilterType::New();
		//	deformableFilter->SetGradient( gradImage );
		//	deformableFilter->SetIntensity( inputImage );
		//	deformableFilter->SetConfidenceThreshold( 1.3 );
		//	deformableFilter->SetAlpha( 0.3 );
		//	deformableFilter->SetBeta( 0.1 );
		//	deformableFilter->SetKappa( 0.3 );
		//	deformableFilter->SetRigidity( 1 );
		//	deformableFilter->SetGamma( 0.15 );
		//	deformableFilter->SetP2p(false);
		//	deformableFilter->SetInput( deformedMesh );
		//	deformableFilter->SetVarianceMap(varianceMap);
		//	deformableFilter->SetIterations( 50 );

		//	unsigned int maxIterations = 0;
		//	double rigidParaDiffTollerance = 0.002;
		//	double shapeParaDiffTollerance = 0.02;
		//	unsigned int iter_ellapsed = 0;

		//	bool fittingByComputing = false;
		//	bool fittingByOptimizing = false;
		//	bool flagLooping = true;

		//	MeshType::Pointer shapeMesh = km::transformMesh<MeshType, CompositeTransformType>( referenceShapeMesh, compositeTransform );
		//	km::loadSimplexMeshGeometryData<MeshType>(geoImage, shapeMesh);
		//	km::ComputeGeometry<MeshType>(shapeMesh, false);

		//	typedef itk::IdentityTransform<double, Dimension> IdentityTransformType;
		//	IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();
		//	MeshType::Pointer bestMesh = km::transformMesh<MeshType, IdentityTransformType>( deformedMesh, identityTransform );
		//	km::loadSimplexMeshGeometryData<MeshType>(geoImage, bestMesh);
		//	km::assigneMesh<MeshType>(bestMesh, 0.0);

		//	while ( flagLooping )
		//	{
		//		std::cout<<"********************************iteration: "<<iter_ellapsed<<"************************************"<<std::endl;

		//		generateBestPointSet<FloatImageType, MeshType, GradientInterpolatorType, IntensityInterpolatorType>
		//			( &ProfileClassifier_Boundary, gradientInterpolator, intensityInterpolator, bestMesh, deformedMesh, varianceMap, 1.5, 8);

		//		deformableFilter->SetTargetMesh( bestMesh );
		//		deformableFilter->SetShapeMesh( shapeMesh );
		//		deformableFilter->Update();

		//		if (fittingByOptimizing && !fittingByComputing)
		//		{
		//			KM_DEBUG_INFO("Fitting model by optimizating.");
		//			//Store rigid parameters before fitting.
		//			RigidTransformType::ParametersType rigidParamPre = rigidTransform->GetParameters();
		//			//Store shape parameters before fitting.
		//			ShapeTransformType::ParametersType shapeParamPre = shapeTransform->GetParameters();

		//			km::transformFitting<MeshType, CompositeTransformType>(
		//				deformedMesh,
		//				referenceShapeMesh,
		//				compositeTransform,
		//				shapeTransform->GetUsedNumberOfCoefficients(),
		//				KdTree,
		//				opt_scales,
		//				200
		//				);

		//			//Store rigid parameters after fitting.
		//			RigidTransformType::ParametersType rigidParamPost = rigidTransform->GetParameters();
		//			//Store shape parameters after fitting.
		//			ShapeTransformType::ParametersType shapeParamPost = shapeTransform->GetParameters();

		//			for (int ttt=0;ttt<shapeTransform->GetUsedNumberOfCoefficients();ttt++)
		//			{
		//				if (shapeParamPost[ttt]<-3.0){
		//					shapeParamPost[ttt] = -3.0;
		//				}else if (shapeParamPost[ttt]>3.0){
		//					shapeParamPost[ttt] = 3.0;
		//				}
		//			}

		//			//Calculate shape parameters variance.
		//			double rigidParamDiff = 0;
		//			for (int p=0;p<rigidTransform->GetNumberOfParameters();p++)
		//			{
		//				rigidParamDiff += std::abs(rigidParamPost[p]-rigidParamPre[p])*(opt_scales[p+shapeTransform->GetNumberOfParameters()]);
		//			}
		//			rigidParamDiff /= rigidTransform->GetNumberOfParameters();
		//			KM_DEBUG_PRINT("Rigid parameters difference", rigidParamDiff);
		//			if (rigidParamDiff < rigidParaDiffTollerance)
		//			{
		//				KM_DEBUG_INFO( "Rigid parameters difference is smaller than tollerance. Switch to fitting by computing now.." );
		//				fittingByComputing = true;
		//			}
		//		}
		//		else
		//		{
		//			RigidTransformType::Pointer inversedRigidTransform = RigidTransformType::New();
		//			rigidTransform->GetInverse( inversedRigidTransform );
		//			MeshType::Pointer displacedShapeMesh = km::transformMesh<MeshType, RigidTransformType>( deformedMesh, inversedRigidTransform );

		//			KM_DEBUG_INFO("Fitting model by computing.");
		//			//Store shape parameters before fitting.
		//			ShapeTransformType::ParametersType shapeParamPre = shapeTransform->GetParameters();
		//			//Store shape parameters after fitting.
		//			ShapeTransformType::ParametersType shapeParamPost = shapeTransform->GetParameters();

		//			StatisticalModelType::VectorType coefficient = model->ComputeCoefficientsForDataset( displacedShapeMesh );
		//			for (int ttt=0;ttt<shapeTransform->GetUsedNumberOfCoefficients();ttt++)
		//			{
		//				if (coefficient[ttt]<-3.0){
		//					shapeParamPost[ttt] = -3.0;
		//				}else if (coefficient[ttt]>3.0){
		//					shapeParamPost[ttt] = 3.0;
		//				}else{
		//					shapeParamPost[ttt] = coefficient[ttt];
		//				}
		//			}
		//			shapeTransform->SetParameters( shapeParamPost );

		//			//Calculate shape parameters variance.
		//			double shapeParamDiff = 0;
		//			for (int p=0;p<shapeTransform->GetUsedNumberOfCoefficients();p++)
		//			{
		//				shapeParamDiff += std::abs(shapeParamPost[p]-shapeParamPre[p]);
		//			}
		//			shapeParamDiff /= shapeTransform->GetUsedNumberOfCoefficients();
		//			KM_DEBUG_PRINT("SSM parameters difference", shapeParamDiff);
		//			if (shapeParamDiff < shapeParaDiffTollerance && iter_ellapsed > 0)
		//			{
		//				KM_DEBUG_INFO( "SSM parameters difference is smaller than tollerance. Stop fitting now.." );
		//				flagLooping = false;
		//			}
		//			else if (iter_ellapsed >= maxIterations)
		//			{
		//				KM_DEBUG_INFO( "Exceed maximum iteration number. Stop fitting now.." );
		//				flagLooping = false;
		//			}
		//		}

		//		std::cout<<"Rigid transform parameters: "<<rigidTransform->GetParameters()<<std::endl;
		//		std::cout<<"Shape transform parameters: "<<shapeTransform->GetParameters()<<std::endl;

		//		km::transformMesh<MeshType, CompositeTransformType>( referenceShapeMesh, shapeMesh, compositeTransform );
		//		km::ComputeGeometry<MeshType>( shapeMesh, false );

		//		if(WRITE_MIDDLE_RESULT)
		//		{
		//			km::writeMesh<MeshType>( outputdir, "BOUNDARY-shapeMesh-iter", iter_ellapsed, ".vtk", shapeMesh );
		//			km::writeMesh<MeshType>( outputdir, "BOUNDARY-bestMesh-iter", iter_ellapsed, ".vtk", bestMesh );
		//			km::writeMesh<MeshType>( outputdir, "BOUNDARY-deformedMesh-iter", iter_ellapsed, ".vtk", deformedMesh );
		//		}

		//		iter_ellapsed++;
		//	}

		//	//Nofify profile based fitting result
		//	notifier->notifyMesh( deformedMesh );
		//}

		//memorymeter.Stop( "fitting" );
		//chronometer.Stop( "fitting" );

		//KM_DEBUG_INFO( "Start to genereate final segmentation image..." );

		//km::loadSimplexMeshGeometryData<MeshType>(geoImage, deformedMesh);
		//MeshType::Pointer triangleLiverMesh = km::simplexMeshToTriangleMesh<MeshType, MeshType>(deformedMesh);

		//km::writeMesh<MeshType>( outputdir, "triangleLiverMesh.vtk", triangleLiverMesh );

		//UCharImageType::Pointer liversegmask = km::generateBinaryFromMesh<MeshType, UCharImageType, ShortImageType>(triangleLiverMesh, inputCached);

		//KM_DEBUG_INFO("Fill slice hole...");
		//km::fillSliceHole<UCharImageType>(liversegmask);

		//notifier->notifyImage( liversegmask );

		////liversegmask->Print(std::cout);

		//km::writeImage<UCharImageType>( outputdir, "finalResult.nii.gz", liversegmask );

		//segmentationResult = liversegmask;

		memorymeter.Stop( "segmentation" );
		chronometer.Stop( "segmentation" );
		chronometer.Report( std::cout );
		memorymeter.Report( std::cout );
	}
}