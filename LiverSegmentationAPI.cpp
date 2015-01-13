#include <fstream>						// file I/O
#include <map>
#include "itkImage.h"
#include "itkMesh.h"
#include "itkSimplexMesh.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"
#include "itkSimpleFilterWatcher.h"
#include "itkDefaultDynamicMeshTraits.h"
#include "Representers/ITK/itkSimplexMeshRepresenter.h"
#include "statismo_ITK/itkStatisticalModel.h"
#include "statismo_ITK/itkStatisticalShapeModelTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkCompositeTransform.h"
#include "itkAffineTransform.h"
#include "itkSimilarity3DTransform.h"
#include "itkSimpleFilterWatcher.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkListSample.h"
#include "itkKdTreeGenerator.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkDeformableSimplexMesh3DWithShapePriorFilter.h"

#include "kmCommon.h"
#include "AdaSegmentAPI.h"

#include "LiverSegmentationAPI.h"

namespace km
{
	const unsigned int Dimension = 3;
	typedef itk::Image<float, Dimension> ShortImageType;
	typedef itk::Image<unsigned char, Dimension> UCharImageType;
	typedef itk::Image<float, Dimension> FloatImageType;
	typedef itk::CovariantVector<double, Dimension> GradientVectorType;
	typedef itk::Image<GradientVectorType, Dimension> GradientImageType;
	typedef ShortImageType::SpacingType SpacingType;

	typedef double MeshPixelType;
	typedef itk::SimplexMeshRepresenter<MeshPixelType, Dimension> RepresenterType;
	typedef itk::StatisticalModel<RepresenterType>                StatisticalModelType;
	typedef RepresenterType::MeshType                             MeshType;
	typedef itk::Mesh<MeshPixelType, 3> TriangleMeshType;

	bool WRITE_MIDDLE_RESULT = true;

	void LiverSeg( 
		km::NotifierBase * notifier,
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
		km::Config::loadConfig(configFile);

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
		downspac.Fill( g_resample_spacing );
		inputImage = km::resampleImage<ShortImageType>( inputImage, downspac, minval );

		if(WRITE_MIDDLE_RESULT){
			km::writeImage<ShortImageType>( outputdir, "inputImage-downsampled.nii.gz", inputImage );
		}

		double zdist = g_resample_spacing * inputImage->GetLargestPossibleRegion().GetSize()[2];
		if ( zdist > 400 ){
			KM_DEBUG_INFO("ROI locating by template matching..");
			inputImage = km::extractRoiByTemplateMatching<ShortImageType, ShortImageType>( inputImage, atlasImage );
		}

		KM_DEBUG_INFO("Histogram matching..");
		inputImage = km::histogramMatch<ShortImageType>( inputImage, atlasImage );

		if(WRITE_MIDDLE_RESULT){
			km::writeImage<ShortImageType>( outputdir, "inputImage.nii.gz", inputImage );
		}

		KM_DEBUG_INFO("Smooth input image..");
		inputImage = km::minMaxSmooth<ShortImageType>( inputImage, 3, 1.0, 1 );

		if(WRITE_MIDDLE_RESULT){
			km::writeImage<ShortImageType>( outputdir, "inputImageSmoothed.nii.gz", inputImage );
		}

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
		if (notifier!=NULL)
		{
			notifier->notify();
		}
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

			km::transformFittingToDistanceMap<FloatImageType, MeshType, CompositeTransformType>(
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
			if (notifier!=NULL)
			{
				notifier->notify();
			}
			if (WRITE_MIDDLE_RESULT)
			{
				writeMesh<MeshType>(outputdir, "initializedMesh.vtk", deformedMesh);
			}
		}

		ComputeGeometry<MeshType>( deformedMesh );

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

			typedef itk::DeformableSimplexMesh3DWithShapePriorFilter<MeshType,
																	 MeshType,
																	 ShortImageType,
																	 StatisticalModelType,
																	 RigidTransformType,
																	 ShapeTransformType> DeformableFilterType;
			DeformableFilterType::Pointer deformaFilter = DeformableFilterType::New();
			deformaFilter->SetAlpha(km::g_alpha);
			deformaFilter->SetKappa(km::g_kappa);
			deformaFilter->SetBeta(km::g_beta);
			deformaFilter->SetGamma(km::g_gamma);
			deformaFilter->SetRigidity(km::g_rigidity);
			deformaFilter->SetIterations(100);
			deformaFilter->SetInput(deformedMesh);
			deformaFilter->SetInputImage(inputImage);
			deformaFilter->SetStatisticalModel(model);
			deformaFilter->SetRigidTransform(rigidTransform);
			deformaFilter->SetShapeTransform(shapeTransform);
			deformaFilter->SetBoundaryClassifier(&ProfileClassifier_Boundary);
			deformaFilter->SetLiverClassifier(&ProfileClassifier_Liver);
			//itk::SimpleFilterWatcher watcher(deformaFilter, "DeformableSimplexMesh3DWithShapePriorFilter");
			deformaFilter->Update();

			if (notifier!=NULL)
			{
				notifier->notify();
			}
			if (WRITE_MIDDLE_RESULT){
				km::writeMesh<MeshType>(outputdir, "deformedMesh.vtk", deformedMesh);
				km::writeMesh<MeshType>(outputdir, "shapeMeshFitted.vtk", deformaFilter->GetShapeMesh());
			}
		}

		memorymeter.Stop( "fitting" );
		chronometer.Stop( "fitting" );

		
		KM_DEBUG_INFO( "Start to genereate final segmentation image..." );

		km::loadSimplexMeshGeometryData<MeshType>(geoImage, deformedMesh);
		MeshType::Pointer triangleLiverMesh = km::simplexMeshToTriangleMesh<MeshType, MeshType>(deformedMesh);
		//km::writeMesh<MeshType>( outputdir, "triangleLiverMesh.vtk", triangleLiverMesh );

		UCharImageType::Pointer liversegmask = km::generateBinaryFromMesh<MeshType, UCharImageType, ShortImageType>(triangleLiverMesh, inputCached);
		KM_DEBUG_INFO("Fill slice hole...");
		//km::fillSliceHole<UCharImageType>(liversegmask);
		km::writeImage<UCharImageType>( outputdir, "finalResult.nii.gz", liversegmask );
		

		memorymeter.Stop( "segmentation" );
		chronometer.Stop( "segmentation" );
		chronometer.Report( std::cout );
		memorymeter.Report( std::cout );
	}
}