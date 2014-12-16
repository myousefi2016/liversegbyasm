#ifndef __kmRegistration_h
#define __kmRegistration_h

#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkImageRegistrationMethod.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include "itkPointSetToPointSetRegistrationMethod.h"
#include "itkEuclideanDistanceKdTreeMultipleValuePointSetMetric.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include "itkPointSetToImageRegistrationMethod.h"
#include "itkMeanSquaresPointSetToImageMetric.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkBSplineTransform.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkLBFGSBOptimizer.h"
#include <itkMeanSquaresHistogramImageToImageMetric.h>
#include "itkNormalizedMutualInformationHistogramImageToImageMetric.h"
#include <itkCorrelationCoefficientHistogramImageToImageMetric.h>
#include <itkKappaStatisticImageToImageMetric.h>

#include <itkScaleSkewVersor3DTransform.h>
#include <itkOnePlusOneEvolutionaryOptimizer.h>
#include <itkNormalVariateGenerator.h>
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkMultiResolutionPyramidImageFilter.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkMattesMutualInformationWithProbabilityImageToImageMetric.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkBSplineTransformInitializer.h"
#include "itkResampleImageFilter.h"
#include "itkBSplineResampleImageFunction.h"
#include "itkPointSetToPointSetRegistrationMethod.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include "itkEuclideanDistanceMultipleValuePointMetric.h"
#include "itkMinValuePointSetToImageMultipleValueMetric.h"
#include "itkPointSetToImageMultipleValueRegistrationMethod.h"

namespace km 
{
	template <typename RegistrationType>
	class RegistrationCommandResolutionUpdate : public itk::Command 
	{
	public:
		typedef  RegistrationCommandResolutionUpdate   Self;
		typedef  itk::Command                   Superclass;
		typedef  itk::SmartPointer<Self>        Pointer;
		itkNewMacro( Self );
	protected:
		RegistrationCommandResolutionUpdate() {};
	public:
		typedef   RegistrationType                           RegistrationType;
		typedef   RegistrationType *                         RegistrationPointer;
		typedef   itk::OnePlusOneEvolutionaryOptimizer       OptimizerType;
		typedef   OptimizerType *                            OptimizerPointer;
		void Execute(itk::Object * object, const itk::EventObject & event)
		{
			if( !(itk::IterationEvent().CheckEvent( &event )) )
			{
				return;
			}
			RegistrationPointer registration = dynamic_cast<RegistrationPointer>( object );
			std::cout << "-------------------------------------" << std::endl;
			std::cout << "MultiResolution Level : " << registration->GetCurrentLevel()  << std::endl;
			std::cout << std::endl;

			OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >( 
				registration->GetOptimizer() );

			if ( registration->GetCurrentLevel() > 0 )
			{
				if(registration->GetCurrentLevel() < registration->GetNumberOfLevels()-1)
				{
					optimizer->SetMaximumIteration( static_cast<int>(optimizer->GetMaximumIteration() * 0.5) );
				}
				else
				{
					optimizer->SetMaximumIteration( 0 );
				}
			}
		}
		void Execute(const itk::Object * , const itk::EventObject & )
		{ return; }
	};

	template<class OptimizerType>
	class OptimizerCommandIterationUpdate : public itk::Command
	{
	public:
		typedef  OptimizerCommandIterationUpdate   Self;
		typedef  itk::Command             Superclass;
		typedef  itk::SmartPointer<Self>  Pointer;
		itkNewMacro( Self );
		typedef const OptimizerType                     *OptimizerPointer;
		void Execute(itk::Object *caller, const itk::EventObject & event)
		{
			Execute( (const itk::Object *)caller, event);
		}
		void Execute(const itk::Object * object, const itk::EventObject & event)
		{
			OptimizerPointer optimizer =
				dynamic_cast< OptimizerPointer >( object );
			if( ! itk::IterationEvent().CheckEvent( &event ) )
			{
				return;
			}
			std::cout << optimizer->GetCurrentIteration() << "   ";
			std::cout << optimizer->GetValue() << "   ";
			//std::cout << optimizer->GetFrobeniusNorm() << "   ";
			std::cout << optimizer->GetCurrentPosition() << std::endl;
		}

	protected:
		OptimizerCommandIterationUpdate() {};
		~OptimizerCommandIterationUpdate() {};
	};

	template < typename ImageType, typename MovingImageType, typename TTransform >
	void
		histogramRigidRegistration( const typename ImageType* fixedImage, 
		const typename MovingImageType* movingImage,
		typename TTransform::Pointer &  transform)
	{
		itkStaticConstMacro(Dimension, unsigned int, ImageType::ImageDimension);
		typedef typename TTransform TransformType;
		typedef typename TransformType::Pointer TransformPointer;
		typedef itk::OnePlusOneEvolutionaryOptimizer OptimizerType;
		typedef typename OptimizerType::Pointer OptimizerPointer;
		typedef itk::MeanSquaresImageToImageMetric< ImageType, MovingImageType> MetricType;
		typedef typename MetricType::Pointer MetricPointer;
		typedef itk::MultiResolutionImageRegistrationMethod< ImageType, MovingImageType > RegistrationType;
		typedef typename RegistrationType::Pointer RegistrationPointer;
		typedef itk::MultiResolutionPyramidImageFilter<ImageType, MovingImageType > FixedImagePyramidType;
		typedef itk::MultiResolutionPyramidImageFilter<ImageType, MovingImageType > MovingImagePyramidType;
		typedef itk::LinearInterpolateImageFunction<MovingImageType> InterpolatorType;
		typedef typename InterpolatorType::Pointer InterpolatorPointer;
		try 
		{ 
			//TransformPointer	      transform     = TransformType::New();
			OptimizerPointer		    optimizer     = OptimizerType::New();
			RegistrationPointer		  registration  = RegistrationType::New();
			MetricPointer		        metric        = MetricType::New();
			InterpolatorPointer     iterpolator   = InterpolatorType::New();

			FixedImagePyramidType::Pointer fixedImagePyramid = FixedImagePyramidType::New();
			MovingImagePyramidType::Pointer movingImagePyramid = MovingImagePyramidType::New();

			registration->SetInterpolator(  iterpolator   );
			registration->SetOptimizer(     optimizer     );
			registration->SetTransform(     transform     );
			registration->SetMetric(        metric  );
			registration->SetFixedImage(    fixedImage   );
			registration->SetMovingImage(   movingImage   );
			registration->SetFixedImagePyramid( fixedImagePyramid );
			registration->SetMovingImagePyramid( movingImagePyramid );

			ImageType::RegionType fixedRegion   = fixedImage->GetLargestPossibleRegion();
			registration->SetFixedImageRegion( fixedRegion );

			typedef itk::Statistics::NormalVariateGenerator  GeneratorType;
			GeneratorType::Pointer generator = GeneratorType::New();
			generator->Initialize(12345);

			optimizer->SetNormalVariateGenerator( generator );
			optimizer->SetInitialRadius( 3.0 );
			optimizer->SetGrowthFactor( 1.1 );
			optimizer->SetEpsilon( 1e-5 );
			optimizer->SetMaximumIteration( 100 );
			//optimizer->MaximizeOn();

			int N = transform->GetNumberOfParameters();
			OptimizerType::ScalesType scales( N );
			scales.Fill( 10.0 );
			scales[2] = 1.0;
			optimizer->SetScales( scales );

			ImageType::PointType fixedCenter = km::getCenter<ImageType>( fixedImage );
			MovingImageType::PointType movingCenter = km::getCenter<MovingImageType>( movingImage );
			transform->Translate( movingCenter - fixedCenter );

			registration->SetInitialTransformParameters(  transform->GetParameters()  );
			registration->SetNumberOfLevels( 3 );

			RegistrationCommandResolutionUpdate<RegistrationType>::Pointer registrationObserver = 
				RegistrationCommandResolutionUpdate<RegistrationType>::New();
			registration->AddObserver( itk::IterationEvent(), registrationObserver );

			OptimizerCommandIterationUpdate<OptimizerType>::Pointer optimizerObserver = 
				OptimizerCommandIterationUpdate<OptimizerType>::New();
			optimizer->AddObserver( itk::IterationEvent(), optimizerObserver );

			registration->Update(); 
			std::cout<<optimizer->GetStopConditionDescription()<<std::endl;
		}
		catch( itk::ExceptionObject & err ) 
		{ 
			std::cout << "Exception thrown ! " << std::endl;
			std::cout << "An error ocurred during Optimization" << std::endl;
			std::cout << "Location    = " << err.GetLocation()    << std::endl;
			std::cout << "Description = " << err.GetDescription() << std::endl;
		}
	}

	template < typename FixedImageType, typename MovingImageType, typename TransformType >
	void
		ProcrustesAlignment( 
		      typename FixedImageType*   fixedImage, 
			  typename MovingImageType*  movingImage,
		      typename TransformType::Pointer & transform,
			  const std::vector<double> & _scales)
	{
		typedef typename TransformType                         TransformType;
		typedef typename TransformType::Pointer                TransformPointer;
		typedef itk::RegularStepGradientDescentOptimizer       OptimizerType;
		typedef typename OptimizerType::Pointer                OptimizerPointer;
		typedef itk::MeanSquaresImageToImageMetric< 
			FixedImageType, 
			MovingImageType>                                   MetricType;
		typedef typename MetricType::Pointer                   MetricPointer;
		typedef itk::MultiResolutionImageRegistrationMethod< 
			FixedImageType, 
			MovingImageType >                                  RegistrationType;
		typedef typename RegistrationType::Pointer             RegistrationPointer;
		typedef itk::MultiResolutionPyramidImageFilter<
			FixedImageType,
			MovingImageType >                                    FixedImagePyramidType;
		typedef itk::MultiResolutionPyramidImageFilter<
			FixedImageType,
			MovingImageType >                                    MovingImagePyramidType;
		typedef itk::NearestNeighborInterpolateImageFunction<
			MovingImageType>                                     InterpolatorType;
		typedef typename InterpolatorType::Pointer               InterpolatorPointer;

		RegistrationType::Pointer registration = RegistrationType::New();
		MetricType::Pointer       metric       = MetricType::New();
		OptimizerType::Pointer    optimizer    = OptimizerType::New();
		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		FixedImagePyramidType::Pointer fixedImagePyramid = FixedImagePyramidType::New();
		MovingImagePyramidType::Pointer movingImagePyramid = MovingImagePyramidType::New();

		FixedImageType::RegionType fixedRegion = fixedImage->GetLargestPossibleRegion();

		//Metric
		const unsigned int numberOfSamples = static_cast<unsigned int>( fixedRegion.GetNumberOfPixels() * 0.05 );
		metric->SetNumberOfSpatialSamples( numberOfSamples );

		//Optimizer
		typedef typename OptimizerType::ScalesType       OptimizerScalesType;
		OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );
		
		for (int i=0;i<transform->GetNumberOfParameters();i++){
			optimizerScales[i] = _scales[i];
		}
		optimizer->SetScales( optimizerScales );
		optimizer->SetMaximumStepLength( 0.15 );
		optimizer->SetMinimumStepLength( 0.01 );
		optimizer->SetRelaxationFactor( 0.9 );
		optimizer->SetNumberOfIterations( 200 );
		optimizer->SetGradientMagnitudeTolerance( 1e-4 );

		//Registration
		registration->SetFixedImageRegion( fixedRegion );
		registration->SetMetric( metric );
		registration->SetOptimizer( optimizer );
		registration->SetTransform( transform );
		registration->SetFixedImagePyramid( fixedImagePyramid );
		registration->SetMovingImagePyramid( movingImagePyramid );
		registration->SetNumberOfLevels( 2 );
		registration->SetFixedImage( fixedImage );
		registration->SetMovingImage( movingImage );
		registration->SetInterpolator( interpolator );
		registration->SetInitialTransformParameters(  transform->GetParameters()  );

		OptimizerCommandIterationUpdate<OptimizerType>::Pointer observer = OptimizerCommandIterationUpdate<OptimizerType>::New();
		optimizer->AddObserver( itk::IterationEvent(), observer );
		try { 
			std::cout<<"Start to euler registration.." <<std::endl;
			registration->Update(); 
			std::cout<<optimizer->GetStopConditionDescription()<<std::endl;
		}catch( itk::ExceptionObject & err ) { 
			std::cout << "Exception thrown ! " << std::endl;
			std::cout << "An error ocurred during Optimization" << std::endl;
			std::cout << "Location    = " << err.GetLocation()    << std::endl;
			std::cout << "Description = " << err.GetDescription() << std::endl;
		}
	}

	template< typename ImageType, typename MeshType, typename TransformType>
	void
		elasticRegistrationImageToImage( 
		typename ImageType::Pointer & fixedImage, 
		typename ImageType::Pointer & movingImage,
		typename TransformType::Pointer & transform,
		unsigned int bsplineMeshSize = 5,
		const char* probabilityImageFile = NULL)
	{
		itkStaticConstMacro(Dimension, unsigned int, ImageType::ImageDimension);
		const unsigned int SplineOrder = 3;

		typedef itk::MultiResolutionImageRegistrationMethod< ImageType, ImageType > RegistrationType;

		typedef itk::MultiResolutionPyramidImageFilter<ImageType, ImageType > ImagePyramidType;
		typedef itk::MattesMutualInformationWithProbabilityImageToImageMetric< ImageType,ImageType> MetricType;
		typedef itk::LBFGSBOptimizer OptimizerType;
		typedef itk::LinearInterpolateImageFunction<ImageType> InterpolatorType;
		typedef typename TransformType::ParametersType ParametersType;

		RegistrationType::Pointer registration = RegistrationType::New();
		MetricType::Pointer       metric       = MetricType::New();
		OptimizerType::Pointer    optimizer    = OptimizerType::New();
		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		ImagePyramidType::Pointer fixedImagePyramid = ImagePyramidType::New();
		ImagePyramidType::Pointer movingImagePyramid = ImagePyramidType::New();

		ImageType::RegionType fixedRegion = fixedImage->GetLargestPossibleRegion();
		registration->SetFixedImageRegion( fixedRegion );

		//Tranform
		TransformType::PhysicalDimensionsType   fixedPhysicalDimensions;
		TransformType::MeshSizeType             meshSize;
		TransformType::OriginType               fixedOrigin;
		fixedImage->TransformIndexToPhysicalPoint( fixedRegion.GetIndex(), fixedOrigin );
		for( unsigned int i=0; i< ImageType::ImageDimension; i++ )
		{
			fixedPhysicalDimensions[i] = movingImage->GetSpacing()[i] * static_cast<double>( fixedRegion.GetSize()[i] - 1 );
		}
		int splineSize = bsplineMeshSize;
		meshSize.Fill( splineSize - SplineOrder );
		transform->SetTransformDomainOrigin( fixedOrigin );
		transform->SetTransformDomainPhysicalDimensions( fixedPhysicalDimensions );
		transform->SetTransformDomainMeshSize( meshSize );
		transform->SetTransformDomainDirection( fixedImage->GetDirection() );

		//Metric
		const unsigned int numberOfSamples = static_cast<unsigned int>( fixedRegion.GetNumberOfPixels() * 0.01 );
		metric->SetNumberOfSpatialSamples( numberOfSamples );
		metric->SetNumberOfHistogramBins( 24 );
		metric->ReinitializeSeed( 7654321 );

		if (probabilityImageFile){
			typedef MetricType::ProbabilityImageType ProbabilityImageType;
			ProbabilityImageType::Pointer probimage = km::readImage<ProbabilityImageType>( probabilityImageFile );
			metric->SetProbabilityImage( probimage );
		}

		//Optimizer
		OptimizerType::BoundSelectionType boundSelect( transform->GetNumberOfParameters() );
		OptimizerType::BoundValueType upperBound( transform->GetNumberOfParameters() );
		OptimizerType::BoundValueType lowerBound( transform->GetNumberOfParameters() );
		boundSelect.Fill( 2 );
		upperBound.Fill( 200 );
		lowerBound.Fill( -200 );
		optimizer->SetBoundSelection( boundSelect );
		optimizer->SetUpperBound( upperBound );
		optimizer->SetLowerBound( lowerBound );
		optimizer->SetCostFunctionConvergenceFactor( 1.e7 );
		optimizer->SetProjectedGradientTolerance( 1e-5);
		optimizer->SetMaximumNumberOfIterations( 100 );
		optimizer->SetMaximumNumberOfEvaluations( 200 );
		optimizer->SetMaximumNumberOfCorrections( 7 );
		optimizer->TraceOn();

		//Registration
		registration->SetMetric(        metric        );
		registration->SetOptimizer(     optimizer     );
		registration->SetTransform(     transform     );
		registration->SetInterpolator(  interpolator );
		registration->SetFixedImagePyramid( fixedImagePyramid );
		registration->SetMovingImagePyramid( movingImagePyramid );
		registration->SetFixedImage(    fixedImage   );
		registration->SetMovingImage(   movingImage   );
		registration->SetInitialTransformParameters( transform->GetParameters() );
		registration->SetNumberOfLevels( 3 );
		try {
			std::cout<<"Start to elastic registration( coarse ).." <<std::endl;
			registration->Update(); 
			std::cout<<"optimization stoped because: "<<optimizer->GetStopConditionDescription()<<std::endl;
		} catch( itk::ExceptionObject & err ) { 
			std::cout << "ExceptionObject caught !" << std::endl; 
			std::cout << err << std::endl; 
		}
		ParametersType finalParams = registration->GetLastTransformParameters();
		double max = -999;
		double min = 999;
		for (int i=0;i<transform->GetNumberOfParameters();i++){
			if (finalParams[i] > max){
				max = finalParams[i];
			}
			if (finalParams[i] < min){
				min = finalParams[i];
			}
		}
		std::cout<<"result*********** max:"<<max<<", min:"<<min<<std::endl;
	}

	template<class MeshType, class TransformType>
	void
		alignMesh(
		const typename MeshType* fixedMesh, 
		const typename MeshType* movingMesh,
		typename TransformType::Pointer & transform,
		const std::vector<double> & _scales)
	{
		typedef itk::PointSetToPointSetRegistrationMethod<MeshType, MeshType>           RegistrationMethodType;
		typedef itk::EuclideanDistanceMultipleValuePointMetric<MeshType>                P2PMetricType;
		typedef itk::LevenbergMarquardtOptimizer                                        OptimizerType;

		OptimizerType::Pointer          optimizer = OptimizerType::New();
		P2PMetricType::Pointer          p2pmetric = P2PMetricType::New();
		RegistrationMethodType::Pointer registration = RegistrationMethodType::New();

		//Optimizer
		int N = transform->GetNumberOfParameters();
		OptimizerType::ScalesType scales( N );
		for (int t=0;t<N;t++)
		{
			scales[t] = _scales[t];
		}
		//std::cout<<"Optimize scale: "<<scales<<std::endl;
		optimizer->SetScales( scales );
		optimizer->SetNumberOfIterations( 500 );
		optimizer->SetUseCostFunctionGradient( false );
		optimizer->SetGradientTolerance( 1e-5 );
		optimizer->SetValueTolerance( 1e-6 );
		optimizer->SetEpsilonFunction( 1e-11);

		//Registration
		registration->SetMetric( p2pmetric );
		registration->SetOptimizer( optimizer );
		registration->SetTransform( transform );
		registration->SetInitialTransformParameters( transform->GetParameters() );
		registration->SetFixedPointSet(  fixedMesh  );
		registration->SetMovingPointSet( movingMesh );

		try 
		{
			std::cout << "starting model fitting" << std::endl;
			registration->Update();
			std::cout << optimizer->GetStopConditionDescription() << std::endl;
		}
		catch ( itk::ExceptionObject& o ) 
		{
			std::cout << "caught exception " << o << std::endl;
		}
	}
	
	template<class DistanceMapType, class MeshType, class TransformType>
	void
		transformFittingToDistanceMap(
		const typename DistanceMapType * fixedDistMap,
		const typename MeshType * movingMesh,
		typename TransformType::Pointer & transform,
		unsigned int numberOfShapeParameters,
		const std::vector<double> & _scales,
		unsigned int maxiterations = 500)
	{
		TransformType::ParametersType preParam = transform->GetParameters();
		typedef itk::PointSetToImageMultipleValueRegistrationMethod<MeshType, DistanceMapType> RegistrationMethodType;
		typedef itk::MinValuePointSetToImageMultipleValueMetric<MeshType, DistanceMapType>     MetricType;
		typedef itk::LevenbergMarquardtOptimizer                                  OptimizerType;
		typedef itk::LinearInterpolateImageFunction<DistanceMapType, double>      InterpolatorType;
		OptimizerType::Pointer          optimizer = OptimizerType::New();
		MetricType::Pointer             metric = MetricType::New();
		RegistrationMethodType::Pointer registration = RegistrationMethodType::New();
		InterpolatorType::Pointer       interpolator = InterpolatorType::New();

		//Metric
		metric->SetNumberOfShapeParameters( numberOfShapeParameters );
		metric->SetWeightForShapePenalty( km::g_shape_penalty ); //Typically we want this to stay between 0.05 and 0.3; 

		//Optimizer
		int N = transform->GetNumberOfParameters();
		OptimizerType::ScalesType scales( N );
		for ( int i=0;i<N;i++ )
		{
			scales[i] = _scales[i];
		}
		//std::cout<<"Optimize scale: "<<scales<<std::endl;
		//std::cout<<maxiterations<<std::endl;
		optimizer->SetScales( scales );
		optimizer->SetNumberOfIterations( maxiterations );
		optimizer->SetUseCostFunctionGradient( false );
		optimizer->SetGradientTolerance( 1e-8 );
		optimizer->SetValueTolerance( 1e-8 );
		optimizer->SetEpsilonFunction( 1e-8);
		//optimizer->AddObserver( itk::IterationEvent(), observer );

		//Registration
		registration->SetOptimizer( optimizer );
		registration->SetInterpolator( interpolator );
		registration->SetMetric(    metric );

		registration->SetTransform( transform );
		registration->SetInitialTransformParameters( transform->GetParameters() );
		registration->SetFixedPointSet(  movingMesh  );
		registration->SetMovingImage( fixedDistMap );

		try {
			std::cout << "starting model fitting" << std::endl;
			registration->Update();
			std::cout << optimizer->GetStopConditionDescription() << std::endl;

		} catch ( itk::ExceptionObject& o ) {
			std::cout << "caught exception " << o << std::endl;
		}
	}
}


#endif