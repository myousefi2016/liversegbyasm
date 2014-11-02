#ifndef __kmModelFitting_h
#define __kmModelFitting_h

#include <map>

#include "itkSimplexMesh.h"
#include "Representers/ITK/itkSimplexMeshRepresenter.h"
#include "statismo_ITK/itkStatisticalModel.h"

#include "itkEuclideanDistanceKdTreeMultipleValuePointSetMetric.h"
#include "itkEuclideanDistanceMultipleValuePointMetric.h"
//#include <itkEuclideanDistancePointMetric.h>
#include "itkLevenbergMarquardtOptimizer.h"
#include "itkOnePlusOneEvolutionaryOptimizer.h"
#include "itkNormalVariateGenerator.h"
#include "itkPointSetToPointSetRegistrationMethod.h"

#include "itkMinValuePointSetToImageMultipleValueMetric.h"
#include "itkPointSetToImageMultipleValueRegistrationMethod.h"
#include "itkPointSetToImageRegistrationMethod.h"

#include "itkPointSetToPointSetSingleValueRegistrationMethod.h"
#include "itkCorrespondingPointsEuclideanDistancePointMetric.h"
#include "itkLBFGSBOptimizer.h"

#include "itkLinearInterpolateImageFunction.h"

#include "itkSimplexMeshGeometry.h"
#include "itkConstNeighborhoodIterator.h"

#include "kmKNNProfileClassifier-FLANN.h"
#include "kmVtkItkUtility.h"
#include "kmGlobal.h"

using namespace km;

namespace km
{
	template<class OptimizerType>
	class IterationStatusObserver : public itk::Command
	{
	public:
		typedef  IterationStatusObserver   Self;
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
			OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >( object );
			if( ! itk::IterationEvent().CheckEvent( &event ) )
			{
				return;
			}
			std::cout << "Iteration: " << ++m_iter_no << " , Metric: " << optimizer->GetValue() << ", Parameters: " << optimizer->GetCachedCurrentPosition() << std::endl;
			//std::cout << "Iteration: " << ++m_iter_no << " , Metric: " << optimizer->GetValue() << std::endl;
		}
	protected:
		IterationStatusObserver():
			 m_iter_no(0)     {};
			 virtual ~IterationStatusObserver(){};
	private:
		int m_iter_no;
	};

	enum KM_MESH_DIST_METRIC
	{
		KdTree = 1,
		P2P
	};

	template<class MeshType, class TransformType>
	void
		transformFitting(
		const typename MeshType * fixedMesh,
		const typename MeshType * movingMesh,
		typename TransformType::Pointer & transform,
		unsigned int numberOfShapeParameters,
		KM_MESH_DIST_METRIC metricType,
		const std::vector<double> & _scales,
		unsigned int maxiterations = 500)
	{
		typedef itk::PointSetToPointSetRegistrationMethod<MeshType, MeshType>           RegistrationMethodType;
		typedef itk::EuclideanDistanceKdTreeMultipleValuePointSetMetric<MeshType>       KdTreeMetricType;
		typedef itk::EuclideanDistanceMultipleValuePointMetric<MeshType>                P2PMetricType;
		
		typedef itk::LevenbergMarquardtOptimizer                                        OptimizerType;

		typedef IterationStatusObserver<OptimizerType>                                  ObserverType;

		OptimizerType::Pointer          optimizer = OptimizerType::New();
		ObserverType::Pointer           observer = ObserverType::New();
		KdTreeMetricType::Pointer       kdtreemetric = KdTreeMetricType::New();
		P2PMetricType::Pointer          p2pmetric = P2PMetricType::New();
		RegistrationMethodType::Pointer registration = RegistrationMethodType::New();

		//Metric
		kdtreemetric->SetNumberOfShapeParameters( numberOfShapeParameters );
		kdtreemetric->SetWeightForShapePenalty( 0 );

		p2pmetric->SetNumberOfShapeParameters( numberOfShapeParameters );
		p2pmetric->SetWeightForShapePenalty( 0 );

		//Optimizer
		int N = transform->GetNumberOfParameters();
		OptimizerType::ScalesType scales( N );
		for (int t=0;t<N;t++)
		{
			scales[t] = _scales[t];
		}
		//std::cout<<"Optimize scale: "<<scales<<std::endl;
		optimizer->SetScales( scales );
		optimizer->SetNumberOfIterations( maxiterations );
		optimizer->SetUseCostFunctionGradient( false );
		optimizer->SetGradientTolerance( 1e-5 );
		optimizer->SetValueTolerance( 1e-6 );
		optimizer->SetEpsilonFunction( 1e-11);
		//optimizer->AddObserver( itk::IterationEvent(), observer );

		//Registration
		if (metricType == KdTree)
		{
			registration->SetMetric( kdtreemetric );
		}
		else if (metricType == P2P)
		{
			registration->SetMetric( p2pmetric );
		}
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
	double
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

		//std::cout<<optimizer->GetNumberOfIterations()<<std::endl;

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

		TransformType::ParametersType postParam = transform->GetParameters();
		//Calculate parameters difference
		double paradiff = 0;
		for (int k=0;k<numberOfShapeParameters;k++)
		{
			double diff = std::abs(preParam[k]-postParam[k]);
			if (diff > paradiff)
			{
				paradiff = diff; 
			}
		}

		return paradiff;
	}
}

#endif