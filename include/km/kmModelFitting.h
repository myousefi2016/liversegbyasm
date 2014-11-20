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
#include "itkCompositeTransform.h"

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

	template<class FloatImageType, class MeshType, class GradientInterpolatorType, class IntensityInterpolatorType>
	void
		generateBestPointSet(
		km::AdaboostProfileClassifier * classifier,
		const typename GradientInterpolatorType* gradientInterpolator, 
		const typename IntensityInterpolatorType* intensityInterpolator,
		typename MeshType* outputMesh,
		const typename MeshType* liverMesh,
		const typename MeshType* varianceMap,
		const double spacingDistForTest = 1.5,
		const int numberOfPointsForTest = 10)
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
					classifier->profile_category);

				double current_confidence = 2*classifier->classify( feature, idx ) - 1.0;

				direction = current_confidence;
			}
			if (classifier->profile_category == LIVER)
			{
				double pointVariance = 1.0;
				varianceMap->GetPointData(idx, &pointVariance);

				double pointError = 0.0;

				int tmpNumberOfPointsForTest = numberOfPointsForTest;//static_cast<int>(numberOfPointsForTest*pointVariance+3.0);

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

						confidence = 2*confidence - 1.0;
					}
					else
					{
						i = tmpNumberOfPointsForTest;
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
					else
					{
						i = tmpNumberOfPointsForTest;
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

		if (classifier->profile_category != BOUNDARY)
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

	template<class MatrixType, class RigidTransformType>
	void FillRigidTransform(const MatrixType & mat, RigidTransformType * rigidTransform)
	{
		RigidTransformType::MatrixType rigidMatrix = rigidTransform->GetMatrix();
		RigidTransformType::OffsetType rigidOffset = rigidTransform->GetOffset();

		for (int i=0;i<3;i++)
		{
			for (int j=0;j<3;j++)
			{
				rigidMatrix[i][j] = mat(i, j);
			}
		}
		//std::cout<<rigidMatrix<<std::endl;

		for (int i=0;i<Dimension;i++)
		{
			rigidOffset[i] = mat(i, Dimension);
		}
		//std::cout<<rigidOffset<<std::endl;

		rigidTransform->SetMatrix( rigidMatrix );
		rigidTransform->SetOffset( rigidOffset );
	}

	template<class VectorType, class ShapeTransformType>
	void FillShapeTransform(const VectorType & vec, ShapeTransformType * shapeTransform)
	{
		ShapeTransformType::ParametersType shapeParams = shapeTransform->GetParameters();
		for (int i=0;i<shapeTransform->GetUsedNumberOfCoefficients();i++)
		{
			if (vec[i]<-3.0){
				shapeParams[i] = -3;
			}else if (vec[i]>3.0){
				shapeParams[i] = 3;
			}else{
				shapeParams[i] = vec[i];	
			}
		}

		shapeTransform->SetParameters(shapeParams);
	}

	template<class MatrixType, class MeshType>
	void fillPointSetIntoMatrix(MatrixType& matrix, const MeshType * mesh)
	{
		typedef MeshType::PointType PointType;
		typedef MeshType::PointIdentifier PointIdentifier;
		typedef MeshType::PointsContainerConstIterator PointsContainerConstIterator;

		PointsContainerConstIterator ptIt = mesh->GetPoints()->Begin();
		PointsContainerConstIterator ptItEnd = mesh->GetPoints()->End();
		while(ptIt != ptItEnd)
		{
			PointIdentifier idx = ptIt.Index();
			PointType pt = ptIt.Value();
			for (int d=0;d<Dimension;d++)
			{
				matrix(idx, d) = pt[d];
			}
			ptIt++;
		}
	}

	template<class MeshType, class StatisticalModelType, class RigidTransformType, class ShapeTransformType>
	void
		compositeTransformFitting( const MeshType* targetMesh, 
		const StatisticalModelType* model, 
		RigidTransformType * rigidTransform,
		ShapeTransformType * shapeTransform)
	{
		unsigned int Dimension = MeshType::PointType::PointDimension;

		if (shapeTransform->GetUsedNumberOfCoefficients() > shapeTransform->GetNumberOfParameters())
		{
			shapeTransform->SetUsedNumberOfCoefficients(shapeTransform->GetNumberOfParameters());
		}

		MeshType::Pointer referenceShapeMesh = model->GetRepresenter()->GetReference();
		MeshType::Pointer meanShapeMesh = model->DrawMean();
		unsigned int numberOfLandmarks = referenceShapeMesh->GetNumberOfPoints();

		typedef statismo::MatrixType StatismoMatrixType;
		typedef statismo::VectorType StatismoVectorType;
		typedef StatisticalModelType::VectorType VnlVectorType;

		StatismoMatrixType matrixToBeFitted;
		matrixToBeFitted.setConstant(numberOfLandmarks, Dimension+1, 1.0);

		StatismoMatrixType matrixTargetForRigid;
		matrixTargetForRigid.setConstant(numberOfLandmarks, Dimension+1, 1.0);

		StatismoMatrixType matrixTargetForShape;
		matrixTargetForShape.setConstant(numberOfLandmarks, Dimension+1, 1.0);

		StatismoVectorType I = StatismoVectorType::Ones(matrixToBeFitted.cols());

		unsigned int iterNum = 0;
		while(iterNum < 1)
		{
			std::cout<<"*****************Fit rigid parameters " << iterNum << " ****************"<<std::endl;

			MeshType::ConstPointer shapeFittedMesh = km::transformMesh<MeshType, ShapeTransformType>(referenceShapeMesh, shapeTransform);
			MeshType::ConstPointer targetMeshForRigid = targetMesh;

			fillPointSetIntoMatrix<StatismoMatrixType, MeshType>(matrixToBeFitted, shapeFittedMesh);
			fillPointSetIntoMatrix<StatismoMatrixType, MeshType>(matrixTargetForRigid, targetMesh);

			StatismoMatrixType Mmatrix = matrixToBeFitted.transpose() * matrixToBeFitted;
			Mmatrix.diagonal() += I;

			StatismoMatrixType MInverseMatrix = Mmatrix.inverse();

			const StatismoMatrixType& WT = matrixToBeFitted.transpose();

			StatismoMatrixType coeffsRigid = MInverseMatrix * (WT * matrixTargetForRigid);

			FillRigidTransform<StatismoMatrixType, RigidTransformType>(coeffsRigid.transpose(), rigidTransform);

			/*****************Fit shape parameters****************/
			RigidTransformType::Pointer inversedRigidTransform = RigidTransformType::New();
			rigidTransform->GetInverse( inversedRigidTransform );

			MeshType::Pointer targetMeshForShape = km::transformMesh<MeshType, RigidTransformType>(targetMesh, inversedRigidTransform);
			VnlVectorType coeffsShape = model->ComputeCoefficientsForDataset( targetMeshForShape );

			FillShapeTransform<VnlVectorType, ShapeTransformType>(coeffsShape, shapeTransform);

			iterNum++;
		}
	}
}

#endif