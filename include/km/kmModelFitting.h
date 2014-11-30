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
#include "itkCovariantVector.h"

#include "kmProfileClassifier.h"
#include "kmVtkItkUtility.h"
#include "kmGlobal.h"
#include "kmUtility.h"

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

	//template<class ClassifiedPointsKdTreeType, class MeshType>
	//void
	//	generateBestPointSet2(
	//	typename ClassifiedPointsKdTreeType* classifiedKdTree, 
	//	typename MeshType* outputMesh,
	//	const typename MeshType* liverMesh)
	//{
	//	typedef MeshType::GeometryMapType GeometryMapType;
	//	typedef GeometryMapType::Iterator GeometryMapIterator;
	//	typedef MeshType::PointType PointType;
	//	typedef PointType::VectorType VectorType;

	//	km::assigneMesh<MeshType>(outputMesh, 0.0);

	//	GeometryMapIterator geoIt = liverMesh->GetGeometryData()->Begin();
	//	GeometryMapIterator geoItEnd = liverMesh->GetGeometryData()->End();

	//	itk::SimplexMeshGeometry *geodata;
	//	while (geoIt!=geoItEnd)
	//	{
	//		int idx = geoIt.Index();
	//		geodata = geoIt.Value();

	//		VectorType normal;
	//		normal.Set_vnl_vector(geodata->normal.Get_vnl_vector());

	//		PointType curPoint = liverMesh->GetPoint(idx);
	//		PointType closedBoundaryPoint;

	//		//std::cout<<"Find closest bounday point for: "<<idx<<std::endl;
	//		bool found = classifiedKdTree->FindClosetBoundaryPoint(closedBoundaryPoint, curPoint, idx);

	//		double offset = 0.0;
	//		if (found)
	//		{
	//			VectorType shiftVec = closedBoundaryPoint - curPoint;
	//			offset = shiftVec*normal;
	//		}

	//		if (isAbnormal(idx))
	//		{
	//			offset = offset>0?-1.0:1.0;
	//		}

	//		outputMesh->SetPointData(idx, offset);

	//		geoIt++;
	//	}

	//	//Remove noise point.
	//	km::smoothMeshData<MeshType>(outputMesh, 0);

	//	geoIt = liverMesh->GetGeometryData()->Begin();
	//	geoItEnd = liverMesh->GetGeometryData()->End();

	//	while (geoIt!=geoItEnd)
	//	{
	//		unsigned int idx = geoIt.Index();
	//		geodata = geoIt.Value();

	//		VectorType normal;
	//		normal.Set_vnl_vector(geodata->normal.Get_vnl_vector());

	//		double offsetVal = 0.0;
	//		outputMesh->GetPointData(idx, &offsetVal);

	//		PointType oldPos = liverMesh->GetPoint(idx);
	//		outputMesh->SetPoint(idx, oldPos + normal*offsetVal);

	//		geoIt++;
	//	}

	//}

	template<class ClassifiedPointsKdTreeType, class ProfileExtractorType, class FloatImageType, class MeshType>
	void
		deformByLiverClassification(
		km::ProfileClassifier * classifier,
		ClassifiedPointsKdTreeType * classifiedKdTree,
		ProfileExtractorType * profileExtractor,
		typename MeshType* outputMesh,
		const typename MeshType* liverMesh,
		double searchStep = 1.5,
		unsigned int maxSearchPoints = 20)
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
			if (isAbnormal(idx))
			{
				//std::cout<<"Find closest bounday point for: "<<idx<<std::endl;
				PointType closedBoundaryPoint;
				bool found = classifiedKdTree->FindClosetBoundaryPoint(closedBoundaryPoint, mpoint, idx);
				if (found)
				{
					VectorType shiftVec = closedBoundaryPoint - mpoint;
					best_offset = shiftVec*normal;
				}
			}
			else
			{
				double maxDist = g_varianceMap[idx];
				double direction = 1.0;
				unsigned int searchedPoints = 1;
				while(searchedPoints<=maxSearchPoints && std::abs(best_offset)<=maxDist)
				{
					PointType pttest = mpoint + normal*best_offset;
					double tmpDirection = -1.0;
					if (profileExtractor->isInsideBuffer(pttest))
					{
						std::vector<FLOATTYPE> feature;
						profileExtractor->extractFeatureSet(feature, classifier->profile_category, geodata, pttest);

						tmpDirection = 2.0*classifier->classify(feature,idx) - 1.0;
					}
					best_offset += tmpDirection * searchStep;
					if (searchedPoints == 1)
					{
						direction = tmpDirection>0?1.0:-1.0;
					}
					else if (tmpDirection*direction<0) //Change direction. Break from here.
					{
						break;
					}
					searchedPoints++;
				}
			}

			if (isAbnormal(idx))
			{
				best_offset = best_offset>0?-1.0*searchStep:searchStep;
			}

			outputMesh->SetPointData(idx, best_offset);

			geoIt++;
		}

		//Remove noise point.
		km::smoothMeshData<MeshType>(outputMesh, 0);

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

	enum FittingType
	{
		Rigid = 0,
		Shape
	};

	template<class MatrixType, class TransformType>
	void FillRigidTransform(const MatrixType & mat, TransformType * transform)
	{
		TransformType::MatrixType rigidMatrix = transform->GetMatrix();
		TransformType::OffsetType rigidOffset = transform->GetOffset();
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
		transform->SetMatrix( rigidMatrix );
		transform->SetOffset( rigidOffset );
	}

	template<class MatrixType, class TransformType>
	void FillShapeTransform(const MatrixType & mat, TransformType * transform)
	{
		TransformType::ParametersType shapeParams = transform->GetParameters();
		for (int i=0;i<transform->GetUsedNumberOfCoefficients();i++)
		{
			if (mat[i]<-3.0){
				shapeParams[i] = -3;
			}else if (mat[i]>3.0){
				shapeParams[i] = 3;
			}else{
				shapeParams[i] = mat[i];	
			}
		}
		transform->SetParameters(shapeParams);
		KM_DEBUG_PRINT("Shape paramters after fitting", shapeParams);
	}

	void initAbnormalMap(unsigned int N, bool reset = false)
	{
		if (g_fittingErrorMap.size()==N && !reset)
		{
			return;
		}
		else
		{
			for (int idx=0;idx<N;idx++)
			{
				g_fittingErrorMap[idx] = 0.0;
			}
		}
	}

	bool isAbnormal(int idx)
	{
		return g_fittingErrorMap[idx]>=g_fitting_error_threshold?true:false;
	}

	unsigned int countAbnormal()
	{
		unsigned int count = 0;
		for (int idx=0;idx<g_fittingErrorMap.size();idx++)
		{
			if(isAbnormal(idx)) count++;
		}
		return count;
	}

	template<class MeshType, class MatrixType>
	void fillPointSetIntoMatrix(const MeshType* mesh, MatrixType& matrix, FittingType rigidOrShape)
	{
		//KM_DEBUG_INFO("fillPointSetIntoMatrix");
		unsigned int Dimension = MeshType::PointDimension;
		unsigned int numberOfLandmarks = mesh->GetNumberOfPoints();

		if (rigidOrShape == Rigid)
		{
			matrix.setConstant(numberOfLandmarks, Dimension+1, 1.0);
			unsigned long rowcnt = 0;
			for (int idx=0;idx<numberOfLandmarks;idx++)
			{
				MeshType::PointType pt = mesh->GetPoint(idx);
				for (int d=0;d<Dimension;d++)
				{
					matrix(rowcnt, d) = pt[d];
				}
				rowcnt++;
			}
		}
		else
		{
			matrix.setZero( numberOfLandmarks*Dimension, 1);
			unsigned long rowcnt = 0;
			for (int idx=0;idx<numberOfLandmarks;idx++)
			{
				MeshType::PointType pt = mesh->GetPoint(idx);
				for (int d=0;d<Dimension;d++)
				{
					matrix(rowcnt, 0) = pt[d];
					rowcnt++;
				}
			}
		}
	}

	template<class MatrixType>
	void
		removeAbnormal(const typename MatrixType & matrix, typename MatrixType& matrixUpdated, FittingType rigidOrShape)
	{
		if (countAbnormal() == 0)
		{
			matrixUpdated = matrix;
			return;
		}

		unsigned int Dimension = 3;
		if (rigidOrShape == Rigid)
		{
			unsigned int numberOfLandmarks = matrix.rows();
			unsigned int numberOfAbnormal = countAbnormal();
			matrixUpdated.setZero(numberOfLandmarks-numberOfAbnormal, matrix.cols());
			unsigned long ptcount = 0;
			for (int idx=0;idx<numberOfLandmarks;idx++)
			{
				if (!isAbnormal(idx))
				{
					for (int c=0;c<matrix.cols();c++)
					{
						matrixUpdated(ptcount, c) = matrix(idx, c);
					}
					ptcount++;
				}
			}
			KM_ASSERT(ptcount == matrixUpdated.rows());
		}
		else
		{
			unsigned int numberOfLandmarks = matrix.rows()/Dimension;
			unsigned int numberOfAbnormal = countAbnormal();
			matrixUpdated.setZero((numberOfLandmarks-numberOfAbnormal)*Dimension, matrix.cols());
			unsigned long ptcount = 0;
			for (int idx=0;idx<numberOfLandmarks;idx++)
			{
				if (!isAbnormal(idx))
				{
					for (int d=0;d<Dimension;d++)
					{
						for (int c=0;c<matrix.cols();c++)
						{
							matrixUpdated(ptcount*Dimension+d, c) = matrix(idx*Dimension+d, c);
						}
					}
					ptcount++;
				}
			}
			KM_ASSERT(ptcount*Dimension == matrixUpdated.rows());
		}

		//KM_DEBUG_PRINT("Matrix rows before removing ", matrix.rows());
		//KM_DEBUG_PRINT("Matrix rows after removing ", matrixUpdated.rows());
	}

	template<class MeshType>
	void
		updateFittingErrorMap(
		const typename MeshType * tagetMesh,
		const typename MeshType * unfittedMesh,
		const typename MeshType * fittedMesh)
	{
		KM_ASSERT(tagetMesh->GetNumberOfPoints() == unfittedMesh->GetNumberOfPoints());

		typedef MeshType::PointType PointType;
		typedef PointType::VectorType VectorType;
		for (int idx=0;idx<tagetMesh->GetNumberOfPoints();idx++)
		{			
			PointType targetPoint = tagetMesh->GetPoint(idx);
			PointType unfittedPoint = unfittedMesh->GetPoint(idx);
			PointType fittedPoint = fittedMesh->GetPoint(idx);

			double errorUpdated = 0.0;
			double unfittedDist = unfittedPoint.EuclideanDistanceTo(targetPoint);

			if (isAbnormal(idx))
			{
				errorUpdated = 0.0;
			}
			else if (unfittedDist < 5.0)
			{
				errorUpdated = 0.0;
			}
			else
			{
				double fittedDist = fittedPoint.EuclideanDistanceTo(targetPoint);
				errorUpdated = fittedDist / unfittedDist;
			}

			g_fittingErrorMap[idx] = g_fittingErrorMap[idx]*0.8 + errorUpdated*0.2;
		}

		unsigned int abnormalCount = countAbnormal();
		if (abnormalCount > g_fittingErrorMap.size() * 0.25)
		{
			KM_DEBUG_ERROR("Abnormal points has exceed 25%");
			//initAbnormalMap(g_fittingErrorMap.size(), true);
		}
	}

	template<class MeshType, class StatisticalModelType, class RigidTransformType, class ShapeTransformType>
	void
		compositeTransformFitting( 
		const MeshType* targetMesh, 
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

		unsigned int numberOfLandmarks = referenceShapeMesh->GetNumberOfPoints();
		km::initAbnormalMap(numberOfLandmarks);

		typedef StatisticalModelType::ImplType StatisticalModelImplType;
		typedef statismo::MatrixType StatismoMatrixType;
		typedef statismo::VectorType StatismoVectorType;

		{
			MeshType::ConstPointer sourceMeshForRigid = km::transformMesh<MeshType, ShapeTransformType>(referenceShapeMesh, shapeTransform);
			MeshType::ConstPointer targetMeshForRigid = targetMesh;

			//StatismoMatrixType sourceMatrixForRigidTmp, targetMatrixForRigidTmp;
			StatismoMatrixType sourceMatrixForRigid, targetMatrixForRigid;
			fillPointSetIntoMatrix<MeshType, StatismoMatrixType>(sourceMeshForRigid, sourceMatrixForRigid, Rigid);
			fillPointSetIntoMatrix<MeshType, StatismoMatrixType>(targetMeshForRigid, targetMatrixForRigid, Rigid);

			//removeAbnormal<StatismoMatrixType>(sourceMatrixForRigidTmp, sourceMatrixForRigid, Rigid);
			//removeAbnormal<StatismoMatrixType>(targetMatrixForRigidTmp, targetMatrixForRigid, Rigid);

			StatismoVectorType I = StatismoVectorType::Ones(sourceMatrixForRigid.cols());
			StatismoMatrixType Mmatrix = sourceMatrixForRigid.transpose() * sourceMatrixForRigid;
			Mmatrix.diagonal() += I;
			StatismoMatrixType MInverseMatrix = Mmatrix.inverse();
			const StatismoMatrixType& WT = sourceMatrixForRigid.transpose();
			StatismoMatrixType coeffsRigid = MInverseMatrix * (WT * targetMatrixForRigid);

			FillRigidTransform<StatismoMatrixType, RigidTransformType>(coeffsRigid.transpose(), rigidTransform);
		}

		{
			/*****************Fit shape parameters****************/

			RigidTransformType::Pointer inversedRigidTransform = RigidTransformType::New();
			rigidTransform->GetInverse( inversedRigidTransform );
			MeshType::Pointer targetMeshForShape = km::transformMesh<MeshType, RigidTransformType>(targetMesh, inversedRigidTransform);
			MeshType::Pointer sourceMeshForShape = model->DrawMean();

			//StatismoVectorType sourceMatrixForShapeTmp, targetMatrixForShapeTmp;
			StatismoVectorType sourceMatrixForShape, targetMatrixForShape;
			//const StatismoMatrixType & basisMatrixForShapeTmp = model->GetstatismoImplObj()->GetPCABasisMatrix();
			//StatismoMatrixType basisMatrixForShape;
			const StatismoMatrixType & basisMatrixForShape = model->GetstatismoImplObj()->GetPCABasisMatrix();
			fillPointSetIntoMatrix<MeshType, StatismoVectorType>(sourceMeshForShape, sourceMatrixForShape, Shape);
			fillPointSetIntoMatrix<MeshType, StatismoVectorType>(targetMeshForShape, targetMatrixForShape, Shape);

			//removeAbnormal<StatismoVectorType>(sourceMatrixForShapeTmp, sourceMatrixForShape, Shape);
			//removeAbnormal<StatismoVectorType>(targetMatrixForShapeTmp, targetMatrixForShape, Shape);
			//removeAbnormal<StatismoMatrixType>(basisMatrixForShapeTmp, basisMatrixForShape, Shape);

			StatismoVectorType I = StatismoVectorType::Ones(basisMatrixForShape.cols());
			StatismoMatrixType Mmatrix = basisMatrixForShape.transpose() * basisMatrixForShape;
			Mmatrix.diagonal() += I;
			StatismoMatrixType MInverseMatrix = Mmatrix.inverse();
			const StatismoMatrixType& WT = basisMatrixForShape.transpose();
			StatismoVectorType coeffsShape = MInverseMatrix * (WT * (targetMatrixForShape-sourceMatrixForShape));

			MeshType::Pointer shapeUnfittedMesh = km::transformMesh<MeshType, ShapeTransformType>(referenceShapeMesh, shapeTransform);

			FillShapeTransform<StatismoVectorType, ShapeTransformType>(coeffsShape, shapeTransform);

			MeshType::Pointer shapeFittedMesh = km::transformMesh<MeshType, ShapeTransformType>(referenceShapeMesh, shapeTransform);

			if(g_disable_abnormal)
			{
				km::updateFittingErrorMap<MeshType>(targetMeshForShape, shapeUnfittedMesh, shapeFittedMesh);
			}
		}
	}

	template<class StatisticalModelType, class MeshType>
	void
		initModelVariance(const typename StatisticalModelType * model, const typename MeshType* mesh, int numberOfComponents = 10)
	{
		typedef itk::SimplexMeshGeometry SimplexMeshGeometryType;
		typedef SimplexMeshGeometryType::CovariantVectorType CovariantVectorType;

		SimplexMeshGeometryType *geodata;

		StatisticalModelType::MatrixType basisMatrix = model->GetPCABasisMatrix();
		for (int id=0;id<mesh->GetNumberOfPoints();id++)
		{
			double variance = 0.0;
			geodata = mesh->GetGeometryData()->GetElement(id);
			for (unsigned p=0;p<numberOfComponents;p++)
			{
				double varOnNormal = 0.0;
				for (unsigned d=0; d<Dimension; d++) 
				{
					unsigned idx = model->GetRepresenter()->MapPointIdToInternalIdx(id, d);
					varOnNormal += geodata->normal[d] * basisMatrix[idx][p];
				}
				variance += std::abs(varOnNormal);
			}
			g_varianceMap[id] = variance;
		}
	}
}

#endif