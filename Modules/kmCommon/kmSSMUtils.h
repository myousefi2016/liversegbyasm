#ifndef __kmSSMUtils_h
#define __kmSSMUtils_h

#include <itkVariableLengthVector.h>
#include "itkListSample.h"
#include "itkKdTree.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkMinimumDecisionRule.h"
#include "itkSampleClassifierFilter.h"
#include "itkKdTreeBasedKmeansEstimator.h"

#include "vtkTimerLog.h"
#include "vtkNew.h"

#include <fstream>

#include <boost/random.hpp>

#include "statismo/CommonTypes.h"
#include "kmUtility.h"
#include "kmGlobal.h"

namespace km
{
	enum LandmarkStatus
	{
		Normal = 0,
		Relaxed,
		Abnormal,
		Leaking,
		Optimized
	};

	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	class SSMUtils
	{
	public:
		typedef typename TMesh                              MeshType;
		typedef typename MeshType::PointType                PointType;
		typedef typename PointType::VectorType              VectorType;
		typedef typename TRigidTransform                    RigidTransformType;
		typedef typename TShapeTransform                    ShapeTransformType;
		typedef typename RigidTransformType::ParametersType RigidParametersType;
		typedef typename ShapeTransformType::ParametersType ShapeParametersType;
		typedef typename TStatisticalModel                  StatisticalModelType;
		typedef typename StatisticalModelType::ImplType     StatisticalModelImplType;
		typedef statismo::MatrixType                        StatismoMatrixType;
		typedef statismo::VectorType                        StatismoVectorType;
		itkStaticConstMacro(Dimension, unsigned int, MeshType::PointDimension);

		class ShapeClusterItem
		{
		public:
			int clusterId;
			std::vector<int> pointIds;
			double shapeProbability;
			double meanError;
			StatismoVectorType sourceMatrixForShape, targetMatrixForShape;
			StatismoMatrixType basisMatrix, MInverseMatrix, WT;
			typename ShapeTransformType::Pointer shapeTransform;
			typename ShapeParametersType shapeParametersPre;

			StatismoMatrixType sourceMatrixForRigid, targetMatrixForRigid;
			typename RigidTransformType::Pointer rigidTransform;
			typename RigidParametersType rigidParamtersPre;

			typedef itk::Statistics::ListSample<PointType> SampleType;
			typedef itk::Statistics::WeightedCentroidKdTreeGenerator<SampleType> TreeGeneratorType;

			typename TreeGeneratorType::Pointer kdTreeGenerator;
			typename SampleType::Pointer meanPoints;

			ShapeClusterItem(int id)
			{
				this->clusterId = id;
				shapeProbability = 0.0;
				meanError = 0.0;

				this->meanPoints = SampleType::New();
				this->meanPoints->SetMeasurementVectorSize( 3 );

				this->kdTreeGenerator = TreeGeneratorType::New();
				this->kdTreeGenerator->SetBucketSize( 4 );
			}

			PointType findClosestPoint(PointType & pt)
			{
				TreeGeneratorType::KdTreeType::InstanceIdentifierVectorType neighbors;
				kdTreeGenerator->GetOutput()->Search( pt, static_cast<unsigned>(1), neighbors );
				PointType closestPt = kdTreeGenerator->GetOutput()->GetMeasurementVector( neighbors[0] );
				return closestPt;
			}

			int getNumberOfPoints()
			{
				return pointIds.size();
			}
		};
		
		SSMUtils()
		{
			m_NumberOfCoefficients = 3;
			m_Iterations = 0;
			m_NumberOfClusters = 1;
		}

		~SSMUtils()
		{

		}
		
		void SetSSM(const StatisticalModelType* model)
		{
			m_SSM = const_cast<StatisticalModelType*>( model );
		}
		
		void SetRigidTransform(RigidTransformType* rtfm)
		{
			m_RigidTranform = rtfm;
		}
		
		void SetShapeTransform(ShapeTransformType* stfm)
		{
			m_ShapeTransform = stfm;
			m_NumberOfCoefficients = m_ShapeTransform->GetUsedNumberOfCoefficients();
			if (m_NumberOfCoefficients > m_ShapeTransform->GetNumberOfParameters()){
				m_NumberOfCoefficients = m_ShapeTransform->GetNumberOfParameters();
			}
			std::cout<<"Number Of Coefficients: "<<m_NumberOfCoefficients<<std::endl;
		}

		void SetNumberOfClusters(int numberOfClusters)
		{
			m_NumberOfClusters = numberOfClusters;
		}

		void Initialize();
		
		void Update(const MeshType * targetMesh, MeshType * outputMesh);
		
		void PrintTransform();
		
		double CalShapeParaDiff();
		
		ShapeClusterItem* GetClusterByClusterId(int clusterId)
		{
			ShapeClusterItem* clusterItem = m_ShapeClusterInstances[clusterId];
			if (clusterItem==NULL){
				std::cerr<<"Cannot find cluster Id: "<<clusterId<<std::endl;
				clusterItem = new ShapeClusterItem(0);
			}

			return clusterItem;
		}
		
		ShapeClusterItem* GetClusterByPointId(int pointId)
		{
			int clusterId = m_ShapeClusterMap[pointId];
			return GetClusterByClusterId(clusterId);
		}

		void SetLandmarkStatus(int pointId, LandmarkStatus status)
		{
			m_LandmarkStatus[pointId] = status;
		}

		LandmarkStatus GetLandmarkStatus(int pointId)
		{
			return m_LandmarkStatus[pointId];
		}

	private:
		vtkNew<vtkTimerLog> timer;
		typename StatisticalModelType* m_SSM;
		typename RigidTransformType::Pointer m_RigidTranform;
		typename ShapeTransformType::Pointer m_ShapeTransform;
		typename ShapeParametersType m_ShapeParametersPre;
		typename MeshType::Pointer m_ReferenceShapeMesh;
		typename MeshType::Pointer m_MeanShapeMesh;
		typename MeshType::Pointer m_UpdateShapeMesh;
		typename MeshType::Pointer m_InputMesh;
		typename MeshType::Pointer m_OutputMesh;
		std::map<int, int> m_ShapeClusterMap; //<pointId, clusterId>
		std::map<int, ShapeClusterItem*> m_ShapeClusterInstances; //<clusterId, clusterItem>
		StatismoMatrixType m_ClusterWeights; //row-col: pointId-clusterWeight
		StatismoVectorType m_ShapeVarianceWeights;
		unsigned m_NumberOfCoefficients;
		std::map<int, LandmarkStatus> m_LandmarkStatus;
		unsigned m_Iterations;
		int m_NumberOfClusters;

		void Cluster(int numberOfClusters);

		void RigidTransformFitting(const MeshType* targetMesh);
		
		void ShapeTransformFitting(const MeshType* targetMesh);
		
		void ClusteredShapeTransformFitting(const MeshType* targetMesh);
		
		void UpdateRigidTransform(const StatismoMatrixType & mat, RigidTransformType* rigidTransform);

		void UpdateShapeTransform(const StatismoVectorType & vec, ShapeTransformType* shapeTransform);

		void AddCluster(int pointId, int clusterId);

		void ClearClusters();

		void AllocateClusters();

		void CalClusterWeights();

		void UpdateShape();

		void UpdateRigid();
	};
}

#include "kmSSMUtils.hxx"

#endif