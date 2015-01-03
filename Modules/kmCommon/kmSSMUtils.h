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
			PointType clusterCentroid;
			double shapeProbability;
			bool enabled;
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
				clusterCentroid.Fill(0.0);
				shapeProbability = 0.0;
				enabled = false;

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
		};
		
		SSMUtils()
		{
			m_ReferenceShapeMesh = NULL;
			m_MeanShapeMesh = NULL;
			m_NumberOfCoefficients = 3;
			m_Iterations = 0;
		}

		~SSMUtils()
		{

		}
		
		void SetSSM(const StatisticalModelType* model)
		{
			m_SSM = const_cast<StatisticalModelType*>( model );
			m_ReferenceShapeMesh = m_SSM->GetRepresenter()->GetReference();
			m_MeanShapeMesh = m_SSM->DrawMean();
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
		
		void initialize();
		
		void update(const MeshType * targetMesh, MeshType * outputMesh);
		
		void printTransform();
		
		double calShapeParaDiff();
		
		double getShapeProbability(int pointId);
		
		void cluster(int numberOfClusters);
		
		ShapeClusterItem* getClusterByClusterId(int clusterId)
		{
			ShapeClusterItem* clusterItem = m_ShapeClusterInstances[clusterId];
			if (clusterItem==NULL){
				std::cerr<<"Cannot find cluster Id: "<<clusterId<<std::endl;
				clusterItem = new ShapeClusterItem(0);
			}

			return clusterItem;
		}
		
		ShapeClusterItem* getClusterByPointId(int pointId)
		{
			int clusterId = m_ShapeClusterMap[pointId];
			return getClusterByClusterId(clusterId);
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
		std::map<int, int> m_ShapeClusterMap; //<pointId, clusterId>
		std::map<int, ShapeClusterItem*> m_ShapeClusterInstances; //<clusterId, clusterItem>
		StatismoMatrixType m_ClusterWeights; //row-col: pointId-clusterWeight
		StatismoVectorType m_ShapeVarianceWeights;
		int m_NumberOfCoefficients;
		std::map<int, double> m_Confidences; //<pointId, shapeProbability>
		unsigned m_Iterations;

		void rigidTransformFitting(const MeshType* targetMesh);
		
		void clusteredRigidTransformFitting(const MeshType* targetMesh);
		
		void shapeTransformFitting(const MeshType* targetMesh);
		
		void clusteredShapeTransformFitting(const MeshType* targetMesh);
		
		void updateRigidTransform(const StatismoMatrixType & mat, RigidTransformType* rigidTransform);

		//return shape probability
		double updateShapeTransform(const StatismoVectorType & vec, ShapeTransformType* shapeTransform);

		void addCluster(int pointId, int clusterId);

		void clearClusters();

		void allocateClusters();

		void calculateClusterWeights();

		void updateShape(bool clustered = true);

		void updateRigid(MeshType * outputMesh, bool clustered = true);
	};
}

#include "kmSSMUtils.hxx"

#endif