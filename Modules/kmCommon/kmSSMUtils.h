#ifndef __kmSSMUtils_h
#define __kmSSMUtils_h

#include "statismo/CommonTypes.h"
//#include "kmUtility.h"
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
		itkStaticConstMacro(PointDimension, unsigned int, MeshType::PointDimension);
		
		SSMUtils()
		{
			m_ReferenceShapeMesh = NULL;
		}
		
		~SSMUtils()
		{
		}
		
		void SetSSM(const StatisticalModelType* model)
		{
			m_SSM = const_cast<StatisticalModelType*>( model );
			m_ReferenceShapeMesh = m_SSM->GetRepresenter()->GetReference();
		}
		void SetRigidTransform(RigidTransformType* rtfm)
		{
			m_RigidTranform = rtfm;
			m_RigidParametersPre = m_RigidTranform->GetParameters();
		}
		void SetShapeTransform(ShapeTransformType* stfm)
		{
			m_ShapeTransform = stfm;
			m_ShapeParametersPre = m_ShapeTransform->GetParameters();
		}
		
		void compositeTransformFitting(const MeshType* targetMesh)
		{
			if (m_SSM == NULL)
			{
				std::cout<<"SSM cannot be null."<<std::endl;
				return;
			}

			if (m_ShapeTransform == NULL)
			{
				m_ShapeTransform = ShapeTransformType::New();
				m_ShapeTransform->SetStatisticalModel( m_SSM );
				m_ShapeTransform->SetIdentity();
			}

			if (m_RigidTranform == NULL)
			{
				m_RigidTranform = RigidTransformType::New();
				m_RigidTranform->SetIdentity();
			}

			if (m_ShapeTransform->GetUsedNumberOfCoefficients() > m_ShapeTransform->GetNumberOfParameters())
			{
				m_ShapeTransform->SetUsedNumberOfCoefficients(m_ShapeTransform->GetNumberOfParameters());
			}
			
			unsigned int numberOfLandmarks = m_ReferenceShapeMesh->GetNumberOfPoints();
			{
				/*****************Fit rigid parameters****************/
				MeshType::ConstPointer sourceMeshForRigid = km::transformMesh<MeshType, ShapeTransformType>(m_ReferenceShapeMesh, m_ShapeTransform);
				MeshType::ConstPointer targetMeshForRigid = targetMesh;

				StatismoMatrixType sourceMatrixForRigid, targetMatrixForRigid;
				fillPointSetIntoStatismoMatrix(sourceMeshForRigid, sourceMatrixForRigid);
				fillPointSetIntoStatismoMatrix(targetMeshForRigid, targetMatrixForRigid);

				StatismoVectorType I = StatismoVectorType::Ones(sourceMatrixForRigid.cols());
				StatismoMatrixType Mmatrix = sourceMatrixForRigid.transpose() * sourceMatrixForRigid;
				Mmatrix.diagonal() += I;
				StatismoMatrixType MInverseMatrix = Mmatrix.inverse();
				const StatismoMatrixType& WT = sourceMatrixForRigid.transpose();
				StatismoMatrixType coeffsRigid = MInverseMatrix * (WT * targetMatrixForRigid);

				UpdateRigidTransform(coeffsRigid.transpose());
			}

			{
				/*****************Fit shape parameters****************/
				RigidTransformType::Pointer inversedRigidTransform = RigidTransformType::New();
				m_RigidTranform->GetInverse( inversedRigidTransform );
				MeshType::Pointer targetMeshForShape = km::transformMesh<MeshType, RigidTransformType>(targetMesh, inversedRigidTransform);
				MeshType::Pointer sourceMeshForShape = this->m_SSM->DrawMean();

				StatismoVectorType sourceMatrixForShape, targetMatrixForShape;
				const StatismoMatrixType & basisMatrixForShape = this->m_SSM->GetstatismoImplObj()->GetPCABasisMatrix();
				fillPointSetIntoStatismoVector(sourceMeshForShape, sourceMatrixForShape);
				fillPointSetIntoStatismoVector(targetMeshForShape, targetMatrixForShape);

				StatismoVectorType I = StatismoVectorType::Ones(basisMatrixForShape.cols());
				StatismoMatrixType Mmatrix = basisMatrixForShape.transpose() * basisMatrixForShape;
				Mmatrix.diagonal() += I;
				StatismoMatrixType MInverseMatrix = Mmatrix.inverse();
				const StatismoMatrixType& WT = basisMatrixForShape.transpose();
				StatismoVectorType coeffsShape = MInverseMatrix * (WT * (targetMatrixForShape-sourceMatrixForShape));

				UpdateShapeTransform(coeffsShape);
			}
		}

		double calShapeParaDiff()
		{
			double shapeParamDiff = 0;
			for (int p=0;p<m_ShapeTransform->GetUsedNumberOfCoefficients();p++)
			{
				shapeParamDiff += std::abs(m_ShapeParametersPost[p]-m_ShapeParametersPre[p]);
			}
			shapeParamDiff /= m_ShapeTransform->GetUsedNumberOfCoefficients();
			return shapeParamDiff;
		}

		void initModelVariance()
		{
			//typedef itk::SimplexMeshGeometry SimplexMeshGeometryType;
			//typedef SimplexMeshGeometryType::CovariantVectorType CovariantVectorType;

			//SimplexMeshGeometryType *geodata;

			//StatisticalModelType::MatrixType basisMatrix = model->GetPCABasisMatrix();
			//for (int id=0;id<mesh->GetNumberOfPoints();id++)
			//{
			//	double variance = 0.0;
			//	geodata = mesh->GetGeometryData()->GetElement(id);
			//	for (unsigned p=0;p<numberOfComponents;p++)
			//	{
			//		double varOnNormal = 0.0;
			//		for (unsigned d=0; d<Dimension; d++) 
			//		{
			//			unsigned idx = model->GetRepresenter()->MapPointIdToInternalIdx(id, d);
			//			varOnNormal += geodata->normal[d] * basisMatrix[idx][p];
			//		}
			//		variance += std::abs(varOnNormal);
			//	}
			//	g_varianceMap[id] = variance;
			//}
		}
	
	private:
		const typename StatisticalModelType* m_SSM;
		typename RigidTransformType* m_RigidTranform;
		typename ShapeTransformType* m_ShapeTransform;
		
		typename RigidParametersType m_RigidParametersPre;
		typename RigidParametersType m_RigidParametersPost;
		typename ShapeParametersType m_ShapeParametersPre;
		typename ShapeParametersType m_ShapeParametersPost;
		
		typename MeshType* m_ReferenceShapeMesh;
		
		void fillPointSetIntoStatismoMatrix(const MeshType* mesh, StatismoMatrixType& matrix)
		{
			unsigned int numberOfLandmarks = mesh->GetNumberOfPoints();
			matrix.setConstant(numberOfLandmarks, PointDimension+1, 1.0);
			unsigned long rowcnt = 0;
			for (int idx=0;idx<numberOfLandmarks;idx++)
			{
				PointType pt = mesh->GetPoint(idx);
				for (int d=0;d<PointDimension;d++)
				{
					matrix(rowcnt, d) = pt[d];
				}
				rowcnt++;
			}
		}
		
		void fillPointSetIntoStatismoVector(const MeshType* mesh, StatismoVectorType& matrix)
		{
			unsigned int numberOfLandmarks = mesh->GetNumberOfPoints();
			matrix.setZero( numberOfLandmarks*PointDimension, 1);
			unsigned long rowcnt = 0;
			for (int idx=0;idx<numberOfLandmarks;idx++)
			{
				PointType pt = mesh->GetPoint(idx);
				for (int d=0;d<PointDimension;d++)
				{
					matrix(rowcnt, 0) = pt[d];
					rowcnt++;
				}
			}
		}
		
		void UpdateRigidTransform(const StatismoMatrixType & mat)
		{
			m_RigidParametersPre = m_RigidTranform->GetParameters();
			RigidTransformType::MatrixType rigidMatrix = m_RigidTranform->GetMatrix();
			RigidTransformType::OffsetType rigidOffset = m_RigidTranform->GetOffset();
			for (int i=0;i<3;i++)
			{
				for (int j=0;j<3;j++)
				{
					rigidMatrix[i][j] = mat(i, j);
				}
			}
			for (int i=0;i<Dimension;i++)
			{
				rigidOffset[i] = mat(i, Dimension);
			}
			m_RigidTranform->SetMatrix( rigidMatrix );
			m_RigidTranform->SetOffset( rigidOffset );
			m_RigidParametersPost = m_RigidTranform->GetParameters();
		}
		
		void UpdateShapeTransform(const StatismoVectorType & vec)
		{
			m_ShapeParametersPre  = m_ShapeTransform->GetParameters();
			m_ShapeParametersPost = m_ShapeTransform->GetParameters();
			for (int i=0;i<m_ShapeTransform->GetUsedNumberOfCoefficients();i++)
			{
				if (vec[i]<-3.0){
					m_ShapeParametersPost[i] = -3;
				}else if (vec[i]>3.0){
					m_ShapeParametersPost[i] = 3;
				}else{
					m_ShapeParametersPost[i] = vec[i];	
				}
			}
			m_ShapeTransform->SetParameters(m_ShapeParametersPost);
		}
	};
}

#endif