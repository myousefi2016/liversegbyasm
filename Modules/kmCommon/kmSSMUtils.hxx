#ifndef __kmSSMUtils_hxx
#define __kmSSMUtils_hxx

#include "kmSSMUtils.h"

namespace km
{	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::Initialize()
	{
		if (m_SSM == NULL){
			std::cout<<"SSM cannot be null."<<std::endl;
			return;
		}

		if (m_ReferenceShapeMesh.IsNull())
		{
			m_ReferenceShapeMesh = m_SSM->GetRepresenter()->GetReference();
		}

		if (m_MeanShapeMesh.IsNull())
		{
			m_MeanShapeMesh = m_SSM->DrawMean();
		}

		if (m_RigidTranform.IsNull()){
			std::cout<<"Rigid transform cannot be null."<<std::endl;
			return;
		}

		if (m_ShapeTransform.IsNull()){
			std::cout<<"Shape transform cannot be null."<<std::endl;
			return;
		}else{
			m_ShapeTransform->SetStatisticalModel(m_SSM);
		}

		if (m_UpdateShapeMesh.IsNull()){
			m_UpdateShapeMesh = km::transformMesh<MeshType, ShapeTransformType>(m_ReferenceShapeMesh, m_ShapeTransform);
		}

		if (m_InputMesh.IsNull())
		{
			m_InputMesh = km::cloneMesh<MeshType, MeshType>(m_MeanShapeMesh);
		}

		if (m_OutputMesh.IsNull())
		{
			m_OutputMesh = km::cloneMesh<MeshType, MeshType>(m_MeanShapeMesh);
		}

		m_ShapeTransform->SetUsedNumberOfCoefficients(m_NumberOfCoefficients);
		m_ShapeParametersPre = m_ShapeTransform->GetParameters();

		StatisticalModelType::VectorType pcaVariance = m_SSM->GetPCAVarianceVector();
		double totalVariances = 0.0;
		for (int i=0;i<m_NumberOfCoefficients;i++){
			totalVariances += pcaVariance[i];
		}

		m_ShapeVarianceWeights.setZero(m_SSM->GetNumberOfPrincipalComponents(), 1);
		for (int i=0;i<m_SSM->GetNumberOfPrincipalComponents();i++){
			m_ShapeVarianceWeights[i] = pcaVariance[i]/totalVariances;
		}

		this->Cluster(m_NumberOfClusters);

		for (int i=0;i<m_MeanShapeMesh->GetNumberOfPoints();i++)
		{
			this->SetLandmarkStatus(i, Normal);
		}
	}
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::Update(const MeshType * targetMesh, MeshType * outputMesh)
	{
		km::copyMeshToMeshPoints<MeshType, MeshType>(targetMesh, m_InputMesh);

		for (int i=0;i<1;i++)
		{
			std::cout<<"-----------------------------------------"<<std::endl;
			std::cout<<"Start rigid transform fitting.."<<std::endl;
			this->RigidTransformFitting(m_InputMesh);

			std::cout<<"Start shape transform fitting.."<<std::endl;
			this->ShapeTransformFitting(m_InputMesh);

			std::cout<<"Start clustered shape transform fitting.."<<std::endl;
			this->ClusteredShapeTransformFitting(m_InputMesh);
		}

		std::cout<<"Apply shape transform.."<<std::endl;
		this->UpdateShape();

		std::cout<<"Apply rigid transform.."<<std::endl;
		this->UpdateRigid();

		km::copyMeshToMeshPoints<MeshType, MeshType>(m_OutputMesh, outputMesh);

		m_Iterations++;
	}
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::PrintTransform()
	{
		std::cout<<"rigid parameters: "<<m_RigidTranform->GetParameters()<<std::endl;
		for(std::map<int, ShapeClusterItem*>::iterator it=m_ShapeClusterInstances.begin();it!=m_ShapeClusterInstances.end();it++)
		{
			ShapeClusterItem* clusterItem = it->second;
			std::cout<<"Cluster "<<clusterItem->clusterId<<", shape probability: "<<clusterItem->shapeProbability<<std::endl;
			std::cout<<"Cluster "<<clusterItem->clusterId<<", shape parameters: "<<clusterItem->shapeTransform->GetParameters()<<std::endl;
		}
	}
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	double
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::CalShapeParaDiff()
	{
		double shapeParamDiff = 0;
		for(std::map<int, ShapeClusterItem*>::iterator it=m_ShapeClusterInstances.begin();it!=m_ShapeClusterInstances.end();it++)
		{
			ShapeClusterItem* clusterItem = it->second;

			ShapeTransformType* shapeTransform = clusterItem->shapeTransform;
			ShapeParametersType shapeParametersPre = clusterItem->shapeParametersPre;
			ShapeParametersType shapeParametersPost = shapeTransform->GetParameters();
			double diff = 0.0;
			for (int p=0;p<m_NumberOfCoefficients;p++)
			{
				diff += std::abs(shapeParametersPost[p]-shapeParametersPre[p])*m_ShapeVarianceWeights[p];
			}
			shapeParamDiff += diff;
		}
		
		shapeParamDiff /= m_ShapeClusterInstances.size();
		return shapeParamDiff;
	}


	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::Cluster(int numberOfClusters)
	{
		timer->StartTimer();

		this->ClearClusters();

		int numberOfModesUsed = km::g_number_principle_components;//this->m_NumberOfCoefficients;
		double multifier = 3.0;
		
		const StatismoMatrixType & basisMatrixForShape = this->m_SSM->GetstatismoImplObj()->GetPCABasisMatrix();
		const unsigned int measureLength = numberOfModesUsed*Dimension;
		const unsigned int measureCount = basisMatrixForShape.rows()/Dimension;

		std::cout<<"number of shape coefficients: "<<m_NumberOfCoefficients<<std::endl;
		std::cout<<"number of modes used: "<<numberOfModesUsed<<std::endl;
		std::cout<<"measure length: "<<measureLength<<", count: "<<measureCount<<std::endl;

		if (numberOfClusters <= 1)
		{
			for (int pointId=0;pointId<m_ReferenceShapeMesh->GetNumberOfPoints();pointId++){
				this->AddCluster(pointId, 0);
			}
		}
		else
		{
			typedef itk::VariableLengthVector< double> MeasurementVectorType;
			typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
			SampleType::Pointer sample = SampleType::New();
			sample->SetMeasurementVectorSize( measureLength );
			MeasurementVectorType mv;
			mv.SetSize( measureLength );
			mv.Fill( 0 );

			std::ofstream clusterSampleFile;
			clusterSampleFile.open ("clusterSamples.txt", std::ofstream::out | std::ofstream::trunc);
			for (int i=0;i<measureCount;i++)
			{
				int offset = 0;
				for(int j=0;j<measureLength/Dimension;j++)
				{
					double shapeVariance = 0.0;
					for(int d=0;d<Dimension;d++)
					{
						//shapeVariance += std::pow(basisMatrixForShape(i*Dimension+d), 2);
						mv[offset++] = multifier*basisMatrixForShape(i*Dimension+d, j);
					}
					//mv[offset++] = shapeVariance;
				}
				sample->PushBack( mv );
				//std::cout<<mv<<std::endl;
				clusterSampleFile << mv <<std::endl;
			}
			clusterSampleFile.close();

			typedef itk::Statistics::WeightedCentroidKdTreeGenerator< SampleType > TreeGeneratorType;
			TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
			treeGenerator->SetSample( sample );
			treeGenerator->SetBucketSize( 16 );
			treeGenerator->Update();

			//static boost::minstd_rand randgen(static_cast<unsigned>(time(0)));
			//static boost::normal_distribution<> dist(0, 1);
			//static boost::variate_generator<boost::minstd_rand, boost::normal_distribution<> > r(randgen, dist);

			typedef TreeGeneratorType::KdTreeType TreeType;
			typedef itk::Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;
			EstimatorType::Pointer estimator = EstimatorType::New();
			EstimatorType::ParametersType initialMeans( measureLength * m_NumberOfClusters );
			initialMeans.Fill( 1.0 );
			estimator->SetParameters( initialMeans );
			estimator->SetKdTree( treeGenerator->GetOutput() );
			estimator->SetMaximumIteration( 200 );
			estimator->SetCentroidPositionChangesThreshold(0.01);
			estimator->StartOptimization();

			EstimatorType::ParametersType estimatedMeans = estimator->GetParameters();
			typedef itk::Statistics::DistanceToCentroidMembershipFunction< MeasurementVectorType > MembershipFunctionType;
			typedef itk::Statistics::MinimumDecisionRule DecisionRuleType;
			DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();
			typedef itk::Statistics::SampleClassifierFilter< SampleType > ClassifierType;
			ClassifierType::Pointer classifier = ClassifierType::New();
			classifier->SetDecisionRule( decisionRule );
			classifier->SetInput( sample );
			classifier->SetNumberOfClasses( numberOfClusters );

			typedef ClassifierType::ClassLabelVectorObjectType
				ClassLabelVectorObjectType;
			typedef ClassifierType::ClassLabelVectorType ClassLabelVectorType;
			typedef ClassifierType::ClassLabelType ClassLabelType;
			ClassLabelVectorObjectType::Pointer classLabelsObject =
				ClassLabelVectorObjectType::New();
			ClassLabelVectorType& classLabelsVector = classLabelsObject->Get();
			for (int k=0;k<m_NumberOfClusters;k++)
			{
				classLabelsVector.push_back( k );
			}
			classifier->SetClassLabels( classLabelsObject );

			typedef ClassifierType::MembershipFunctionVectorObjectType MembershipFunctionVectorObjectType;
			typedef ClassifierType::MembershipFunctionVectorType MembershipFunctionVectorType;
			MembershipFunctionVectorObjectType::Pointer membershipFunctionVectorObject =
				MembershipFunctionVectorObjectType::New();
			MembershipFunctionVectorType& membershipFunctionVector =
				membershipFunctionVectorObject->Get();

			int index = 0;
			for ( unsigned int i = 0 ; i < m_NumberOfClusters ; i++ ){
				MembershipFunctionType::Pointer membershipFunction = MembershipFunctionType::New();
				MembershipFunctionType::CentroidType centroid( sample->GetMeasurementVectorSize() );
				for ( unsigned int j = 0 ; j < sample->GetMeasurementVectorSize(); j++ ){
					centroid[j] = estimatedMeans[index++];
				}
				membershipFunction->SetCentroid( centroid );
				membershipFunctionVector.push_back( membershipFunction.GetPointer() );
			}
			classifier->SetMembershipFunctions( membershipFunctionVectorObject );
			classifier->Update();

			const ClassifierType::MembershipSampleType* membershipSample = classifier->GetOutput();
			ClassifierType::MembershipSampleType::ConstIterator labelIt = membershipSample->Begin();
			ClassifierType::MembershipSampleType::ConstIterator labelItEnd = membershipSample->End();
			int pointId = 0;
			while(labelIt!=labelItEnd)
			{
				int clusterId = static_cast<int>(labelIt.GetClassLabel());
				this->AddCluster(pointId, clusterId);

				++pointId;
				++labelIt;
			}
		}

		timer->StopTimer();
		std::cout<<"Cluster done. "<<timer->GetElapsedTime() << " s." <<endl;
		std::cout<<"Actual cluster number: "<<this->m_ShapeClusterInstances.size()<<std::endl;

		this->AllocateClusters();

		this->CalClusterWeights();
	}
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::RigidTransformFitting(const MeshType* targetMesh)
	{
		//timer->StartTimer();

		unsigned int numberOfLandmarks = m_ReferenceShapeMesh->GetNumberOfPoints();
		unsigned int numberOfLandmarksNormal = 0;
		for (int i=0;i<numberOfLandmarks;i++)
		{
			if (this->GetLandmarkStatus(i) == Normal) numberOfLandmarksNormal++;
		}

		/*****************Fit rigid parameters****************/
		MeshType::ConstPointer sourceMeshForRigid = m_UpdateShapeMesh;//km::transformMesh<MeshType, ShapeTransformType>(m_ReferenceShapeMesh, m_ShapeTransform);;
		MeshType::ConstPointer targetMeshForRigid = targetMesh;

		StatismoMatrixType sourceMatrixForRigid = StatismoMatrixType::Ones(numberOfLandmarksNormal, Dimension+1);
		StatismoMatrixType targetMatrixForRigid = StatismoMatrixType::Ones(numberOfLandmarksNormal, Dimension+1);
		//StatismoMatrixType sourceMatrixForRigid = StatismoMatrixType::Ones(numberOfLandmarks, Dimension+1);
		//StatismoMatrixType targetMatrixForRigid = StatismoMatrixType::Ones(numberOfLandmarks, Dimension+1);
		unsigned long rowcnt = 0;
		for (int idx=0;idx<numberOfLandmarks;idx++)
		{
			if (this->GetLandmarkStatus(idx) == Normal)
			{
				PointType sourcePt = sourceMeshForRigid->GetPoint(idx);
				PointType targetPt = targetMeshForRigid->GetPoint(idx);
				for (int d=0;d<Dimension;d++)
				{
					sourceMatrixForRigid(rowcnt, d) = sourcePt[d];
					targetMatrixForRigid(rowcnt, d) = targetPt[d];
				}
				rowcnt++;
			}


			//PointType sourcePt = sourceMeshForRigid->GetPoint(idx);
			//PointType targetPt = targetMeshForRigid->GetPoint(idx);
			//if (!m_Flags[idx])
			//{
			//	targetPt = m_RigidTranform->TransformPoint(sourcePt);
			//}
			//for (int d=0;d<Dimension;d++)
			//{
			//	sourceMatrixForRigid(rowcnt, d) = sourcePt[d];
			//	targetMatrixForRigid(rowcnt, d) = targetPt[d];
			//}
			//rowcnt++;
		}

		StatismoVectorType I = StatismoVectorType::Ones(sourceMatrixForRigid.cols());
		StatismoMatrixType Mmatrix = sourceMatrixForRigid.transpose() * sourceMatrixForRigid;
		Mmatrix.diagonal() += I;
		StatismoMatrixType MInverseMatrix = Mmatrix.inverse();
		const StatismoMatrixType& WT = sourceMatrixForRigid.transpose();
		StatismoMatrixType coeffsRigid = MInverseMatrix * (WT * targetMatrixForRigid);

		UpdateRigidTransform(coeffsRigid.transpose(), this->m_RigidTranform);

		//timer->StopTimer();
		//std::cout<<"Rigid transform fitting done. "<<timer->GetElapsedTime() << " s." <<endl;
	}
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::ShapeTransformFitting(const MeshType* targetMesh)
	{
		//timer->StartTimer();

		unsigned int numberOfLandmarks = m_ReferenceShapeMesh->GetNumberOfPoints();
		unsigned int numberOfLandmarksNormal = 0;
		for (int i=0;i<numberOfLandmarks;i++)
		{
			if (this->GetLandmarkStatus(i) == Normal) numberOfLandmarksNormal++;
		}

		/*****************Fit shape parameters****************/
		RigidTransformType::Pointer inversedRigidTransform = RigidTransformType::New();
		m_RigidTranform->GetInverse( inversedRigidTransform );
		MeshType::Pointer targetMeshForShape = km::transformMesh<MeshType, RigidTransformType>(targetMesh, inversedRigidTransform);
		MeshType::Pointer sourceMeshForShape = this->m_MeanShapeMesh;
		const StatismoMatrixType & originalBasisMatrixForShape = this->m_SSM->GetstatismoImplObj()->GetPCABasisMatrix();

		StatismoMatrixType basisMatrixForShape = StatismoMatrixType::Zero(numberOfLandmarksNormal*Dimension, m_NumberOfCoefficients);
		StatismoVectorType sourceMatrixForShape = StatismoVectorType::Zero(basisMatrixForShape.rows());
		StatismoVectorType targetMatrixForShape = StatismoVectorType::Zero(basisMatrixForShape.rows());
		
		unsigned long rowcnt = 0;
		for (int idx=0;idx<numberOfLandmarks;idx++)
		{
			if (this->GetLandmarkStatus(idx) == Normal)
			{
				PointType sourcePt = sourceMeshForShape->GetPoint(idx);
				PointType targetPt = targetMeshForShape->GetPoint(idx);
				for (int d=0;d<Dimension;d++)
				{
					sourceMatrixForShape[rowcnt] = sourcePt[d];
					targetMatrixForShape[rowcnt] = targetPt[d];
					for (int c=0;c<basisMatrixForShape.cols();c++){
						basisMatrixForShape(rowcnt, c) = originalBasisMatrixForShape(idx*Dimension+d, c);
					}
					rowcnt++;
				}
			}

			//PointType sourcePt = sourceMeshForShape->GetPoint(idx);
			//PointType targetPt = targetMeshForShape->GetPoint(idx);
			//if (!m_Flags[idx])
			//{
			//	targetPt = m_ShapeTransform->TransformPoint(m_ReferenceShapeMesh->GetPoint(idx));
			//}
			//for (int d=0;d<Dimension;d++)
			//{
			//	sourceMatrixForShape[rowcnt] = sourcePt[d];
			//	targetMatrixForShape[rowcnt] = targetPt[d];
			//	//for (int c=0;c<basisMatrixForShape.cols();c++){
			//	//	basisMatrixForShape(rowcnt, c) = originalBasisMatrixForShape(idx*Dimension+d, c);
			//	//}
			//	rowcnt++;
			//}
		}

		StatismoVectorType I = StatismoVectorType::Ones(basisMatrixForShape.cols());
		StatismoMatrixType Mmatrix = basisMatrixForShape.transpose() * basisMatrixForShape;
		Mmatrix.diagonal() += I;
		StatismoMatrixType MInverseMatrix = Mmatrix.inverse();
		const StatismoMatrixType& WT = basisMatrixForShape.transpose();
		StatismoVectorType coeffsShape = MInverseMatrix * (WT * (targetMatrixForShape-sourceMatrixForShape));

		this->m_ShapeParametersPre = m_ShapeTransform->GetParameters();
		UpdateShapeTransform(coeffsShape, m_ShapeTransform);

		//timer->StopTimer();
		//std::cout<<"Shape transform fitting done. "<<timer->GetElapsedTime() << " s." <<endl;
	}
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::ClusteredShapeTransformFitting(const MeshType* targetMesh)
	{
		//timer->StartTimer();

		/*****************Fit shape parameters****************/
		RigidTransformType::Pointer inversedRigidTransform = RigidTransformType::New();
		m_RigidTranform->GetInverse( inversedRigidTransform );
		MeshType::Pointer targetMeshForShape = km::transformMesh<MeshType, RigidTransformType>(targetMesh, inversedRigidTransform);
		MeshType::Pointer sourceMeshForShape = this->m_MeanShapeMesh;

		for(std::map<int, ShapeClusterItem*>::iterator it=m_ShapeClusterInstances.begin();it!=m_ShapeClusterInstances.end();it++)
		{
			ShapeClusterItem* clusterItem = it->second;
			clusterItem->shapeParametersPre = clusterItem->shapeTransform->GetParameters();

			if(false)
			{
				for (int idx=0;idx<clusterItem->pointIds.size();idx++)
				{
					int pointId = clusterItem->pointIds[idx];
					PointType targetPt = targetMeshForShape->GetPoint(pointId);
					PointType sourcePt = sourceMeshForShape->GetPoint(pointId);
					for (int d=0;d<Dimension;d++)
					{
						clusterItem->targetMatrixForShape[idx*Dimension+d] = targetPt[d];
						clusterItem->sourceMatrixForShape[idx*Dimension+d] = sourcePt[d];
					}
				}
				StatismoVectorType coeffs = clusterItem->MInverseMatrix * (clusterItem->WT * (clusterItem->targetMatrixForShape-clusterItem->sourceMatrixForShape));
				this->UpdateShapeTransform(coeffs, clusterItem->shapeTransform);
			}
			else
			{
				clusterItem->shapeTransform->SetParameters(m_ShapeTransform->GetParameters());
			}
		}

		//timer->StopTimer();
		//std::cout<<"Clustered shape transform fitting done. "<<timer->GetElapsedTime() << " s." <<endl;
	}	
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::UpdateRigidTransform(const StatismoMatrixType & mat, RigidTransformType* rigidTransform)
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
		for (int i=0;i<Dimension;i++)
		{
			rigidOffset[i] = mat(i, Dimension);
		}
		rigidTransform->SetMatrix( rigidMatrix );
		rigidTransform->SetOffset( rigidOffset );
	}	
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::UpdateShapeTransform(const StatismoVectorType & vec, ShapeTransformType* shapeTransform)
	{
		ShapeParametersType shapeParameters = shapeTransform->GetParameters();
		for (int i=0;i<m_NumberOfCoefficients;i++)
		{
			shapeParameters[i] = km::Math::setToBetween(vec[i], -2.5, 2.5);
		}
		shapeTransform->SetParameters(shapeParameters);
	}
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::AddCluster(int pointId, int clusterId)
	{
		this->m_ShapeClusterMap[pointId] = clusterId;
		ShapeClusterItem* item = this->m_ShapeClusterInstances[clusterId];
		if (item == NULL)
		{
			item = new ShapeClusterItem(clusterId);
			item->clusterId = clusterId;
			this->m_ShapeClusterInstances[clusterId] = item;
		}
		item->pointIds.push_back(pointId);
	}
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::ClearClusters()
	{
		for(std::map<int, ShapeClusterItem*>::iterator it=m_ShapeClusterInstances.begin();it!=m_ShapeClusterInstances.end();it++)
		{
			ShapeClusterItem* item = it->second;
			delete item;
			item = NULL;
		}
		m_ShapeClusterInstances.clear();
		m_ShapeClusterMap.clear();
	}
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::AllocateClusters()
	{
		timer->StartTimer();

		m_MeanShapeMesh->GetPointData()->Reserve(m_MeanShapeMesh->GetNumberOfPoints());
		for(std::map<int, ShapeClusterItem*>::iterator it=m_ShapeClusterInstances.begin();it!=m_ShapeClusterInstances.end();it++)
		{
			ShapeClusterItem* item = it->second;
			if (item!=NULL)
			{
				item->shapeTransform = ShapeTransformType::New();
				item->shapeTransform->SetStatisticalModel( m_SSM );
				item->shapeTransform->SetUsedNumberOfCoefficients( m_NumberOfCoefficients );
				item->shapeTransform->SetParameters(m_ShapeTransform->GetParameters());
				item->shapeParametersPre = item->shapeTransform->GetParameters();

				item->sourceMatrixForShape.setZero( item->pointIds.size()*Dimension, 1);
				item->targetMatrixForShape.setZero( item->pointIds.size()*Dimension, 1);
				item->basisMatrix.setZero( item->pointIds.size()*Dimension, m_NumberOfCoefficients);
				const StatismoMatrixType & origBasisMatrix = this->m_SSM->GetstatismoImplObj()->GetPCABasisMatrix();

				item->rigidTransform = RigidTransformType::New();
				item->rigidTransform->SetParameters(m_RigidTranform->GetParameters());
				item->rigidTransform->SetFixedParameters(m_RigidTranform->GetFixedParameters());

				item->sourceMatrixForRigid.setConstant( item->pointIds.size(), Dimension+1, 1.0 );
				item->targetMatrixForRigid.setConstant( item->pointIds.size(), Dimension+1, 1.0 );

				for(int idx=0; idx < item->pointIds.size(); idx++)
				{
					int pointId = item->pointIds[idx];
					PointType pt = m_MeanShapeMesh->GetPoint(pointId);
					for (int d=0;d<Dimension;d++)
					{
						for (int c=0;c<item->basisMatrix.cols();c++){
							item->basisMatrix(idx*Dimension+d, c) = origBasisMatrix(pointId*Dimension+d, c);
						}
					}
					item->meanPoints->PushBack(pt);
					m_MeanShapeMesh->SetPointData(pointId, item->clusterId);
				}

				item->kdTreeGenerator->SetSample( item->meanPoints );
				item->kdTreeGenerator->Update();

				StatismoVectorType I = StatismoVectorType::Ones(item->basisMatrix.cols());
				StatismoMatrixType Mmatrix = item->basisMatrix.transpose() * item->basisMatrix;
				Mmatrix.diagonal() += I;
				item->MInverseMatrix = Mmatrix.inverse();
				item->WT = item->basisMatrix.transpose();

				//std::cout<<"**** Cluster "<<item->clusterId<<": "<<std::endl;
				//std::cout<<"       Number of points: "<<item->pointIds.size()<<std::endl;
			}
			else
			{
				std::cerr<<"Null cluster item occurs!!!!"<<std::endl;
			}
		}

		timer->StopTimer();
		std::cout<<"Allocate clusters done. "<<timer->GetElapsedTime() << " s." <<endl;
		km::writeMesh<MeshType>("ClusterLabels.vtk", m_MeanShapeMesh);
	}

	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::CalClusterWeights()
	{
		timer->StartTimer();

		m_MeanShapeMesh->GetPointData()->Reserve(m_MeanShapeMesh->GetNumberOfPoints());
		m_ClusterWeights.setZero(m_MeanShapeMesh->GetNumberOfPoints(), m_ShapeClusterInstances.size());
		for (int id=0;id<m_MeanShapeMesh->GetNumberOfPoints();id++)
		{
			PointType point = m_MeanShapeMesh->GetPoint(id);
			double weightSum = 0.0;
			for(std::map<int, ShapeClusterItem*>::iterator it=m_ShapeClusterInstances.begin();it!=m_ShapeClusterInstances.end();it++)
			{
				int clusterId = it->first;

				ShapeClusterItem* item = it->second;
				PointType closestPt = item->findClosestPoint(point);

				double dist = point.EuclideanDistanceTo(closestPt);
				double weight = 1.0/(km::g_cluster_min_dist+dist);

				//weight = std::sqrt(weight);
				
				m_ClusterWeights(id, clusterId) = weight;

				weightSum += weight;
			}

			for(std::map<int, ShapeClusterItem*>::iterator it=m_ShapeClusterInstances.begin();it!=m_ShapeClusterInstances.end();it++)
			{
				int clusterId = it->first;
				m_ClusterWeights(id, clusterId) /= weightSum;
			}

			m_MeanShapeMesh->SetPointData(id, m_ClusterWeights(id, m_ShapeClusterMap[id]));
		}

		timer->StopTimer();
		std::cout<<"Calculate cluster weights done. "<<timer->GetElapsedTime() << " s." <<endl;
		km::writeMesh<MeshType>("ClusterWeights.vtk", m_MeanShapeMesh);
	}
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::UpdateShape()
	{
		//Store updated shape mesh into m_UpdateShapeMesh.
		ShapeTransformType::Pointer clusteredShapeTransform = ShapeTransformType::New();
		clusteredShapeTransform->SetStatisticalModel(m_SSM);
		clusteredShapeTransform->SetUsedNumberOfCoefficients(m_NumberOfCoefficients);
		clusteredShapeTransform->SetIdentity();

		ShapeParametersType clusteredShapeParameters = clusteredShapeTransform->GetParameters();
		for (int id=0;id<m_ReferenceShapeMesh->GetNumberOfPoints();id++)
		{
			clusteredShapeParameters.Fill(0.0);
			for(std::map<int, ShapeClusterItem*>::iterator it=m_ShapeClusterInstances.begin();it!=m_ShapeClusterInstances.end();it++)
			{
				int clusterId = it->first;
				ShapeClusterItem* item = it->second;
				const ShapeParametersType & shapeParam = item->shapeTransform->GetParameters();
				for (int p=0;p<clusteredShapeTransform->GetUsedNumberOfCoefficients();p++)
				{
					clusteredShapeParameters[p] += shapeParam[p]*m_ClusterWeights(id, clusterId);
				}
			}
			clusteredShapeTransform->SetParameters(clusteredShapeParameters);
			m_UpdateShapeMesh->SetPoint(id, clusteredShapeTransform->TransformPoint(m_ReferenceShapeMesh->GetPoint(id)));
		}
	}
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::UpdateRigid()
	{
		//Store updated rigid mesh into m_OutputMesh.
		for (int id=0;id<m_UpdateShapeMesh->GetNumberOfPoints();id++)
		{
			m_OutputMesh->SetPoint(id, m_RigidTranform->TransformPoint(m_UpdateShapeMesh->GetPoint(id)));
		}
	}
}

#endif