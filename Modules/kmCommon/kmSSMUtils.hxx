#ifndef __kmSSMUtils_hxx
#define __kmSSMUtils_hxx

#include "kmSSMUtils.h"

namespace km
{	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::initialize()
	{
		if (m_SSM == NULL){
			std::cout<<"SSM cannot be null."<<std::endl;
			return;
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
	}
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::update(const MeshType * targetMesh, MeshType * outputMesh)
	{
		if (m_Iterations == 0)
		{
			this->initialize();
		}

		if (m_Iterations>=20)
		{
			this->getClusterByClusterId(4)->enabled = true;
		}

		this->rigidTransformFitting(targetMesh);
		this->shapeTransformFitting(targetMesh);
		this->clusteredShapeTransformFitting(targetMesh);
		this->updateShape(true);
		
		this->updateRigid(outputMesh, false);

		m_Iterations++;
	}
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::printTransform()
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
	::calShapeParaDiff()
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
	double
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::getShapeProbability(int pointId)
	{
		return m_ShapeClusterInstances[m_ShapeClusterMap[pointId]]->shapeProbability;
	}

	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::cluster(int numberOfClusters)
	{
		timer->StartTimer();

		this->clearClusters();

		int numberOfModesUsed = km::g_number_principle_components;//this->m_NumberOfCoefficients;
		double multifier = 3.0;
		
		const StatismoMatrixType & basisMatrixForShape = this->m_SSM->GetstatismoImplObj()->GetPCABasisMatrix();
		const unsigned int measureLength = numberOfModesUsed*Dimension;
		const unsigned int measureCount = basisMatrixForShape.rows()/Dimension;

		typedef itk::VariableLengthVector< double> MeasurementVectorType;
		typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
		SampleType::Pointer sample = SampleType::New();
		sample->SetMeasurementVectorSize( measureLength );
		MeasurementVectorType mv;
		mv.SetSize( measureLength );
		mv.Fill( 0 );

		std::cout<<"number of shape coefficients: "<<m_NumberOfCoefficients<<std::endl;
		std::cout<<"number of modes used: "<<numberOfModesUsed<<std::endl;
		std::cout<<"measure length: "<<measureLength<<", count: "<<measureCount<<std::endl;

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
		EstimatorType::ParametersType initialMeans( measureLength * numberOfClusters );
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
		for (int k=0;k<numberOfClusters;k++)
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
		for ( unsigned int i = 0 ; i < numberOfClusters ; i++ ){
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
			this->addCluster(pointId, clusterId);

			++pointId;
			++labelIt;
		}

		this->allocateClusters();

		this->calculateClusterWeights();

		std::cout<<"Actual cluster number: "<<this->m_ShapeClusterInstances.size()<<std::endl;
		timer->StopTimer();
		std::cout<<"Cluster done. "<<timer->GetElapsedTime() << " s." <<endl;
	}
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::rigidTransformFitting(const MeshType* targetMesh)
	{
		//timer->StartTimer();

		unsigned int numberOfLandmarks = m_ReferenceShapeMesh->GetNumberOfPoints();

		/*****************Fit rigid parameters****************/
		MeshType::ConstPointer sourceMeshForRigid = m_UpdateShapeMesh;//km::transformMesh<MeshType, ShapeTransformType>(m_ReferenceShapeMesh, m_ShapeTransform);;
		MeshType::ConstPointer targetMeshForRigid = targetMesh;

		StatismoMatrixType sourceMatrixForRigid, targetMatrixForRigid;
		sourceMatrixForRigid.setConstant(numberOfLandmarks, Dimension+1, 1.0);
		targetMatrixForRigid.setConstant(numberOfLandmarks, Dimension+1, 1.0);
		unsigned long rowcnt = 0;
		for (int idx=0;idx<numberOfLandmarks;idx++)
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

		StatismoVectorType I = StatismoVectorType::Ones(sourceMatrixForRigid.cols());
		StatismoMatrixType Mmatrix = sourceMatrixForRigid.transpose() * sourceMatrixForRigid;
		Mmatrix.diagonal() += I;
		StatismoMatrixType MInverseMatrix = Mmatrix.inverse();
		const StatismoMatrixType& WT = sourceMatrixForRigid.transpose();
		StatismoMatrixType coeffsRigid = MInverseMatrix * (WT * targetMatrixForRigid);

		updateRigidTransform(coeffsRigid.transpose(), this->m_RigidTranform);

		//timer->StopTimer();
		//std::cout<<"Rigid transform fitting done. "<<timer->GetElapsedTime() << " s." <<endl;
	}
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::clusteredRigidTransformFitting(const MeshType* targetMesh)
	{
		//timer->StartTimer();

		unsigned int numberOfLandmarks = m_ReferenceShapeMesh->GetNumberOfPoints();

		/*****************Fit rigid parameters****************/
		MeshType::ConstPointer sourceMeshForRigid = m_UpdateShapeMesh;//km::transformMesh<MeshType, ShapeTransformType>(m_ReferenceShapeMesh, m_ShapeTransform);;
		MeshType::ConstPointer targetMeshForRigid = targetMesh;

		for(std::map<int, ShapeClusterItem*>::iterator it=m_ShapeClusterInstances.begin();it!=m_ShapeClusterInstances.end();it++)
		{
			ShapeClusterItem* clusterItem = it->second;
			int rowcnt = 0;
			for (int idx=0;idx<clusterItem->pointIds.size();idx++)
			{
				int pointId = clusterItem->pointIds[idx];
				PointType targetPt = targetMeshForRigid->GetPoint(pointId);
				PointType sourcePt = sourceMeshForRigid->GetPoint(pointId);
				for (int d=0;d<Dimension;d++)
				{
					clusterItem->sourceMatrixForRigid(rowcnt, d) = sourcePt[d];
					clusterItem->targetMatrixForRigid(rowcnt, d) = targetPt[d];
				}
				rowcnt++;
			}

			StatismoVectorType I = StatismoVectorType::Ones(clusterItem->sourceMatrixForRigid.cols());
			StatismoMatrixType Mmatrix = clusterItem->sourceMatrixForRigid.transpose() * clusterItem->sourceMatrixForRigid;
			Mmatrix.diagonal() += I;
			StatismoMatrixType MInverseMatrix = Mmatrix.inverse();
			const StatismoMatrixType& WT = clusterItem->sourceMatrixForRigid.transpose();
			StatismoMatrixType coeffsRigid = MInverseMatrix * (WT * clusterItem->targetMatrixForRigid);

			clusterItem->rigidParamtersPre = clusterItem->rigidTransform->GetParameters();
			updateRigidTransform(coeffsRigid.transpose(), clusterItem->rigidTransform);
		}

		//timer->StopTimer();
		//std::cout<<"Rigid transform fitting done. "<<timer->GetElapsedTime() << " s." <<endl;
	}
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::shapeTransformFitting(const MeshType* targetMesh)
	{
		//timer->StartTimer();

		unsigned int numberOfLandmarks = m_ReferenceShapeMesh->GetNumberOfPoints();

		/*****************Fit shape parameters****************/
		RigidTransformType::Pointer inversedRigidTransform = RigidTransformType::New();
		m_RigidTranform->GetInverse( inversedRigidTransform );
		MeshType::Pointer targetMeshForShape = km::transformMesh<MeshType, RigidTransformType>(targetMesh, inversedRigidTransform);
		MeshType::Pointer sourceMeshForShape = this->m_MeanShapeMesh;

		const StatismoMatrixType & basisMatrixForShape = this->m_SSM->GetstatismoImplObj()->GetPCABasisMatrix();
		StatismoVectorType sourceMatrixForShape = StatismoVectorType::Zero(basisMatrixForShape.rows());
		StatismoVectorType targetMatrixForShape = StatismoVectorType::Zero(basisMatrixForShape.rows());
		unsigned long rowcnt = 0;
		for (int idx=0;idx<numberOfLandmarks;idx++)
		{
			PointType sourcePt = sourceMeshForShape->GetPoint(idx);
			PointType targetPt = targetMeshForShape->GetPoint(idx);
			for (int d=0;d<Dimension;d++)
			{
				sourceMatrixForShape[rowcnt] = sourcePt[d];
				targetMatrixForShape[rowcnt] = targetPt[d];
				rowcnt++;
			}
		}

		StatismoVectorType I = StatismoVectorType::Ones(basisMatrixForShape.cols());
		StatismoMatrixType Mmatrix = basisMatrixForShape.transpose() * basisMatrixForShape;
		Mmatrix.diagonal() += I;
		StatismoMatrixType MInverseMatrix = Mmatrix.inverse();
		const StatismoMatrixType& WT = basisMatrixForShape.transpose();
		StatismoVectorType coeffsShape = MInverseMatrix * (WT * (targetMatrixForShape-sourceMatrixForShape));

		this->m_ShapeParametersPre = m_ShapeTransform->GetParameters();
		updateShapeTransform(coeffsShape, m_ShapeTransform);

		//timer->StopTimer();
		//std::cout<<"Shape transform fitting done. "<<timer->GetElapsedTime() << " s." <<endl;
	}
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::clusteredShapeTransformFitting(const MeshType* targetMesh)
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
			if (!clusterItem->enabled)
			{
				clusterItem->shapeParametersPre = clusterItem->shapeTransform->GetParameters();
				clusterItem->shapeTransform->SetParameters(m_ShapeTransform->GetParameters());
				continue;
			}

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
			clusterItem->shapeParametersPre  = clusterItem->shapeTransform->GetParameters();
			clusterItem->shapeProbability = this->updateShapeTransform(coeffs, clusterItem->shapeTransform);
		}

		//this->updateShape();

		//timer->StopTimer();
		//std::cout<<"Clustered shape transform fitting done. "<<timer->GetElapsedTime() << " s." <<endl;
	}	
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::updateRigidTransform(const StatismoMatrixType & mat, RigidTransformType* rigidTransform)
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
	double
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::updateShapeTransform(const StatismoVectorType & vec, ShapeTransformType* shapeTransform)
	{
		ShapeParametersType shapeParameters = shapeTransform->GetParameters();
		double shapeProbability = 0.0;
		bool abnormalFlag = false;
		for (int i=0;i<m_NumberOfCoefficients;i++)
		{
			shapeProbability += km::Math::cdf_outside(vec[i])*m_ShapeVarianceWeights[i];
			shapeParameters[i] = km::Math::setToBetween(vec[i], -2.5, 2.5);
			//shapeProbability += km::Math::cdf_outside(shapeParameters[i])*m_ShapeVarianceWeights[i];
		}
		shapeTransform->SetParameters(shapeParameters);

		return shapeProbability;
	}
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::addCluster(int pointId, int clusterId)
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
	::clearClusters()
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
	::allocateClusters()
	{
		timer->StartTimer();

		m_MeanShapeMesh->GetPointData()->Reserve(m_MeanShapeMesh->GetNumberOfPoints());
		for(std::map<int, ShapeClusterItem*>::iterator it=m_ShapeClusterInstances.begin();it!=m_ShapeClusterInstances.end();it++)
		{
			ShapeClusterItem* item = it->second;
			if (item!=NULL)
			{
				item->enabled = false;
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

				item->clusterCentroid.Fill(0);
				for(int idx=0; idx < item->pointIds.size(); idx++)
				{
					int pointId = item->pointIds[idx];
					PointType pt = m_MeanShapeMesh->GetPoint(pointId);
					for (int d=0;d<Dimension;d++)
					{
						item->clusterCentroid[d] += pt[d];
						for (int c=0;c<item->basisMatrix.cols();c++){
							item->basisMatrix(idx*Dimension+d, c) = origBasisMatrix(pointId*Dimension+d, c);
						}
					}
					item->meanPoints->PushBack(pt);
					m_MeanShapeMesh->SetPointData(pointId, item->clusterId);
				}
				for (int d=0;d<Dimension;d++){
					item->clusterCentroid[d] /= item->pointIds.size();
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
				//std::cout<<"       Centroid: "<<item->clusterCentroid<<std::endl;
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
	::calculateClusterWeights()
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
	::updateShape(bool clustered)
	{
		if (!clustered)
		{
			km::transformMesh<MeshType, ShapeTransformType>(m_ReferenceShapeMesh, m_UpdateShapeMesh, m_ShapeTransform);
		}
		else
		{
			ShapeTransformType::Pointer clusteredShapeTransform = ShapeTransformType::New();
			clusteredShapeTransform->SetStatisticalModel(m_SSM);
			clusteredShapeTransform->SetUsedNumberOfCoefficients(m_NumberOfCoefficients);
			clusteredShapeTransform->SetIdentity();

			//ShapeParametersType shapeParameters = m_ShapeTransform->GetParameters();

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
	}
	
	template< class TMesh, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
	SSMUtils< TMesh, TStatisticalModel, TRigidTransform, TShapeTransform>
	::updateRigid(MeshType * outputMesh, bool clustered)
	{
		if (!clustered)
		{
			km::transformMesh<MeshType, RigidTransformType>(m_UpdateShapeMesh, outputMesh, m_RigidTranform);
		}
		else
		{
			RigidTransformType::Pointer clusteredRigidTransform = RigidTransformType::New();
			clusteredRigidTransform->SetFixedParameters(m_RigidTranform->GetFixedParameters());
			clusteredRigidTransform->SetIdentity();

			RigidParametersType clusteredRigidParameters = clusteredRigidTransform->GetParameters();
			for (int id=0;id<m_ReferenceShapeMesh->GetNumberOfPoints();id++)
			{
				clusteredRigidParameters.Fill(0.0);
				for(std::map<int, ShapeClusterItem*>::iterator it=m_ShapeClusterInstances.begin();it!=m_ShapeClusterInstances.end();it++)
				{
					int clusterId = it->first;
					ShapeClusterItem* item = it->second;
					const ShapeParametersType & rigidParam = item->rigidTransform->GetParameters();
					for (int p=0;p<m_RigidTranform->GetNumberOfParameters();p++)
					{
						clusteredRigidParameters[p] += rigidParam[p]*m_ClusterWeights(id, clusterId);
					}
				}
				clusteredRigidTransform->SetParameters(clusteredRigidParameters);
				outputMesh->SetPoint(id, clusteredRigidTransform->TransformPoint(m_UpdateShapeMesh->GetPoint(id)));
			}
		}
	}
}

#endif