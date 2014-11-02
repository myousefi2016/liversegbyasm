#ifndef __kmKmeansClassifier_h
#define __kmKmeansClassifier_h

#include "itkMesh.h"
#include "itkVector.h"
#include "itkListSample.h"
#include "itkKdTree.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkMinimumDecisionRule.h"
#include "itkImageToListSampleAdaptor.h"
#include "itkSampleClassifierFilter.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkImageKmeansModelEstimator.h"
#include "itkImageRegionIterator.h"

#define MAX_MESH_NUM (200)
#define MV_DIM (PROFILE_DIM*CLASS_DIM*MAX_MESH_NUM)

unsigned int K = 5;

namespace km
{
	template< class MeshType, class ProfileType >
	void
		kmeansClassify( typename MeshType::Pointer & mesh, std::vector<ProfileType>& profiles, const unsigned int numberOfMeshes)
	{
		if (K <= 1)
		{
			km::assigneMesh<MeshType>( mesh, 0 );
			return;
		}

		unsigned int numberOfPoints = mesh->GetNumberOfPoints();

		unsigned int numberOfProfiles = profiles.size();

		typedef itk::Vector< double, MV_DIM > MeasurementVectorType;
		typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
		SampleType::Pointer sample = SampleType::New();
		sample->SetMeasurementVectorSize( MV_DIM );

		MeasurementVectorType mv;
		mv.Fill( 0 );
		for (int i=0;i<numberOfPoints;i++)
		{
			for (int j=0;j<numberOfMeshes;j++)
			{
				for (int k=0;k<CLASS_DIM;k++)
				{
					ProfileType p = profiles[ j*numberOfPoints*CLASS_DIM + i*CLASS_DIM + k ];

					for (int t=0;t<PROFILE_DIM;t++)
					{
						mv[ j*CLASS_DIM*PROFILE_DIM + k*PROFILE_DIM + t ] = p.value[t];
					}
				}
			}

			sample->PushBack( mv );
		}

		typedef itk::Statistics::WeightedCentroidKdTreeGenerator< SampleType > TreeGeneratorType;
		TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();

		treeGenerator->SetSample( sample );
		treeGenerator->SetBucketSize( 16 );
		treeGenerator->Update();

		typedef TreeGeneratorType::KdTreeType TreeType;
		typedef itk::Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;
		EstimatorType::Pointer estimator = EstimatorType::New();

		EstimatorType::ParametersType initialMeans( MV_DIM * K );
		initialMeans.Fill(0.0);

		estimator->SetParameters( initialMeans );
		estimator->SetKdTree( treeGenerator->GetOutput() );
		estimator->SetMaximumIteration( 200 );
		estimator->SetCentroidPositionChangesThreshold(0.0);
		estimator->StartOptimization();

		EstimatorType::ParametersType estimatedMeans = estimator->GetParameters();
		//for ( unsigned int i = 0 ; i < K ; ++i )
		//{
		//	std::cout << "cluster[" << i << "] " << std::endl;
		//	std::cout << "    estimated mean : " << estimatedMeans[i] << std::endl;
		//}
		//std::cout << "estimatedMeans: " << estimatedMeans << std::endl;

		typedef itk::Statistics::DistanceToCentroidMembershipFunction< MeasurementVectorType > MembershipFunctionType;
		typedef itk::Statistics::MinimumDecisionRule DecisionRuleType;
		DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();

		typedef itk::Statistics::SampleClassifierFilter< SampleType > ClassifierType;
		ClassifierType::Pointer classifier = ClassifierType::New();

		classifier->SetDecisionRule( decisionRule );
		classifier->SetInput( sample );
		classifier->SetNumberOfClasses( K );

		typedef ClassifierType::ClassLabelVectorObjectType
			ClassLabelVectorObjectType;
		typedef ClassifierType::ClassLabelVectorType ClassLabelVectorType;
		typedef ClassifierType::ClassLabelType ClassLabelType;

		ClassLabelVectorObjectType::Pointer classLabelsObject =
			ClassLabelVectorObjectType::New();
		ClassLabelVectorType& classLabelsVector = classLabelsObject->Get();
		for (int k=0;k<K;k++)
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
		for ( unsigned int i = 0 ; i < K ; i++ )
		{
			MembershipFunctionType::Pointer membershipFunction = MembershipFunctionType::New();
			MembershipFunctionType::CentroidType centroid( sample->GetMeasurementVectorSize() );
			for ( unsigned int j = 0 ; j < sample->GetMeasurementVectorSize(); j++ )
			{
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

		//while ( labelIt != labelItEnd )
		//{
		//	std::cout << "measurement vector = " << labelIt.GetMeasurementVector()
		//		<< " class label = " << labelIt.GetClassLabel()
		//		<< std::endl;
		//	++labelIt;
		//}

		typedef MeshType::PointDataContainer MeshPointDataContainer;
		MeshPointDataContainer::Pointer pointdata = mesh->GetPointData();
		pointdata->Reserve( numberOfPoints );
		MeshPointDataContainer::Iterator pointdataIt = pointdata->Begin();
		MeshPointDataContainer::Iterator pointdataItEnd = pointdata->End();

		std::vector<unsigned int> kmclasses;
		labelIt = membershipSample->Begin();
		while(pointdataIt!=pointdataItEnd && labelIt!=labelItEnd)
		{
			unsigned int labelVal = labelIt.GetClassLabel();
			pointdataIt.Value() = labelVal;

			kmclasses.push_back( labelVal );

			++labelIt;
			++pointdataIt;
		}

		for (int i=0;i<profiles.size();i++)
		{
			ProfileType * p = &profiles[i];
			p->kmclass = kmclasses[p->index];
		}
	}

	template< class MeshType, class ProfilePackageType >
	void
		kmeansClassifyWithWeight( typename MeshType::Pointer & mesh, ProfilePackageType & profilePackage, const unsigned int numberOfMeshes)
	{
		if (K <= 1)
		{
			km::assigneMesh<MeshType>( mesh, 0 );
			return;
		}

		typedef itk::Vector< double, MV_DIM > MeasurementVectorType;
		typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
		SampleType::Pointer sample = SampleType::New();
		sample->SetMeasurementVectorSize( MV_DIM );

		ProfilePackageType::ProfileContainerType & profiles = profilePackage.profiles;

		unsigned int numberOfPoints = mesh->GetNumberOfPoints();

		unsigned int numberOfProfiles = profiles.size();

		MeasurementVectorType mv;
		mv.Fill( 0 );
		for (int i=0;i<numberOfPoints;i++)
		{
			for (int j=0;j<numberOfMeshes;j++)
			{
				for (int k=0;k<CLASS_DIM;k++)
				{
					ProfileType p = profiles[ j*numberOfPoints*CLASS_DIM + i*CLASS_DIM + k ];

					for (int t=0;t<PROFILE_DIM;t++)
					{
						mv[ j*CLASS_DIM*PROFILE_DIM + k*PROFILE_DIM + t ] = p.value[t];
					}
				}
			}

			sample->PushBack( mv );
		}

		typedef itk::Statistics::WeightedCentroidKdTreeGenerator< SampleType > TreeGeneratorType;
		TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();

		treeGenerator->SetSample( sample );
		treeGenerator->SetBucketSize( 16 );
		treeGenerator->Update();

		typedef TreeGeneratorType::KdTreeType TreeType;
		typedef itk::Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;
		EstimatorType::Pointer estimator = EstimatorType::New();

		EstimatorType::ParametersType initialMeans( MV_DIM * K );
		initialMeans.Fill(0.0);

		estimator->SetParameters( initialMeans );
		estimator->SetKdTree( treeGenerator->GetOutput() );
		estimator->SetMaximumIteration( 200 );
		estimator->SetCentroidPositionChangesThreshold(0.0);
		estimator->StartOptimization();

		EstimatorType::ParametersType estimatedMeans = estimator->GetParameters();
		//for ( unsigned int i = 0 ; i < K ; ++i )
		//{
		//	std::cout << "cluster[" << i << "] " << std::endl;
		//	std::cout << "    estimated mean : " << estimatedMeans[i] << std::endl;
		//}
		//std::cout << "estimatedMeans: " << estimatedMeans << std::endl;

		typedef itk::Statistics::DistanceToCentroidMembershipFunction< MeasurementVectorType > MembershipFunctionType;
		typedef itk::Statistics::MinimumDecisionRule DecisionRuleType;
		DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();

		typedef itk::Statistics::SampleClassifierFilter< SampleType > ClassifierType;
		ClassifierType::Pointer classifier = ClassifierType::New();

		classifier->SetDecisionRule( decisionRule );
		classifier->SetInput( sample );
		classifier->SetNumberOfClasses( K );

		typedef ClassifierType::ClassLabelVectorObjectType
			ClassLabelVectorObjectType;
		typedef ClassifierType::ClassLabelVectorType ClassLabelVectorType;
		typedef ClassifierType::ClassLabelType ClassLabelType;

		ClassLabelVectorObjectType::Pointer classLabelsObject =
			ClassLabelVectorObjectType::New();
		ClassLabelVectorType& classLabelsVector = classLabelsObject->Get();
		for (int k=0;k<K;k++)
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
		for ( unsigned int i = 0 ; i < K ; i++ )
		{
			MembershipFunctionType::Pointer membershipFunction = MembershipFunctionType::New();
			MembershipFunctionType::CentroidType centroid( sample->GetMeasurementVectorSize() );
			for ( unsigned int j = 0 ; j < sample->GetMeasurementVectorSize(); j++ )
			{
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

		//while ( labelIt != labelItEnd )
		//{
		//	std::cout << "measurement vector = " << labelIt.GetMeasurementVector()
		//		<< " class label = " << labelIt.GetClassLabel()
		//		<< std::endl;
		//	++labelIt;
		//}

		typedef MeshType::PointDataContainer MeshPointDataContainer;
		MeshPointDataContainer::Pointer pointdata = mesh->GetPointData();
		pointdata->Reserve( numberOfPoints );
		MeshPointDataContainer::Iterator pointdataIt = pointdata->Begin();
		MeshPointDataContainer::Iterator pointdataItEnd = pointdata->End();

		std::vector<unsigned int> kmclasses;
		labelIt = membershipSample->Begin();
		while(pointdataIt!=pointdataItEnd && labelIt!=labelItEnd)
		{
			unsigned int labelVal = labelIt.GetClassLabel();
			pointdataIt.Value() = labelVal;

			kmclasses.push_back( labelVal );

			++labelIt;
			++pointdataIt;
		}

		for (int i=0;i<profiles.size();i++)
		{
			ProfileType * p = &profiles[i];
			p->kmclass = kmclasses[p->index];
		}

		//Test KNN to get weight for each classifier.
		testKmeansWeight( profilePackage );
	}

	template< class ProfilePackageType >
	void
		testKmeansWeight( ProfilePackageType & profilePackage )
	{
		typedef ProfilePackageType::ProfileType ProfileType;
		typedef ProfilePackageType::ProfileContainerType ProfileContainerType;
		typedef km::KNNProfileClassifier<ProfileType> KNNProfileClassifierType;

		ProfileContainerType allProfiles = profilePackage.profiles;

		double testRatio = 0.20;

		unsigned int numberOfSamples = allProfiles.size();
		KM_DEBUG_PRINT( "Number of all samples ", numberOfSamples );

		unsigned int numberOfTest = static_cast<unsigned int>(numberOfSamples*testRatio);
		KM_DEBUG_PRINT( "Number of test samples ", numberOfTest );

		unsigned int startIndexOfTestSample = numberOfSamples-numberOfTest;

		if (false)
		{
			ProfileContainerType profilesForTraining;

			ProfileContainerType profilesForTestig;

			for (int i=0;i<startIndexOfTestSample-1;i++)
			{
				profilesForTraining.push_back( allProfiles[i] );
			}

			KM_DEBUG_PRINT( "Number of training samples ", profilesForTraining.size() );

			for (int i=startIndexOfTestSample-1;i<numberOfSamples;i++)
			{
				profilesForTestig.push_back( allProfiles[i] );
			}

			KM_DEBUG_PRINT( "Number of testing samples ", profilesForTestig.size() );

			KNNProfileClassifierType classifier;
			classifier.buildFromPoints( profilesForTraining );
			classifier.testFromPoints( profilesForTestig );
		}
		else
		{
			std::vector<ProfileContainerType*> allTrainingProfilesVec;
			std::vector<ProfileContainerType*> allTestingProfilesVec;
			std::vector<KNNProfileClassifierType*> allClassifiersVeco;

			for (int i=0;i<K;i++)
			{
				ProfileContainerType * pctraining = new ProfileContainerType;
				allTrainingProfilesVec.push_back(pctraining);

				ProfileContainerType * pctesting = new ProfileContainerType;
				allTestingProfilesVec.push_back(pctesting);

				KNNProfileClassifierType * classifier = new KNNProfileClassifierType;
				allClassifiersVeco.push_back(classifier);
			}

			for (int i=0;i<startIndexOfTestSample;i++)
			{
				unsigned kidx = allProfiles[i].kmclass;
				if (kidx>=K || kidx<0)
				{
					KM_DEBUG_ERROR( "Number of kmeans class is not 5" );
				}
				ProfileContainerType * pctmp = allTrainingProfilesVec[kidx];
				pctmp->push_back(allProfiles[i]);
			}

			for (int i=startIndexOfTestSample;i<numberOfSamples;i++)
			{
				unsigned kidx = allProfiles[i].kmclass;
				if (kidx>=K || kidx<0)
				{
					KM_DEBUG_ERROR( "Number of kmeans class is not 5" );
				}
				ProfileContainerType * pctmp = allTestingProfilesVec[kidx];
				pctmp->push_back(allProfiles[i]);

				//std::cout<<kidx<<", "<<i<<std::endl;
			}

			for (int i=0;i<K;i++)
			{
				KM_DEBUG_PRINT( "Test classifier: ", i );
				allClassifiersVeco[i]->buildFromPoints( *allTrainingProfilesVec[i] );
				double successRate = allClassifiersVeco[i]->testFromPoints( *allTestingProfilesVec[i] );

				profilePackage.AddKmeansPair( i, successRate );
			}
		}
	}

}

#endif