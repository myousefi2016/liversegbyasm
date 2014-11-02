
#include <iostream>
#include <fstream>

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include "kmUtility.h"
#include "kmKNNProfileClassifier.h"

using namespace km;
using namespace std;


int main(int argc, char * argv[] )
{
	if(argc<3)
	{
		std::cout<<"Usage:"<<std::endl;
		std::cout<<"SamplePts TestRatio";
		return EXIT_FAILURE;
	}

	itk::TimeProbesCollectorBase chronometer;
	itk::MemoryProbesCollectorBase memorymeter;
	memorymeter.Start( "starting" );
	chronometer.Start( "starting" );

	const char* sampleFile = argv[1];
	const double testRaio = atof(argv[2]);

	//Read profile
	typedef std::vector<Profile> ProfileContainerType;

	ProfileContainerType allProfiles;
	readProfilesFromFile( sampleFile, allProfiles );

	unsigned int numberOfSamples = allProfiles.size();
	KM_DEBUG_PRINT( "Number of all samples ", numberOfSamples );

	unsigned int numberOfTest = static_cast<unsigned int>(numberOfSamples*testRaio);
	KM_DEBUG_PRINT( "Number of test samples ", numberOfSamples );

	unsigned int startIndexOfTestSample = numberOfSamples-numberOfTest;

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


	KNNProfileClassifier classifier;
	classifier.buildFromPoints( profilesForTraining );
	classifier.testFromPoints( profilesForTestig );

	return EXIT_SUCCESS;

	//const char* testSampleFile = argv[2];

	//int numberOfClassifier = 1;
	//if (argc>3)
	//{
	//	numberOfClassifier = atoi(argv[3]);
	//}

	//ProfileContainerType allQueryProfiles;
	//readProfilesFromFile( testSampleFile, allQueryProfiles );

	//if (numberOfClassifier==1)
	//{
	//	KNNProfileClassifier classifier;
	//	classifier.buildFromPoints( allProfiles );

	//	classifier.testFromPoints( allQueryProfiles );

	//	return EXIT_SUCCESS;
	//}
	//

	//typedef std::vector<ProfileContainerType> ProfileContainerVectorType;
	//ProfileContainerVectorType separatedProfiles;

	////Initial classifiers
	//typedef std::vector<KNNProfileClassifier*> ClassifierContainerType;
	//ClassifierContainerType classifiers;

	//for (int i=0;i<numberOfClassifier;i++)
	//{
	//	KNNProfileClassifier * c = new KNNProfileClassifier;
	//	classifiers.push_back(c);

	//	ProfileContainerType profilesForClassifier;
	//	separatedProfiles.push_back( profilesForClassifier );
	//}


	////Separate profiles into different container according to the number of classifiers.
	//for (int i=0;i<allProfiles.size();i++)
	//{
	//	Profile profilesample = allProfiles[i];

	//	if (profilesample.index<0 || profilesample.index>numberOfClassifier-1)
	//	{
	//		continue;
	//	}

	//	separatedProfiles[profilesample.index].push_back( profilesample );
	//}


	////Start to classify
	//for (int i=0;i<numberOfClassifier;i++)
	//{
	//	classifiers[i]->buildFromPoints( separatedProfiles[i] );
	//}

	//unsigned int numberOfPositiveCases = 0;
	//unsigned int numberOfNegativeCases = 0;
	//unsigned int numberOfTestCases = 0;

	//for (int i=0;i<allQueryProfiles.size();i++)
	//{
	//	Profile testsample = allQueryProfiles[i];

	//	if (testsample.index<0 || testsample.index>numberOfClassifier-1)
	//	{
	//		continue;
	//	}

	//	KNNProfileClassifier* classifier = classifiers[testsample.index];

	//	ClassType maximumPossiblerType;
	//	double maximumProbablity = 0.0;
	//	
	//	classifier->classify( testsample, maximumPossiblerType, maximumProbablity );

	//	if (maximumPossiblerType == testsample.type){
	//		numberOfPositiveCases++;
	//	}else{
	//		numberOfNegativeCases++;
	//	}

	//	numberOfTestCases++;
	//}

	//std::cout<<"Total test cases: "<<numberOfTestCases<<std::endl;
	//std::cout<<"Positive test cases: "<<numberOfPositiveCases<<std::endl;
	//std::cout<<"Negative test cases: "<<numberOfNegativeCases<<std::endl;
	//std::cout<<"Success positive rate: "<<static_cast<double>(numberOfPositiveCases)/numberOfTestCases<<std::endl;

	chronometer.Stop( "starting" );
	memorymeter.Stop( "starting" );
	chronometer.Report( std::cout );
	memorymeter.Report( std::cout );
}