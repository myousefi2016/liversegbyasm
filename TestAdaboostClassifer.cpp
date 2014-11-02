
#include <iostream>
#include <fstream>
#include <string>

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include "flann/flann.hpp"
#include "flann/io/hdf5.h"

#include "kmKNNProfileClassifier-FLANN.h"
#include "kmUtility.h"

using namespace km;
using namespace std;
using namespace flann;

int main(int argc, char * argv[] )
{
	if(argc<4)
	{
		std::cout<<"Usage:"<<std::endl;
		std::cout<<"SamplePts AdaboostFile Category";
		return EXIT_FAILURE;
	}

	const char* samplefile = argv[1];
	const char* adaboostfile = argv[2];
	const char* categorystr = argv[3];

	PROFILE_CATEGORY category;
	if (strcmp(categorystr,"BOUNDARY")==0)
	{
		category = BOUNDARY;
	}
	else if (strcmp(categorystr,"LIVER")==0)
	{
		category = LIVER;
	}
	else
	{
		category = DEFAULT;
	}

	typedef KNNProfileClassifier KNNProfileClassifierType;
	
	KNNProfileClassifierType knnclassifier(category);
	knnclassifier.load( samplefile );
	
	typedef AdaboostProfileClassifier AdaboostProfileClassifierType;
	AdaboostProfileClassifierType adaboostclassifier(knnclassifier.profile_category);

	std::cout<<adaboostclassifier.profile_category<<std::endl;

	////Adaboost Training
	//itk::TimeProbesCollectorBase chronometer1;
	//itk::MemoryProbesCollectorBase memorymeter1;
	//memorymeter1.Start( "Adaboost_Training" );
	//chronometer1.Start( "Adaboost_Training" );

	//adaboostclassifier.train(knnclassifier, outputdir);
	//
	//chronometer1.Stop( "Adaboost_Training" );
	//memorymeter1.Stop( "Adaboost_Training" );
	//chronometer1.Report( std::cout );
	//memorymeter1.Report( std::cout );

	adaboostclassifier.load(adaboostfile);

	//Adaboost Testing
	itk::TimeProbesCollectorBase chronometer2;
	itk::MemoryProbesCollectorBase memorymeter2;
	memorymeter2.Start( "Adaboost_Testing" );
	chronometer2.Start( "Adaboost_Testing" );

	adaboostclassifier.test(knnclassifier, 0.5);
	
	chronometer2.Stop( "Adaboost_Testing" );
	memorymeter2.Stop( "Adaboost_Testing" );
	chronometer2.Report( std::cout );
	memorymeter2.Report( std::cout );

	return EXIT_SUCCESS;
}