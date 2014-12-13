
#include <iostream>
#include <fstream>
#include <string>

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include "kmProfileContainer.h"
#include "kmProfileClassifier.h"
#include "kmGlobal.h"
#include "kmUtility.h"

using namespace km;
using namespace std;

int main(int argc, char * argv[] )
{
	
	itk::TimeProbesCollectorBase chronometer2;
	itk::MemoryProbesCollectorBase memorymeter2;
	memorymeter2.Start( "Adaboost_Testing" );
	chronometer2.Start( "Adaboost_Testing" );

	const char* profileTxtFile = "D:\\Workspace\\ASM\\projects\\LiverSegbyASM\\experiments\\training\\appearance\\output_20140815\\samples_PLAIN.txt";
	const char* profilesHDF5File = "samples.h5";
	km::ProfileContainer profileContainer;
	profileContainer.setShapeNumber(10);
	profileContainer.setShapePointsNumber(3996);
	profileContainer.loadSamples(profileTxtFile);
	profileContainer.print();
	//profileContainer.cluster(9);
	profileContainer.save(profilesHDF5File);

	std::cout<<"Finish stage 1."<<std::endl;

	const char* classifierHDF5File = "profileClassifier.h5";
	km::ProfileClassifier profileClassifier;
	profileClassifier.train(profileContainer, PLAIN);
	profileClassifier.save(classifierHDF5File);
	profileClassifier.test(profileContainer);
	profileClassifier.print();

	std::cout<<"Finish stage 2."<<std::endl;

	const char* classifierHDF5File2 = "profileClassifier2.h5";
	km::ProfileClassifier profileClassifier2;
	profileClassifier2.load(classifierHDF5File);
	profileClassifier2.test(profileContainer);
	profileClassifier2.save(classifierHDF5File2);

	std::cout<<"Finish stage 3."<<std::endl;
	
	chronometer2.Stop( "Adaboost_Testing" );
	memorymeter2.Stop( "Adaboost_Testing" );
	chronometer2.Report( std::cout );
	memorymeter2.Report( std::cout );

	return EXIT_SUCCESS;
}