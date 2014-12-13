
#include <iostream>
#include <fstream>

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"
#include "itkMesh.h"

#include "kmProfileClassifier.h"
#include "kmGlobal.h"
#include "kmUtility.h"

using namespace km;
using namespace std;

int main(int argc, char * argv[] )
{
	if(argc<7)
	{
		std::cout<<"Usage:"<<std::endl;
		std::cout<<"samplePts labeledMeshInput labeledMeshOutput numberOfData numberOfClusters outputHDF5File";
		return EXIT_FAILURE;
	}

	itk::TimeProbesCollectorBase chronometer;
	itk::MemoryProbesCollectorBase memorymeter;
	memorymeter.Start( "starting" );
	chronometer.Start( "starting" );

	const char* sampleFile = argv[1];
	const char* labeledMeshInputFile = argv[2];
	const char* labeledMeshOutputFile = argv[3];
	const int numberOfData = atoi(argv[4]);
	const int numberOfClusters = atoi(argv[5]);
	const char* outputHDF5File = argv[6];

	typedef itk::Mesh<float, 3> MeshType;
	MeshType::Pointer labeledMesh = km::readMesh<MeshType>(labeledMeshInputFile);
	const int numberOfPoints = labeledMesh->GetNumberOfPoints();

	km::ProfileContainer profileContainer;
	profileContainer.setShapeNumber(numberOfData);
	profileContainer.setShapePointsNumber(numberOfPoints);
	profileContainer.loadSamples( sampleFile );
	profileContainer.cluster(numberOfClusters);
	profileContainer.save( outputHDF5File );

	km::ProfileContainerUtils<MeshType>::assignClusterLabels(profileContainer, labeledMesh);
	km::writeMesh<MeshType>(labeledMeshOutputFile, labeledMesh);

	chronometer.Stop( "starting" );
	memorymeter.Stop( "starting" );
	chronometer.Report( std::cout );
	memorymeter.Report( std::cout );
	
	return EXIT_SUCCESS;
}