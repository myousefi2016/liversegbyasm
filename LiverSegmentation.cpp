#include <iostream>

#include "LiverSegmentationAPI.h"

namespace km
{
	class NotifierImp:public Notifier
	{
	public:
		void notifyMesh( MeshType* mesh )
		{
			//km::writeMesh<MeshType>( "tempMesh.vtk", mesh );
			std::cout<<"tempMesh.vtk"<<std::endl;
		}

		void notifyImage( UCharImageType* image )
		{
			//km::writeImage<UCharImageType>( "tempImage.nii.gz", image );
			std::cout<<"tempImage.nii.gz"<<std::endl;
		}
		//virtual notifyFinish();
	};
}

using namespace km;

int main(int argc, char* argv[])
{
	if(argc<10)
	{
		std::cerr<<"Usage: "<<std::endl;
		std::cout<<"       inputImage ssmFile outputDir boundaryClassifierFile liverClassifierFile adaboostSegmentFile geoFile atlasImage configFile"<<std::endl;
		return -1;
	}

	NotifierImp* notifier = new NotifierImp;

	int paramidx = 1;

	const char* inputimagefile = argv[paramidx++];
	const char* ssmfile = argv[paramidx++];
	outputdir = argv[paramidx++];
	const char* boundaryClassifierFile = argv[paramidx++];
	const char* liverClassifierFile = argv[paramidx++];
	const char* adaboostSegmentFile = argv[paramidx++];
	const char* geofile = argv[paramidx++];
	const char* atlasimagefile = argv[paramidx++];
	const char* configFile = argv[paramidx++];

	const KM_STRATEGY strategy = static_cast<KM_STRATEGY>(1);

	std::cout<<"input image: " << inputimagefile << std::endl;
	std::cout<<"SSM file: " << ssmfile << std::endl;
	std::cout<<"output dir: " << outputdir << std::endl;
	std::cout<<"boundary classifier file: " << boundaryClassifierFile << std::endl;
	std::cout<<"liver classifier file: " << liverClassifierFile << std::endl;
	std::cout<<"adaboost segmentation file: " << adaboostSegmentFile << std::endl;
	std::cout<<"geometry file: " << geofile << std::endl;
	std::cout<<"atlas image: " << atlasimagefile << std::endl;
	std::cout<<"config file: "<<configFile<<std::endl;
	std::cout<<"strategy: " << strategy << std::endl;

	LiverSeg( notifier, inputimagefile, ssmfile, boundaryClassifierFile, liverClassifierFile, adaboostSegmentFile, geofile, atlasimagefile, configFile, strategy );

	return 0;
}