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
	if(argc<9)
	{
		std::cerr<<"Usage: "<<std::endl;
		std::cout<<" inputImage"
		         <<" SSMFile"
                 <<" outputDir"
				 <<" boundaryClassifierFile"
				 <<" liverClassifierFile"
				 <<" adaboostSegmentFile"
				 <<" geoFile"
				 <<" atlasImageFile"
				 <<" configFile"
				 <<" varianceMapFile"
				 <<std::endl;
		return -1;
	}

	NotifierImp* notifier = new NotifierImp;

	int paramidx = 1;

	const char* inputImageFile = argv[paramidx++];
	const char* SSMFile = argv[paramidx++];
	outputdir = argv[paramidx++];
	const char* boundaryClassifierFile = argv[paramidx++];
	const char* liverClassifierFile = argv[paramidx++];
	const char* adaboostSegmentFile = argv[paramidx++];
	const char* geoFile = argv[paramidx++];
	const char* atlasImageFile = argv[paramidx++];
	const char* configFile = argv[paramidx++];
	const char* varianceMapFile = argv[paramidx++];

	std::cout<<"** input image                 : " << inputImageFile << std::endl;
	std::cout<<"** SSM file                    : " << SSMFile << std::endl;
	std::cout<<"** output dir                  : " << outputdir << std::endl;
	std::cout<<"** boundary classifier file    : " << boundaryClassifierFile << std::endl;
	std::cout<<"** liver classifier file       : " << liverClassifierFile << std::endl;
	std::cout<<"** adaboost segmentation file  : " << adaboostSegmentFile << std::endl;
	std::cout<<"** geometry file               : " << geoFile << std::endl;
	std::cout<<"** atlas image                 : " << atlasImageFile << std::endl;
	std::cout<<"** config file                 : " <<configFile<<std::endl;
	std::cout<<"** variance map                : " <<varianceMapFile<<std::endl;

	LiverSeg( notifier,
	          inputImageFile,
			  SSMFile,
			  boundaryClassifierFile,
			  liverClassifierFile,
			  adaboostSegmentFile,
			  geoFile,
			  atlasImageFile,
			  configFile,
			  varianceMapFile);

	return 0;
}