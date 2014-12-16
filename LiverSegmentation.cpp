#include <iostream>
#include <sstream>

#include "LiverSegmentationAPI.h"

using namespace km;

int main(int argc, char* argv[])
{
	if(argc<9)
	{
		std::cerr <<"Usage: "<<std::endl;
		std::cout <<" outputDir"
				  <<" inputImage"
		          <<" SSMFile"
				  <<" boundaryClassifierFile"
				  <<" liverClassifierFile"
				  <<" adaboostSegmentFile"
				  <<" geoFile"
				  <<" atlasImageFile"
				  <<" configFile"
				  <<std::endl;
		return -1;
	}

	int paramidx = 1;

	const char* outputdir = argv[paramidx++];
	const char* inputImageFile = argv[paramidx++];
	const char* SSMFile = argv[paramidx++];
	const char* boundaryClassifierFile = argv[paramidx++];
	const char* liverClassifierFile = argv[paramidx++];
	const char* adaboostSegmentFile = argv[paramidx++];
	const char* geoFile = argv[paramidx++];
	const char* atlasImageFile = argv[paramidx++];
	const char* configFile = argv[paramidx++];

	std::cout<<"** output dir                  : " << outputdir << std::endl;
	std::cout<<"** input image                 : " << inputImageFile << std::endl;
	std::cout<<"** SSM file                    : " << SSMFile << std::endl;
	std::cout<<"** boundary classifier file    : " << boundaryClassifierFile << std::endl;
	std::cout<<"** liver classifier file       : " << liverClassifierFile << std::endl;
	std::cout<<"** adaboost segmentation file  : " << adaboostSegmentFile << std::endl;
	std::cout<<"** geometry file               : " << geoFile << std::endl;
	std::cout<<"** atlas image                 : " << atlasImageFile << std::endl;
	std::cout<<"** config file                 : " <<configFile<<std::endl;

	std::cout<<"Start to call LiverSeg(..)"<<std::endl;

	LiverSeg( NULL,
			  outputdir,
	          inputImageFile,
			  SSMFile,
			  boundaryClassifierFile,
			  liverClassifierFile,
			  adaboostSegmentFile,
			  geoFile,
			  atlasImageFile,
			  configFile);

	return 0;
}