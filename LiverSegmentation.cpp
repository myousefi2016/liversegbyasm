#include <iostream>

#include "LiverSegmentationAPI.h"
#include "kmNotifierBase.h"
#include "kmUtility.h"

namespace km
{
	template<class TImage, class TPointSet>
	class NotifierImp:public NotifierBase<TImage, TPointSet>
	{
	public:
		char* outputdir;

		NotifierImp():
		{
			outputdir = NULL;
		}

		~NotifierImp()
		{

		}

		void notifyMesh( TPointSet* mesh, const char* filename )
		{
			if (outputdir != NULL)
			{
				km::writeMesh<TPointSet>(outputdir, filename, mesh);
			}
			else
			{
				km::writeMesh<TPointSet>(filename, mesh);
			}
		}

		void notifyImage( TImage* image, const char* filename )
		{
			if (outputdir != NULL)
			{
				km::writeImage<TImage>(outputdir, filename, image);
			}
			else
			{
				km::writeImage<TImage>(filename, image);
			}
		}
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

	NotifierImp* notifier = new NotifierImp;
	notifier->outputdir = outputdir;

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