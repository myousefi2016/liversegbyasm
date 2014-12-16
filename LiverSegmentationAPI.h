#ifndef _LiverSegmentation_H
#define _LiverSegmentation_H

#include "LiverSegmentationAPI_Export.h"
#include "kmNotifierBase.h"

using namespace std;
namespace km
{
	LiverSegmentationAPI_EXPORT 
		void LiverSeg( 
		km::NotifierBase * notifier,
		const char* outputdir,
		const char* inputImageFile,
		const char* SSMFile,
		const char* boundaryClassifierFile,
		const char* liverClassifierFile,
		const char* adaboostSegmentFile,
		const char* geoFile,
		const char* atlasImageFile,
		const char* configFile);

}

#endif //_LiverSegmentation_H