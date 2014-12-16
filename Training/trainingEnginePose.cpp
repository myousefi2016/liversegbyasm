
#include <sstream>
#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include "kmPoseTraining.h"

/************************************************************************
����: ѵ������

����:
1: �����б��ı�. �ı�������ʽ: xxx/xxx/xxx1.mhd xxx/xxx/xxx2.mhd
...
2: ���ڶ���Ĳο������������б��е�λ�� ��1��ʼ,����xxx1.mhaΪ�ο�����,�����Ϊ1
3: ����ļ���
4: ��������ļ���ǰ׺,��������Ϊ"Liver", �������������ΪLiver.1.vtk, Liver.2.vtk...

************************************************************************/

int main(int argc, char* argv[])
{
	if( argc < 5 )
	{
		std::cerr<<"Usage:   origDataList segDataList referenceOrig referenceSeg OutputDir"<<std::endl;
		return EXIT_FAILURE;
	}
	itk::TimeProbesCollectorBase chronometer;
	chronometer.Start( "shapetraining" );

	const char* origlistfile = argv[1];
	const char* seglistfile  = argv[2];
	//int refindex = atoi( argv[3] );
	const char* referenceorigfile = argv[3];
	const char* referencesegfile  = argv[4];
	const char* outputdir = argv[5];

	int resultcode = km::trainPose(origlistfile, seglistfile, referenceorigfile, referencesegfile, outputdir);
	
	if(resultcode==EXIT_SUCCESS)
	{
		KM_DEBUG_INFO( "Train pose done without error!" );
	}
	else
	{
		KM_DEBUG_INFO( "Train pose done with error!" );
	}

	chronometer.Stop( "shapetraining" );
	chronometer.Report( std::cout );

	return resultcode;
}

