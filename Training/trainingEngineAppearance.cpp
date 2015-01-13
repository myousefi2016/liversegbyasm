
#include <sstream>

#include "kmAppearanceTraining.h"

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
		std::cerr<<"Usage:   origDataList meshList referenceGeometry configFile OutputDir"<<std::endl;
		return EXIT_FAILURE;
	}

	const char* origlistfile = argv[1];
	const char* meshlistfile  = argv[2];
	const char* referencegeometry = argv[3];
	const char* configfile = argv[4];
	const char* outputdir = argv[5];

	itk::TimeProbesCollectorBase chronometer;
	itk::MemoryProbesCollectorBase memorymeter;
	chronometer.Start( "appearancetraining" );
	memorymeter.Start( "appearancetraining" );

	KM_DEBUG_INFO("Load config file...");
	km::Config::loadConfig(configfile);

	//For liver tissue classification
	int resultcode = km::trainAppearance(
		origlistfile 
		,meshlistfile
		,referencegeometry
		,outputdir);

	if(resultcode==EXIT_SUCCESS)
	{
		KM_DEBUG_INFO( "Train appearance done successfully!" );
	}
	else
	{
		KM_DEBUG_INFO( "Train appearance done with error!" );
	}

	chronometer.Stop( "appearancetraining" );
	memorymeter.Stop( "appearancetraining" );
	chronometer.Report( std::cout );
	memorymeter.Report( std::cout );

	system("pause");

	return resultcode;
}

