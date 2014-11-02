
#include <sstream>
#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include "kmShapeTraining.h"

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
	if( argc < 8 )
	{
		std::cerr<<"Usage:   origDataList segDataList referenceOrig referenceSeg referenceShape referenceProbability outputDir startDataIndex"<<std::endl;
		return EXIT_FAILURE;
	}
	itk::TimeProbesCollectorBase chronometer;
	chronometer.Start( "shapetraining" );

	const char* origlistfile = argv[1];
	const char* seglistfile  = argv[2];
	const char* referenceorigfile = argv[3];
	const char* referencesegfile  = argv[4];
	const char* referenceshapefile = argv[5];
	const char* referenceprobabilityfile = argv[6];
	const char* outputdir = argv[7];

	unsigned int startdataindex = 1;
	if (argc>8)
	{
		startdataindex = atoi(argv[8]);
	}

	startdataindex -= 1; //�ĳ�Cϰ�ߣ���0��ʼ

	int resultcode = km::trainShape(origlistfile, 
																	seglistfile, 
																	referenceorigfile, 
																	referencesegfile, 
																	referenceshapefile,
																	referenceprobabilityfile, 
																	outputdir,
																	startdataindex);

	if(resultcode==EXIT_SUCCESS)
	{
		KM_DEBUG_INFO( "Train shape done without error!" );
	}
	else
	{
		KM_DEBUG_INFO( "Train shape done with error!" );
	}

	chronometer.Stop( "shapetraining" );
	chronometer.Report( std::cout );

	return resultcode;
}

