
#include <sstream>
#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include "kmShapeAlignmentTraining.h"

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
		std::cerr<<"Usage:   origDataList segDataList shapeDataList referenceShape refGeometryData outputDir startDataIndex"<<std::endl;
		return EXIT_FAILURE;
	}
	itk::TimeProbesCollectorBase chronometer;
	chronometer.Start( "shapetraining" );

	const char* origlistfile = argv[1];
	const char* seglistfile  = argv[2];
	const char* shapelistfile  = argv[3];
	const char* referenceshapefile = argv[4];
	const char* refgeodatafile = argv[5];
	const char* outputdir = argv[6];

	unsigned int startdataindex = 1;
	if (argc>8)
	{
		startdataindex = atoi(argv[7]);
	}

	startdataindex -= 1; //�ĳ�Cϰ�ߣ���0��ʼ

	int resultcode = km::alignData(
		origlistfile, 
		seglistfile, 
		shapelistfile,
		referenceshapefile,
		refgeodatafile, 
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

