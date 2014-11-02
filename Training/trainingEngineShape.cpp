
#include <sstream>
#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include "kmShapeTraining.h"

/************************************************************************
作用: 训练数据

输入:
1: 数据列表文本. 文本内容形式: xxx/xxx/xxx1.mhd xxx/xxx/xxx2.mhd
...
2: 用于对齐的参考数据在数据列表中的位置 从1开始,比如xxx1.mha为参考数据,则参数为1
3: 输出文件夹
4: 输出数据文件名前缀,比如输入为"Liver", 则最终输出数据为Liver.1.vtk, Liver.2.vtk...

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

	startdataindex -= 1; //改成C习惯，从0开始

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

