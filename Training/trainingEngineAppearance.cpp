
#include <sstream>

#include "kmAppearanceTraining.h"

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
	if( argc < 5 )
	{
		std::cerr<<"Usage:   origDataList meshList referenceGeometry OutputDir [Profile Type: INSIDE/BOUNDARY]"<<std::endl;
		return EXIT_FAILURE;
	}

	const char* origlistfile = argv[1];
	const char* meshlistfile  = argv[2];
	const char* referencegeometry = argv[3];
	const char* outputdir = argv[4];

	itk::TimeProbesCollectorBase chronometer;
	itk::MemoryProbesCollectorBase memorymeter;
	chronometer.Start( "appearancetraining" );
	memorymeter.Start( "appearancetraining" );

	//For liver tissue classification
	int resultcode = km::trainAppearance(
		origlistfile 
		,meshlistfile
		,referencegeometry
		,outputdir
		/*,category*/);

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

	return resultcode;
}

