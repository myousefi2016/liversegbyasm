#include <iostream>
#include <string>

#include "itkDefaultDynamicMeshTraits.h"
#include "itkSimplexMesh.h"
#include "Representers/ITK/itkSimplexMeshRepresenter.h"
#include "statismo_ITK/itkStatisticalModel.h"

#include "kmUtility.h"

int main(int argc, char* argv[])
{
	//const char* SSMFile = argv[1];
	const char* SSMFile = "D:\\Workspace\\ASM\\projects\\LiverSegbyASM\\experiments\\buidModel\\shape\\output_20140815\\StatisticalShapeModel_20140815.h5";

	std::cout<<"test"+3<<std::endl;
	
	std::string str = argv[1];

	const int Dimension = 3;
	typedef double MeshPixelType;
	typedef itk::SimplexMeshRepresenter<MeshPixelType, Dimension> RepresenterType;
	typedef itk::StatisticalModel<RepresenterType>                StatisticalModelType;
	typedef RepresenterType::MeshType                             MeshType;
	
	StatisticalModelType::Pointer model = StatisticalModelType::New();
	model->Load( SSMFile );
	
	return 0;
}