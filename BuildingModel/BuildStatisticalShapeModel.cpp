
//#include "Representers/ITK/itkMeshRepresenter.h"
#include "Representers/ITK/itkSimplexMeshRepresenter.h"
#include "statismo_ITK/itkStatisticalModel.h"
#include "statismo_ITK/itkPCAModelBuilder.h"
#include "statismo_ITK/itkDataManager.h"
#include "itkDirectory.h"
#include "itkMesh.h"
#include "itkMeshFileWriter.h"
#include "itkMeshFileReader.h"
#include <sys/types.h>
#include <errno.h>
#include <iostream>
#include <string>

#include "kmCommon.h"

/*
* This example shows the ITK Wrapping of statismo can be used to build a shape model.
*/

void buildShapeModel(const char* datalist, const int referenceindex, const char* modelname) 
{
	const unsigned Dimensions = 3;
	typedef float                                              PixelType;
	typedef itk::SimplexMeshRepresenter<PixelType, Dimensions> RepresenterType;
	typedef RepresenterType::MeshType                          MeshType;
	typedef itk::PCAModelBuilder<RepresenterType>              ModelBuilderType;
	typedef itk::StatisticalModel<RepresenterType>             StatisticalModelType;
	typedef std::vector<std::string>                           StringVectorType;
	typedef itk::DataManager<RepresenterType>                  DataManagerType;
	typedef itk::MeshFileReader<MeshType>                      MeshReaderType;

	//////////////////////////////////////////////////////////////////////////
	//导入训练数据文件
	KM_DEBUG_INFO("Import Datalist");
	typedef std::vector<std::string>    StringVectorType;
	StringVectorType filenames;
	int numberOfshapes = km::getDataList(datalist, filenames);
	for (int i=0;i<numberOfshapes;i++){
		std::cout<<filenames[i]<<std::endl;
	}

	KM_DEBUG_PRINT( "Number of shapes: ", numberOfshapes );
	RepresenterType::Pointer representer = RepresenterType::New();
	MeshReaderType::Pointer refReader = MeshReaderType::New();
	refReader->SetFileName( filenames[referenceindex] );
	refReader->Update();
	representer->SetReference(refReader->GetOutput());
	DataManagerType::Pointer dataManager = DataManagerType::New();
	dataManager->SetRepresenter(representer);

	for (StringVectorType::const_iterator it = filenames.begin(); it != filenames.end(); it++) 
	{
		const std::string fullpath = *it;//(std::string(dir) + "/") + *it;
		MeshReaderType::Pointer reader = MeshReaderType::New();
		reader->SetFileName(fullpath.c_str());
		reader->Update();
		MeshType::Pointer mesh = reader->GetOutput();
		dataManager->AddDataset(mesh, fullpath.c_str());
	}

	ModelBuilderType::Pointer pcaModelBuilder = ModelBuilderType::New();
	StatisticalModelType::Pointer model = pcaModelBuilder->BuildNewModel(dataManager->GetSampleDataStructure(), 0);
	model->Save(modelname);
}

int main(int argc, char* argv[]) {

	if (argc < 4) {
		std::cout << "usage " << argv[0] << " shapelist referenceindex modelname" << std::endl;
		exit(-1);
	}
	const char* shapelist = argv[1];
	int referenceindex = atoi(argv[2]);
	const char* modelname = argv[3];
	referenceindex -= 1;
	buildShapeModel(shapelist, referenceindex, modelname);
	std::cout << "Model building is completed successfully." << std::endl;
}

