/*
* This file is part of the statismo library.
*
* Author: Marcel Luethi (marcel.luethi@unibas.ch)
*
* Copyright (c) 2011 University of Basel
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*
* Redistributions of source code must retain the above copyright notice,
* this list of conditions and the following disclaimer.
*
* Redistributions in binary form must reproduce the above copyright
* notice, this list of conditions and the following disclaimer in the
* documentation and/or other materials provided with the distribution.
*
* Neither the name of the project's author nor the names of its
* contributors may be used to endorse or promote products derived from
* this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
* TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
* LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*/


#include <itkCovariantVector.h>
#include "Representers/ITK/itkVectorSimplexMeshRepresenter.h"
#include "statismo_ITK/itkStatisticalModel.h"
#include "statismo_ITK/itkPCAModelBuilder.h"
#include "statismo_ITK/itkDataManager.h"
#include "itkDirectory.h"
#include "itkSimplexMesh.h"
#include "itkMeshFileWriter.h"
#include "itkMeshFileReader.h"
#include <sys/types.h>
#include <errno.h>
#include <iostream>

#include "kmUtility.h"

/*
* This example shows the ITK Wrapping of statismo can be used to build a shape model.
*/

void buildAppearanceModel(const char* datalist, const int referenceindex, const char* modelname) 
{
	const unsigned Dimensions = 3;
	typedef itk::CovariantVector< double, 3 >                        PixelType;
	typedef itk::VectorSimplexMeshRepresenter<PixelType, Dimensions> RepresenterType;
	typedef RepresenterType::MeshType                                MeshType;
	typedef itk::PCAModelBuilder<RepresenterType>                    ModelBuilderType;
	typedef itk::StatisticalModel<RepresenterType>                   StatisticalModelType;
	typedef std::vector<std::string>                                 StringVectorType;
	typedef itk::DataManager<RepresenterType>                        DataManagerType;
	typedef itk::MeshFileReader<MeshType>                            MeshReaderType;

	//////////////////////////////////////////////////////////////////////////
	//导入训练数据文件
	KM_DEBUG_INFO("Import Datalist");
	typedef std::vector<std::string>    StringVectorType;
	StringVectorType filenames;
	int numberOfshapes = km::getDataList(datalist, filenames);

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

	MeshType::Pointer meanAppearance = model->DrawMean();
	km::writeMesh<MeshType>( "meanAppearance.vtk", meanAppearance );
}

int main(int argc, char* argv[]) {

	if (argc < 4) {
		std::cout << "usage " << argv[0] << " meshlist referenceindex modelname" << std::endl;
		exit(-1);
	}

	const char* meshlist = argv[1];
	int referenceindex = atoi(argv[2]);
	const char* modelname = argv[3];

	referenceindex -= 1;

	buildAppearanceModel(meshlist, referenceindex, modelname);

	std::cout << "Model building is completed successfully." << std::endl;
}

