#ifndef __kmShapeTraining_h
#define __kmShapeTraining_h

#include <iostream>
#include <fstream>
#include  <direct.h> 

#include "itkImage.h"
#include "itkEuler3DTransform.h"

#include <itkScaleVersor3DTransform.h>
#include "itkSimilarity3DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkBSplineTransform.h"

#include "itkDefaultDynamicMeshTraits.h"
#include "itkMesh.h"
#include "itkSimplexMesh.h"

#include "kmUtility.h"
#include "kmRegistration.h"
#include "kmVtkItkUtility.h"

using namespace km;

namespace km
{


	/************************************************************************/
	//                                                                       
	//形状训练
	//步骤：


	//(part 1: 预处理，对齐)
	//
	/************************************************************************/

	int alignData( const char* origlistfile,
		const char* seglistfile,
		const char* shapelistfile,
		const char* referenceshapefile,
		const char* refgeodatafile,
		const char* outputdir,
		const unsigned int startdataindex = 0)
	{

		KM_DEBUG_INFO( "======================================" );
		KM_DEBUG_PRINT( "original image list", origlistfile );
		KM_DEBUG_PRINT( "liver mask list", seglistfile );
		KM_DEBUG_PRINT( "liver mesh list", shapelistfile );
		KM_DEBUG_PRINT( "reference triangle mesh", referenceshapefile );
		KM_DEBUG_PRINT( "reference geometry data", refgeodatafile );
		KM_DEBUG_PRINT( "Output directory", outputdir );

		const unsigned int Dimension = 3;
		typedef itk::Image<signed short, Dimension>    ShortImageType;
		typedef itk::Image<unsigned char, Dimension>   UCharImageType;
		typedef itk::Image<float, Dimension>           FloatImageType;
		typedef itk::Image<signed short, Dimension-1>  ShortSliceType;
		typedef itk::Image<unsigned char, Dimension-1> UCharSliceType;
		typedef double MeshPixelType;
		typedef itk::DefaultDynamicMeshTraits<MeshPixelType, Dimension, Dimension,MeshPixelType,MeshPixelType> MeshTraitsType;
		typedef itk::SimplexMesh<MeshPixelType,Dimension, MeshTraitsType>                                      SimplexMeshType;

		//////////////////////////////////////////////////////////////////////////
		//导入训练数据文件
		KM_DEBUG_INFO("Import Datalist");
		typedef std::vector<std::string>    StringVectorType;

		StringVectorType origDataFileNamesList;
		StringVectorType segDataFileNamesList;
		StringVectorType shapeDataFileNameList;

		int numberOfGrayData = km::getDataList(origlistfile, origDataFileNamesList);
		int numberOfBinaryData = km::getDataList(seglistfile, segDataFileNamesList);
		int numberOfShapeData = km::getDataList(shapelistfile, shapeDataFileNameList);

		int numberOfData = std::min(numberOfGrayData, numberOfBinaryData);
		numberOfData = std::min(numberOfShapeData, numberOfData);

		KM_DEBUG_PRINT( "Number Of Data:" , numberOfData);

		std::stringstream alignedShapeListFilename, alignedOrigListFilename, alignedSegListFilename, geometryFilename;
		alignedShapeListFilename << outputdir << "\\aligned-shapesList.txt";
		alignedOrigListFilename << outputdir << "\\aligned-origList.txt";
		alignedSegListFilename << outputdir << "\\aligned-segList.txt";
		geometryFilename  << outputdir << "\\geoImage.mha";

		KM_DEBUG_INFO(alignedShapeListFilename.str());
		KM_DEBUG_INFO(alignedOrigListFilename.str());
		KM_DEBUG_INFO(alignedSegListFilename.str());
		KM_DEBUG_INFO(geometryFilename.str());

		//参考图象对应稀疏三角面片数据
		SimplexMeshType::Pointer refShapeMesh = km::readMesh<SimplexMeshType>( referenceshapefile );
		
		KM_DEBUG_INFO("Copy geometry data.");
		km::GeometryImageType::Pointer geoImage = km::readImage<km::GeometryImageType>( refgeodatafile );
		km::writeImage<km::GeometryImageType>( geometryFilename.str().c_str(), geoImage );

		ofstream outputShapeListStream;
		outputShapeListStream.open( alignedShapeListFilename.str().c_str() );

		ofstream outputOrigListStream;
		outputOrigListStream.open( alignedOrigListFilename.str().c_str() );

		ofstream outputSegListStream;
		outputSegListStream.open( alignedSegListFilename.str().c_str() );

		char absolutePath[255];
		getcwd(absolutePath, 255);

		typedef itk::Similarity3DTransform<double> RigidTransformType;
		std::vector<double> opt_scales;
		opt_scales.resize( RigidTransformType::ParametersDimension );
		opt_scales[0]   =  1.0;
		opt_scales[1]   =  1.0;
		opt_scales[2]   =  1.0;
		opt_scales[3]   =  1e-3;
		opt_scales[4]   =  1e-3;
		opt_scales[5]   =  1e-3;
		opt_scales[6]   =  1.0;

		for (unsigned int i=startdataindex;i<numberOfData;i++)
		{
			KM_DEBUG_PRINT("Start to align data: ", i+1);
			//读入图像
			KM_DEBUG_INFO("Reading images..");
			std::string origimagefile = origDataFileNamesList[i];
			std::string livermaskfile = segDataFileNamesList[i];
			std::string livermeshfile = shapeDataFileNameList[i];
			KM_DEBUG_INFO(origimagefile);
			KM_DEBUG_INFO(livermaskfile);
			ShortImageType::Pointer origimage = km::readImage<ShortImageType>(origimagefile);
			UCharImageType::Pointer livermask = km::readImage<UCharImageType>(livermaskfile);
			SimplexMeshType::Pointer livermesh = km::readMesh<SimplexMeshType>(livermeshfile);

			//统计最小灰度值
			double standardMinimum = -1024;
			double minValue, maxValue;
			km::calculateMinAndMax<ShortImageType>(origimage, minValue, maxValue);
			standardMinimum = minValue;

			RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
			SimplexMeshType::PointType refCentroid = km::getMeshCentroid<SimplexMeshType>( refShapeMesh );
			rigidTransform->SetIdentity();
			rigidTransform->SetCenter( refCentroid );

			km::alignMesh<SimplexMeshType, RigidTransformType>( refShapeMesh, livermesh, rigidTransform, opt_scales );
			livermesh = km::transformMesh<SimplexMeshType, RigidTransformType>( livermesh, rigidTransform );

			RigidTransformType::Pointer inversedRigidTransform = RigidTransformType::New();
			rigidTransform->GetInverse( inversedRigidTransform );

			origimage = km::transformImage<ShortImageType, RigidTransformType>( origimage, origimage, inversedRigidTransform, standardMinimum );
			livermask = km::transformImage<UCharImageType, RigidTransformType>( livermask, livermask, inversedRigidTransform, 0 );

			km::writeMesh<SimplexMeshType>( outputdir, "liverMesh-rigidAligned", i, ".vtk", livermesh );
			km::writeImage<ShortImageType>( outputdir, "origImage-rigidAligned", i, ".nii.gz", origimage );
			km::writeImage<UCharImageType>( outputdir, "liverMask-rigidAligned", i, ".nii.gz", livermask );

			outputShapeListStream << absolutePath<< "\\" << outputdir << "\\liverMesh-rigidAligned" << "." << i+1 << ".vtk" << std::endl;
			outputOrigListStream << absolutePath<< "\\" << outputdir << "\\origImage-rigidAligned" << "." << i+1 << ".nii.gz" << std::endl;
			outputSegListStream << absolutePath<< "\\" << outputdir << "\\liverMask-rigidAligned" << "." << i+1 << ".nii.gz" << std::endl;

			KM_DEBUG_PRINT( "Shape generating done: ", i+1 );
		}

		return EXIT_SUCCESS;
	}

}

#endif