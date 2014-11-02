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
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"

#include "kmUtility.h"
#include "kmProcessing.h"
#include "kmRegistration.h"
#include "kmVtkItkUtility.h"
#include "kmSegmentation.h"

using namespace km;

namespace km
{


	/************************************************************************/
	//                                                                       
	//形状训练
	//步骤：


	//(part 1: 预处理，对齐)
	//		1. 读取训练集中的灰度图像Iorig及对应分割mask图像Sliver
	//    2. 提取灰度图像Iorig人体部分Ibody，生成人体mask图像Sbody
	//    3. 横断面旋转角度纠正
	//    4. 根据body横断面平均面积，重新设置训练图像的spacing
	//    5. 对齐训练图像肝脏（依据肝脏最大面积切片位置）


	//(part 2: 生成点对应肝脏形状mesh)
	//    6. 将训练图像的参考图像与其他训练图像进行B样条配准
	//    7. 通过B样条变换，warp参考图像对应肝脏mesh，再根据warp后的mesh在训练图像上找到对应实际的肝脏表面点，refine mesh
	//    8. 输出
	//
	/************************************************************************/

	int trainShape( const char* origlistfile,
		const char* seglistfile,
		const char* referenceorigfile,
		const char* referencesegfile,
		const char* referenceshapefile,
		const char* referenceprobabilityfile,
		const char* outputdir,
		const unsigned int startdataindex = 0)
	{

		KM_DEBUG_INFO( "======================================" );
		KM_DEBUG_PRINT( "original image list", origlistfile );
		KM_DEBUG_PRINT( "liver mask list", seglistfile );
		KM_DEBUG_PRINT( "Reference original image", referenceorigfile );
		KM_DEBUG_PRINT( "Reference liver mask", referencesegfile );
		KM_DEBUG_PRINT( "Reference triangle mesh", referenceshapefile );
		KM_DEBUG_PRINT( "Reference probability image", referenceprobabilityfile );
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

		int numberOfGrayData = km::getDataList(origlistfile, origDataFileNamesList);
		int numberOfBinaryData = km::getDataList(seglistfile, segDataFileNamesList);

		if(numberOfGrayData != numberOfBinaryData)
		{
			std::cout<<"number of gray data "<<numberOfGrayData<<" does not match the number of segmented data "<<numberOfBinaryData<< "!"<<std::endl;
			return EXIT_FAILURE;
		}
		KM_DEBUG_PRINT( "Number Of Data:" , numberOfGrayData);

		ShortImageType::Pointer     refOrigImage = km::readImage<ShortImageType>( referenceorigfile );
		UCharImageType::Pointer     refSegImage  = km::readImage<UCharImageType>( referencesegfile );
		ShortImageType::Pointer     refLiverImage = km::maskImage<ShortImageType, UCharImageType>( refOrigImage, refSegImage );
		UCharImageType::PointType   refCentroid  = km::getCentroid<UCharImageType>( refSegImage,1,1 );
		ShortImageType::PointType   refOrgin     = refOrigImage->GetOrigin();
		ShortImageType::SpacingType refSpacing   = refOrigImage->GetSpacing();
		ShortImageType::SpacingType downSpacing;
		downSpacing.Fill( 2.0 );
		ShortImageType::Pointer coarseRefLiver = km::resampleImage<ShortImageType>( 
			refLiverImage, 
			downSpacing, 
			0 );

		std::stringstream shapeListFilename, geometryFilename, reflivermeshFilename;
		shapeListFilename << outputdir << "\\corresponded-shapesList.txt";
		geometryFilename  << outputdir << "\\geoImage.mha";
		reflivermeshFilename << outputdir << "\\refSimplexMesh.vtk";

		KM_DEBUG_INFO(shapeListFilename.str());
		KM_DEBUG_INFO(reflivermeshFilename.str());
		KM_DEBUG_INFO(geometryFilename.str());

		//参考图象对应稀疏三角面片数据
		SimplexMeshType::Pointer refTriangleMesh = km::readMesh<SimplexMeshType>( referenceshapefile );

		KM_DEBUG_INFO("Read reference mesh done.");

		//生成Simplex Mesh之后，用于建立形状模型
		SimplexMeshType::Pointer refSimplexMesh
			= km::triangleMeshToSimplexMesh<SimplexMeshType, SimplexMeshType>( refTriangleMesh );

		KM_DEBUG_INFO("Generate reference simplex mesh done.");

		km::writeMesh<SimplexMeshType>( reflivermeshFilename.str().c_str(), refSimplexMesh );

		GeometryImageType::Pointer geoImage = km::generateGeoImage<SimplexMeshType>( refSimplexMesh );
		km::writeImage<GeometryImageType>( geometryFilename.str().c_str(), geoImage );
		//km::writeSimplexMeshGeometryData<SimplexMeshType>( geometryFilename.str().c_str(), refSimplexMesh );

		//用以向corresponded-shapesList.txt写入形状文件名
		ofstream outputShapeListStream;
		outputShapeListStream.open( shapeListFilename.str().c_str() );

		char absolutePath[255];
		getcwd(absolutePath, 255);

		for (unsigned int i=startdataindex;i<numberOfGrayData;i++)
		{
			KM_DEBUG_PRINT("Start to align data: ", i+1);
			//读入图像
			KM_DEBUG_INFO("Reading images..");
			std::string origimagefile = origDataFileNamesList[i];
			std::string livermaskfile = segDataFileNamesList[i];
			KM_DEBUG_INFO(origimagefile);
			KM_DEBUG_INFO(livermaskfile);
			ShortImageType::Pointer origimage = km::readImage<ShortImageType>(origimagefile);
			UCharImageType::Pointer livermask = km::readImage<UCharImageType>(livermaskfile);
			livermask = km::binaryThresholdImage<UCharImageType, UCharImageType>( livermask, 1, 255, 1, 0 );
			ShortImageType::Pointer liverimage =km::maskImage<ShortImageType, UCharImageType>( origimage, livermask );
			liverimage = km::thresholdImage<ShortImageType>( liverimage, 1, 300, 0 );

			//km::writeImage<ShortImageType>( outputdir, "origimage", i, ".nii.gz", origimage );
			//km::writeImage<UCharImageType>( outputdir, "livermask", i, ".nii.gz", livermask );
			//km::writeImage<ShortImageType>( outputdir, "liverimage", i, ".nii.gz", liverimage );

			//统计最小灰度值
			double standardMinimum = -1024;
			double minValue, maxValue;
			km::calculateMinAndMax<ShortImageType>(origimage, minValue, maxValue);
			standardMinimum = minValue;

			ShortImageType::Pointer coarseLiver = km::resampleImage<ShortImageType>( 
				liverimage, 
				downSpacing, 
				0 );

			KM_DEBUG_PRINT( "Elastic registration: ", i+1 );

			typedef itk::BSplineTransform<double, 3, 3> BSplineTransformType;
			BSplineTransformType::Pointer bsplineTransform = BSplineTransformType::New();

			km::elasticRegistrationImageToImage<ShortImageType, ShortImageType, BSplineTransformType>(
				coarseRefLiver,
				coarseLiver,
				bsplineTransform,
				8,
				referenceprobabilityfile
				);
			KM_DEBUG_INFO( "Elastic registration done!" );

			SimplexMeshType::Pointer liverMesh = km::transformMesh<SimplexMeshType, BSplineTransformType>(
				refSimplexMesh,
				bsplineTransform
				);

			vtkSmartPointer<vtkPolyData> polydata = km::mesh2PolyData<SimplexMeshType>(liverMesh);
			polydata = km::smoothPolyData( polydata, 100 );
			km::copyPointsFromPolydataToMesh<SimplexMeshType>( polydata, liverMesh );
			km::loadSimplexMeshGeometryData<SimplexMeshType>(geoImage, liverMesh);

			km::writeMesh<SimplexMeshType>( outputdir, "elasticWarpedMesh", i, ".vtk", liverMesh );

			SimplexMeshType::Pointer accurateMesh = km::generateMeshFromBinary<UCharImageType, SimplexMeshType>( livermask );
			polydata = km::mesh2PolyData<SimplexMeshType>(accurateMesh);
			polydata = km::smoothPolyData( polydata, 100 );
			polydata = km::decimatePolydata(polydata, 5000.0/accurateMesh->GetNumberOfPoints());
			polydata = km::smoothPolyData( polydata, 100 );
			accurateMesh = km::polyData2Mesh<SimplexMeshType>( polydata );

			km::writeMesh<SimplexMeshType>( outputdir, "accurateMesh", i, ".vtk", accurateMesh );

			KM_DEBUG_INFO( "Deformation fitting..." );
			km::deformSegSimplexMeshByICP<SimplexMeshType>(
				accurateMesh, 
				liverMesh);

			//Smooth and automatic adaption.
			km::adaptMesh<SimplexMeshType>( liverMesh, 0.05, 50 );

			km::writeMesh<SimplexMeshType>( outputdir, "livermesh", i, ".vtk", liverMesh );

			outputShapeListStream << absolutePath<< "\\" << outputdir << "\\livermesh" << "." << i+1 << ".vtk" << std::endl;

			KM_DEBUG_PRINT( "Shape generating done: ", i+1 );
		}

		return EXIT_SUCCESS;
	}

}

#endif