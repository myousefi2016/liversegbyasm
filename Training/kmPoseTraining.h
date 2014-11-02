#ifndef __kmPoseTraining_h
#define __kmPoseTraining_h

#include <iostream>
#include <fstream>
#include  <direct.h> 

#include "itkMatrix.h"

#include "itkImage.h"
#include "itkEuler3DTransform.h"

#include "itkEuler3DTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkScaleVersor3DTransform.h"
#include "itkSimilarity3DTransform.h"
#include "itkCenteredTransformInitializer.h"

#include "itkDefaultDynamicMeshTraits.h"
#include "itkMesh.h"
#include "itkSimplexMesh.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"

#include "kmUtility.h"
#include "kmProcessing.h"
#include "kmVtkItkUtility.h"
#include "kmSegmentation.h"
#include "kmRegistration.h"

#define  PI 3.14159265359

using namespace std;

namespace km
{
	int trainPose( const char* origlistdata,
		const char* seglistdata,
		const char* referenceorigfile,
		const char* referencesegfile,
		const char* outputdir)
	{
		const char* origDataListFile = origlistdata;
		const char* segDataListFile = seglistdata;
		//int referenceDataIndex = refindex;
		const char* outputDataDir = outputdir;

		//referenceDataIndex -= 1; //改成C习惯的索引
		const unsigned int Dimension = 3;
		typedef itk::Image<signed short, Dimension>    ShortImageType;
		typedef itk::Image<unsigned char, Dimension>   UCharImageType;
		typedef itk::Image<float, Dimension>           FloatImageType;
		typedef itk::SimplexMesh<float, Dimension>     MeshType;

		//////////////////////////////////////////////////////////////////////////
		//导入训练数据文件
		KM_DEBUG_INFO("Import Datalist");
		typedef std::vector<std::string>    StringVectorType;

		StringVectorType origDataFileNamesList;
		StringVectorType segDataFileNamesList;

		int numberOfGrayData = km::getDataList(origDataListFile, origDataFileNamesList);
		int numberOfBinaryData = km::getDataList(segDataListFile, segDataFileNamesList);

		if(numberOfGrayData != numberOfBinaryData)
		{
			std::cout<<"number of gray data "<<numberOfGrayData<<" does not match the number of segmented data "<<numberOfBinaryData<< "!"<<std::endl;
			return EXIT_FAILURE;
		}
		KM_DEBUG_PRINT( "Number Of Data:" , numberOfGrayData);

		ShortImageType::Pointer     refOrigImage  = km::readImage<ShortImageType>( referenceorigfile );
		UCharImageType::Pointer     refLiverMask   = km::readImage<UCharImageType>( referencesegfile );

		UCharImageType::SpacingType downspac;
		downspac.Fill( 2 );

		UCharImageType::Pointer downRefLiverMask = km::resampleImage<UCharImageType>( refLiverMask, downspac, 0, 1 );
		downRefLiverMask = km::binaryThresholdImage<UCharImageType, UCharImageType>( downRefLiverMask, 1, 255, 100, 0 );

		//MeshType::Pointer refMesh = km::generateMeshFromBinary<UCharImageType, MeshType>( refLiverMask );
		//refMesh = km::decimateMesh<MeshType>( refMesh );

		////MeshType::PointType refCentroid = km::getCentroid<MeshType>( refMesh );
		//UCharImageType::PointType refCentroid = km::getCentroid<UCharImageType>(refLiverMask);

		std::stringstream origListFilename, segListFilename, rigidTransformListFilename;
		origListFilename << outputDataDir << "\\aligned-origList.txt";
		segListFilename  << outputDataDir << "\\aligned-segList.txt";
		rigidTransformListFilename << outputDataDir << "\\aligned-transformList.txt";

		ofstream outputOrigListStream, outputSegListStream, outputTransformListStream;
		outputOrigListStream.open( origListFilename.str().c_str() );
		outputSegListStream.open( segListFilename.str().c_str() );
		outputTransformListStream.open( rigidTransformListFilename.str().c_str() );

		char absolutePath[255];
		getcwd(absolutePath, 255);

		//typedef itk::VersorRigid3DTransform<double> RigidTransformType;
		typedef itk::Similarity3DTransform<double> RigidTransformType;
		std::vector<double> opt_scales;
		opt_scales.resize( RigidTransformType::ParametersDimension );

		opt_scales[0]   =  1e3;
		opt_scales[1]   =  1e3;
		opt_scales[2]   =  1.0;
		opt_scales[3]   =  1e-3;
		opt_scales[4]   =  1e-3;
		opt_scales[5]   =  1e-3;
		opt_scales[6]   =  1.0;

		/*******************************************
		//循环处理每个训练数据
		*******************************************/
		for (unsigned int i=0;i<numberOfGrayData;i++)
		{
			try
			{
				KM_DEBUG_PRINT("Start to align data: ", i+1);
				//读入图像
				KM_DEBUG_INFO("Reading images..");
				std::string origimagefile = origDataFileNamesList[i];
				std::string livermaskfile = segDataFileNamesList[i];

				KM_DEBUG_PRINT("Original Image: ", origimagefile);
				KM_DEBUG_PRINT("Liver Mask: ", livermaskfile);

				ShortImageType::Pointer origimage = km::readImage<ShortImageType>(origimagefile);
				UCharImageType::Pointer livermask = km::readImage<UCharImageType>(livermaskfile);
				livermask = km::binaryThresholdImage<UCharImageType, UCharImageType>( livermask, 1, 255, 1, 0 );

				KM_DEBUG_INFO("Read original image & liver mask done!");

				KM_DEBUG_INFO("Down-sample images to speed up..");
				UCharImageType::Pointer downLiverMask = km::resampleImage<UCharImageType>( livermask, downspac, 0, 1 );
				downLiverMask    = km::binaryThresholdImage<UCharImageType, UCharImageType>( downLiverMask, 1, 255, 100, 0 );

				//计算最低灰度
				double minValue, maxValue;
				km::calculateMinAndMax<ShortImageType>(origimage, minValue, maxValue);

				RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
				//rigidTransform->SetIdentity();
				//rigidTransform->SetCenter( refCentroid );
				//rigidTransform->SetOffset( centroidOffset );
				typedef itk::CenteredTransformInitializer<RigidTransformType, UCharImageType, UCharImageType> TransformInitializerType;
				TransformInitializerType::Pointer transformInitializer = TransformInitializerType::New();
				transformInitializer->SetFixedImage(downRefLiverMask);
				transformInitializer->SetMovingImage(downLiverMask);
				transformInitializer->SetTransform(rigidTransform);
				transformInitializer->MomentsOn();
				transformInitializer->InitializeTransform();

				std::cout<<"[BEFORE]"<<rigidTransform->GetParameters()<<std::endl;
				std::cout<<"[BEFORE]"<<rigidTransform->GetFixedParameters()<<std::endl;

				//km::eulerRegistrationMeshToMesh<MeshType, RigidTransformType>(
				//	refMesh,
				//	livermesh,
				//	rigidTransform,
				//	opt_scales);

				KM_DEBUG_INFO("Start to procrustes aligment..");
				km::ProcrustesAlignment<UCharImageType, UCharImageType, RigidTransformType>( 
					downRefLiverMask,
					downLiverMask,
					rigidTransform,
					opt_scales);

				RigidTransformType::Pointer rigidTransformInversed = RigidTransformType::New();
				rigidTransform->GetInverse( rigidTransformInversed );

				origimage = km::transformImageByReference<ShortImageType, ShortImageType, RigidTransformType>(
					origimage,
					refOrigImage,
					rigidTransform,
					minValue);
				livermask = km::transformImageByReference<UCharImageType, ShortImageType, RigidTransformType>(
					livermask,
					refOrigImage,
					rigidTransform,
					0,
					1);

				std::cout<<"[AFTER]"<<rigidTransform->GetParameters()<<std::endl;
				std::cout<<"[AFTER]"<<rigidTransform->GetFixedParameters()<<std::endl;

				km::writeImage<ShortImageType>( outputDataDir, "origImage-rigidAligned", i, ".nii.gz", origimage );
				km::writeImage<UCharImageType>( outputDataDir, "liverMask-rigidAligned", i, ".nii.gz", livermask );
				km::writeTransform<RigidTransformType>( outputDataDir, "rigidTransform", i, ".tfm", rigidTransform );
				outputOrigListStream << absolutePath << "\\" << outputDataDir << "\\origImage-rigidAligned" << "." << i+1 << ".nii.gz" << std::endl;
				outputSegListStream  << absolutePath << "\\" << outputDataDir << "\\liverMask-rigidAligned" << "." << i+1 << ".nii.gz" << std::endl;
				outputTransformListStream << absolutePath << "\\" << outputDataDir << "\\rigidTransform"    << "." << i+1 << ".tfm" << std::endl;

				KM_DEBUG_PRINT( "Rigid alignment done: ", i+1 );
			}
			catch(itk::ExceptionObject & e)
			{
				std::cerr<<"Exception thrown! "<<e<<std::endl;
			}
			catch(...)
			{
				std::cerr<<"Exception thrown! Unknown exception!"<<std::endl;
			}

		}

		outputOrigListStream.close();
		outputSegListStream.close();
		outputTransformListStream.close();

		return EXIT_SUCCESS;
	}
}

#endif