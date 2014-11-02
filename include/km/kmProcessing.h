#ifndef __kmProcessing_h
#define __kmProcessing_h

#include <vector>
#include <iostream>
#include <fstream>

#include "itkImage.h"
#include <itkFixedArray.h>
#include <itkVector.h>
#include <itkFlipImageFilter.h>
#include "itkChangeInformationImageFilter.h"
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageDuplicator.h>
#include <itkCenteredRigid2DTransform.h>
#include <itkTranslationTransform.h>
#include "itkMinimumMaximumImageCalculator.h"
#include "itkTranslationTransform.h"
#include "itkNormalizeImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkCropImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include <itkConnectedComponentImageFilter.h>
#include <itkImageToHistogramFilter.h>

#include "itkTranslationTransform.h"

#include <itkConfidenceConnectedImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>

#include <itkBinaryFillholeImageFilter.h>
#include <itkNeighborhoodConnectedImageFilter.h>

#include <itkBinaryThresholdImageFunction.h>
#include <itkFloodFilledImageFunctionConditionalIterator.h>
#include <itkShapedFloodFilledImageFunctionConditionalIterator.h>

#include <itkBinaryImageToLabelMapFilter.h>
#include "itkLabelMapToLabelImageFilter.h"

#include <itkScaleVersor3DTransform.h>

#include "itkOrImageFilter.h"

#include "kmUtility.h"

namespace km
{
	/************************************************************************/
	//(新)预处理API                                                            
	/************************************************************************/


	char* dir(char* dir)
	{
		char* resultStr;
		std::string s(dir);
		int dirLen = s.length();
		if(s[dirLen-1] != '\\'){
			s.push_back('\\');
		}
		resultStr = (char*)malloc(sizeof(char)*s.length());
		return strcpy(resultStr, s.c_str());
	}

	template<class ImageType>
	void
		shiftMinimum(typename ImageType::Pointer & image, double targetMin = -1024)
	{
		double minval, maxval;
		km::calculateMinAndMax<ShortImageType>( image, minval, maxval );

		image = km::shiftScale<ImageType>( image, targetMin - minval, 1.0 );
	}

	template<class InputImageType, class AtlasImageType>
	typename InputImageType::Pointer
		extractRoiByTemplateMatching(typename InputImageType* inputImage, typename AtlasImageType* atlasImage)
	{
		InputImageType::Pointer inputTmp = km::resampleImage<InputImageType>( inputImage, 3.0, -1024.0 );
		AtlasImageType::Pointer atlasTmp = km::resampleImage<AtlasImageType>( atlasImage, 3.0, -1024.0 );
		
		inputTmp = km::thresholdImage<InputImageType>( inputTmp, -1024, 400, 400 );
		atlasTmp = km::thresholdImage<AtlasImageType>( atlasTmp, -1024, 400, 400 );

		typedef itk::TranslationTransform<double, 3> TransformType;
		TransformType::Pointer transformInversed = TransformType::New();
		transformInversed->SetIdentity();
		km::histogramRigidRegistration<InputImageType, AtlasImageType, TransformType>( atlasTmp, inputTmp, transformInversed );

		TransformType::Pointer transform = TransformType::New();
		transformInversed->GetInverse( transform );

		AtlasImageType::Pointer alignedAtlas = km::transformImageByReference<AtlasImageType, InputImageType, TransformType>( atlasImage, inputImage, transform, -1024, 0 );
		//km::writeImage<AtlasImageType>( "alignedAtlas.nii.gz", alignedAtlas );

		TransformType::OutputVectorType offset = transform->GetOffset();

		InputImageType::IndexType atlasStartIdx, atlasUpperIdx, inputStartIdx, inputUpperIdx;
		InputImageType::PointType atlasStartPt, atlasUpperPt, inputStartPt, inputUpperPt;
		atlasStartIdx = atlasImage->GetLargestPossibleRegion().GetIndex();
		atlasUpperIdx = atlasImage->GetLargestPossibleRegion().GetUpperIndex();
		inputStartIdx = inputImage->GetLargestPossibleRegion().GetIndex();
		inputUpperIdx = inputImage->GetLargestPossibleRegion().GetUpperIndex();

		atlasImage->TransformIndexToPhysicalPoint( atlasStartIdx, atlasStartPt );
		atlasImage->TransformIndexToPhysicalPoint( atlasUpperIdx, atlasUpperPt );

		inputImage->TransformIndexToPhysicalPoint( inputStartIdx, inputStartPt );
		inputImage->TransformIndexToPhysicalPoint( inputUpperIdx, inputUpperPt );

		double paddedAtlasStartZPos = atlasStartPt[2]-offset[2];//-80;
		double paddedAtlasUpperZPos = atlasUpperPt[2]-offset[2];//+80;

		inputStartPt[2] = std::max(paddedAtlasStartZPos, inputStartPt[2]);
		inputUpperPt[2] = std::min(paddedAtlasUpperZPos, inputUpperPt[2]);

		inputImage->TransformPhysicalPointToIndex( inputStartPt, inputStartIdx );
		inputImage->TransformPhysicalPointToIndex( inputUpperPt, inputUpperIdx );

		InputImageType::RegionType inputRoi = inputImage->GetLargestPossibleRegion();
		inputRoi.SetIndex( inputStartIdx );
		inputRoi.SetUpperIndex( inputUpperIdx );

		InputImageType::Pointer extractedRoi = km::extractByBoundRegion<InputImageType>(inputImage, inputRoi);

		//If ROI is too small, we will pad the roi.
		if (inputUpperPt[2] < paddedAtlasUpperZPos || inputStartPt[2] > paddedAtlasStartZPos)
		{
			double minval, maxval;
			km::calculateMinAndMax<InputImageType>( extractedRoi, minval, maxval );

			InputImageType::SizeType padUpperBound, padLowerBound;
			padUpperBound.Fill(0);
			padLowerBound.Fill(0);

			padUpperBound[2] = std::max( paddedAtlasUpperZPos - inputUpperPt[2], 0.0 ) / (extractedRoi->GetSpacing())[2];
			//padLowerBound[2] = std::max( inputStartPt[2] - paddedAtlasStartZPos, 0.0 ) / (extractedRoi->GetSpacing())[2];

			extractedRoi = km::padImage<InputImageType>( extractedRoi, padLowerBound, padUpperBound, minval );
		}

		return extractedRoi;
	}

	//Extract single connected component from a binary image.
	template<class ImageType>
	typename ImageType::Pointer
		extractMaxConnectedComponent(typename ImageType* image)
	{
		typedef itk::BinaryImageToLabelMapFilter<ImageType> BinaryImageToLabelMapFilterType;
		BinaryImageToLabelMapFilterType::Pointer binaryImageToLabelMapFilter = BinaryImageToLabelMapFilterType::New();
		binaryImageToLabelMapFilter->SetInput(image);
		binaryImageToLabelMapFilter->SetFullyConnected(false);
		binaryImageToLabelMapFilter->SetInputForegroundValue(1);
		binaryImageToLabelMapFilter->SetOutputBackgroundValue(0);
		binaryImageToLabelMapFilter->Update();

		std::cout<<"Number Of Label Objects: "<<binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects()<<std::endl;

		itk::SizeValueType maxSize = 0;
		unsigned int maxLabel = 0;
		for(unsigned int i = 0; i < binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++)
		{
			// Get the ith region
			BinaryImageToLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToLabelMapFilter->GetOutput()->GetNthLabelObject(i);

			if (labelObject->Size() > maxSize)
			{
				maxLabel = labelObject->GetLabel();
				maxSize = labelObject->Size();
			}
		}

		std::cout<<"Max size: "<<maxSize<<", Max label: "<<maxLabel<<std::endl;

		typedef itk::LabelMapToLabelImageFilter<BinaryImageToLabelMapFilterType::OutputImageType, ImageType> LabelMapToLabelImageFilterType;
		LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
		labelMapToLabelImageFilter->SetInput(binaryImageToLabelMapFilter->GetOutput());
		labelMapToLabelImageFilter->Update();

		//km::writeImage<ImageType>( "label2.nii.gz", labelMapToLabelImageFilter->GetOutput() );

		ImageType::Pointer output = km::binaryThresholdImage<ImageType, ImageType>( labelMapToLabelImageFilter->GetOutput(), maxLabel, maxLabel, 1, 0 );
		return output;
	}


	template<class ImageType>
	void
		locateLiverFromProbablityMap(typename ImageType* image, typename ImageType::PointType& centroid)
	{
		double spacingZ = (image->GetSpacing())[2];

		int searchBandHeight = std::min( static_cast<int>(80.0/spacingZ), static_cast<int>((image->GetLargestPossibleRegion().GetSize())[2]) );
		ImageType::IndexType lowerIndex = image->GetLargestPossibleRegion().GetIndex();
		ImageType::SizeType searchSize  = image->GetLargestPossibleRegion().GetSize();
		searchSize[0] = 0.6*searchSize[0];
		searchSize[2] = searchBandHeight;

		ImageType::RegionType searchRegion( lowerIndex, searchSize );

		ImageType::RegionType liverCenterRegion( lowerIndex, searchSize );

		int maxLiverPixelsCount = 0;
		while( image->GetLargestPossibleRegion().IsInside( lowerIndex ) )
		{
			searchRegion.SetIndex( lowerIndex );

			int tmpPixelsCount = km::countPixels<ImageType>( image, searchRegion, 0.5, 1.0 );

			std::cout<<lowerIndex<<", "<<searchSize<< ", " << tmpPixelsCount<< std::endl;

			if (tmpPixelsCount > maxLiverPixelsCount)
			{
				liverCenterRegion.SetIndex( lowerIndex );

				maxLiverPixelsCount = tmpPixelsCount;
			}

			lowerIndex[2] = lowerIndex[2] + searchBandHeight/2;
		}

		liverCenterRegion.Crop( image->GetLargestPossibleRegion() );

		centroid = getCentroid<ImageType>(image, liverCenterRegion, 0.2, 1.0);
	}

	template<class ImageType, class MeshType>
	void
		detectLiverIntensityRangeWithoutTumor(typename ImageType* image, typename MeshType* mesh, std::vector<std::pair<double, double>>& ranges)
	{
		g_liverThresholds.clear();

		double liverThresholdLow, liverThresholdHigh;
		liverThresholdLow = liverThresholdHigh = 0;

		MeshType::Pointer triangleMesh = km::simplexMeshToTriangleMesh<MeshType, MeshType>( mesh );
		//km::writeMesh<MeshType>( "triangleMesh.vtk", triangleMesh );

		typedef itk::Image<unsigned char, ImageType::ImageDimension> MaskImageType;
		MaskImageType::Pointer mask = km::generateBinaryFromMesh<MeshType, MaskImageType, ImageType>( triangleMesh, image );
		//km::writeImage<MaskImageType>( "mask.nii.gz", mask );

		ImageType::Pointer roi = km::maskImage<ImageType, MaskImageType>( image, mask );
		//km::writeImage<ImageType>( "roi.nii.gz", roi );

		const double statisticMin = 15.0;
		const double statisticMax = 300.0;
		const double binWidth = 3.0;

		const double thresholdSpan = 30;

		Histogram histogram(statisticMin, statisticMax, binWidth);
		km::calculateHistogram<ImageType>(roi, &histogram);

		double minLiverTissueRatio = 0.60;

		do 
		{
			double thtmp1, thtmp2;
			histogram.calcThresholdOfRatio(minLiverTissueRatio + 0.05, thtmp1, thtmp2);

			liverThresholdLow = thtmp1;
			liverThresholdHigh = thtmp2;

			std::cout<<"<" <<thtmp1 << ", " <<thtmp2 <<  ">: "<< minLiverTissueRatio + 0.05 <<std::endl;

			double thresholdDiff = thtmp2 - thtmp1;
			if (thresholdDiff<thresholdSpan)
			{
				liverThresholdLow = thtmp1;
				liverThresholdHigh = thtmp2;
				minLiverTissueRatio += 0.05;
			}
			else
			{
				break;
			}
		} while (minLiverTissueRatio<0.95);

		ranges.push_back(std::pair<double, double>( liverThresholdLow, liverThresholdHigh)); 
	}

	template<class ImageType, class MeshType>
	void
		detectLiverIntensityRangeIncludingTumor(typename ImageType* image, typename MeshType* mesh, std::vector<std::pair<double, double>>& ranges)
	{
		g_liverThresholds.clear();

		double liverThresholdLow, liverThresholdHigh;
		liverThresholdLow = liverThresholdHigh = 0;

		MeshType::Pointer triangleMesh = km::simplexMeshToTriangleMesh<MeshType, MeshType>( mesh );
		//km::writeMesh<MeshType>( "triangleMesh.vtk", triangleMesh );

		typedef itk::Image<unsigned char, ImageType::ImageDimension> MaskImageType;
		MaskImageType::Pointer mask = km::generateBinaryFromMesh<MeshType, MaskImageType, ImageType>( triangleMesh, image );
		//km::writeImage<MaskImageType>( "mask.nii.gz", mask );

		ImageType::Pointer roi = km::maskImage<ImageType, MaskImageType>( image, mask );
		//km::writeImage<ImageType>( "roi.nii.gz", roi );

		const double statisticMin = 15.0;
		const double statisticMax = 300.0;
		const double binWidth = 3.0;

		const double thresholdSpan = 40;

		Histogram histogram(statisticMin, statisticMax, binWidth);
		km::calculateHistogram<ImageType>(roi, &histogram);

		double minLiverTissueRatio = 0.60;

		do 
		{
			double thtmp1, thtmp2;
			histogram.calcThresholdOfRatio(minLiverTissueRatio + 0.05, thtmp1, thtmp2);

			liverThresholdLow = thtmp1;
			liverThresholdHigh = thtmp2;

			std::cout<<"<" <<thtmp1 << ", " <<thtmp2 <<  ">: "<< minLiverTissueRatio + 0.05 <<std::endl;

			double thresholdDiff = thtmp2 - thtmp1;
			if (thresholdDiff<thresholdSpan)
			{
				minLiverTissueRatio += 0.05;
			}
			else
			{
				break;
			}
		} while (minLiverTissueRatio<0.95);

		std::cout<<"Minimum liver tissue ratio: "<<minLiverTissueRatio<<std::endl;

		ranges.push_back(std::pair<double, double>( liverThresholdLow, liverThresholdHigh )); 

		//Start to check possible tumors
		if ( minLiverTissueRatio < 0.85 )
		{
			//km::thresholdOutImage<ImageType>( roi, liverThresholdLow, liverThresholdHigh, 0 );
			//km::calculateHistogram<ImageType>(roi, &histogram);

			const double minTumorPercentage = 0.05;
			const double thresholdSpan2 = 30.0;
			double tmpThresholdLower, tmpThresholdUpper;
			
			tmpThresholdLower = histogram.statisticMin;
			tmpThresholdUpper = tmpThresholdLower + thresholdSpan2;
			while(tmpThresholdUpper < liverThresholdLow - 15.0 )
			{
				double tmpRatio = histogram.calcRatioOfThreshold(tmpThresholdLower, tmpThresholdUpper);
				std::cout<<"<" <<tmpThresholdLower << ", " <<tmpThresholdUpper <<  ">: "<< tmpRatio <<std::endl;

				if (tmpRatio>=minTumorPercentage)
				{
					ranges.push_back(std::pair<double,double>(tmpThresholdLower, tmpThresholdUpper));
					tmpThresholdLower += binWidth;
					tmpThresholdUpper += binWidth;
				}
				else
				{
					tmpThresholdLower += binWidth;
					tmpThresholdUpper += binWidth;
				}
			}

			tmpThresholdLower = liverThresholdHigh + 30.0;
			tmpThresholdUpper = tmpThresholdLower + thresholdSpan2;
			while(tmpThresholdUpper<statisticMax)
			{
				double tmpRatio = histogram.calcRatioOfThreshold(tmpThresholdLower, tmpThresholdUpper);
				std::cout<<"<" <<tmpThresholdLower << ", " <<tmpThresholdUpper <<  ">: "<< tmpRatio <<std::endl;

				if (tmpRatio>=minTumorPercentage)
				{
					ranges.push_back(std::pair<double,double>(tmpThresholdLower, tmpThresholdUpper));
					tmpThresholdLower += 0.33*thresholdSpan2;
					tmpThresholdUpper += 0.33*thresholdSpan2;
				}
				else
				{
					tmpThresholdLower += binWidth;
					tmpThresholdUpper += binWidth;
				}

				tmpThresholdLower += binWidth;
				tmpThresholdUpper += binWidth;
			}

			if (g_liverThresholds.size() == 1 && (g_liverThresholds[0].second - g_liverThresholds[0].first)<thresholdSpan )
			{
				//std::cout<<"***"<<g_liverThresholds[0].first<<","<<g_liverThresholds[0].second<<"***, "<<g_liverThresholds[0].second - g_liverThresholds[0].first<<"***"<<thresholdSpan<<std::endl;
				g_liverThresholds[0].first -= 10;
				g_liverThresholds[0].second += 10;
			}
		}
	}


	template<class ImageType>
	typename ImageType::PointType
		getBinaryTopestPoint( const ImageType* image, double backgroundValue )
	{
		typedef itk::ImageRegionConstIterator<ImageType> IteratorType;
		IteratorType it( image, image->GetLargestPossibleRegion() );
		it.GoToBegin();
		ImageType::IndexType topestIndex;
		topestIndex.Fill( -9999 );
		while( !it.IsAtEnd() )
		{
			if( it.Get() != backgroundValue )
			{
				if( (it.GetIndex())[2] > topestIndex[2] )
				{
					topestIndex = it.GetIndex();
				}
			}

			++it;
		}

		ImageType::PointType pt;
		image->TransformIndexToPhysicalPoint( topestIndex, pt );

		return pt;
	}

	template<typename TImageType>
	typename TImageType::Pointer
		changeImage(const typename TImageType* inputImage, 
								const typename TImageType::SpacingType newSpacing,
								const typename TImageType::PointType   newOrigin)
	{
		typedef itk::ChangeInformationImageFilter<TImageType> ChangeInformationImageFilterType;
		ChangeInformationImageFilterType::Pointer changer = ChangeInformationImageFilterType::New();
		changer->SetInput( inputImage );
		changer->SetOutputSpacing( newSpacing );
		changer->ChangeSpacingOn();
		changer->SetOutputOrigin( newOrigin );
		changer->ChangeOriginOn();
		changer->Update();

		return changer->GetOutput();
	}


	////提取人体部分API
	////其中有对输入图像进行重采样，加快速度
	//template<class TImageType, class MaskImageType>
	//void
	//	extractBody( typename TImageType* inputImage,
	//							 const double standardMinmum,
	//							 typename TImageType::Pointer &bodyImage,
	//							 typename MaskImageType::Pointer &bodyMask )
	//{
	//	//重采样，加快速度
	//	typedef TImageType::SpacingType SpacingType;

	//	SpacingType inputspacing = inputImage->GetSpacing();

	//	TImageType::Pointer downsampledImage = TImageType::New();

	//	if (inputspacing[0]<3.0)
	//	{
	//		SpacingType downsampleSpacing;
	//		downsampleSpacing.Fill(4.0);
	//		TImageType::Pointer downsampledImage = km::resampleImage<TImageType>( inputImage,  downsampleSpacing, standardMinmum);
	//	}
	//	else
	//	{
	//		downsampledImage = inputImage;
	//	}

	//	//除去空气部分
	//	downsampledImage = km::removeAir<TImageType, MaskImageType>( downsampledImage, bodyMask, standardMinmum );


	//	//将bodyMask重采样回测试图像的尺寸
	//	bodyMask = km::resampleImageByReference<MaskImageType, TImageType>( bodyMask, inputImage, 0 );
	//	bodyImage = km::maskImage<TImageType, MaskImageType>( inputImage, bodyMask );
	//}

	//移除空气的API
	//参数：
	//	     airValue 输入图像空气灰度值
	template<class TImageType, class MaskImageType> 
	typename TImageType::Pointer 
		removeAir(typename TImageType* inputImage, typename MaskImageType::Pointer &maskImage, double airValue, bool reversed = true)
	{
		//通过区域增长，提取人体部分
		kmStaticImageMacro(TImageType);
		const unsigned int Dimension = TImageType::ImageDimension;

		//为了保证所有空气都被连通，首先将图像向外pad 5个像素
		TImageType::SizeType padLowerBound, padUpperBound;
		padLowerBound.Fill(10);
		padUpperBound.Fill(10);

		if ( Dimension == 3 )
		{
			padUpperBound[2] = 0;
		}

		TImageType::Pointer paddedImage = km::padImage<TImageType>( inputImage, padLowerBound, padUpperBound, airValue );

		typedef itk::ConnectedThresholdImageFilter< TImageType, MaskImageType > ConnectedFilterType;

		//首先从空气部分开始增长，将空气与人体隔开
		RegionType region = paddedImage->GetRequestedRegion();
		SizeType   size   = region.GetSize();

		//选取横断面左上角、冠状面最上面的点作为种子点，因为这个点肯定是空气，灰度值为最低
		IndexType seed = region.GetIndex();
		seed[0] += 5;
		seed[1] += 5;

		if ( Dimension == 3 )
		{
			seed[2] += size[2]-5;
		}

		std::cout<<seed<<std::endl;

		ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
		connectedThreshold->SetInput( paddedImage );
		connectedThreshold->SetLower(  static_cast<float>(airValue) );
		connectedThreshold->SetUpper(  static_cast<float>(airValue+500) );
		connectedThreshold->SetReplaceValue( 1 );
		//connectedThreshold->SetRadius( radius );
		connectedThreshold->SetSeed( seed );
		try
		{
			connectedThreshold->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception caught !" << std::endl;
			std::cerr << excep << std::endl;
		}
		maskImage = connectedThreshold->GetOutput();
		MaskImageType::SizeType radius;
		radius[0] = 3;
		radius[1] = 3;
		if ( Dimension == 3 )
		{
			radius[2] = 1;
		}
		
		maskImage = km::binaryMedianSmooth<MaskImageType>( maskImage, radius );
		
		if (reversed)
		{
			//灰度值交换，交换后，人体部分灰度值为1，空气部分灰度值为0
			maskImage = km::shiftScale<MaskImageType>( maskImage, -1, -1 );
		}

		//由于先前有pad，所以需要重新采样，将pad区域去掉
		maskImage = resampleByCopying<MaskImageType,TImageType>( maskImage, inputImage );

		//将输入图像与人体mask图像相乘，得到去掉空气后的人体灰度图像
		TImageType::Pointer bodyImage = km::maskImage<TImageType, MaskImageType>(inputImage, maskImage);

		return bodyImage;
	}


	//旋转纠正API

	//旋转纠正的原理：
	//	人体横断面，如果是正常角度，在X方向应该是基本对称的，尤其是肺部区域。假设有旋转角度alpha,那么将图像按X轴方向翻转，
	//	然后将原图像与翻转后的图像做旋转配准，得到的旋转角度除以2，则可以将图像旋转成在X方向基本对称的图像。

	//旋转纠正的过程：

	//	1. 从输入图像中，搜索平均灰度值最低的2维切片作为候选切片lung slice
	//	2. 将lung slice沿着x轴方向翻转，得到翻转后切片fliped lung slice
	//	3. 将两张lung slice进行旋转配准，得到2维刚体变换
	//	4. 将2维刚体变换，转换为3维刚体变换
	//	5. 图像变换，返回

	template<class TImageType, class Rigid3DTransformType>
	void
		rotationCorrection(const typename TImageType* inputImage, 
											 typename Rigid3DTransformType::Pointer &rigid3DTransform,
											 typename TImageType::PointType& lungCentroidPoint)
	{

		typedef itk::Image<TImageType::PixelType, 2> SliceType;
		typedef itk::Image<unsigned char, 2> MaskSliceType;

		//搜索肺部切片
		int indexOfCoreLungSlice = km::getMaximumSliceIndex<TImageType>( inputImage, -1024, -500 );

		std::cout<<"index of core lung slice: "<<indexOfCoreLungSlice<<std::endl;

		SliceType::Pointer lungCoreSlice = extractSlice<TImageType, SliceType>( inputImage, indexOfCoreLungSlice );

		//km::writeImage<SliceType>( "lungCoreSlice.nii.gz", lungCoreSlice );

		MaskSliceType::Pointer lungSliceBinary = km::binaryThresholdImage<SliceType, MaskSliceType>( lungCoreSlice, -1024, -500, 1, 0 );
		MaskSliceType::SizeType radius;
		radius[0] = 1;
		radius[1] = 3;
		lungSliceBinary = km::binaryMedianSmooth<MaskSliceType>( lungSliceBinary, radius );

		//标记每个连通区域
		unsigned int numberOfLabels = 0;
		MaskSliceType::Pointer labels = km::generateLabelImage<MaskSliceType, MaskSliceType>( lungSliceBinary,  numberOfLabels);
		//km::writeImage<MaskSliceType>("lungLabel.nii.gz", labels);

		int* sizeArray = new int[numberOfLabels];
		int* labelArray = new int[numberOfLabels];
		
		//统计每个label像素数目，并排序。排序后的value数组和size数组（记录label像素数目）是一一对应的
		km::calculateAndSortLabelComponent<MaskSliceType>( labels, labelArray, sizeArray, numberOfLabels );

		MaskSliceType::PointType centroid_pt1 = getCentroid<MaskSliceType>( labels, labelArray[0], labelArray[0] );
		MaskSliceType::PointType centroid_pt2 = getCentroid<MaskSliceType>( labels, labelArray[1], labelArray[1] );

		MaskSliceType::PointType leftLungCentroid, rightLungCetroid;
		double leftLungLabelValue, rightLungLabelValue;
		if(centroid_pt1[0] < centroid_pt2[0] )
		{
			leftLungLabelValue = labelArray[0];
			rightLungLabelValue = labelArray[1];
		}
		else
		{
			leftLungLabelValue = labelArray[1];
			rightLungLabelValue = labelArray[0];
			std::swap( centroid_pt1[0], centroid_pt2[0] );
			std::swap( centroid_pt1[1], centroid_pt2[1] );
		}

		//std::cout<<"left lung slice centroid: "<<centroid_pt1<<std::endl;
		//std::cout<<"leftLungMaskValue: "<< leftLungLabelValue<<"  rightLungMaskValue:"<< rightLungLabelValue<<std::endl;
		//std::cout<<"lung centroid point: "<<lungCentroidPoint<<std::endl;

		typedef itk::Vector<double, 2> VectorType;
		VectorType vector1;
		vector1[0] = centroid_pt2[0] - centroid_pt1[0];
		vector1[1] = centroid_pt2[1] - centroid_pt1[1];

		VectorType vector2;
		vector2[0] = 1.0;
		vector2[1] = 0.0;

		double dotp = vector1*vector2;
		double cosn = dotp / (vector1.GetNorm() * vector2.GetNorm() );
		double angle = std::acos( cosn );

		if( vector2[2]<0 )
		{
			angle *= -1;
		}

		angle += 0.1;//左右肺本来就不对称，顺时针旋转0.1(15度)，消除一点偏差

		TImageType::PointType bodyCentroid = getCenter<TImageType>( inputImage );
		Rigid3DTransformType::ParametersType fixedParameters = rigid3DTransform->GetFixedParameters();
		Rigid3DTransformType::ParametersType rigidParameters = rigid3DTransform->GetParameters();
		for (int i=0;i<3;i++)
		{
			fixedParameters[i] = bodyCentroid[i];
		}

		//角度比较小的情况下，不纠正
		if(std::abs(angle) < 0.35)
		{
			rigidParameters[2] = 0;
		}
		else if(angle > 0.35)
		{
			rigidParameters[2] = angle;
		}

		std::cout<<"angle: "<<angle<<std::endl;

		rigid3DTransform->SetFixedParameters( fixedParameters );
		rigid3DTransform->SetParameters( rigidParameters );


		//get lung centroid

		MaskSliceType::IndexType leftLungSliceCentroidIndex;
		labels->TransformPhysicalPointToIndex( centroid_pt1, leftLungSliceCentroidIndex );

		TImageType::IndexType lungCentroidIndex;
		lungCentroidIndex[0] = leftLungSliceCentroidIndex[0];
		lungCentroidIndex[1] = leftLungSliceCentroidIndex[1];
		lungCentroidIndex[2] = indexOfCoreLungSlice;

		inputImage->TransformIndexToPhysicalPoint( lungCentroidIndex, lungCentroidPoint );
		lungCentroidPoint = rigid3DTransform->TransformPoint( lungCentroidPoint );

		rigid3DTransform->GetInverse(rigid3DTransform);

		//TImageType::Pointer rotatedImage = km::transformImage<TImageType, Rigid3DTransformType>( inputImage, inputImage, rigid3DTransform );

		//return rotatedImage;
	}

	template<class PointType>
	double
		calcuateAngle(typename PointType pt1, typename PointType pt2)
	{
		typedef itk::Vector<double, PointType::Dimension> VectorType;
		VectorType vector1;
		for (int i=0;i<PointType::Dimension;i++)
		{
			vector1[i] = pt2[i]-pt1[i];
		}

		VectorType vector2;
		vector2.Fill( 1 );
		vector2[ PointType::Dimension-1] = 0;

		double dotp = vector1*vector2;
		double cosn = dotp / (vector1.GetNorm() * vector2.GetNorm() );
		double angle = std::acos( cosn );

		if (vector1[PointType::Dimension-1]<0)
		{
			angle *= -1.0;
		}

		return angle;
	}

	template<typename TMeshType>
	typename TMeshType::PointType
		findPointOfLiverLeftLopTip(const TMeshType * inputMesh)
	{
		TMeshType::PointType tipPt;
		tipPt.Fill( 0 );
		bool foundFirst = false;

		MeshType::PointsContainerConstIterator pointIt = inputMesh->GetPoints()->Begin();
		MeshType::PointsContainerConstIterator pointItEnd = inputMesh->GetPoints()->End();
		while(pointIt!=pointItEnd)
		{
			TMeshType::PointType pt = pointIt->Value();

			if (!foundFirst)
			{
				tipPt = pt;
			}
			else
			{
				foundFirst = true;
				if (pt[0] > tipPt[0])
				{
					tipPt = pt;
				}
			}

			++pointIt;
		}

		if (!foundFirst)
		{
			std::cout<<"Warning: no tip of liver left lop found!"<<std::endl;
		}

		return tipPt;
	}

	//寻找肝左叶最靠左的体素点坐标
	template<typename TImageType>
	typename TImageType::IndexType
		findIndexOfLiverLeftLopTip(const typename TImageType* input3DImage, double lowerthreshold, double upperthreshold)
	{
		TImageType::IndexType tipIdx;
		tipIdx.Fill( 0 );
		bool foundFirst = false;

		itk::ImageRegionConstIteratorWithIndex< TImageType > it(input3DImage, input3DImage->GetLargestPossibleRegion());
		it.GoToBegin();
		while(!it.IsAtEnd())
		{
			TImageType::PixelType val = it.Get();
			TImageType::IndexType idx = it.GetIndex();

			if (val>=lowerthreshold && val<=upperthreshold)
			{
				if (!foundFirst)
				{
					tipIdx = idx;
				}
				else
				{
					foundFirst = true;
					if ( idx[0] > tipIdx[0] )
					{
						tipIdx = idx;
					}
				}
			}

			++it;
		}

		if (!foundFirst)
		{
			std::cout<<"Warning: no tip of liver left lop found!"<<std::endl;
		}

		return tipIdx;
	}

	//寻找肺部切片API
	template<typename TImageType>
	int
		getMaximumSliceIndex(const typename TImageType* input3DImage, double lowerthreshold, double upperthreshold)
	{
		unsigned int maximumPixelCounts = 0;
		int searchSliceIndex = 0;

		typedef typename TImageType::IndexType IndexType;
		typedef typename TImageType::PixelType PixelType;

		typedef itk::ImageSliceConstIteratorWithIndex< TImageType > Input3DConstIteratorType;
		Input3DConstIteratorType it( input3DImage, input3DImage->GetLargestPossibleRegion() );
		it.GoToBegin();
		it.SetFirstDirection( 0 );  // 0=x, 1=y, 2=z
		it.SetSecondDirection( 1 ); // 0=x, 1=y, 2=z

		int sliceIndex = 0;
		while( !it.IsAtEnd() )
		{
			unsigned long pixelCounted = 0;
			while( !it.IsAtEndOfSlice() )
			{
				while( !it.IsAtEndOfLine() )
				{
					IndexType ind = it.GetIndex();
					PixelType val = it.Get();

					//std::cout<<val<<",";

					if( val>=lowerthreshold && val<=upperthreshold)
					{
						pixelCounted++;
					}
					++it;
				}
				it.NextLine();
			}

			if(pixelCounted > maximumPixelCounts)
			{
				searchSliceIndex = sliceIndex;
				maximumPixelCounts = pixelCounted;
			}

			++sliceIndex;
			it.NextSlice();

		} //end of while 

		//std::cout<<"lung slice: "<<sliceIndex<<", mean slice value:"<<minMeanValue<<std::endl;
		return searchSliceIndex;

	}//end of method

	template<class TImageType>
	typename TImageType::Pointer
		alignedCentroid(typename TImageType* fixedImage, typename TImageType* movingImage, double defaultValue)
	{
		TImageType::PointType fixedCentroid = km::getCentroid<TImageType>(fixedImage);
		TImageType::PointType movingCentroid = km::getCentroid<TImageType>(movingImage);
		typedef itk::TranslationTransform<double, TImageType::ImageDimension> TransformType;
		TransformType::Pointer transform = TransformType::New();
		itk::Vector<double, TImageType::ImageDimension> offset;
		for(int i=0;i<TImageType::ImageDimension;i++)
		{
			offset[i] = movingCentroid[i] - fixedCentroid[i];
		}
		transform->SetOffset( offset );

		TImageType::Pointer alignedMoving = transformImage<TImageType, TransformType>( movingImage, fixedImage, transform, defaultValue );

		return alignedMoving;
	}

	template<class TImageType>
	void
		findTopOfLiver(const typename TImageType* image, 
															 typename TImageType::PointType lungCentroidPoint,
															 typename TImageType::PointType& liverTopPoint)
	{
		//std::cout<<liverTopPoint<<std::endl;
		int lowestNonLungSliceIndex = 9999;
		TImageType::IndexType tempIndex;
		image->TransformPhysicalPointToIndex( lungCentroidPoint, tempIndex );

		TImageType::SpacingType spacing = image->GetSpacing();
		TImageType::SizeType searchsize;
		searchsize.Fill(1);
		searchsize[2] = 150/spacing[2];
		TImageType::IndexType searchindex = tempIndex;
		searchindex[2] -= 100/spacing[2];
		TImageType::RegionType searchRegion( searchindex, searchsize );
		std::cout<<searchindex<<std::endl;
		std::cout<<searchsize<<std::endl;
		searchRegion.Crop( image->GetLargestPossibleRegion() );

		typedef itk::ImageRegionConstIterator<TImageType> IteratorType;
		IteratorType it( image, searchRegion );
		it.GoToBegin();
		while(!it.IsAtEnd())
		{
			double val = static_cast<double>( it.Get() );
			if (val>-1024 && val<-500)
			{
				if( (it.GetIndex())[2] < lowestNonLungSliceIndex )
				{
					lowestNonLungSliceIndex = (it.GetIndex())[2];
				}
			}
			//std::cout<<it.GetIndex()<<" "<<val<<std::endl;
			++it;
		}

		if(lowestNonLungSliceIndex != 9999)
		{
			tempIndex[2] = lowestNonLungSliceIndex;
		}
		
		image->TransformIndexToPhysicalPoint( tempIndex, liverTopPoint );
	}

	//填充图像在所有横断面切片上的空洞
	template<typename ImageType>
	void
		fillSliceHole( typename ImageType::Pointer & inputImage )
	{
		ImageType::Pointer inputImageBackup = ImageType::New();
		inputImageBackup->CopyInformation( inputImage );

		//ImageType::SpacingType spacing;
		//spacing.Fill(2);
		//ImageType::Pointer tmpImage = km::resampleImage<ImageType>(inputImage, spacing, 0, 1);
		ImageType::Pointer tmpImage = km::castImage<ImageType, ImageType>(inputImage);

		ImageType::RegionType roiRegion;
		km::getBoundRegion<ImageType>(tmpImage, roiRegion, 0, 5);

		typedef BinaryThresholdImageFunction< ImageType, double > FunctionType;
		FunctionType::Pointer function = FunctionType::New();
		function->SetInputImage (tmpImage);
		function->ThresholdBetween (0, 0);

		ImageType::RegionType sliceRegion = roiRegion;
		ImageType::SizeType sliceSize = sliceRegion.GetSize();
		sliceSize[2] = 1;
		ImageType::IndexType sliceIndex = sliceRegion.GetIndex();

		sliceRegion.SetSize(sliceSize);
		sliceRegion.SetIndex(sliceIndex);

		ImageType::Pointer sliceImage = ImageType::New();
		sliceImage->SetRegions(sliceRegion);

		for (unsigned i=0;i<roiRegion.GetSize()[2];i++)
		{
			sliceRegion.SetIndex(sliceIndex);
			sliceImage->SetBufferedRegion(sliceRegion);

			ImageType::IndexType seed = sliceIndex;
			seed[0]+=1;
			seed[1]+=1;

			typedef ShapedFloodFilledImageFunctionConditionalIterator< ImageType, FunctionType > IteratorType;
			IteratorType it (sliceImage, function, seed);
			it.FullyConnectedOn();
			it.GoToBegin();
			while ( !it.IsAtEnd() )
			{
				tmpImage->SetPixel(it.GetIndex(), 2);
				++it;
			}

			sliceIndex[2] += 1;
		}
		tmpImage = km::extractByBoundRegion<ImageType>(tmpImage, roiRegion);
		tmpImage = km::binaryThresholdImage<ImageType, ImageType>(tmpImage, 0, 1, 1, 0);
		
		inputImage = km::resampleImageByReference<ImageType, ImageType>(tmpImage, inputImageBackup, 0, 1);
	}

	template<class ImageType, class MeshType>
	void
		evaluateIntensity(
												const typename ImageType* inputImage,
												const typename MeshType* mesh,
												double staMin,
												double staMax,
												double & thresholdLower,
												double & thresholdUpper,
												double minRatio = 0.6)
	{
		MeshType::Pointer triangleMesh = km::simplexMeshToTriangleMesh<MeshType, MeshType>( mesh );

		typedef itk::Image<unsigned char, ImageType::ImageDimension> MaskImageType;
		MaskImageType::Pointer mask = km::generateBinaryFromMesh<MeshType, MaskImageType, ImageType>( triangleMesh, inputImage );
		ImageType::Pointer roi = km::maskImage<ImageType, MaskImageType>( inputImage, mask );

		evaluateIntensity<ImageType>( roi, staMin, staMax, thresholdLower, thresholdUpper, minRatio );
	}

	template<class ImageType>
	void
		evaluateIntensity(
											const typename ImageType* inputImage, 
											double staMin, 
											double staMax, 
											double & thresholdLower, 
											double & thresholdUpper, 
											double minRatio = 0.6)
	{
		Histogram histogram(staMin, staMax, 1.0);
		calculateHistogram<ImageType>( inputImage, &histogram );

		histogram.calcThresholdOfRatio( minRatio, thresholdLower, thresholdUpper );
	}

	template<typename LabelImageType>
	void
		calculateAndSortLabelComponent( const typename LabelImageType* labels, int* labelArray, int* sizeArray, int numberOfLabels )
	{
		//记录下每个连同区域的label值和像素点
		for(int i = 0; i < numberOfLabels; i++)
		{
			labelArray[i] = i+1;
			sizeArray[i] = 0;
		}

		typedef itk::ImageRegionConstIterator<LabelImageType> IteratorType;
		IteratorType it(labels, labels->GetLargestPossibleRegion());
		it.GoToBegin();
		while (!it.IsAtEnd())
		{
			int labelValue = static_cast<int>( it.Get() );
			if(labelValue!=0)
			{
				sizeArray[labelValue-1]++;
			} 
			it++;
		}

		km::sort_insert<int, int>( sizeArray, labelArray, numberOfLabels, false );
	}

	template<class InputImageType, class OutputImageType>
	typename OutputImageType::Pointer
		regionGrow(typename InputImageType* inputImage, typename InputImageType::IndexType seed, double thLower, double thUpper, double replaceVal )
	{
		typedef itk::ConnectedThresholdImageFilter< InputImageType, OutputImageType > ConnectedFilterType;
		ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
		connectedThreshold->SetInput( inputImage );
		connectedThreshold->SetLower(  thLower );
		connectedThreshold->SetUpper(  thUpper );
		connectedThreshold->SetReplaceValue( replaceVal );
		//connectedThreshold->SetRadius( radius );
		connectedThreshold->SetSeed( seed );

		try
		{
			connectedThreshold->Update();
		}
		catch (itk::ExceptionObject* e)
		{
			std::cout<<e<<std::endl;
		}

		return connectedThreshold->GetOutput();
	}

	template<class InputImageType, class OutputImageType>
	typename OutputImageType::Pointer
		binaryThresholdImage(const typename InputImageType* inputImage, const std::vector<std::pair<double, double>> & thresholds)
	{
		OutputImageType::Pointer mask = OutputImageType::New();
		mask->SetRegions( inputImage->GetLargestPossibleRegion() );
		mask->SetSpacing( inputImage->GetSpacing() );
		mask->SetOrigin( inputImage->GetOrigin() );
		mask->SetDirection( inputImage->GetDirection() );
		mask->Allocate();
		mask->FillBuffer(0);

		typedef itk::OrImageFilter <OutputImageType> OrImageFilterType;

		for(int i=0;i<thresholds.size();i++)
		{
			OutputImageType::Pointer masktmp = km::binaryThresholdImage<InputImageType, OutputImageType>(inputImage, thresholds[i].first, thresholds[i].second, 1, 0);
			
			OrImageFilterType::Pointer orFilter = OrImageFilterType::New();
			orFilter->SetInput(0, mask);
			orFilter->SetInput(1, masktmp);
			orFilter->Update();

			mask = orFilter->GetOutput();
		}

		return mask;
	}
}

#endif

