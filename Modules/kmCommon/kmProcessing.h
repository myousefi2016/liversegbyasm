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

#include "kmCommon.h"

namespace km
{
	/************************************************************************/
	//(新)预处理API                                                            
	/************************************************************************/

	/*
	inline char* dir(char* dir)
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
	*/

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
		if (inputUpperPt[2] < paddedAtlasUpperZPos || inputStartPt[2] > paddedAtlasStartZPos){
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
		for(unsigned int i = 0; i < binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++){
			// Get the ith region
			BinaryImageToLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToLabelMapFilter->GetOutput()->GetNthLabelObject(i);
			if (labelObject->Size() > maxSize){
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
		typedef itk::Image<unsigned char, ImageType::ImageDimension> MaskImageType;
		MaskImageType::Pointer mask = km::generateBinaryFromMesh<MeshType, MaskImageType, ImageType>( triangleMesh, image );
		ImageType::Pointer roi = km::maskImage<ImageType, MaskImageType>( image, mask );
		const double statisticMin = 15.0;
		const double statisticMax = 300.0;
		const double binWidth = 3.0;
		const double thresholdSpan = 30;
		Histogram histogram(statisticMin, statisticMax, binWidth);
		km::calculateHistogram<ImageType>(roi, &histogram);
		double minLiverTissueRatio = 0.60;
		do {
			double thtmp1, thtmp2;
			histogram.calcThresholdOfRatio(minLiverTissueRatio + 0.05, thtmp1, thtmp2);
			liverThresholdLow = thtmp1;
			liverThresholdHigh = thtmp2;
			std::cout<<"<" <<thtmp1 << ", " <<thtmp2 <<  ">: "<< minLiverTissueRatio + 0.05 <<std::endl;
			double thresholdDiff = thtmp2 - thtmp1;
			if (thresholdDiff<thresholdSpan){
				liverThresholdLow = thtmp1;
				liverThresholdHigh = thtmp2;
				minLiverTissueRatio += 0.05;
			}else{
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
		typedef itk::Image<unsigned char, ImageType::ImageDimension> MaskImageType;
		MaskImageType::Pointer mask = km::generateBinaryFromMesh<MeshType, MaskImageType, ImageType>( triangleMesh, image );
		ImageType::Pointer roi = km::maskImage<ImageType, MaskImageType>( image, mask );
		const double statisticMin = 15.0;
		const double statisticMax = 300.0;
		const double binWidth = 3.0;
		const double thresholdSpan = 40;
		Histogram histogram(statisticMin, statisticMax, binWidth);
		km::calculateHistogram<ImageType>(roi, &histogram);
		double minLiverTissueRatio = 0.60;
		do {
			double thtmp1, thtmp2;
			histogram.calcThresholdOfRatio(minLiverTissueRatio + 0.05, thtmp1, thtmp2);
			liverThresholdLow = thtmp1;
			liverThresholdHigh = thtmp2;
			std::cout<<"<" <<thtmp1 << ", " <<thtmp2 <<  ">: "<< minLiverTissueRatio + 0.05 <<std::endl;
			double thresholdDiff = thtmp2 - thtmp1;
			if (thresholdDiff<thresholdSpan){
				minLiverTissueRatio += 0.05;
			}else{
				break;
			}
		} while (minLiverTissueRatio<0.95);

		std::cout<<"Minimum liver tissue ratio: "<<minLiverTissueRatio<<std::endl;
		ranges.push_back(std::pair<double, double>( liverThresholdLow, liverThresholdHigh )); 
		//Start to check possible tumors
		if ( minLiverTissueRatio < 0.85 ){
			const double minTumorPercentage = 0.05;
			const double thresholdSpan2 = 30.0;
			double tmpThresholdLower, tmpThresholdUpper;
			tmpThresholdLower = histogram.statisticMin;
			tmpThresholdUpper = tmpThresholdLower + thresholdSpan2;
			while(tmpThresholdUpper < liverThresholdLow - 15.0 ){
				double tmpRatio = histogram.calcRatioOfThreshold(tmpThresholdLower, tmpThresholdUpper);
				std::cout<<"<" <<tmpThresholdLower << ", " <<tmpThresholdUpper <<  ">: "<< tmpRatio <<std::endl;
				if (tmpRatio>=minTumorPercentage){
					ranges.push_back(std::pair<double,double>(tmpThresholdLower, tmpThresholdUpper));
					tmpThresholdLower += binWidth;
					tmpThresholdUpper += binWidth;
				}else{
					tmpThresholdLower += binWidth;
					tmpThresholdUpper += binWidth;
				}
			}
			tmpThresholdLower = liverThresholdHigh + 30.0;
			tmpThresholdUpper = tmpThresholdLower + thresholdSpan2;
			while(tmpThresholdUpper<statisticMax){
				double tmpRatio = histogram.calcRatioOfThreshold(tmpThresholdLower, tmpThresholdUpper);
				std::cout<<"<" <<tmpThresholdLower << ", " <<tmpThresholdUpper <<  ">: "<< tmpRatio <<std::endl;
				if (tmpRatio>=minTumorPercentage){
					ranges.push_back(std::pair<double,double>(tmpThresholdLower, tmpThresholdUpper));
					tmpThresholdLower += 0.33*thresholdSpan2;
					tmpThresholdUpper += 0.33*thresholdSpan2;
				}else{
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

		if ( Dimension == 3 ){
			seed[2] += size[2]-5;
		}

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

	//填充图像在所有横断面切片上的空洞
	template<typename ImageType>
	void
		fillSliceHole( typename ImageType::Pointer & inputImage )
	{
		ImageType::Pointer inputImageBackup = ImageType::New();
		inputImageBackup->CopyInformation( inputImage );

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
}

#endif

