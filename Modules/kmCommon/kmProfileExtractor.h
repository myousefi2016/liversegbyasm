#ifndef __kmProfileExtractor_h
#define __kmProfileExtractor_h

#include <vector>
#include <map>
#include <memory>
#include "itkMacro.h"
#include "itkObject.h"
#include "itkSimplexMeshGeometry.h"
#include "itkImage.h"
#include "itkCovariantVector.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMeanImageFunction.h"
#include "itkMedianImageFunction.h"
#include "itkVarianceImageFunction.h"
#include "itkGaussianDerivativeImageFunction.h"
#include "itkCentralDifferenceImageFunction.h"
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkGradientRecursiveGaussianImageFilter.h>

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

#include "kmGlobal.h"

namespace km
{
	template<typename TImage>
	class ProfileExtractor
	{
	public:
		typedef typename TImage                                                             ImageType;
		typedef typename ImageType::PixelType                                               PixelType;
		itkStaticConstMacro(ImageDimension, unsigned int, ImageType::ImageDimension);
		typedef typename ImageType::IndexType                                               IndexType;
		typedef typename ImageType::PointType                                               PointType;
		typedef typename PointType::VectorType                                              VectorType;
		typedef itk::Image<float, ImageDimension>                                           GradientImageType;
		typedef typename itk::LinearInterpolateImageFunction<ImageType>                     LinearInterpolateImageFunctionType;
		typedef typename itk::MedianImageFunction<ImageType>                                MedianImageFunctionType;
		typedef typename itk::MeanImageFunction<ImageType>                                  MeanImageFunctionType;
		typedef typename itk::LinearInterpolateImageFunction<GradientImageType>             GradientImageFunctionType;
		typedef typename itk::VarianceImageFunction<ImageType, double>                      VarianceImageFunctionType;
		typedef std::vector<double>                                                         FeatureSetType;

		typedef typename ImageType::OffsetValueType                                         OffsetValueType;
		typedef std::map<OffsetValueType, FeatureSetType*>                                  FeatureSetMap;
		typedef std::map<PROFILE_CATEGORY, FeatureSetMap*>                                  CategorizedFeatureSetMap;
		
		ProfileExtractor()
		{
			this->inputImage = NULL;
			normal.Fill(0);

			linerInterpolatorFunc = LinearInterpolateImageFunctionType::New();
			medianFunc = MedianImageFunctionType::New();
			meanFunc = MeanImageFunctionType::New();
			gradientFunc = GradientImageFunctionType::New();
			varianceFunc = VarianceImageFunctionType::New();

			enableCacheFlag = true;
		}

		~ProfileExtractor()
		{
		}

		void enableCache(bool flag)
		{
			this->enableCacheFlag = flag;
		}
		
		void clear()
		{
			pointSet.clear();
		}

		bool isInsideBuffer(typename PointType & pt)
		{
			return this->linerInterpolatorFunc->IsInsideBuffer(pt);
		}

		void setImage(const ImageType* img)
		{
			this->inputImage = img;

			linerInterpolatorFunc->SetInputImage(inputImage);

			medianFunc->SetInputImage(inputImage);
			medianFunc->SetNeighborhoodRadius(2);

			meanFunc->SetInputImage(inputImage);
			meanFunc->SetNeighborhoodRadius(2);

			varianceFunc->SetInputImage(inputImage);
			varianceFunc->SetNeighborhoodRadius(2);

			typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType, GradientImageType> GradientMagnitudeRecursiveGaussianImageFilterType;
			GradientMagnitudeRecursiveGaussianImageFilterType::Pointer gradientMagFilter = GradientMagnitudeRecursiveGaussianImageFilterType::New();
			gradientMagFilter->SetInput( inputImage );
			gradientMagFilter->SetSigma( SIGMA );
			gradientMagFilter->Update();

			gradientFunc->SetInputImage(gradientMagFilter->GetOutput());

			//std::cout<<g_liverThresholds[0].first<<", "<<g_liverThresholds[0].second<<std::endl;
		}

		void extractFeatureSet( FeatureSetType & featureSet, PROFILE_CATEGORY profileCategory, const itk::SimplexMeshGeometry * geoData, const PointType & cur_pos )
		{
			featureSet.clear();
			pointSet.clear();

			this->currentPoint = cur_pos;
			this->offsetVector = cur_pos - geoData->pos;
			this->geoData = const_cast<itk::SimplexMeshGeometry*>(geoData);
			this->normal.Set_vnl_vector(geoData->normal.Get_vnl_vector());

			OffsetValueType offset = itk::NumericTraits<OffsetValueType>::max();
			IndexType cur_idx;
			bool isInside = this->inputImage->TransformPhysicalPointToIndex(cur_pos, cur_idx);
			if (isInside)
			{
				offset = this->inputImage->ComputeOffset(cur_idx);
			}

			if (profileCategory==PLAIN){
				this->addPoints(profileCategory, PROFILE_DIM, PROFILE_SPACING);
				this->extractIntensitySet(featureSet);
			}else if (profileCategory==LIVER)
			{
				if (enableCacheFlag){
					getCachedFeatures(featureSet, profileCategory, offset);
				}

				if (featureSet.size() == 0){
					this->addPoints(profileCategory, PROFILE_DIM, PROFILE_SPACING);
					this->extractIntensityStatistics(featureSet);
					this->extractGradientSet(featureSet);
					//this->extractVarianceSet(featureSet);
					if (enableCacheFlag){
						setCachedFeatures(featureSet, profileCategory, offset);
					}
				}
			}
			else if (profileCategory==BOUNDARY)
			{
				if (enableCacheFlag){
					getCachedFeatures(featureSet, profileCategory, offset);
				}

				if (featureSet.size() == 0){
					this->addPoints(profileCategory, PROFILE_DIM, PROFILE_SPACING);
					this->extractIntensityStatistics(featureSet);
					this->extractTangentFeatureSet(featureSet);
					if (enableCacheFlag){
						setCachedFeatures(featureSet, profileCategory, offset);
					}
				}
			}
			else if(profileCategory==COORDINATE){
				this->extractCoordinate(featureSet);
			}else{
				std::cout<<"Unknown profile category: "<<ProfileCategoryUtils::category2string(profileCategory)<<std::endl;
			}
		}
	private:
		typename ImageType::ConstPointer inputImage;
		typename LinearInterpolateImageFunctionType::Pointer linerInterpolatorFunc;
		typename MedianImageFunctionType::Pointer medianFunc;
		typename MeanImageFunctionType::Pointer meanFunc;
		typename GradientImageFunctionType::Pointer gradientFunc;
		typename VarianceImageFunctionType::Pointer varianceFunc;

		VectorType normal;
		std::vector<PointType> pointSet;

		bool enableCacheFlag;
		CategorizedFeatureSetMap cachedFeaturesSetMap;

		PointType currentPoint;
		VectorType offsetVector;
		itk::SimplexMeshGeometry * geoData;

		void addPoints(PROFILE_CATEGORY profileCategory, int profile_dimension, double profile_spacing)
		{
			int shiftDim;
			if (profileCategory == PLAIN){
				shiftDim = profile_dimension-1;
			}else if (profileCategory == BOUNDARY){
				shiftDim = profile_dimension-1;
			}else if (profileCategory == LIVER){
				shiftDim = profile_dimension-1;
			}else{
				shiftDim = profile_dimension-1;
			}

			PointType startPt = this->currentPoint - normal*shiftDim*profile_spacing;
			for ( int i=0;i<profile_dimension;i++ ){
				pointSet.push_back(startPt + normal*profile_spacing*i);
			}
		}

		double mappingItensity( double val )
		{
			double intensity = val;
			intensity = intensity - 0.5*( g_liverThresholds[0].first + g_liverThresholds[0].second );

			return intensity;
		}

		double extractIntensity( PointType & pt )
		{
			double grayValue = -1024;
			if (this->linerInterpolatorFunc->IsInsideBuffer( pt ))
			{
				grayValue = this->linerInterpolatorFunc->Evaluate( pt );
				grayValue = mappingItensity(grayValue);
			}
			return grayValue;
		}

		void extractIntensitySet( FeatureSetType & featureSet )
		{
			for (int i=0;i<pointSet.size();i++)
			{
				double greyValue = extractIntensity(pointSet[i]);
				featureSet.push_back(greyValue);
			}
		}

		void extractIntensityStatistics( FeatureSetType & featureSet )
		{
			std::vector<double> gray_list;

			double gray_mean, gray_sd, gray_contrast, norm = 0;;
			gray_mean = gray_sd = gray_contrast = norm = 0;
			for (int i=0;i<pointSet.size();i++)
			{
				double grayValue = extractIntensity(pointSet[i]);

				gray_list.push_back(grayValue);
				featureSet.push_back(grayValue);

				norm += grayValue*grayValue;
				gray_mean += grayValue;
				if (i<pointSet.size()/2){
					gray_contrast += grayValue;
				}else{
					gray_contrast -= grayValue;
				}
			}

			//featureSet.push_back(std::sqrt(norm));

			gray_mean /= gray_list.size();
			featureSet.push_back(gray_mean);
			featureSet.push_back(gray_contrast);

			for (int i=0;i<gray_list.size();i++)
			{
				gray_sd += (gray_list[i]-gray_mean)*(gray_list[i]-gray_mean);
			}
			gray_sd = std::sqrt( gray_sd / gray_list.size() );
			featureSet.push_back(gray_sd);
		}

		void extractMean( FeatureSetType & featureSet )
		{
			double meanValue = -1024;
			if (this->meanFunc->IsInsideBuffer(pointSet[i]))
			{
				meanValue = this->meanFunc->Evaluate(pointSet[i]);
				meanValue = mappingItensity(meanValue);
			}
			featureSet.push_back(meanValue);
		}

		void extractMedian( FeatureSetType & featureSet )
		{
			double medianValue = -1024;
			if (this->medianFunc->IsInsideBuffer(pointSet[i]))
			{
				medianValue = this->medianFunc->Evaluate(pointSet[i]);
				medianValue = mappingItensity(medianValue);
			}
			featureSet.push_back(medianValue);
		}

		typename GradientImageType::PixelType extractGradient( PointType & pt )
		{
			GradientImageType::PixelType gradientValue = 0.0;
			if (this->gradientFunc->IsInsideBuffer( pt ))
			{
				gradientValue = this->gradientFunc->Evaluate( pt );
			}
			return gradientValue;
		}

		void extractGradientSet( FeatureSetType & featureSet )
		{
			for (int i=0;i<pointSet.size();i++)
			{
				GradientImageType::PixelType gradientValue = extractGradient(pointSet[i]);
				featureSet.push_back(gradientValue);
			}
		}

		void extractGradientMoment( FeatureSetType & featureSet )
		{
			for (int i=0;i<pointSet.size();i++)
			{
				GradientImageType::PixelType gradientValue = extractGradient(pointSet[i]);
				featureSet.push_back(std::pow(2.0, i)*gradientValue);
			}
		}

		void extractTangentFeatureSet( FeatureSetType & featureSet )
		{
			featureSet.push_back(extractGradient(geoData->pos+offsetVector));
			featureSet.push_back(extractGradient(geoData->neighbors[0]+offsetVector));
			featureSet.push_back(extractGradient(geoData->neighbors[1]+offsetVector));
			featureSet.push_back(extractGradient(geoData->neighbors[2]+offsetVector));
		}

		double extractVariance( PointType & pt )
		{
			double varianceValue = 0.0;
			if (this->varianceFunc->IsInsideBuffer( pt ))
			{
				varianceValue = this->varianceFunc->Evaluate( pt );
			}
			return varianceValue;
		}

		void extractVarianceSet( FeatureSetType & featureSet )
		{
			featureSet.push_back(extractVariance(currentPoint));
		}

		void extractCoordinate( FeatureSetType & featureSet )
		{
			VectorType coordiOffset = this->currentPoint - g_liverCentroid;
			for (int i=0;i<ImageDimension;i++)
			{
				featureSet.push_back(coordiOffset[i]);
			}

			VectorType coordiNormal = normal*coordiOffset.GetNorm();
			for (int i=0;i<ImageDimension;i++)
			{
				featureSet.push_back(coordiNormal[i]);
			}
		}

		void getCachedFeatures( FeatureSetType & featureSet, PROFILE_CATEGORY profileCategory, const OffsetValueType & offset )
		{
			CategorizedFeatureSetMap::iterator cachedFeaturesSetMapIt = this->cachedFeaturesSetMap.find(profileCategory);
			CategorizedFeatureSetMap::iterator cachedFeaturesSetMapItEnd = this->cachedFeaturesSetMap.end();

			if (cachedFeaturesSetMapIt == cachedFeaturesSetMapItEnd)
			{
				FeatureSetMap * newFeatureSetMap = new FeatureSetMap;
				cachedFeaturesSetMap[profileCategory] = newFeatureSetMap;

				FeatureSetType * newFeatureSet = new FeatureSetType;
				(*newFeatureSetMap)[offset] = newFeatureSet;
				return;
			}

			FeatureSetMap & featureSetMap = *(cachedFeaturesSetMapIt->second);

			FeatureSetMap::iterator featureSetIt, featureSetItEnd;
			featureSetIt = featureSetMap.find(offset);
			featureSetItEnd = featureSetMap.end();

			if (featureSetIt == featureSetItEnd)
			{
				FeatureSetType * newFeatureSet = new FeatureSetType;
				featureSetMap[offset] = newFeatureSet;
				return;
			}
			else
			{
				FeatureSetType& cachedFeatureSet = *(featureSetIt->second);
				for (int t=0;t<cachedFeatureSet.size();t++)
				{
					featureSet.push_back(cachedFeatureSet[t]);
				}
			}
		}

		void setCachedFeatures( FeatureSetType & featureSet, PROFILE_CATEGORY profileCategory, const OffsetValueType & offset )
		{
			FeatureSetMap & featureSetMap = *(cachedFeaturesSetMap[profileCategory]);
			FeatureSetType & cachedFeatures = *(featureSetMap[offset]);
			for (int t=0;t<featureSet.size();t++)
			{
				cachedFeatures.push_back(featureSet[t]);
			}
		}
	};
}

#endif