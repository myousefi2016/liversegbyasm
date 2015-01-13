#ifndef __kmProfileExtractor_h
#define __kmProfileExtractor_h

#include <vector>
#include <map>
#include <memory>
#include "itkMacro.h"
#include "itkObject.h"
#include "itkNumericTraits.h"
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

#include "kmCommon.h"

namespace km
{
	template<typename TImage, typename TMesh>
	class ProfileExtractor
	{
	public:
		typedef typename TImage                                                             ImageType;
		typedef typename ImageType::PixelType                                               PixelType;
		itkStaticConstMacro(ImageDimension, unsigned int, ImageType::ImageDimension);
		typedef typename ImageType::IndexType                                               IndexType;

		typedef typename TMesh                                                              MeshType;
		typedef typename MeshType::PointType                                                PointType;
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
			gradientMagFilter->SetSigma( g_sigma );
			gradientMagFilter->Update();

			gradientFunc->SetInputImage(gradientMagFilter->GetOutput());

			//std::cout<<g_liverThresholds[0].first<<", "<<g_liverThresholds[0].second<<std::endl;
		}

		void extractFeatureSet( FeatureSetType & featureSet, PROFILE_CATEGORY profileCategory, const itk::SimplexMeshGeometry * geoData, const PointType & cur_pos, int pointId = GLOBAL_POINT_ID )
		{
			featureSet.clear();
			pointSet.clear();

			this->localIntenDist = intensityDistributions[pointId];
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
				this->addPoints(profileCategory, g_profile_dim, g_profile_spacing);
				this->extractIntensitySet(featureSet);
			}else if (profileCategory==LIVER)
			{
				if (enableCacheFlag){
					getCachedFeatures(featureSet, profileCategory, offset);
				}

				if (featureSet.size() == 0){
					this->addPoints(profileCategory, g_profile_dim, g_profile_spacing);
					this->extractBinaryIntensitySet(featureSet);
					this->extractIntensityStatistics(featureSet);
					//this->extractGradientSet(featureSet);
					this->extractGradientMoments(featureSet);
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
					this->addPoints(profileCategory, g_profile_dim, g_profile_spacing);
					this->extractBinaryIntensitySet(featureSet);
					this->extractIntensityStatistics(featureSet);
					//this->extractGradientSet(featureSet);
					this->extractGradientMoments(featureSet);
					this->extractTangentFeatureSet(featureSet);
					if (enableCacheFlag){
						setCachedFeatures(featureSet, profileCategory, offset);
					}
				}
			}else{
				std::cout<<"Unknown profile category: "<<ProfileCategoryUtils::category2string(profileCategory)<<std::endl;
			}
		}

		void extractLocalIntensityRange(const MeshType* mesh)
		{
			for (int i=0;i<mesh->GetNumberOfPoints();i++)
			{
				this->extractLocalIntensityRange(mesh, i);
			}
		}

		void extractLocalIntensityRange(const MeshType* mesh, int pointId)
		{
			std::vector<double> intensities;
			double mu, sigma;
			mu = sigma = 0.0;

			int localRadius = 2;
			MeshType::NeighborListType * localIndexList = mesh->GetNeighbors(pointId, localRadius);
			localIndexList->push_back(pointId);

			MeshType::NeighborListType::iterator localIt = localIndexList->begin();
			while ( localIt != localIndexList->end() )
			{
				int tmpPointId = *localIt++;
				const itk::SimplexMeshGeometry * tmpGeoData = mesh->GetGeometryData()->GetElement(tmpPointId);
				if (tmpGeoData != NULL)
				{
					VectorType tmpNormal;
					tmpNormal.Set_vnl_vector(tmpGeoData->normal.Get_vnl_vector());

					PointType pos = mesh->GetPoint(tmpPointId);
					int localDepth = 10;
					for (int i=1;i<=localDepth;i++){
						PointType samplePos = pos - tmpNormal*i;
						double grayValue = 0;
						if (this->linerInterpolatorFunc->IsInsideBuffer( samplePos )){
							grayValue = this->linerInterpolatorFunc->Evaluate( samplePos );
							grayValue = std::max(grayValue, 0.0);
						}
						intensities.push_back(grayValue);
					}
				}
				else
				{
					std::cerr<<"[extractLocalIntensityRange] Cannot find geometry data for point: "<<pointId<<std::endl;
				}	
			}

			km::Math::calculateMeanAndSigma(intensities, mu, sigma);

			NormalDistribution newDistribution;
			newDistribution.set_mu(mu);
			newDistribution.set_sigma(sigma);
			intensityDistributions[pointId] = newDistribution;
		}

		void extractGlobalIntensityRange(MeshType* mesh)
		{
			//std::vector<double> intensities;
			//double mu, sigma;
			//mu = sigma = 0.0;

			////for (int pointId=0;pointId<mesh->GetNumberOfPoints();pointId++)
			////{
			////	const itk::SimplexMeshGeometry * tmpGeoData = mesh->GetGeometryData()->GetElement(pointId);
			////	if (tmpGeoData != NULL)
			////	{
			////		VectorType tmpNormal;
			////		tmpNormal.Set_vnl_vector(tmpGeoData->normal.Get_vnl_vector());

			////		PointType pos = mesh->GetPoint(pointId);
			////		int N = 10;
			////		for (int i=1;i<=N;i++){
			////			PointType samplePos = pos - tmpNormal*i;
			////			double grayValue = 0;
			////			if (this->linerInterpolatorFunc->IsInsideBuffer( samplePos )){
			////				grayValue = this->linerInterpolatorFunc->Evaluate( samplePos );
			////				grayValue = std::max(grayValue, 0.0);
			////			}
			////			intensities.push_back(grayValue);
			////		}
			////	}
			////	else
			////	{
			////		std::cerr<<"[extractGlobalIntensityRange] Cannot find geometry data for point: "<<pointId<<std::endl;
			////	}
			////}

			//MeshType::Pointer triangleMesh = km::simplexMeshToTriangleMesh<MeshType, MeshType>( mesh );
			//typedef itk::Image<unsigned char, ImageType::ImageDimension> MaskImageType;
			//MaskImageType::Pointer mask = km::generateBinaryFromMesh<MeshType, MaskImageType, ImageType>( triangleMesh, this->inputImage );
			//itk::ImageRegionConstIterator<MaskImageType> it( mask, mask->GetLargestPossibleRegion() );
			//it.GoToBegin();
			//while(!it.IsAtEnd())
			//{
			//	unsigned char val = it.Get();
			//	IndexType idx = it.GetIndex();

			//	if (val == 1){
			//		intensities.push_back(this->inputImage->GetPixel(idx));
			//	}

			//	it++;
			//}

			//km::Math::calculateMeanAndSigma(intensities, mu, sigma);

			//NormalDistribution newDistribution;
			//newDistribution.set_mu(mu);
			//newDistribution.set_sigma(sigma);
			//intensityDistributions[GLOBAL_POINT_ID] = newDistribution;
			//newDistribution.print();

			//globalIntenDist = newDistribution;

			//for (int i=0;i<mesh->GetNumberOfPoints();i++)
			//{
			//	intensityDistributions[i] = newDistribution; //By default, all landmarks use global intensity threshold.
			//}

			double liverThresholdLow, liverThresholdHigh;
			liverThresholdLow = liverThresholdHigh = 0;
			MeshType::Pointer triangleMesh = km::simplexMeshToTriangleMesh<MeshType, MeshType>( mesh );
			typedef itk::Image<unsigned char, ImageType::ImageDimension> MaskImageType;
			MaskImageType::Pointer mask = km::generateBinaryFromMesh<MeshType, MaskImageType, ImageType>( triangleMesh, this->inputImage );
			ImageType::Pointer roi = km::maskImage<ImageType, MaskImageType>( this->inputImage, mask );
			const double statisticMin = 15.0;
			const double statisticMax = 300.0;
			const double binWidth = 3.0;
			const double thresholdSpan = 30;
			Histogram histogram(statisticMin, statisticMax, binWidth);
			km::calculateHistogram<ImageType>(roi, &histogram);
			double minLiverTissueRatio = 0.85;
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

			NormalDistribution newRange;
			newRange.mu = (liverThresholdLow+liverThresholdHigh)/2;
			newRange.sigma = ((liverThresholdHigh-liverThresholdLow)/2)/2.5;
			globalIntenDist = newRange;

			intensityDistributions[GLOBAL_POINT_ID] = newRange;
			for (int i=0;i<mesh->GetNumberOfPoints();i++)
			{
				intensityDistributions[i] = newRange; //By default, all landmarks use global intensity threshold.
			}

			globalIntenDist.print();
		}

		km::NormalDistribution& getIntensityRange(int pointId)
		{
			return intensityDistributions[pointId];
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

		std::map<int, km::NormalDistribution> intensityDistributions;
		NormalDistribution localIntenDist;
		NormalDistribution globalIntenDist;
		static const int GLOBAL_POINT_ID = -1;

		void addPoints(PROFILE_CATEGORY profileCategory, int profile_dimension, double profile_spacing)
		{
			//int shiftDim;
			//if (profileCategory == PLAIN){
			//	shiftDim = profile_dimension/2;
			//}else if (profileCategory == BOUNDARY){
			//	shiftDim = profile_dimension/2;
			//}else if (profileCategory == LIVER){
			//	shiftDim = profile_dimension/2;
			//}else{
			//	shiftDim = profile_dimension-1;
			//}

			PointType startPt = this->currentPoint;
			for ( int i=0;i<profile_dimension;i++ ){
				pointSet.push_back(startPt - normal*profile_spacing*i);
			}
		}

		double extractIntensity( PointType & pt )
		{
			double grayValue = -1024;
			if (this->linerInterpolatorFunc->IsInsideBuffer( pt ))
			{
				grayValue = this->linerInterpolatorFunc->Evaluate( pt );
			}
			return grayValue;
		}

		void extractIntensitySet( FeatureSetType & featureSet )
		{
			for (int i=0;i<pointSet.size();i++)
			{
				double grayValue = extractIntensity(pointSet[i]);
				featureSet.push_back(grayValue - localIntenDist.mu);
			}
		}

		void extractNormalizedIntensitySet( FeatureSetType & featureSet )
		{
			for (int i=0;i<pointSet.size();i++)
			{
				double grayValue = extractIntensity(pointSet[i]);
				featureSet.push_back( (grayValue - localIntenDist.mu)/globalIntenDist.sigma );
			}
		}

		void extractBinaryIntensitySet( FeatureSetType & featureSet )
		{
			for (int i=0;i<pointSet.size();i++)
			{
				double grayValue = extractIntensity(pointSet[i]);
				if (grayValue<localIntenDist.mu-2.5*globalIntenDist.sigma || grayValue>localIntenDist.mu+2.5*globalIntenDist.sigma){
					featureSet.push_back(0);
				}else{
					featureSet.push_back(1);
				}
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
				grayValue -= localIntenDist.mu;

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

		void extractGradient( FeatureSetType & featureSet )
		{
			GradientImageType::PixelType gradientValue = 0.0;
			if (this->gradientFunc->IsInsideBuffer( this->currentPoint ))
			{
				gradientValue = this->gradientFunc->Evaluate( this->currentPoint );
			}
			featureSet.push_back(gradientValue);
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

		void extractGradientMoments( FeatureSetType & featureSet )
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