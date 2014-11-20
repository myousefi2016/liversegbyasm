#ifndef __kmProfileClassifier_hxx
#define __kmProfileClassifier_hxx

#include "kmProfileClassifier.h"

namespace km
{
	template<class MatrixType>
	void fillFlannMatrix( typename MatrixType& matrix, typename MatrixType::type val )
	{
		for (int i=0;i<matrix.rows;i++)
		{
			for (int j=0;j<matrix.cols;j++)
			{
				matrix[i][j] = val;
			}
		}
	}
	
	double mappingItensity( double val )
	{
		bool belongToLiver = false;
		int intensity = val;
		for (int i=0;i<g_liverThresholds.size();i++)
		{
			if (g_liverThresholds[i].first <= val && g_liverThresholds[i].second >= val)
			{
				belongToLiver = true;
				intensity = intensity - 0.5*( g_liverThresholds[i].first + g_liverThresholds[i].second );
				break;
			}
		}

		if (!belongToLiver)
		{
			intensity = intensity - 0.5*( g_liverThresholds[0].first + g_liverThresholds[0].second );
		}

		return intensity;
	}

	template<class T>
	void
		calculateMeanAndStd(std::vector<T> & list, typename T& mean, typename T& sd)
	{
		if (list.size()==0)
		{
			return;
		}

		mean = sd = 0.0;

		double powered_sum = 0;
		for (int i=0;i<list.size();i++)
		{
			powered_sum += list[i]*list[i];
			mean += list[i];
		}

		mean /= list.size();

		if (powered_sum <= 1e-6)
		{
			return;
		}

		for (int i=0;i<list.size();i++)
		{
			sd += (list[i]-mean)*(list[i]-mean);
			//list[i] = std::sqrt( (list[i]*list[i])/powered_sum);
		}

		sd = std::sqrt( sd / list.size() );
	}

	template<class T>
	void
		normalizeList(std::vector<T> & list, typename T& mean, typename T& sd)
	{
		if (list.size()==0)
		{
			return;
		}

		mean = sd = 0.0;

		double powered_sum = 0;
		for (int i=0;i<list.size();i++)
		{
			powered_sum += list[i]*list[i];
			mean += list[i];
		}

		mean /= list.size();

		if (powered_sum <= 1e-6)
		{
			return;
		}

		for (int i=0;i<list.size();i++)
		{
			sd += (list[i]-mean)*(list[i]-mean);
			list[i] = std::sqrt( (list[i]*list[i])/powered_sum);
		}

		sd = std::sqrt( sd / list.size() );
	}

	template<class T>
	void
		normalizeList(std::vector<T> & list)
	{
		T mean, sd;
		mean = sd = 0.0;
		normalizeList<T>( list, mean, sd );
	}

	template<class GradientInterpolatorType, class IntensityInterpolatorType>
	void 
		extractFeature(
			const typename GradientInterpolatorType* gradientInterpolator, 
			const typename IntensityInterpolatorType* intensityInterpolator,
			const itk::SimplexMeshGeometry* geoData,
			const itk::SimplexMeshGeometry::PointType & ipoint,
			std::vector<double> & feature,
			PROFILE_CATEGORY profile_category)
	{
		typedef IntensityInterpolatorType::InputImageType IntensityImageType;
		typedef IntensityImageType::PixelType IntensityPixelType;
		typedef GradientInterpolatorType::InputImageType GradientImageType;
		typedef GradientImageType::PixelType GradientPixelType;

		typedef itk::SimplexMeshGeometry::PointType PointType;
		typedef PointType::VectorType VectorType;

		VectorType normal;
		normal.Set_vnl_vector(geoData->normal.Get_vnl_vector());
		
		std::vector<double> gray_list;
		std::vector<double> gradient_list;

		double defaultIntensity = -1024.0;
		int N = PROFILE_DIM;
		int K = N-1;
		double spa = PROFILE_SPACING;
		PointType startPt = ipoint - normal*(K)*spa;
		PointType profilePos = startPt;

		for ( int i=0;i<N;i++ )
		{
			double val = -1024;
			if (intensityInterpolator->IsInsideBuffer(profilePos))
			{
				val = intensityInterpolator->Evaluate(profilePos);
				val = mappingItensity(val);
			}
			gray_list.push_back(val);

			double gradient = 100.0;
			if (gradientInterpolator->IsInsideBuffer(profilePos))
			{
				GradientPixelType gradientvec = gradientInterpolator->Evaluate(profilePos);
				gradient = dot_product( gradientvec.GetVnlVector(), normal.GetVnlVector());
			}
			gradient_list.push_back(gradient);

			profilePos += normal*spa;
		}

		if (profile_category == PLAIN)
		{
			for (int i=0;i<gray_list.size();i++)
			{
				feature.push_back(gray_list[i]);
			}
		}
		else if (profile_category == LIVER)
		{
			//For distinguish liver tissue and plain background.
			double gray_mean, gray_sd, gray_contrast;
			gray_mean = gray_sd = gray_contrast = 0;
			//Calculate gray statistics
			{
				double powered_sum = 0;
				for (int i=0;i<gray_list.size();i++)
				{
					powered_sum += gray_list[i]*gray_list[i];
					gray_mean += gray_list[i];
					if (i<gray_list.size()/2)
					{
						gray_contrast += gray_list[i];
					}
					else
					{
						gray_contrast -= gray_list[i];
					}
				}

				gray_mean /= gray_list.size();

				for (int i=0;i<gray_list.size();i++)
				{
					gray_sd += (gray_list[i]-gray_mean)*(gray_list[i]-gray_mean);
				}

				gray_sd = std::sqrt( gray_sd / gray_list.size() );
			}
			feature.push_back(gray_mean);
			feature.push_back(gray_sd);
			feature.push_back(gray_contrast);
			for (int i=0;i<gray_list.size();i++)
			{
				feature.push_back(gray_list[i]);
			}

			km::normalizeList<double>(gray_list);
			for (int i=0;i<gradient_list.size();i++)
			{
				feature.push_back(gradient_list[i]);
			}

			//VectorType coordiOffset = ipoint-g_liverCentroid;
			//feature.push_back(coordiOffset.GetNorm());
		}
		else if (profile_category == COORDINATE)
		{
			VectorType coordiOffset = ipoint-g_liverCentroid;//normal*dot_product( (ipoint-g_liverCentroid).GetVnlVector(), normal.GetVnlVector());
			for (int i=0;i<coordiOffset.Dimension;i++)
			{
				feature.push_back(coordiOffset[i]);
			}

			VectorType coordiNormal = normal*coordiOffset.GetNorm();
			for (int i=0;i<coordiNormal.Dimension;i++)
			{
				feature.push_back(coordiNormal[i]);
			}
		}
	}
}
#endif