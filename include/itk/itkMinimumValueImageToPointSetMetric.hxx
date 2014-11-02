#ifndef __itkGradientMagnitudeImageToPointSetMetric_hxx
#define __itkGradientMagnitudeImageToPointSetMetric_hxx

#include "itkMinimumValueImageToPointSetMetric.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{
	/**
	* Constructor
	*/
	template< class TFixedImage, class TMovingPointSet >
	MinimumValueImageToPointSetMetric<TFixedImage, TMovingPointSet>
		::MinimumValueImageToPointSetMetric()
	{
		m_SubtractMean = false;
		m_Interpolator = InterpolatorType::New();
	}

	template< class TFixedImage, class TMovingPointSet >
	void
		MinimumValueImageToPointSetMetric<TFixedImage, TMovingPointSet>
		::Initialize(void)
	{
		Superclass::Initialize();
	}

	/**
	* Get the match Measure
	*/
	template< class TFixedImage, class TMovingPointSet >
	typename MinimumValueImageToPointSetMetric<TFixedImage, TMovingPointSet>::MeasureType
		MinimumValueImageToPointSetMetric<TFixedImage, TMovingPointSet>
		::GetValue(const TransformParametersType & parameters) const
	{
		MeasureType measure;
		DerivativeType derivative;

		GetValueAndDerivative( parameters, measure, derivative );

		return measure;
	}

	/**
	* Get the Derivative Measure
	*/
	template< class TFixedImage, class TMovingPointSet >
	void
		MinimumValueImageToPointSetMetric<TFixedImage, TMovingPointSet>
		::GetDerivative(const TransformParametersType & parameters,
		DerivativeType & derivative) const
	{
		MeasureType measure;

		GetValueAndDerivative( parameters, measure, derivative );
	}

	/*
	* Get both the match Measure and theDerivative Measure
	*/
	template< class TFixedImage, class TMovingPointSet >
	void
		MinimumValueImageToPointSetMetric<TFixedImage, TMovingPointSet>
		::GetValueAndDerivative(const TransformParametersType & parameters,
		MeasureType & value, DerivativeType  & derivative) const
	{
		MovingPointSetConstPointer movingPointSet = this->GetMovingPointSet();
		if( !movingPointSet )
		{
			itkExceptionMacro(<< "Fixed point set has not been assigned");
		}
		PointIterator pointItr = movingPointSet->GetPoints()->Begin();
		PointIterator pointEnd = movingPointSet->GetPoints()->End();
		value = 0;
		this->m_NumberOfPixelsCounted = 0;
		this->SetTransformParameters(parameters);

		while( pointItr != pointEnd )
		{
			typename Superclass::InputPointType inputPoint;
			inputPoint.CastFrom( pointItr.Value() );
			typename Superclass::OutputPointType transformedPoint = this->m_Transform->TransformPoint(inputPoint);
			double minimumDistance = NumericTraits< double >::max();
			typename FixedImageType::IndexType index;

			if ( m_FixedImage->TransformPhysicalPointToIndex(transformedPoint, index) )
			{
				minimumDistance = m_FixedImage->GetPixel(index);
				if ( minimumDistance < 0.0 )
				{
					minimumDistance = -minimumDistance;
				}
			}

			value += minimumDistance;
			m_NumberOfPixelsCounted++;
			++pointItr;
		}
		value /= m_NumberOfPixelsCounted;
	}

	template <class TFixedImage, class TMovingImage>
	void
		MinimumValueImageToPointSetMetric<TFixedImage, TMovingImage>
		::PrintSelf(std::ostream & os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
		os << indent << "SubtractMean: " << m_SubtractMean << std::endl;
	}

} // end namespace itk

#endif
