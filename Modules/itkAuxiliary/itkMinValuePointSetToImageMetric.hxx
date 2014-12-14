/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkMinValuePointSetToImageMetric_hxx
#define __itkMinValuePointSetToImageMetric_hxx

#include "itkMinValuePointSetToImageMetric.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "kmMath.h"

namespace itk
{
/** Constructor */
template< class TFixedPointSet, class TMovingPointSet >
MinValuePointSetToImageMetric< TFixedPointSet, TMovingPointSet >
::MinValuePointSetToImageMetric()
{
  // when set to true it will be a bit faster, but it will result in minimizing
  // the sum of distances^4 instead of the sum of distances^2
  m_ComputeSquaredDistance = false;

  m_NumberOfShapeParameters = 5;
  m_WeightForShapePenalty = 10.0;
}

template< class TFixedPointSet, class TMovingPointSet >
MinValuePointSetToImageMetric< TFixedPointSet, TMovingPointSet >
::~MinValuePointSetToImageMetric()
{

}

template< class TFixedPointSet, class TMovingPointSet >
void
MinValuePointSetToImageMetric< TFixedPointSet, TMovingPointSet >
::Initialize(void) throw ( ExceptionObject )
{
	Superclass::Initialize();

	const MovingImageType * movingImage = this->GetMovingImage();
	MovingImageType::RegionType region = movingImage->GetLargestPossibleRegion();
	MovingImageType::IndexType  origin = region.GetIndex();
	MovingImageType::SizeType   size   = region.GetSize();
	MovingImageType::IndexType  centerIndex;
	for(int i=0;i<MovingImageType::ImageDimension;i++)
	{
		centerIndex[i] = origin[i] + size[i]/2;
	}

	movingImage->TransformIndexToPhysicalPoint(centerIndex, m_MovingImageCenter);
}

/** Get the match Measure */
template< class TFixedPointSet, class TMovingPointSet >
typename MinValuePointSetToImageMetric< TFixedPointSet, TMovingPointSet >::MeasureType
MinValuePointSetToImageMetric< TFixedPointSet, TMovingPointSet >
::GetValue(const TransformParametersType & parameters) const
{
 MeasureType value;
 DerivativeType derivative;

 this->GetValueAndDerivative( parameters, value, derivative );

 return value;
}

/** Get the Derivative Measure */
template< class TFixedPointSet, class TMovingPointSet >
void
MinValuePointSetToImageMetric< TFixedPointSet, TMovingPointSet >
::GetDerivative( const TransformParametersType & parameters,
                 DerivativeType & derivative ) const
{
	MeasureType value;
	this->GetValueAndDerivative( parameters, value, derivative );
}

/** Get both the match Measure and theDerivative Measure  */
template< class TFixedPointSet, class TMovingPointSet >
void
MinValuePointSetToImageMetric< TFixedPointSet, TMovingPointSet >
::GetValueAndDerivative(const TransformParametersType & parameters,
                        MeasureType & value, DerivativeType  & derivative) const
{
	FixedPointSetConstPointer fixedPointSet = this->GetFixedPointSet();

	if ( !fixedPointSet )
	{
		itkExceptionMacro(<< "Fixed point set has not been assigned");
	}

	MovingImageConstPointer movingImage = this->GetMovingImage();

	if ( !movingImage )
	{
		itkExceptionMacro(<< "Moving image has not been assigned");
	}

	PointIterator pointItr = fixedPointSet->GetPoints()->Begin();
	PointIterator pointEnd = fixedPointSet->GetPoints()->End();

	this->m_Transform->SetParameters(parameters);

	double shapePenalty = 0.0;
	double shapeProbability = 0.0;
	if (m_WeightForShapePenalty > 0)
	{
		double maxParam = 0.0;
		for (int i=0;i<m_NumberOfShapeParameters;i++)
		{
			if (std::abs(parameters[i]) > maxParam)
			{
				maxParam = std::abs(parameters[i]);
			}
		}
		//1e-6 is make sure the log take positive input.
		shapeProbability = log( km::cdf_outside(maxParam) + 1e-10 );
		shapePenalty =  -10.0*shapeProbability;
	}

	value = 0;

	unsigned int identifier = 0;
	while ( pointItr != pointEnd  )
	{
		typename Superclass::InputPointType inputPoint;
		inputPoint.CastFrom( pointItr.Value() );
		typename Superclass::OutputPointType transformedPoint =
			this->m_Transform->TransformPoint(inputPoint);

		double minval = 1024;//NumericTraits< double >::max();
		
		if (m_Interpolator->IsInsideBuffer(transformedPoint))
		{
			minval = m_Interpolator->Evaluate( transformedPoint );
		}
		else
		{
			minval = transformedPoint.EuclideanDistanceTo(m_MovingImageCenter);
		}

		value += minval+m_WeightForShapePenalty*shapePenalty;

		++pointItr;
		++identifier;
	}

	if (identifier>0)
	{
		value /= identifier;
	}
}

/** PrintSelf method */
template< class TFixedPointSet, class TMovingPointSet >
void
MinValuePointSetToImageMetric< TFixedPointSet, TMovingPointSet >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  if ( m_ComputeSquaredDistance )
    {
    os << indent << "m_ComputeSquaredDistance: True" << std::endl;
    }
  else
    {
    os << indent << "m_ComputeSquaredDistance: False" << std::endl;
    }
}
} // end namespace itk

#endif
