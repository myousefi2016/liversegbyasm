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
#ifndef __itkMinValuePointSetToImageMultipleValueMetric_hxx
#define __itkMinValuePointSetToImageMultipleValueMetric_hxx

#include "itkMinValuePointSetToImageMultipleValueMetric.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "kmUtility.h"

namespace itk
{
/** Constructor */
template< class TFixedPointSet, class TMovingPointSet >
MinValuePointSetToImageMultipleValueMetric< TFixedPointSet, TMovingPointSet >
::MinValuePointSetToImageMultipleValueMetric()
{
  // when set to true it will be a bit faster, but it will result in minimizing
  // the sum of distances^4 instead of the sum of distances^2
  m_ComputeSquaredDistance = false;

  m_NumberOfShapeParameters = 5;
  m_WeightForShapePenalty = 10.0;

  m_MeanVal = 0;
  m_Iterations = 0;
}

template< class TFixedPointSet, class TMovingPointSet >
MinValuePointSetToImageMultipleValueMetric< TFixedPointSet, TMovingPointSet >
::~MinValuePointSetToImageMultipleValueMetric()
{

}

template< class TFixedPointSet, class TMovingPointSet >
void
MinValuePointSetToImageMultipleValueMetric< TFixedPointSet, TMovingPointSet >
::Initialize(void) throw ( ExceptionObject )
{
	Superclass::Initialize();
}

/** Return the number of values, i.e the number of points in the moving set */
template< class TFixedPointSet, class TMovingPointSet >
unsigned int
MinValuePointSetToImageMultipleValueMetric< TFixedPointSet, TMovingPointSet >
::GetNumberOfValues() const
{
  FixedPointSetConstPointer fixedPointSet = this->GetFixedPointSet();

  if ( !fixedPointSet )
    {
    itkExceptionMacro(<< "Moving point set has not been assigned");
    }

  return fixedPointSet->GetPoints()->Size();
}

/** Get the match Measure */
template< class TFixedPointSet, class TMovingPointSet >
typename MinValuePointSetToImageMultipleValueMetric< TFixedPointSet, TMovingPointSet >::MeasureType
MinValuePointSetToImageMultipleValueMetric< TFixedPointSet, TMovingPointSet >
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
MinValuePointSetToImageMultipleValueMetric< TFixedPointSet, TMovingPointSet >
::GetDerivative( const TransformParametersType & parameters,
                 DerivativeType & derivative ) const
{
	MeasureType value;
	this->GetValueAndDerivative( parameters, value, derivative );
}

/** Get both the match Measure and theDerivative Measure  */
template< class TFixedPointSet, class TMovingPointSet >
void
MinValuePointSetToImageMultipleValueMetric< TFixedPointSet, TMovingPointSet >
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
	value.set_size( this->GetNumberOfValues() );
	derivative.set_size( this->GetNumberOfValues(), m_Transform->GetNumberOfParameters()  );

	//value.Fill(0);
	derivative.Fill( 0 );

	double shapeProbability = 0.0;
	km::NormalDistribution shapeDistribution(0,1);
	if (m_WeightForShapePenalty > 0)
	{
		double maxParam = 0.0;
		for (int i=0;i<m_NumberOfShapeParameters;i++)
		{
			//if (std::abs(parameters[i]) > maxParam)
			//{
			//	maxParam = std::abs(parameters[i]);
			//}
			shapeProbability += log( shapeDistribution.cdf_outside(std::abs(parameters[i])) + 1e-10 );
		}
		shapeProbability /= m_NumberOfShapeParameters;
		//1e-6 is make sure the log take positive input.
	}

	if (m_Iterations == 0)
	{
		double sumVal = 0.0;
		while ( pointItr != pointEnd  )
		{
			typename Superclass::InputPointType inputPoint;
			inputPoint.CastFrom( pointItr.Value() );
			typename Superclass::OutputPointType transformedPoint =
				this->m_Transform->TransformPoint(inputPoint);

			double minval = m_MeanVal;//NumericTraits< double >::max();

			if (m_Interpolator->IsInsideBuffer(transformedPoint))
			{
				minval = m_Interpolator->Evaluate( transformedPoint );
			}

			sumVal += minval;

			++pointItr;
		}
		m_MeanVal = sumVal/this->GetNumberOfValues();
	}

	unsigned int identifier = 0;
	pointItr = fixedPointSet->GetPoints()->Begin();
	while ( pointItr != pointEnd  )
	{
		typename Superclass::InputPointType inputPoint;
		inputPoint.CastFrom( pointItr.Value() );
		typename Superclass::OutputPointType transformedPoint =
			this->m_Transform->TransformPoint(inputPoint);

		double minval = m_MeanVal;//NumericTraits< double >::max();
		
		if (m_Interpolator->IsInsideBuffer(transformedPoint))
		{
			minval = m_Interpolator->Evaluate( transformedPoint );
		}

		if (m_WeightForShapePenalty > 0.0)
		{
			minval -= m_WeightForShapePenalty*m_MeanVal*shapeProbability;
		}

		value.put(identifier, minval);

		++pointItr;
		++identifier;
	}

	m_Iterations++;
}

/** PrintSelf method */
template< class TFixedPointSet, class TMovingPointSet >
void
MinValuePointSetToImageMultipleValueMetric< TFixedPointSet, TMovingPointSet >
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
