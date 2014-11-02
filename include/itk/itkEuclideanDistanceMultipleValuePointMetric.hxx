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
#ifndef __itkEuclideanDistanceMultipleValuePointMetric_hxx
#define __itkEuclideanDistanceMultipleValuePointMetric_hxx

#include "itkEuclideanDistanceMultipleValuePointMetric.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "kmMath.h"

namespace itk
{
/** Constructor */
template< class TFixedPointSet, class TMovingPointSet >
EuclideanDistanceMultipleValuePointMetric< TFixedPointSet, TMovingPointSet >
::EuclideanDistanceMultipleValuePointMetric()
{
  // when set to true it will be a bit faster, but it will result in minimizing
  // the sum of distances^4 instead of the sum of distances^2
  m_ComputeSquaredDistance = false;

  m_NumberOfShapeParameters = 5;
  m_WeightForShapePenalty = 0.0;
}

template< class TFixedPointSet, class TMovingPointSet >
EuclideanDistanceMultipleValuePointMetric< TFixedPointSet, TMovingPointSet >
::~EuclideanDistanceMultipleValuePointMetric()
{

}

template< class TFixedPointSet, class TMovingPointSet >
void
EuclideanDistanceMultipleValuePointMetric< TFixedPointSet, TMovingPointSet >
::Initialize(void) throw ( ExceptionObject )
{
	Superclass::Initialize();
}

/** Return the number of values, i.e the number of points in the moving set */
template< class TFixedPointSet, class TMovingPointSet >
unsigned int
EuclideanDistanceMultipleValuePointMetric< TFixedPointSet, TMovingPointSet >
::GetNumberOfValues() const
{
  MovingPointSetConstPointer movingPointSet = this->GetMovingPointSet();

  if ( !movingPointSet )
    {
    itkExceptionMacro(<< "Moving point set has not been assigned");
    }

  return movingPointSet->GetPoints()->Size();
}

/** Get the match Measure */
template< class TFixedPointSet, class TMovingPointSet >
typename EuclideanDistanceMultipleValuePointMetric< TFixedPointSet, TMovingPointSet >::MeasureType
EuclideanDistanceMultipleValuePointMetric< TFixedPointSet, TMovingPointSet >
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
EuclideanDistanceMultipleValuePointMetric< TFixedPointSet, TMovingPointSet >
::GetDerivative( const TransformParametersType & parameters,
                 DerivativeType & derivative ) const
{
	MeasureType value;
	this->GetValueAndDerivative( parameters, value, derivative );
}

/** Get both the match Measure and theDerivative Measure  */
template< class TFixedPointSet, class TMovingPointSet >
void
EuclideanDistanceMultipleValuePointMetric< TFixedPointSet, TMovingPointSet >
::GetValueAndDerivative(const TransformParametersType & parameters,
                        MeasureType & value, DerivativeType  & derivative) const
{
	FixedPointSetConstPointer fixedPointSet = this->GetFixedPointSet();

	if ( !fixedPointSet )
	{
		itkExceptionMacro(<< "Fixed point set has not been assigned");
	}

	MovingPointSetConstPointer movingPointSet = this->GetMovingPointSet();

	if ( !movingPointSet )
	{
		itkExceptionMacro(<< "Moving point set has not been assigned");
	}

	PointIterator pointItr = movingPointSet->GetPoints()->Begin();
	PointIterator pointEnd = movingPointSet->GetPoints()->End();

	this->SetTransformParameters(parameters);
	value.set_size( movingPointSet->GetPoints()->Size() );
	derivative.set_size( movingPointSet->GetNumberOfPoints(), m_Transform->GetNumberOfParameters()  );

	double shapePenalty = 0.0;
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
		shapePenalty =  -10.0*log( km::cdf_outside(maxParam) + 1e-10 );
	}

	unsigned int identifier = 0;
	while ( pointItr != pointEnd  )
	{
		typename Superclass::InputPointType inputPoint;
		inputPoint.CastFrom( pointItr.Value() );
		typename Superclass::OutputPointType transformedPoint =
			this->m_Transform->TransformPoint(inputPoint);

		double minimumDistance = 0.0;
		bool   closestPoint = false;

		typename OutputPointType nearestPt = fixedPointSet->GetPoint(pointItr.Index());
		
		minimumDistance = nearestPt.EuclideanDistanceTo( transformedPoint );

		minimumDistance += m_WeightForShapePenalty*shapePenalty;

		value.put(identifier, minimumDistance);

		++pointItr;
		++identifier;
	}
}

/** PrintSelf method */
template< class TFixedPointSet, class TMovingPointSet >
void
EuclideanDistanceMultipleValuePointMetric< TFixedPointSet, TMovingPointSet >
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
