/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: EuclideanDistanceKdTreeMultipleValuePointSetMetric.hxx,v $
Language:  C++
Date:      $Date: 2008/12/06 05:19:54 $
Version:   $Revision: 1.5 $

Copyright (c) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __EuclideanDistanceKdTreeMultipleValuePointSetMetric_hxx
#define __EuclideanDistanceKdTreeMultipleValuePointSetMetric_hxx

#include "itkEuclideanDistanceKdTreeMultipleValuePointSetMetric.h"
#include "kmUtility.h"

namespace itk {

	template <class TPointSet>
	EuclideanDistanceKdTreeMultipleValuePointSetMetric<TPointSet>
		::EuclideanDistanceKdTreeMultipleValuePointSetMetric()
	{
		this->m_UseWithRespectToTheMovingPointSet = true;

		m_NumberOfShapeParameters = 5;
		m_WeightForShapePenalty = 0.1;

		this->m_KdTreeGenerator = NULL;
		this->m_SamplePoints = NULL;

		typename DefaultTransformType::Pointer transform
			= DefaultTransformType::New();
		transform->SetIdentity();

		Superclass::SetTransform( transform );
	}

	/** Initialize the metric */
	template <class TPointSet>
	void
		EuclideanDistanceKdTreeMultipleValuePointSetMetric<TPointSet>
		::Initialize( void ) throw ( ExceptionObject )
	{
		Superclass::Initialize();

		/**
		* Generate KdTrees for the opposite point set
		*/
		this->m_SamplePoints = SampleType::New();
		this->m_SamplePoints->SetMeasurementVectorSize( PointDimension );

		typename PointSetType::PointsContainerConstIterator It
			= this->m_FixedPointSet->GetPoints()->Begin();
		while( It != this->m_FixedPointSet->GetPoints()->End() )
		{
			PointType point = It.Value();
			this->m_SamplePoints->PushBack( point );
			++It;
		}

		this->m_KdTreeGenerator = TreeGeneratorType::New();
		this->m_KdTreeGenerator->SetSample( this->m_SamplePoints );
		this->m_KdTreeGenerator->SetBucketSize( 4 );
		this->m_KdTreeGenerator->Update();
	}

	/** Return the number of values, i.e the number of points in the moving set */
	template <class TPointSet>
	unsigned int
		EuclideanDistanceKdTreeMultipleValuePointSetMetric<TPointSet>
		::GetNumberOfValues() const
	{
		if( this->m_MovingPointSet )
		{
			return this->m_MovingPointSet->GetPoints()->Size();
		}

		return 0;
	}


	/** Get the match Measure */
	template <class TPointSet>
	typename EuclideanDistanceKdTreeMultipleValuePointSetMetric
		<TPointSet>::MeasureType
		EuclideanDistanceKdTreeMultipleValuePointSetMetric<TPointSet>
		::GetValue( const TransformParametersType & parameters ) const
	{
		this->m_Transform->SetParameters( parameters );

		PointSetConstPointer fixedPointSet = this->GetFixedPointSet();

		if ( !fixedPointSet )
		{
			itkExceptionMacro(<< "Fixed point set has not been assigned");
		}

		PointSetConstPointer movingPointSet = this->GetMovingPointSet();

		if ( !movingPointSet )
		{
			itkExceptionMacro(<< "Moving point set has not been assigned");
		}
		
		MeasureType measure;
		measure.set_size( this->GetNumberOfValues() );

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

		PointIterator pointItr = movingPointSet->GetPoints()->Begin();
		PointIterator pointEnd = movingPointSet->GetPoints()->End();

		unsigned int identifier = 0;
		while( pointItr != pointEnd )
		{
			PointType point = pointItr.Value();
			PointType transformedPoint = this->m_Transform->TransformPoint(point);

			typename TreeGeneratorType::KdTreeType::InstanceIdentifierVectorType neighbors;
			this->m_KdTreeGenerator->GetOutput()->Search( transformedPoint, 1u, neighbors );

			PointType neighbor =
				this->m_KdTreeGenerator->GetOutput()->GetMeasurementVector( neighbors[0] );

			RealType dist = 0.0;
			for( unsigned int d = 0; d < PointDimension; d++ )
			{
				dist += vnl_math_sqr( neighbor[d] - transformedPoint[d] );
			}

			dist = dist + m_WeightForShapePenalty*shapePenalty;

			measure.put(identifier, dist);

			++identifier;
			++pointItr;
		}

		return measure;
	}

	/** Get the Derivative Measure */
	template <class TPointSet>
	void
		EuclideanDistanceKdTreeMultipleValuePointSetMetric<TPointSet>
		::GetDerivative( const TransformParametersType & parameters,
		DerivativeType & derivative ) const
	{
		/*std::cout<<"GetDerivative()!"<<std::endl;

		PointSetPointer points[2];

		if( this->m_UseWithRespectToTheMovingPointSet )
		{
			points[0] = const_cast<PointSetType *>(
				static_cast<const PointSetType *>( this->m_FixedPointSet ) );
			points[1] = const_cast<PointSetType *>(
				static_cast<const PointSetType *>( this->m_MovingPointSet ) );
		}
		else
		{
			points[1] = const_cast<PointSetType *>(
				static_cast<const PointSetType *>( this->m_FixedPointSet ) );
			points[0] = const_cast<PointSetType *>(
				static_cast<const PointSetType *>( this->m_MovingPointSet ) );
		}
		derivative.SetSize( points[1]->GetPoints()->Size(), PointDimension );
		derivative.Fill( 0 );

		unsigned long count = 0;
		typename PointSetType::PointsContainerConstIterator It
			= points[1]->GetPoints()->Begin();
		while( It != points[1]->GetPoints()->End() )
		{
			PointType point = It.Value();

			MeasurementVectorType queryPoint;
			for( unsigned int d = 0; d < PointDimension; d++ )
			{
				queryPoint[d] = point[d];
			}
			typename TreeGeneratorType::KdTreeType
				::InstanceIdentifierVectorType neighbors;
			this->m_KdTreeGenerator->GetOutput()->Search( queryPoint, 1u, neighbors );

			MeasurementVectorType neighbor =
				this->m_KdTreeGenerator->GetOutput()->GetMeasurementVector( neighbors[0] );

			for( unsigned int d = 0; d < PointDimension; d++ )
			{
				derivative( count, d ) = -( neighbor[d] - queryPoint[d] );
			}
			count++;
			++It;
		}*/
	}

	/** Get both the match Measure and theDerivative Measure  */
	template <class TPointSet>
	void
		EuclideanDistanceKdTreeMultipleValuePointSetMetric<TPointSet>
		::GetValueAndDerivative( const TransformParametersType & parameters,
		MeasureType & value, DerivativeType  & derivative ) const
	{
		/*std::cout<<"GetValueAndDerivative()!"<<std::endl;

		PointSetPointer points[2];

		if( this->m_UseWithRespectToTheMovingPointSet )
		{
			points[0] = const_cast<PointSetType *>(
				static_cast<const PointSetType *>( this->m_FixedPointSet ) );
			points[1] = const_cast<PointSetType *>(
				static_cast<const PointSetType *>( this->m_MovingPointSet ) );
		}
		else
		{
			points[1] = const_cast<PointSetType *>(
				static_cast<const PointSetType *>( this->m_FixedPointSet ) );
			points[0] = const_cast<PointSetType *>(
				static_cast<const PointSetType *>( this->m_MovingPointSet ) );
		}
		derivative.SetSize( points[1]->GetPoints()->Size(), PointDimension );
		derivative.Fill( 0.0 );

		value.SetSize( 1 );
		value.Fill( 0.0 );

		unsigned long count = 0;
		typename PointSetType::PointsContainerConstIterator It
			= points[1]->GetPoints()->Begin();
		while( It != points[1]->GetPoints()->End() )
		{
			PointType point = It.Value();

			MeasurementVectorType queryPoint;
			for( unsigned int d = 0; d < PointDimension; d++ )
			{
				queryPoint[d] = point[d];
			}
			typename TreeGeneratorType::KdTreeType
				::InstanceIdentifierVectorType neighbors;
			this->m_KdTreeGenerator->GetOutput()->Search( queryPoint, 1u, neighbors );

			MeasurementVectorType neighbor =
				this->m_KdTreeGenerator->GetOutput()->GetMeasurementVector( neighbors[0] );

			RealType sum = 0.0;
			for( unsigned int d = 0; d < PointDimension; d++ )
			{
				derivative( count, d ) = -( neighbor[d] - queryPoint[d] );
				sum += vnl_math_sqr( neighbor[d] - queryPoint[d] );
			}
			value[0] += vcl_sqrt( sum );

			count++;
			++It;
		}

		value[0] /= static_cast<RealType>( points[1]->GetNumberOfPoints() );*/
	}

	template <class TPointSet>
	void
		EuclideanDistanceKdTreeMultipleValuePointSetMetric<TPointSet>
		::PrintSelf( std::ostream& os, Indent indent ) const
	{
		Superclass::PrintSelf( os, indent );

		os << indent << "Use with respect to the moving point set: "
			<< this->m_UseWithRespectToTheMovingPointSet << std::endl;
	}

} // end namespace itk


#endif
