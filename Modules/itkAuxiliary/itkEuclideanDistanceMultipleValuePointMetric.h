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
#ifndef __itkEuclideanDistanceMultipleValuePointMetric_h
#define __itkEuclideanDistanceMultipleValuePointMetric_h

#include "itkPointSetToPointSetMetric.h"
#include "itkCovariantVector.h"
#include "itkPointSet.h"
#include "itkImage.h"

namespace itk
{
	/** \class EuclideanDistanceMultipleValuePointMetric
	* \brief Computes the minimum distance between a moving point-set
	*  and a fixed point-set. A vector of minimum closest point distance is
	*  created for each point in the moving point-set.
	*  No correspondance is needed.
	*  For speed consideration, the point-set with the minimum number of points
	*  should be used as the moving point-set.
	*  If the number of points is high, the possibility of setting a distance map
	*  should improve the speed of the closest point computation.
	*
	*  Reference: "A Method for Registration of 3-D Shapes",
	*             IEEE PAMI, Vol 14, No. 2, February 1992
	*
	* \ingroup RegistrationMetrics
	* \ingroup ITKRegistrationCommon
	*/
	template< class TFixedPointSet, class TMovingPointSet = TFixedPointSet>
	class ITK_EXPORT EuclideanDistanceMultipleValuePointMetric:
		public PointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
	{
	public:

		/** Standard class typedefs. */
		typedef EuclideanDistanceMultipleValuePointMetric                   Self;
		typedef PointSetToPointSetMetric< TFixedPointSet, TMovingPointSet > Superclass;

		typedef SmartPointer< Self >       Pointer;
		typedef SmartPointer< const Self > ConstPointer;

		/** Method for creation through the object factory. */
		itkNewMacro(Self);

		/** Run-time type information (and related methods). */
		itkTypeMacro(EuclideanDistanceMultipleValuePointMetric, Object);

		itkStaticConstMacro( PointDimension, unsigned int,
			TFixedPointSet::PointDimension );

		/** Types transferred from the base class */
		typedef typename Superclass::TransformType           TransformType;
		typedef typename Superclass::TransformPointer        TransformPointer;
		typedef typename Superclass::TransformParametersType TransformParametersType;
		typedef typename Superclass::TransformJacobianType   TransformJacobianType;

		typedef typename Superclass::MeasureType                MeasureType;
		typedef typename Superclass::DerivativeType             DerivativeType;
		typedef typename Superclass::FixedPointSetType          FixedPointSetType;
		typedef typename Superclass::MovingPointSetType         MovingPointSetType;
		typedef typename Superclass::FixedPointSetConstPointer  FixedPointSetConstPointer;
		typedef typename Superclass::MovingPointSetConstPointer MovingPointSetConstPointer;

		typedef typename Superclass::PointIterator     PointIterator;
		typedef typename Superclass::PointDataIterator PointDataIterator;

		typedef itk::VariableLengthVector<double> VariableVectorType;

		/** Get the number of values */
		unsigned int GetNumberOfValues() const;

		/** Get the derivatives of the match measure. */
		void GetDerivative(const TransformParametersType & parameters,
			DerivativeType & Derivative) const;

		/**  Get the value for single valued optimizers. */
		MeasureType GetValue(const TransformParametersType & parameters) const;

		/**  Get value and derivatives for multiple valued optimizers. */
		void GetValueAndDerivative(const TransformParametersType & parameters,
			MeasureType & Value, DerivativeType & Derivative) const;

		/** Set/Get if the distance should be squared. Default is true for computation
		speed */
		itkSetMacro(ComputeSquaredDistance, bool);
		itkGetConstMacro(ComputeSquaredDistance, bool);
		itkBooleanMacro(ComputeSquaredDistance);

		void SetNumberOfShapeParameters(unsigned int num)
		{
			m_NumberOfShapeParameters = num;
		}

		void SetWeightForShapePenalty(double weight)
		{
			m_WeightForShapePenalty = weight;
		}

		virtual void Initialize(void) throw ( ExceptionObject );
	protected:
		EuclideanDistanceMultipleValuePointMetric();
		virtual ~EuclideanDistanceMultipleValuePointMetric();

		/** PrintSelf function */
		void PrintSelf(std::ostream & os, Indent indent) const;

	private:
		EuclideanDistanceMultipleValuePointMetric(const Self &); //purposely not implemented
		void operator=(const Self &);               //purposely not implemented

		bool               m_ComputeSquaredDistance;

		unsigned int       m_NumberOfShapeParameters;
		double             m_WeightForShapePenalty;
	};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkEuclideanDistanceMultipleValuePointMetric.hxx"
#endif

#endif
