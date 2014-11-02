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
#ifndef __itkMinValuePointSetToImageMultipleValueMetric_h
#define __itkMinValuePointSetToImageMultipleValueMetric_h

#include "itkPointSetToImageMultipleValueMetric.h"
#include "itkCovariantVector.h"
#include "itkPointSet.h"
#include "itkImage.h"
#include "itkVariableLengthVector.h"

namespace itk
{
	/** \class MinValuePointSetToImageMultipleValueMetric
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
	template< class TFixedPointSet, class TMovingImage>
	class ITK_EXPORT MinValuePointSetToImageMultipleValueMetric:
		public PointSetToImageMultipleValueMetric< TFixedPointSet, TMovingImage >
	{
	public:

		/** Standard class typedefs. */
		typedef MinValuePointSetToImageMultipleValueMetric                   Self;
		typedef PointSetToImageMultipleValueMetric< TFixedPointSet, TMovingImage > Superclass;

		typedef SmartPointer< Self >       Pointer;
		typedef SmartPointer< const Self > ConstPointer;

		/** Method for creation through the object factory. */
		itkNewMacro(Self);

		/** Run-time type information (and related methods). */
		itkTypeMacro(MinValuePointSetToImageMultipleValueMetric, Object);

		/** Types transferred from the base class */
		typedef typename Superclass::TransformType           TransformType;
		typedef typename Superclass::TransformPointer        TransformPointer;
		typedef typename Superclass::TransformParametersType TransformParametersType;
		typedef typename Superclass::TransformJacobianType   TransformJacobianType;

		typedef typename Superclass::MeasureType                MeasureType;
		typedef typename Superclass::DerivativeType             DerivativeType;
		typedef typename Superclass::FixedPointSetType          FixedPointSetType;
		typedef typename Superclass::MovingImageType            MovingImageType;
		typedef typename Superclass::FixedPointSetConstPointer  FixedPointSetConstPointer;
		typedef typename Superclass::MovingImageConstPointer    MovingImageConstPointer;

		typedef typename Superclass::PointIterator     PointIterator;
		typedef typename Superclass::PointDataIterator PointDataIterator;

		typedef itk::VariableLengthVector<double> VariableVectorType;

		typedef typename FixedPointSetType::PointType PointType;

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

		itkSetMacro(Weight, bool);
		itkGetConstMacro(Weight, bool);
		itkBooleanMacro(Weight);

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
		MinValuePointSetToImageMultipleValueMetric();
		virtual ~MinValuePointSetToImageMultipleValueMetric();

		/** PrintSelf function */
		void PrintSelf(std::ostream & os, Indent indent) const;

	private:
		MinValuePointSetToImageMultipleValueMetric(const Self &); //purposely not implemented
		void operator=(const Self &);               //purposely not implemented

		bool               m_ComputeSquaredDistance;
		bool               m_Weight;

		unsigned int       m_NumberOfShapeParameters;
		double             m_WeightForShapePenalty;

		mutable double     m_MeanVal;
		mutable int        m_Iterations;
	};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMinValuePointSetToImageMultipleValueMetric.hxx"
#endif

#endif
