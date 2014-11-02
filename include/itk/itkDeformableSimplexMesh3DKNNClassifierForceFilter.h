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
/*=========================================================================
*
*  Portions of this file are subject to the VTK Toolkit Version 3 copyright.
*
*  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
*
*  For complete copyright, license and disclaimer of warranty information
*  please refer to the NOTICE file at the top of the ITK source tree.
*
*=========================================================================*/
#ifndef __itkDeformableSimplexMesh3DKNNClassifierForceFilter_h
#define __itkDeformableSimplexMesh3DKNNClassifierForceFilter_h

//#include "itkDeformableSimplexMesh3DFilter.h"
#include "itkDeformableSimplexMesh3DFilterWithAdaptedConfidence.h"
#include "itkMesh.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkCovariantVector.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkListSample.h"
#include "itkKdTreeGenerator.h"
#include "itkWeightedCentroidKdTreeGenerator.h"

#include "kmKNNProfileClassifier-FLANN.h"

#include <set>
#include <map>

namespace itk
{
	/** \class DeformableSimplexMesh3DKNNClassifierForceFilter
	* \brief
	* Additional to its superclass this model adds an balloon force component to the
	* internal forces.
	*
	* The balloon force can be scaled, by setting the parameter kappa.
	*
	* \author Thomas Boettger. Division Medical and Biological Informatics, German Cancer Research Center, Heidelberg.
	*
	* \ingroup ITKDeformableMesh
	*/
	template< class TInputMesh, class TOutputMesh, class TKNNClassifier >
	class ITK_EXPORT DeformableSimplexMesh3DKNNClassifierForceFilter
		:public DeformableSimplexMesh3DFilterWithAdaptedConfidence< TInputMesh, TOutputMesh >
	{
	public:
		/** Standard "Self" typedef. */
		typedef DeformableSimplexMesh3DKNNClassifierForceFilter Self;

		/** Standard "Superclass" typedef. */
		typedef DeformableSimplexMesh3DFilterWithAdaptedConfidence< TInputMesh, TOutputMesh > Superclass;

		/** Smart pointer typedef support */
		typedef SmartPointer< Self >       Pointer;
		typedef SmartPointer< const Self > ConstPointer;

		/** Method of creation through the object factory. */
		itkNewMacro(Self);

		/** Run-time type information (and related methods). */
		itkTypeMacro(DeformableSimplexMesh3DKNNClassifierForceFilter, DeformableSimplexMesh3DFilterWithAdaptedConfidence);

		/** Some typedefs. */
		typedef TInputMesh                                  InputMeshType;
		typedef TOutputMesh                                 OutputMeshType;
		typedef typename Superclass::PointType              PointType;
		typedef typename Superclass::GradientIndexType      GradientIndexType;
		typedef typename Superclass::GradientIndexValueType GradientIndexValueType;
		typedef typename Superclass::GradientPixelType      GradientPixelType;

		typedef itk::Image<float, GradientImageType::ImageDimension> IntensityImageType;
		typedef typename IntensityImageType::Pointer                 IntensityImagePointer;
		typedef itk::LinearInterpolateImageFunction<IntensityImageType, double> IntensityInterpolatorType;
		typedef typename IntensityInterpolatorType::Pointer                     IntensityInterpolatorPointer;

		/* Mesh pointer definition. */
		typedef typename InputMeshType::Pointer  InputMeshPointer;
		typedef typename OutputMeshType::Pointer OutputMeshPointer;

		typedef typename Statistics::ListSample<PointType>                       SampleType;
		typedef typename Statistics::WeightedCentroidKdTreeGenerator<SampleType> TreeGeneratorType;
		typedef typename TreeGeneratorType::KdTreeType                           KdTreeType;
		typedef typename KdTreeType::InstanceIdentifierVectorType                NeighborhoodIdentifierType;

		typedef typename InputMeshType::PixelType PixelType;

		typedef SimplexMeshGeometry::CovariantVectorType                       NormalVectorType;

		typedef typename TKNNClassifier KNNProfileClassifierType;
		typedef km::Profile<PROFILE_DIM> ProfileType;

		itkStaticConstMacro(Dimension,   unsigned int, InputMeshType::PointDimension);

		itkSetMacro(Kappa, double);
		itkGetConstMacro(Kappa, double);

		virtual void Initialize();

		void SetItensityClassifier(  KNNProfileClassifierType * classifier )
		{ 
			m_IntensityClassifier = classifier; 
		}

		void SetGradientClassifier(  KNNProfileClassifierType * classifier )
		{ 
			m_GradientClassifier = classifier; 
		}

		void SetIntensity( IntensityImageType* intensityImage )
		{
			m_IntensityImage = intensityImage;
		}

		IntensityInterpolatorPointer GetIntensity()
		{
			return m_IntensityInterpolator;
		}

		void SetEdgeGradientImage( GradientImageType * edgeGradientImage )
		{
			m_EdgeGradientImage = edgeGradientImage;	
		}

		GradientImagePointer GetEdgeGradientImage()
		{
			return m_EdgeGradientImage;
		}

		void SetInnerMesh( InputMeshType* innerMesh/*, InputMeshType* outerMesh*/ )
		{
			m_InnerMesh = innerMesh;
			//m_OuterMesh = outerMesh;
		}

	protected:
		DeformableSimplexMesh3DKNNClassifierForceFilter();
		~DeformableSimplexMesh3DKNNClassifierForceFilter();
		DeformableSimplexMesh3DKNNClassifierForceFilter(const Self &)
		{}

		void operator=(const Self &)
		{}

		void PrintSelf(std::ostream & os, Indent indent) const;

		/**
		* Compute the external force component
		*/
		virtual void ComputeExternalForce(SimplexMeshGeometry *data, const GradientImageType *gradient, unsigned int idx, double & confidence);

		/** Parameters definitions. */

		/**
		* scalar for balloon force
		*/
		double m_Kappa;

		//GradientInterpolatorPointer m_GradientInterpolator;

		KNNProfileClassifierType * m_IntensityClassifier;
		KNNProfileClassifierType * m_GradientClassifier;

		IntensityImagePointer m_IntensityImage;
		IntensityInterpolatorPointer m_IntensityInterpolator;

		GradientImagePointer m_EdgeGradientImage;
		GradientInterpolatorPointer m_EdgeGradientInterpolator;

		InputMeshPointer m_InnerMesh;
		//InputMeshPointer m_OuterMesh;

		typename TreeGeneratorType::Pointer      m_InnerKdTreeGenerator;
		typename SampleType::Pointer             m_InnerSamplePoints;

		//typename TreeGeneratorType::Pointer      m_OuterKdTreeGenerator;
		//typename SampleType::Pointer             m_OuterSamplePoints;
	}; // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDeformableSimplexMesh3DKNNClassifierForceFilter.hxx"
#endif

#endif //__itkDeformableSimplexMesh3DKNNClassifierForceFilter_H
