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
#ifndef __itkDeformableSimplexMesh3DAdaboostClassifierForceFilter_h
#define __itkDeformableSimplexMesh3DAdaboostClassifierForceFilter_h

//#include "itkDeformableSimplexMesh3DFilter.h"
#include "itkDeformableSimplexMesh3DFilterWithAdaptedConfidence.h"
#include "itkMesh.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkCovariantVector.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkListSample.h"
#include "itkKdTreeGenerator.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkDifferenceOfGaussiansGradientImageFilter.h"

#include <set>
#include <map>

namespace itk
{
	/** \class DeformableSimplexMesh3DAdaboostClassifierForceFilter
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
	template< class TInputMesh, class TOutputMesh>
	class ITK_EXPORT DeformableSimplexMesh3DAdaboostClassifierForceFilter
		:public DeformableSimplexMesh3DFilterWithAdaptedConfidence< TInputMesh, TOutputMesh >
	{
	public:
		/** Standard "Self" typedef. */
		typedef DeformableSimplexMesh3DAdaboostClassifierForceFilter Self;

		/** Standard "Superclass" typedef. */
		typedef DeformableSimplexMesh3DFilterWithAdaptedConfidence< TInputMesh, TOutputMesh > Superclass;

		/** Smart pointer typedef support */
		typedef SmartPointer< Self >       Pointer;
		typedef SmartPointer< const Self > ConstPointer;

		/** Method of creation through the object factory. */
		itkNewMacro(Self);

		/** Run-time type information (and related methods). */
		itkTypeMacro(DeformableSimplexMesh3DAdaboostClassifierForceFilter, DeformableSimplexMesh3DFilterWithAdaptedConfidence);

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

		typedef itk::Image<float, InputMeshType::PointDimension> FloatImageType;
		typedef typename FloatImageType::Pointer FloatImagePointer;

		typedef itk::LinearInterpolateImageFunction<FloatImageType, double>     FloatInterpolatorType;
		typedef typename FloatInterpolatorType::Pointer                         FloatInterpolatorPointer;

		itkStaticConstMacro(Dimension,   unsigned int, InputMeshType::PointDimension);

		itkSetMacro(Kappa, double);
		itkGetConstMacro(Kappa, double);

		virtual void Initialize();

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

		void SetProbabilityMap(FloatImageType* p)
		{
			m_ProbabilityMap = p;
		}

		FloatImageType* GetProbabilityMap()
		{
			return m_ProbabilityMap;
		}

		void SetBestMesh( InputMeshType* bestMesh )
		{
			m_BestMesh = bestMesh;
		}

		void SetBoundaryMesh( InputMeshType* boundaryMesh )
		{
			m_BoundaryMesh = boundaryMesh;
		}

		void SetShapeMesh( InputMeshType* shapeMesh )
		{
			m_ShapeMesh = shapeMesh;
		}

		void SetP2p(bool flag)
		{
			m_P2p = flag;
		}

	protected:
		DeformableSimplexMesh3DAdaboostClassifierForceFilter();
		~DeformableSimplexMesh3DAdaboostClassifierForceFilter();
		DeformableSimplexMesh3DAdaboostClassifierForceFilter(const Self &)
		{}

		void operator=(const Self &)
		{}

		void PrintSelf(std::ostream & os, Indent indent) const;

		/**
		* Compute the external force component
		*/
		virtual void ComputeExternalForce(SimplexMeshGeometry *data, const GradientImageType *gradient, unsigned int idx, double & confidence);

		virtual void Intervene( );

		/** Parameters definitions. */

		/**
		* scalar for balloon force
		*/
		double m_Kappa;

		//GradientInterpolatorPointer m_GradientInterpolator;

		IntensityImagePointer m_IntensityImage;
		IntensityInterpolatorPointer m_IntensityInterpolator;

		GradientImagePointer m_EdgeGradientImage;
		GradientInterpolatorPointer m_EdgeGradientInterpolator;

		FloatImagePointer m_ProbabilityMap;
		FloatInterpolatorPointer m_ProbabilityMapInterpolator;

		InputMeshPointer m_BestMesh;
		typename SampleType::Pointer m_BestSamplePoints;
		typename TreeGeneratorType::Pointer m_BestKdTreeGenerator;

		InputMeshPointer m_BoundaryMesh;
		typename SampleType::Pointer m_BoundarySamplePoints;
		typename TreeGeneratorType::Pointer m_BoundaryKdTreeGenerator;

		InputMeshPointer m_ShapeMesh;
		typename SampleType::Pointer m_ShapeSamplePoints;
		typename TreeGeneratorType::Pointer m_ShapeKdTreeGenerator;

		bool m_P2p;
	}; // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDeformableSimplexMesh3DAdaboostClassifierForceFilter.hxx"
#endif

#endif //__itkDeformableSimplexMesh3DAdaboostClassifierForceFilter_H
