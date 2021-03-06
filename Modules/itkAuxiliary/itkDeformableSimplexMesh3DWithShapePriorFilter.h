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
#ifndef __itkDeformableSimplexMesh3DWithShapePriorFilter_h
#define __itkDeformableSimplexMesh3DWithShapePriorFilter_h

#include "itkDeformableSimplexMesh3DFilter.h"
#include "itkMesh.h"
#include "itkImage.h"
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

#include "itkVector.h"
#include "itkCompositeTransform.h"
#include "itkAffineTransform.h"
#include "statismo_ITK/itkStatisticalShapeModelTransform.h"
#include "kmCommon.h"

using namespace km;

namespace itk
{
	/** \class DeformableSimplexMesh3DWithShapePriorFilter
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

	enum Phase
	{
		Lv1,
		Lv2,
		Lv3
	};

	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	class ITK_EXPORT DeformableSimplexMesh3DWithShapePriorFilter
		:public DeformableSimplexMesh3DFilter< TInputMesh, TOutputMesh >
	{
	public:
		/** Standard "Self" typedef. */
		typedef DeformableSimplexMesh3DWithShapePriorFilter Self;

		/** Standard "Superclass" typedef. */
		typedef DeformableSimplexMesh3DFilter< TInputMesh, TOutputMesh > Superclass;

		/** Smart pointer typedef support */
		typedef SmartPointer< Self >       Pointer;
		typedef SmartPointer< const Self > ConstPointer;

		/** Method of creation through the object factory. */
		itkNewMacro(Self);

		/** Run-time type information (and related methods). */
		itkTypeMacro(DeformableSimplexMesh3DWithShapePriorFilter, DeformableSimplexMesh3DFilter);

		/** Some typedefs. */
		typedef TInputMesh                                  InputMeshType;
		typedef TOutputMesh                                 OutputMeshType;
		typedef typename Superclass::PointType              PointType;
		typedef typename Superclass::GradientIndexType      GradientIndexType;
		typedef typename Superclass::GradientIndexValueType GradientIndexValueType;
		typedef typename Superclass::GradientPixelType      GradientPixelType;

		/* Mesh pointer definition. */
		typedef typename InputMeshType::Pointer  InputMeshPointer;
		typedef typename OutputMeshType::Pointer OutputMeshPointer;

		typedef typename InputMeshType::PointDataContainer             InputPointDataContainer;
		typedef typename InputMeshType::PointDataContainerIterator     InputPointDataContainerIterator;

		typedef typename Statistics::ListSample<PointType>                       SampleType;
		typedef typename Statistics::WeightedCentroidKdTreeGenerator<SampleType> TreeGeneratorType;
		typedef typename TreeGeneratorType::KdTreeType                           KdTreeType;
		typedef typename KdTreeType::InstanceIdentifierVectorType                NeighborhoodIdentifierType;

		typedef typename InputMeshType::PixelType PixelType;
		typedef SimplexMeshGeometry::CovariantVectorType NormalVectorType;

		itkStaticConstMacro(Dimension, unsigned int, InputMeshType::PointDimension);

		typedef typename TInputImage                  InputImageType;
		typedef typename InputImageType::Pointer      InputImagePointer;
		typedef typename InputImageType::ConstPointer InputImageConstPointer;

		typedef typename TStatisticalModel                  StatisticalModelType;
		typedef typename TRigidTransform                    RigidTransformType;
		typedef typename RigidTransformType::ParametersType RigidTransformParametersType;
		typedef typename RigidTransformType::Pointer        RigidTransformPointer;
		typedef typename RigidTransformType::ConstPointer   RigidTransformConstPointer;
		typedef typename TShapeTransform                    ShapeTransformType;
		typedef typename ShapeTransformType::ParametersType ShapeTransformParametersType;
		typedef typename ShapeTransformType::Pointer        ShapeTransformPointer;
		typedef typename ShapeTransformType::ConstPointer   ShapeTransformConstPointer;
		typedef CompositeTransform<double, Dimension>       CompositeTransformType;
		typedef typename CompositeTransformType::Pointer    CompositeTransformPointer;

		typedef km::ProfileExtractor<InputImageType, InputMeshType> ProfileExtractorType;
		typedef km::SSMUtils<InputMeshType, 
							 StatisticalModelType, 
							 RigidTransformType, 
							 ShapeTransformType>     SSMUtilsType;
		typedef km::ProfileClassifier                ProfileClassifierType;
		typedef km::ClassifierUtils<ProfileExtractorType> ClassifierUtilsType;

		itkSetMacro(Kappa, double);
		itkGetConstMacro(Kappa, double);

		itkSetConstObjectMacro( InputImage, InputImageType );
		itkGetConstObjectMacro( InputImage, InputImageType );

		itkSetObjectMacro( RigidTransform, RigidTransformType );
		itkGetObjectMacro( RigidTransform, RigidTransformType );

		itkSetObjectMacro( ShapeTransform, ShapeTransformType );
		itkGetObjectMacro( ShapeTransform, ShapeTransformType );

		void SetStatisticalModel(const StatisticalModelType* model)
		{
			this->m_Model = model;
		}
		const StatisticalModelType* GetStatisticalModel()
		{
			return this->m_Model;
		}

		void SetBoundaryClassifier(const ProfileClassifierType* c)
		{
			this->m_BoundaryClassifier = const_cast<ProfileClassifierType*>(c);
		}
		ProfileClassifierType* GetBoundaryClassifier()
		{
			return this->m_BoundaryClassifier;
		}

		void SetLiverClassifier(const ProfileClassifierType* c)
		{
			this->m_LiverClassifier = const_cast<ProfileClassifierType*>(c);
		}
		ProfileClassifierType* GetLiverClassifier()
		{
			return this->m_LiverClassifier;
		}

		InputMeshPointer GetShapeMesh()
		{
			return this->m_UpdatedShapeMesh;
		}

	protected:
		DeformableSimplexMesh3DWithShapePriorFilter();
		~DeformableSimplexMesh3DWithShapePriorFilter();
		DeformableSimplexMesh3DWithShapePriorFilter(const Self &)
		{}

		void operator=(const Self &)
		{}

		void PrintSelf(std::ostream & os, Indent indent) const;
		
		virtual void Initialize();

		virtual void ComputeDisplacement();
		
		virtual void GenerateData();
		
		virtual void IntervenePre();

		virtual void IntervenePost();

		virtual void UpdateShape();

		virtual void UpdateLandmarkStatus();

		/**
		* Compute the external force component
		*/
		virtual void ComputeExternalForce(SimplexMeshGeometry *data, unsigned int indx);

		/**
		* scalar for balloon force
		*/
		double m_Kappa;

		InputPointDataContainer* m_Pointdata;
		
		const StatisticalModelType*  m_Model;
		InputImageConstPointer       m_InputImage;
		RigidTransformPointer        m_RigidTransform;
		ShapeTransformPointer        m_ShapeTransform;
		ProfileExtractorType         m_ProfileExtractor;
		SSMUtilsType                 m_SSMUtils;
		ProfileClassifierType*       m_BoundaryClassifier;
		ProfileClassifierType*       m_LiverClassifier;
		CompositeTransformPointer    m_CompositeTransform;

		ClassifierUtilsType          m_ClassifierUtils;

		InputMeshPointer             m_ReferenceShapeMesh;
		InputMeshPointer             m_UpdatedShapeMesh;

		//cluster related
		int                          m_NumberOfShapeClusters;

		Phase                        m_Phase;
		std::map<int, double>        m_Forces;

		double                       m_MinShapeDifference;
		unsigned                     m_IterationsLv1;
		
	}; // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDeformableSimplexMesh3DWithShapePriorFilter.hxx"
#endif

#endif //__itkDeformableSimplexMesh3DWithShapePriorFilter_H
