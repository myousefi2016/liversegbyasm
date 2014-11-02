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
#ifndef __itkDeformableSimplexMesh3DFilterWithAdaptedConfidence_hxx
#define __itkDeformableSimplexMesh3DFilterWithAdaptedConfidence_hxx

#include "itkDeformableSimplexMesh3DFilterWithAdaptedConfidence.h"
#include "itkNumericTraits.h"

#include <set>
#include <math.h>

#include "vxl_version.h"
#if VXL_VERSION_DATE_FULL > 20040406
#include "vnl/vnl_cross.h"
#define itk_cross_3d vnl_cross_3d
#else
#define itk_cross_3d cross_3d
#endif

namespace itk
{
	/* Constructor. */
	template< typename TInputMesh, typename TOutputMesh >
	DeformableSimplexMesh3DFilterWithAdaptedConfidence< TInputMesh, TOutputMesh >
		::DeformableSimplexMesh3DFilterWithAdaptedConfidence()
	{
		error_count = 0;
		m_ConfidenceThreshold = 0.8;

		m_GradientInterpolator = GradientInterpolatorType::New();
		m_VarianceMap = NULL;
	}

	template< typename TInputMesh, typename TOutputMesh >
	DeformableSimplexMesh3DFilterWithAdaptedConfidence< TInputMesh, TOutputMesh >
		::~DeformableSimplexMesh3DFilterWithAdaptedConfidence()
	{}

	/* PrintSelf. */
	template< typename TInputMesh, typename TOutputMesh >
	void
		DeformableSimplexMesh3DFilterWithAdaptedConfidence< TInputMesh, TOutputMesh >
		::PrintSelf(std::ostream & os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
		os << indent << "Alpha = " << this->GetAlpha() << std::endl;
		os << indent << "Beta = " << this->GetBeta() << std::endl;
		os << indent << "Gamma = " << this->GetGamma() << std::endl;
		os << indent << "Rigidity = " << this->GetRigidity() << std::endl;
		os << indent << "ConfidenceThreshold = " << this->GetConfidenceThreshold() << std::endl;
		os << indent << "Iterations = " << this->GetIterations() << std::endl;
		os << indent << "Step = " << this->GetStep() << std::endl;
		os << indent << "ImageDepth = " << this->GetImageDepth() << std::endl;

		const GradientImageType * gradientImage = this->GetGradient();

		if ( gradientImage )
		{
			os << indent << "Gradient = " << gradientImage << std::endl;
		}
		else
		{
			os << indent << "Gradient = " << "(None)" << std::endl;
		}
		os << indent << "ImageHeight = " << this->GetImageHeight() << std::endl;
		os << indent << "ImageWidth = " << this->GetImageWidth() << std::endl;
		os << indent << "Damping = " << this->GetDamping() << std::endl;
		if ( this->m_Data.IsNotNull() )
		{
			os << indent << "Data = " << this->GetData() << std::endl;
		}
		else
		{
			os << indent << "Data = " << "(None)" << std::endl;
		}
	} /* End PrintSelf. */

	template< typename TInputMesh, typename TOutputMesh >
	void
		DeformableSimplexMesh3DFilterWithAdaptedConfidence< TInputMesh, TOutputMesh >
		::Initialize()
	{
		//Superclass::Initialize();

		std::cout<<"DeformableSimplexMesh3DFilterWithAdaptedConfidence::Initialize()..."<<std::endl;

		const InputMeshType *             inputMesh = this->GetInput(0);
		const InputPointsContainer *      points    = inputMesh->GetPoints();
		InputPointsContainerConstIterator pointItr  = points->Begin();

		const GradientImageType * gradientImage = this->GetGradient();

		if ( gradientImage )
		{
			GradientImageSizeType imageSize = gradientImage->GetBufferedRegion().GetSize();

			m_GradientInterpolator->SetInputImage( gradientImage );
		}

		this->m_Data = inputMesh->GetGeometryData();
		if ( this->m_Data.IsNull() || this->m_Data->Size()==0)
		{
			std::cout<<"Geometry data is empty!!!"<<std::endl;
			//this->ComputeGeometry();
		}

		while ( pointItr != points->End() )
		{
			SimplexMeshGeometry *data;
			IdentifierType       idx = pointItr.Index();

			data = this->m_Data->GetElement(idx);
			data->pos = pointItr.Value();

			//        InputMeshType::ArrayType neighbors =
			// this->GetInput(0)->GetNeighbors( pointItr.Index() );

			data->neighbors[0] = points->GetElement(data->neighborIndices[0]);
			data->neighbors[1] = points->GetElement(data->neighborIndices[1]);
			data->neighbors[2] = points->GetElement(data->neighborIndices[2]);

			// store neighborset with a specific radius

			unsigned int adaptedRigidity = m_Rigidity;
			
			if (m_VarianceMap)
			{
				double variance = 0.0;
				m_VarianceMap->GetPointData(idx, &variance);
				adaptedRigidity = static_cast<unsigned int>(3.0*(1.0-variance));
			}

			InputNeighbors *       neighborsList = inputMesh->GetNeighbors(pointItr.Index(), adaptedRigidity);
			InputNeighborsIterator neighborIt = neighborsList->begin();

			NeighborSetType *neighborSet = new NeighborSetType();
			while ( neighborIt != neighborsList->end() )
			{
				neighborSet->insert(*neighborIt++);
			}
			// garbage collection (from itkSimplexMesh)
			delete neighborsList;
			data->neighborSet =  neighborSet;

			pointItr++;
		}

		InputMeshType* nonconstInputMesh = const_cast< InputMeshType * >( inputMesh );

		InputPointDataContainer* pointdata = nonconstInputMesh->GetPointData();
		if (pointdata->Size() <= 0)
		{
			pointdata->Reserve( nonconstInputMesh->GetNumberOfPoints() );
			for (int i=0;i<nonconstInputMesh->GetNumberOfPoints();i++)
			{
				nonconstInputMesh->SetPointData( i, 0.0 );	
			}
		}

		std::cout<<"DeformableSimplexMesh3DFilterWithAdaptedConfidence::Initialize() end..."<<std::endl;
	}

	template< typename TInputMesh, typename TOutputMesh >
	void
		DeformableSimplexMesh3DFilterWithAdaptedConfidence< TInputMesh, TOutputMesh >
		::GenerateData()
	{
		this->Initialize();

		m_Step = 0;

		while ( m_Step < m_Iterations )
		{
			const float progress = static_cast< float >( m_Step )
				/ static_cast< float >( m_Iterations );

			this->UpdateProgress(progress);

			this->ComputeGeometry();

			if ( m_Step % 10 == 0 && m_Step > 0 )
			{
				this->UpdateReferenceMetrics();
			}

			this->ComputeDisplacement();

			this->Intervene();

			m_Step++;
		}
	}


	template< typename TInputMesh, typename TOutputMesh >
	void
		DeformableSimplexMesh3DFilterWithAdaptedConfidence< TInputMesh, TOutputMesh >
		::ComputeDisplacement()
	{
		const InputMeshType *inputMesh = this->GetInput(0);

		const GradientImageType *gradientImage = this->GetGradient();

		// Filters should not modify their input...
		// There is a design flaw here.
		InputPointsContainer *nonConstPoints =
			const_cast< InputPointsContainer * >( inputMesh->GetPoints() );

		InputPointDataContainer *pointdata = const_cast< InputPointDataContainer * >( inputMesh->GetPointData() );

		typename GeometryMapType::Iterator dataIt = this->m_Data->Begin();
		SimplexMeshGeometry *data;
		unsigned int idx;
		VectorType           displacement;

		while ( dataIt != this->m_Data->End() )
		{
			idx = dataIt.Index();
			data = dataIt.Value();

			this->ComputeInternalForce(data);

			//this->ComputeExternalForce(data,gradientImage);

			double confidence = pointdata->GetElement(idx);

			if (confidence>= this->GetConfidenceThreshold() )
			{
				data->externalForce.Fill(0);
			}
			else
			{
				this->ComputeExternalForce(data,gradientImage, idx, confidence);
			}

			displacement.SetVnlVector( m_Alpha * ( data->internalForce ).GetVnlVector()
				+ ( data->externalForce ).GetVnlVector() );

			data->pos += displacement;

			nonConstPoints->InsertElement(idx, data->pos);

			dataIt++;
		}

		//std::cout<<"Error count: "<<error_count<<std::endl;
		error_count = 0;
	}

	template< typename TInputMesh, typename TOutputMesh >
	void
		DeformableSimplexMesh3DFilterWithAdaptedConfidence< TInputMesh, TOutputMesh >
		::ComputeExternalForce(SimplexMeshGeometry *data,const GradientImageType *gradientImage, unsigned int idx, double & confidence)
	{
		//this->ComputeExternalForce( data, gradientImage );
	}

	template< typename TInputMesh, typename TOutputMesh >
	void
		DeformableSimplexMesh3DFilterWithAdaptedConfidence< TInputMesh, TOutputMesh >
		::Intervene()
	{
		
	}

} /* end namespace itk. */

#endif //__itkDeformableSimplexMesh3DFilterWithAdaptedConfidence_TXX
