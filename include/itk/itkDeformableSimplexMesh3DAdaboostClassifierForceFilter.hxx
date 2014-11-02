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
#ifndef __itkDeformableSimplexMesh3DAdaboostClassifierForceFilter_hxx
#define __itkDeformableSimplexMesh3DAdaboostClassifierForceFilter_hxx

#include "itkDeformableSimplexMesh3DAdaboostClassifierForceFilter.h"
#include "itkNumericTraits.h"
#include <itkNormalizedCorrelationPointSetToImageMetric.h>

#include <set>

namespace itk
{
	/* Constructor. */
	template< typename TInputMesh, typename TOutputMesh>
	DeformableSimplexMesh3DAdaboostClassifierForceFilter< TInputMesh, TOutputMesh>
		::DeformableSimplexMesh3DAdaboostClassifierForceFilter()

	{
		m_Kappa = 0.1;
		m_EdgeGradientImage = NULL;
		m_ProbabilityMap = NULL;

		m_BestMesh = NULL;
		m_BestSamplePoints = NULL;
		m_BestKdTreeGenerator = NULL;

		m_BoundaryMesh = NULL;
		m_BoundarySamplePoints = NULL;
		m_BoundaryKdTreeGenerator = NULL;

		m_ShapeMesh = NULL;
		m_ShapeSamplePoints = NULL;
		m_ShapeKdTreeGenerator = NULL;
	}

	template< typename TInputMesh, typename TOutputMesh>
	DeformableSimplexMesh3DAdaboostClassifierForceFilter< TInputMesh, TOutputMesh>
		::~DeformableSimplexMesh3DAdaboostClassifierForceFilter()
	{}

	/* PrintSelf. */
	template< typename TInputMesh, typename TOutputMesh>
	void
		DeformableSimplexMesh3DAdaboostClassifierForceFilter< TInputMesh, TOutputMesh>
		::PrintSelf(std::ostream & os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
		os << indent << "Kappa = " << m_Kappa << std::endl;
	} /* End PrintSelf. */


	template< typename TInputMesh, typename TOutputMesh>
	void
		DeformableSimplexMesh3DAdaboostClassifierForceFilter< TInputMesh, TOutputMesh>
		::Initialize()
	{
		Superclass::Initialize();

		const GradientImageType* gradimg = this->GetGradient();
		if (!gradimg)
		{
			std::cerr<<"Error! No gradient image!"<<std::endl;
			return;
		}

		if (!m_IntensityImage)
		{
			std::cerr<<"Error! No intensity image!"<<std::endl;
			return;
		}
		else
		{
			m_IntensityInterpolator = IntensityInterpolatorType::New();
			m_IntensityInterpolator->SetInputImage( m_IntensityImage );
		}

		if (m_ProbabilityMap)
		{
			m_ProbabilityMap->FillBuffer( -1 );
			m_ProbabilityMapInterpolator = FloatInterpolatorType::New();
			m_ProbabilityMapInterpolator->SetInputImage( m_ProbabilityMap );
		}
		else
		{
			m_ProbabilityMap = FloatImageType::New();
			m_ProbabilityMap->SetRegions( m_IntensityImage->GetLargestPossibleRegion() );
			m_ProbabilityMap->SetOrigin( m_IntensityImage->GetOrigin() );
			m_ProbabilityMap->SetSpacing( m_IntensityImage->GetSpacing() );
			m_ProbabilityMap->SetDirection( m_IntensityImage->GetDirection() );
			m_ProbabilityMap->Allocate();
			m_ProbabilityMap->FillBuffer( -1 );
			m_ProbabilityMapInterpolator = FloatInterpolatorType::New();
			m_ProbabilityMapInterpolator->SetInputImage( m_ProbabilityMap );
		}

		if (m_EdgeGradientImage)
		{
			m_EdgeGradientInterpolator = GradientInterpolatorType::New();
			m_EdgeGradientInterpolator->SetInputImage( m_EdgeGradientImage );

			int imageWidth  = m_EdgeGradientImage->GetBufferedRegion().GetSize()[0];
			int imageHeight = m_EdgeGradientImage->GetBufferedRegion().GetSize()[1];
			int imageDepth  = m_EdgeGradientImage->GetBufferedRegion().GetSize()[2];

			//std::cout<<imageWidth<<", "<<imageHeight<<", "<<imageDepth<<std::endl;
		}
		else
		{
			std::cout<<"No edge gradient image!"<<std::endl;
		}

		if (m_BestMesh)
		{
			m_BestSamplePoints = SampleType::New();
			m_BestSamplePoints->SetMeasurementVectorSize( 3 );

			typename InputMeshType::PointsContainerConstIterator It = m_BestMesh->GetPoints()->Begin();
			while( It != m_BestMesh->GetPoints()->End() )
			{
				PointType point = It.Value();
				m_BestSamplePoints->PushBack( point );
				++It;
			}

			m_BestKdTreeGenerator = TreeGeneratorType::New();
			m_BestKdTreeGenerator->SetSample( m_BestSamplePoints );
			m_BestKdTreeGenerator->SetBucketSize( 4 );
			m_BestKdTreeGenerator->Update();
		}

		if (m_BoundaryMesh)
		{
			m_BoundarySamplePoints = SampleType::New();
			m_BoundarySamplePoints->SetMeasurementVectorSize( 3 );

			typename InputMeshType::PointsContainerConstIterator It = m_BoundaryMesh->GetPoints()->Begin();
			while( It != m_BoundaryMesh->GetPoints()->End() )
			{
				PointType point = It.Value();
				m_BoundarySamplePoints->PushBack( point );
				++It;
			}

			m_BoundaryKdTreeGenerator = TreeGeneratorType::New();
			m_BoundaryKdTreeGenerator->SetSample( m_BoundarySamplePoints );
			m_BoundaryKdTreeGenerator->SetBucketSize( 4 );
			m_BoundaryKdTreeGenerator->Update();
		}

		if (m_ShapeMesh)
		{
			m_ShapeSamplePoints = SampleType::New();
			m_ShapeSamplePoints->SetMeasurementVectorSize( 3 );

			typename InputMeshType::PointsContainerConstIterator It = m_ShapeMesh->GetPoints()->Begin();
			while( It != m_ShapeMesh->GetPoints()->End() )
			{
				PointType point = It.Value();
				m_ShapeSamplePoints->PushBack( point );
				++It;
			}

			m_ShapeKdTreeGenerator = TreeGeneratorType::New();
			m_ShapeKdTreeGenerator->SetSample( m_ShapeSamplePoints );
			m_ShapeKdTreeGenerator->SetBucketSize( 4 );
			m_ShapeKdTreeGenerator->Update();
		}
	}

	template< typename TInputMesh, typename TOutputMesh>
	void
		DeformableSimplexMesh3DAdaboostClassifierForceFilter< TInputMesh, TOutputMesh>
		::Intervene( )
	{
		//std::cout<<m_Step<<std::endl;
	}

	template< typename TInputMesh, typename TOutputMesh>
	void
		DeformableSimplexMesh3DAdaboostClassifierForceFilter< TInputMesh, TOutputMesh>
		::ComputeExternalForce(SimplexMeshGeometry *data,const GradientImageType *gradientImage, unsigned int idx, double & confidence)
	{
		VectorType         vec_for;
		VectorType         vec_for_gradient;
		VectorType         vec_normal;
		VectorType         vec_normal_no_cross;
		VectorType         vec_tmp;

		vec_normal.SetVnlVector( data->normal.GetVnlVector() );
		vec_for_gradient.Fill( 0 );

		PointType pos_cur;
		pos_cur.CastFrom( data->pos );

		double step_sum = 0.0;
		double step_shape = 0.0;
		double step_best = 0.0;

		double shape_force_factor = 0.6;
		if (m_ShapeMesh )
		{
			typename TreeGeneratorType::KdTreeType::InstanceIdentifierVectorType neighbors;
			this->m_ShapeKdTreeGenerator->GetOutput()->Search( pos_cur, 1u, neighbors );
			PointType closestShapePt = this->m_ShapeKdTreeGenerator->GetOutput()->GetMeasurementVector( neighbors[0] );

			double variance = 0.0;
			this->m_VarianceMap->GetPointData(neighbors[0], &variance);

			shape_force_factor *= (1.0-variance);

			SimplexMeshGeometry* data_shape = this->m_ShapeMesh->GetGeometryData()->GetElement(neighbors[0]);
			double phi_shape = data_shape->phi;

			double phi_cur = data->phi;

			step_shape = (180.0/(20.0*3.14))*(phi_cur - phi_shape); //20 degree as a threshold

			if (step_shape > 1.0)
			{
				step_shape = 1.0;
			}
			else if (step_shape < -1.0)
			{
				step_shape = -1.0;
			}
		}

		if (m_BestMesh )
		{
			PointType closestBestPt;
			if (m_P2p)
			{
				closestBestPt = m_BestMesh->GetPoint( idx );
			}
			else
			{
				typename TreeGeneratorType::KdTreeType::InstanceIdentifierVectorType neighbors;
				this->m_BestKdTreeGenerator->GetOutput()->Search( pos_cur, 1u, neighbors );
				closestBestPt = this->m_BestKdTreeGenerator->GetOutput()->GetMeasurementVector( neighbors[0] );
			}

			double dist = closestBestPt.EuclideanDistanceTo( pos_cur );
			vec_tmp = closestBestPt - pos_cur;

			if (vec_tmp.GetNorm() > 0.0)
			{
				vec_tmp.Normalize();
			}

			step_best = dot_product( data->normal.GetVnlVector(), vec_tmp.GetVnlVector() );
		}

		step_sum = shape_force_factor*step_shape + (1.0-shape_force_factor)*step_best;

		vec_for = vec_normal * step_sum;

		data->externalForce = vec_for * this->GetKappa() + vec_for_gradient * this->GetBeta();
	}

}/* end namespace itk. */

#endif //__itkDeformableSimplexMesh3DAdaboostClassifierForceFilter_hxx
