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
#ifndef __itkDeformableSimplexMesh3DKNNClassifierForceFilter_hxx
#define __itkDeformableSimplexMesh3DKNNClassifierForceFilter_hxx

#include "itkDeformableSimplexMesh3DKNNClassifierForceFilter.h"
#include "itkNumericTraits.h"
#include <itkNormalizedCorrelationPointSetToImageMetric.h>

#include <set>

namespace itk
{
	/* Constructor. */
	template< typename TInputMesh, typename TOutputMesh, typename TKNNClassifier >
	DeformableSimplexMesh3DKNNClassifierForceFilter< TInputMesh, TOutputMesh, TKNNClassifier >
		::DeformableSimplexMesh3DKNNClassifierForceFilter()

	{
		m_Kappa = 0.1;
		m_EdgeGradientImage = NULL;
		
		m_InnerMesh = NULL;
		//m_OuterMesh = NULL;

		this->m_InnerKdTreeGenerator = NULL;
		this->m_InnerSamplePoints = NULL;

		//this->m_OuterKdTreeGenerator = NULL;
		//this->m_OuterSamplePoints = NULL;
	}

	template< typename TInputMesh, typename TOutputMesh, typename TKNNClassifier >
	DeformableSimplexMesh3DKNNClassifierForceFilter< TInputMesh, TOutputMesh, TKNNClassifier >
		::~DeformableSimplexMesh3DKNNClassifierForceFilter()
	{}

	/* PrintSelf. */
	template< typename TInputMesh, typename TOutputMesh, typename TKNNClassifier >
	void
		DeformableSimplexMesh3DKNNClassifierForceFilter< TInputMesh, TOutputMesh, TKNNClassifier >
		::PrintSelf(std::ostream & os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
		os << indent << "Kappa = " << m_Kappa << std::endl;
	} /* End PrintSelf. */


	template< typename TInputMesh, typename TOutputMesh, typename TKNNClassifier >
	void
		DeformableSimplexMesh3DKNNClassifierForceFilter< TInputMesh, TOutputMesh, TKNNClassifier >
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

		if (m_InnerMesh)
		{
			//Inner
			this->m_InnerSamplePoints = SampleType::New();
			this->m_InnerSamplePoints->SetMeasurementVectorSize( 3 );

			typename InputMeshType::PointsContainerConstIterator It = this->m_InnerMesh->GetPoints()->Begin();
			while( It != this->m_InnerMesh->GetPoints()->End() )
			{
				PointType point = It.Value();
				this->m_InnerSamplePoints->PushBack( point );
				++It;
			}

			this->m_InnerKdTreeGenerator = TreeGeneratorType::New();
			this->m_InnerKdTreeGenerator->SetSample( this->m_InnerSamplePoints );
			this->m_InnerKdTreeGenerator->SetBucketSize( 4 );
			this->m_InnerKdTreeGenerator->Update();

			////Outer
			//this->m_OuterSamplePoints = SampleType::New();
			//this->m_OuterSamplePoints->SetMeasurementVectorSize( 3 );

			//It = this->m_OuterMesh->GetPoints()->Begin();
			//while( It != this->m_OuterMesh->GetPoints()->End() )
			//{
			//	PointType point = It.Value();
			//	this->m_OuterSamplePoints->PushBack( point );
			//	++It;
			//}

			//this->m_OuterKdTreeGenerator = TreeGeneratorType::New();
			//this->m_OuterKdTreeGenerator->SetSample( this->m_OuterSamplePoints );
			//this->m_OuterKdTreeGenerator->SetBucketSize( 4 );
			//this->m_OuterKdTreeGenerator->Update();
		}
	}


	template< typename TInputMesh, typename TOutputMesh, typename TKNNClassifier >
	void
		DeformableSimplexMesh3DKNNClassifierForceFilter< TInputMesh, TOutputMesh, TKNNClassifier >
		::ComputeExternalForce(SimplexMeshGeometry *data,const GradientImageType *gradientImage, unsigned int idx, double & confidence)
	{
		VectorType         vec_for;
		VectorType         vec_for_gradient;
		VectorType         vec_normal;
		VectorType         vec_normal_no_cross;

		vec_normal.SetVnlVector( data->normal.GetVnlVector() );
		vec_for_gradient.Fill( 0 );

		MeshPointType pos_cur = data->pos;
		MeshPointType pos_si = data->pos - vec_normal*3.0;  //shift inside
		MeshPointType pos_so = data->pos + vec_normal*3.0;  //shift outside

		ProbabilityContainer probablities_cur;
		ProbabilityContainer probablities_si;
		ProbabilityContainer probablities_so;

		DistanceContainer dists_cur;
		DistanceContainer dists_si;
		DistanceContainer dists_so;

		ProfileType profile;
		profile.spacing = PROFILE_SPACING;
		profile.index = idx;
		BClassType btype_cur;
		BClassType btype_si;
		BClassType btype_so;

		double step_sum = 0.0;

		PointMatrixType pmatrix = PointMatrixType(new float[profile.Dimension], 1, profile.Dimension);

		const GradientInterpolatorType* gradientInterpolator = this->GetGradientIntepolator();

		km::extractAdvancedProfile<GradientInterpolatorType, IntensityInterpolatorType, MeshPointType, NormalVectorType, ProfileType>(
			gradientInterpolator,
			m_IntensityInterpolator,
			pos_cur,
			data->normal,
			profile,
			m_IntensityClassifier->profile_category);
		profile2matrix( profile, pmatrix );
		m_IntensityClassifier->classify( pmatrix,idx,btype_cur,probablities_cur,dists_cur );

		step_sum = probablities_cur[IPClass] - probablities_cur[OPClass];

		if ( step_sum > 0.01 )
		{
			km::extractAdvancedProfile<GradientInterpolatorType, IntensityInterpolatorType, MeshPointType, NormalVectorType, ProfileType>(
				gradientInterpolator,
				m_IntensityInterpolator,
				pos_cur,
				data->normal,
				profile,
				m_GradientClassifier->profile_category);
			profile2matrix( profile, pmatrix );
			m_GradientClassifier->classify( pmatrix,idx,btype_cur,probablities_cur,dists_cur );

			if (btype_cur != IPClass)
			{
				step_sum = probablities_cur[IPClass] - probablities_cur[OPClass];
			}
		}

		InputPointDataContainer *pointdata = const_cast< InputPointDataContainer * >( this->GetInput(0)->GetPointData() );
		pointdata->SetElement( idx, step_sum );

		PointType origpt;
		origpt.CastFrom( data->pos );
		if (m_InnerMesh)
		{
			PointType closestInnerPt = this->m_InnerMesh->GetPoint( idx );
			vec_normal_no_cross = origpt - closestInnerPt;
			if (vec_normal_no_cross.GetSquaredNorm() > 1)
			{
				vec_normal_no_cross.Normalize();
			}
		}
		else
		{
			vec_normal_no_cross = vec_normal;
		}

		vec_for = vec_normal_no_cross * step_sum;

		//Compute gradient force
		if (m_EdgeGradientImage  && step_sum >= 0.01)
		{
			if (m_EdgeGradientInterpolator->IsInsideBuffer( data->pos ))
			{
				GradientType gradient = m_EdgeGradientInterpolator->Evaluate( data->pos );

				double mag = dot_product( vec_normal_no_cross.GetVnlVector(), gradient.GetVnlVector() );

				if ( mag > 0.5 )
				{
					mag = 0.5;
				}
				else if (mag < -0.5)
				{
					mag = -0.5;
				}
				vec_for_gradient = vec_normal_no_cross * mag;
			}
			else
			{
				vec_for_gradient = vec_normal_no_cross * -1.0;
				vec_for.Fill( 0 );
			}
		}
		else
		{
			vec_for_gradient.Fill( 0.0 );
		}

		data->externalForce = vec_for * this->GetKappa() + vec_for_gradient * this->GetBeta();
	}

}/* end namespace itk. */

#endif //__itkDeformableSimplexMesh3DKNNClassifierForceFilter_hxx
