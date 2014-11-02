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
#ifndef __itkDeformableSimplexMesh3DThresholdForceFilter_hxx
#define __itkDeformableSimplexMesh3DThresholdForceFilter_hxx

#include "itkDeformableSimplexMesh3DThresholdForceFilter.h"
#include "itkNumericTraits.h"

#include <set>

namespace itk
{
	/* Constructor. */
	template< typename TInputMesh, typename TOutputMesh >
	DeformableSimplexMesh3DThresholdForceFilter< TInputMesh, TOutputMesh >
		::DeformableSimplexMesh3DThresholdForceFilter()

	{
		m_Kappa = 0.1;
		m_FeatureImage = NULL;
		m_Threshold = 1.0;
	}

	template< typename TInputMesh, typename TOutputMesh >
	DeformableSimplexMesh3DThresholdForceFilter< TInputMesh, TOutputMesh >
		::~DeformableSimplexMesh3DThresholdForceFilter()
	{}

	/* PrintSelf. */
	template< typename TInputMesh, typename TOutputMesh >
	void
		DeformableSimplexMesh3DThresholdForceFilter< TInputMesh, TOutputMesh >
		::PrintSelf(std::ostream & os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
		os << indent << "Kappa = " << m_Kappa << std::endl;
	} /* End PrintSelf. */

	template< typename TInputMesh, typename TOutputMesh >
	void
		DeformableSimplexMesh3DThresholdForceFilter< TInputMesh, TOutputMesh >
		::Initialize()
	{
		Superclass::Initialize();

		m_GradientInterpolator = GradientInterpolatorType::New();
		m_GradientInterpolator->SetInputImage( this->GetGradient() );

		m_FeatureInterpolator = FeatureInterpolatorType::New();
		m_FeatureInterpolator->SetInputImage( this->GetFeatureImage() );
	}

	/** Compute model Displacement according to image gradient forces */
	template< typename TInputMesh, typename TOutputMesh >
	void
		DeformableSimplexMesh3DThresholdForceFilter< TInputMesh, TOutputMesh >
		::ComputeExternalForce(SimplexMeshGeometry *data,const GradientImageType *gradientImage)
	{
		PointType         vec_for, pt_cur;
		GradientPixelType gradient;

		gradient.Fill(0);
		pt_cur.CastFrom(data->pos);

		if ( m_GradientInterpolator->IsInsideBuffer( pt_cur ))
		{
			gradient = m_GradientInterpolator->Evaluate(pt_cur);
			vec_for[0] = gradient[0];
			vec_for[1] = gradient[1];
			vec_for[2] = gradient[2];
		}
		else
		{
			vec_for.Fill(0);
		}

		double mag = dot_product( data->normal.GetVnlVector(), vec_for.GetVnlVector() );

		vec_for[0] = this->GetBeta() * mag * ( data->normal )[0];
		vec_for[1] = this->GetBeta() * mag * ( data->normal )[1];
		vec_for[2] = this->GetBeta() * mag * ( data->normal )[2];

		double thr = 0.0;
		if(m_FeatureInterpolator->IsInsideBuffer(pt_cur))
		{
			double val = m_FeatureInterpolator->Evaluate(pt_cur);
			thr = val - m_Threshold;
		}
		else
		{
			thr = 0.0 - m_Threshold;
		}

		vec_for[0] += m_Kappa * thr * data->normal[0];
		vec_for[1] += m_Kappa * thr * data->normal[1];
		vec_for[2] += m_Kappa * thr * data->normal[2];

		//mag = vec_for.GetVectorFromOrigin().GetNorm();

		//if (mag > 0.5)
		//{
		//  for (int i=0; i<3; i++)
		//    vec_for[i] = (0.5 * vec_for[i])/mag;
		//}
		//
		data->externalForce[0] = vec_for[0];
		data->externalForce[1] = vec_for[1];
		data->externalForce[2] = vec_for[2];
	}
} /* end namespace itk. */

#endif //__itkDeformableSimplexMesh3DThresholdForceFilter_hxx
