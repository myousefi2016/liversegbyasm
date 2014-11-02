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
#ifndef __itkAdaptSimplexMesh3DFilter_hxx
#define __itkAdaptSimplexMesh3DFilter_hxx

#include "itkAdaptSimplexMesh3DFilter.h"
#include "itkNumericTraits.h"

#include <set>

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
	AdaptSimplexMesh3DFilter< TInputMesh, TOutputMesh >
		::AdaptSimplexMesh3DFilter()
	{
		m_Step = 0;
		m_Iterations = 20;
		m_Alpha = 0.2;
		m_Beta  = 0.01;
		m_Gamma = 0.05;
		m_Damping = 0.65;
		m_Rigidity = 1;

		m_TangentFactor = 0.5;
		m_NormalFactor = 0.5;

		m_ImageDepth = 0;
		m_ImageHeight = 0;
		m_ImageWidth = 0;

		this->ProcessObject::SetNumberOfRequiredInputs(1);

		OutputMeshPointer output = OutputMeshType::New();
		this->ProcessObject::SetNumberOfRequiredOutputs(1);
		this->ProcessObject::SetNthOutput( 0, output.GetPointer() );

		this->m_Data = NULL;
	}

	template< typename TInputMesh, typename TOutputMesh >
	AdaptSimplexMesh3DFilter< TInputMesh, TOutputMesh >
		::~AdaptSimplexMesh3DFilter()
	{}

	/* Generate Data */
	template< typename TInputMesh, typename TOutputMesh >
	void
		AdaptSimplexMesh3DFilter< TInputMesh, TOutputMesh >
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
			m_Step++;
		}

		const InputMeshType *             inputMesh = this->GetInput(0);
		const InputPointsContainer *      points = inputMesh->GetPoints();
		InputPointsContainerConstIterator pointItr = points->Begin();

		while ( pointItr != points->End() )
		{
			SimplexMeshGeometry *data;
			IdentifierType       idx = pointItr.Index();
			data = this->m_Data->GetElement(idx);
			delete data->neighborSet;
			data->neighborSet = NULL;
			pointItr++;
		}

		//this->ComputeOutput();
	}

	template< typename TInputMesh, typename TOutputMesh >
	void
		AdaptSimplexMesh3DFilter< TInputMesh, TOutputMesh >
		::ComputeDisplacement()
	{
		const InputMeshType *inputMesh = this->GetInput(0);

		// Filters should not modify their input...
		// There is a design flaw here.
		InputPointsContainer *nonConstPoints =
			const_cast< InputPointsContainer * >( inputMesh->GetPoints() );

		typename GeometryMapType::Iterator dataIt = this->m_Data->Begin();
		SimplexMeshGeometry *data;
		VectorType           displacement;

		while ( dataIt != this->m_Data->End() )
		{
			data = dataIt.Value();

			this->ComputeInternalForce(data);

			displacement.SetVnlVector( m_Alpha * ( data->internalForce ).GetVnlVector() );

			data->pos += displacement;
			nonConstPoints->InsertElement(dataIt.Index(), data->pos);

			dataIt++;
		}
	}

	/* */
	template< typename TInputMesh, typename TOutputMesh >
	void
		AdaptSimplexMesh3DFilter< TInputMesh, TOutputMesh >
		::ComputeInternalForce(SimplexMeshGeometry *data)
	{
		VectorType tangentForce, normalForce;
		double     eps1Diff, eps2Diff, eps3Diff;
		//    double diffAbsSum;
		double           d, phi, r;
		NeighborSetType *neighborSet;
		PointType        xOrig;
		PointType        eps, epsRef;
		double           phiRef;
		PointType        f_int;

		xOrig = data->pos;
		eps = data->eps;
		epsRef = data->referenceMetrics;

		eps1Diff = epsRef[0] - eps[0];
		eps2Diff = epsRef[1] - eps[1];
		eps3Diff = epsRef[2] - eps[2];
		//    diffAbsSum = vcl_abs(eps1Diff)+vcl_abs(eps2Diff)+vcl_abs(eps3Diff);

		tangentForce.SetVnlVector(eps1Diff * ( data->neighbors[0] ).GetVnlVector()
			+ eps2Diff * ( data->neighbors[1] ).GetVnlVector()
			+ eps3Diff * ( data->neighbors[2] ).GetVnlVector()
			);

		r = data->circleRadius;
		d = data->distance;
		phi = data->phi;

		neighborSet = data->neighborSet;

		NeighborSetIterator neighborIt = neighborSet->begin();
		phiRef = 0.0;

		while ( neighborIt != neighborSet->end() )
		{
			phiRef += this->m_Data->GetElement(*neighborIt++)->phi;
		}
		phiRef /= (double)neighborSet->size();

		double L = L_Func(r, d, phi);
		double L_Ref = L_Func(r, d, phiRef);

		normalForce.SetVnlVector( -1.0 * ( L_Ref - L ) * ( data->normal ).GetVnlVector() );

		data->internalForce.Fill(0.0);

		// quick hack fixing for div by zero error
		if ( L_Ref != (double)NumericTraits< IdentifierType >::max() && L != (double)NumericTraits< IdentifierType >::max() )
		{
			data->internalForce += m_TangentFactor*tangentForce + m_NormalFactor*normalForce;
		}
	}

} /* end namespace itk. */

#endif //__itkAdaptSimplexMesh3DFilter_TXX
