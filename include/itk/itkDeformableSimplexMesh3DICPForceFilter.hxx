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
#ifndef __itkDeformableSimplexMesh3DICPForceFilter_hxx
#define __itkDeformableSimplexMesh3DICPForceFilter_hxx

#include "itkDeformableSimplexMesh3DICPForceFilter.h"
#include "itkNumericTraits.h"

#include <set>

namespace itk
{
/* Constructor. */
template< typename TInputMesh, typename TOutputMesh >
DeformableSimplexMesh3DICPForceFilter< TInputMesh, TOutputMesh >
::DeformableSimplexMesh3DICPForceFilter()

{
  m_Kappa = 0.1;

  this->m_KdTreeGenerator = NULL;
  this->m_SamplePoints = NULL;
  this->m_TargetMesh = NULL;
}

template< typename TInputMesh, typename TOutputMesh >
DeformableSimplexMesh3DICPForceFilter< TInputMesh, TOutputMesh >
::~DeformableSimplexMesh3DICPForceFilter()
{}

/* PrintSelf. */
template< typename TInputMesh, typename TOutputMesh >
void
DeformableSimplexMesh3DICPForceFilter< TInputMesh, TOutputMesh >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Kappa = " << m_Kappa << std::endl;
} /* End PrintSelf. */

template< typename TInputMesh, typename TOutputMesh >
void
DeformableSimplexMesh3DICPForceFilter< TInputMesh, TOutputMesh >
::Initialize()
{
	Superclass::Initialize();

	m_GradientInterpolator = GradientInterpolatorType::New();
	m_GradientInterpolator->SetInputImage( this->GetGradient() );

    this->m_SamplePoints = SampleType::New();
    this->m_SamplePoints->SetMeasurementVectorSize( 3 );

    typename InputMeshType::PointsContainerConstIterator It = this->m_TargetMesh->GetPoints()->Begin();
	while( It != this->m_TargetMesh->GetPoints()->End() )
	{
		PointType point = It.Value();
		this->m_SamplePoints->PushBack( point );
		++It;
	}

	this->m_KdTreeGenerator = TreeGeneratorType::New();
	this->m_KdTreeGenerator->SetSample( this->m_SamplePoints );
	this->m_KdTreeGenerator->SetBucketSize( 4 );
	this->m_KdTreeGenerator->Update();
}

/** Compute model Displacement according to image gradient forces */
template< typename TInputMesh, typename TOutputMesh >
void
DeformableSimplexMesh3DICPForceFilter< TInputMesh, TOutputMesh >
::ComputeExternalForce(SimplexMeshGeometry *data,const GradientImageType *gradientImage)
{
  PointType vec_for, vec_tmp1;
  
  PointType origpt;
  origpt.CastFrom( data->pos );

  unsigned int num_neighbor = 1u;
  
  typename TreeGeneratorType::KdTreeType::InstanceIdentifierVectorType neighbors;
  this->m_KdTreeGenerator->GetOutput()->Search( origpt, num_neighbor, neighbors );

  double max_mag = 0.0;
  
  for ( unsigned int i=0;i<num_neighbor;i++ )
  {
    PointType tmppt = this->m_KdTreeGenerator->GetOutput()->GetMeasurementVector( neighbors[i] );
	vec_tmp1[0] = tmppt[0] - origpt[0];
	vec_tmp1[1] = tmppt[1] - origpt[1];
	vec_tmp1[2] = tmppt[2] - origpt[2];

    double tmpmag = dot_product( data->normal.GetVnlVector(), vec_tmp1.GetVnlVector() );

	if ( std::abs(tmpmag) > std::abs(max_mag) )
	{
		max_mag = tmpmag;
	}
  }

  if (max_mag > 0)
  {
	max_mag = std::min(max_mag, 0.3);
  }
  else if(max_mag <= 0)
  {
	max_mag = std::max(max_mag, -0.3);
  }

  max_mag += 0.1;

  vec_for[0] = this->GetKappa() * max_mag * data->normal[0];
  vec_for[1] = this->GetKappa() * max_mag * data->normal[1];
  vec_for[2] = this->GetKappa() * max_mag * data->normal[2];
  
  data->externalForce[0] = vec_for[0];
  data->externalForce[1] = vec_for[1];
  data->externalForce[2] = vec_for[2];
}
} /* end namespace itk. */

#endif //__itkDeformableSimplexMesh3DICPForceFilter_hxx
