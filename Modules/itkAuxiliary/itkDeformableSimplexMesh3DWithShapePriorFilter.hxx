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
#ifndef __itkDeformableSimplexMesh3DWithShapePriorFilter_hxx
#define __itkDeformableSimplexMesh3DWithShapePriorFilter_hxx

#include "itkDeformableSimplexMesh3DWithShapePriorFilter.h"
#include "itkNumericTraits.h"

#include <set>

namespace itk
{
	/* Constructor. */
	template< typename TInputMesh, typename TOutputMesh>
	DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh>
		::DeformableSimplexMesh3DWithShapePriorFilter()

	{
		m_Kappa = 0.1;

		m_TargetMesh = NULL;
		m_TargetSamplePoints = NULL;
		m_TargetKdTreeGenerator = NULL;

		m_ShapeMesh = NULL;
		m_ShapeSamplePoints = NULL;
		m_ShapeKdTreeGenerator = NULL;
		
		m_VarianceMap = NULL;
	}

	template< typename TInputMesh, typename TOutputMesh>
	DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh>
		::~DeformableSimplexMesh3DWithShapePriorFilter()
	{}

	/* PrintSelf. */
	template< typename TInputMesh, typename TOutputMesh>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh>
		::PrintSelf(std::ostream & os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
		os << indent << "Kappa = " << m_Kappa << std::endl;
	} /* End PrintSelf. */
	
	template< typename TInputMesh, typename TOutputMesh>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh>
		::Initialize()
	{
		const InputMeshType *             inputMesh = this->GetInput(0);
		const InputPointsContainer *      points    = inputMesh->GetPoints();
		InputPointsContainerConstIterator pointItr  = points->Begin();

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

		const GradientImageType* gradimg = this->GetGradient();
		if (!gradimg)
		{
			std::cerr<<"Error! No gradient image!"<<std::endl;
			return;
		}

		if (m_TargetMesh)
		{
			m_TargetSamplePoints = SampleType::New();
			m_TargetSamplePoints->SetMeasurementVectorSize( 3 );

			typename InputMeshType::PointsContainerConstIterator It = m_TargetMesh->GetPoints()->Begin();
			while( It != m_TargetMesh->GetPoints()->End() )
			{
				PointType point = It.Value();
				m_TargetSamplePoints->PushBack( point );
				++It;
			}

			m_TargetKdTreeGenerator = TreeGeneratorType::New();
			m_TargetKdTreeGenerator->SetSample( m_TargetSamplePoints );
			m_TargetKdTreeGenerator->SetBucketSize( 4 );
			m_TargetKdTreeGenerator->Update();
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
	
	template< typename TInputMesh, typename TOutputMesh >
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh >
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
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh >
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

			this->ComputeExternalForce(data, idx);

			displacement.SetVnlVector( m_Alpha * ( data->internalForce ).GetVnlVector()
				+ ( data->externalForce ).GetVnlVector() );

			data->pos += displacement;

			nonConstPoints->InsertElement(idx, data->pos);

			dataIt++;
		}
	}

	template< typename TInputMesh, typename TOutputMesh>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh>
		::ComputeExternalForce(SimplexMeshGeometry *data, unsigned int idx)
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
		double step_Target = 0.0;

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

		if (m_TargetMesh )
		{
			PointType closestTargetPt;
			if (m_P2p)
			{
				closestTargetPt = m_TargetMesh->GetPoint( idx );
			}
			else
			{
				typename TreeGeneratorType::KdTreeType::InstanceIdentifierVectorType neighbors;
				this->m_TargetKdTreeGenerator->GetOutput()->Search( pos_cur, 1u, neighbors );
				closestTargetPt = this->m_TargetKdTreeGenerator->GetOutput()->GetMeasurementVector( neighbors[0] );
			}

			double dist = closestTargetPt.EuclideanDistanceTo( pos_cur );
			vec_tmp = closestTargetPt - pos_cur;

			if (vec_tmp.GetNorm() > 0.0)
			{
				vec_tmp.Normalize();
			}

			step_Target = dot_product( data->normal.GetVnlVector(), vec_tmp.GetVnlVector() );
		}

		step_sum = shape_force_factor*step_shape + (1.0-shape_force_factor)*step_Target;

		vec_for = vec_normal * step_sum;

		data->externalForce = vec_for * this->GetKappa() + vec_for_gradient * this->GetBeta();
	}
	
	template< typename TInputMesh, typename TOutputMesh >
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh >
		::Intervene()
	{
		
	}

}/* end namespace itk. */

#endif //__itkDeformableSimplexMesh3DWithShapePriorFilter_hxx
