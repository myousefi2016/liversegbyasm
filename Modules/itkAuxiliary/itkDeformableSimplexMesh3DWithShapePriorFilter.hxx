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
	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::DeformableSimplexMesh3DWithShapePriorFilter()

	{
		m_Kappa = 0.1;
	}

	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::~DeformableSimplexMesh3DWithShapePriorFilter()
	{}

	/* PrintSelf. */
	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::PrintSelf(std::ostream & os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
		os << indent << "Kappa = " << m_Kappa << std::endl;
	} /* End PrintSelf. */
	
	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::Initialize()
	{
		//Check input.
		if (this->GetInputImage() != NULL){
			this->m_ProfileExtractor.setImage(this->GetInputImage());
			this->m_ProfileExtractor.enableCache(false);
		}else{
			std::cerr<<"Input image is NULL."<<std::endl;
			return;
		}

		if (this->GetBoundaryClassifier()==NULL){
			std::cerr<<"Boundary classifier is NULL."<<std::endl;
			return;
		}
		if (this->GetLiverClassifier()==NULL){
			std::cerr<<"Liver classifier is NULL."<<std::endl;
			return;
		}

		if (this->GetStatisticalModel() != NULL){
			this->m_SSMUtils.SetSSM(this->GetStatisticalModel());
			std::cerr<<"SSM is NULL."<<std::endl;
			return;
		}
		if (this->GetRigidTransform() != NULL){
			this->m_SSMUtils.SetRigidTransform(this->GetRigidTransform());
			std::cerr<<"Rigid transform is NULL."<<std::endl;
			return;
		}
		if (this->GetShapeTransform() != NULL){
			this->m_SSMUtils.SetShapeTransform(this->GetShapeTransform());
			std::cerr<<"Shape transform is NULL."<<std::endl;
			return;
		}
		//Check input done.

		const InputMeshType *             inputMesh = this->GetInput(0);
		const InputPointsContainer *      points    = inputMesh->GetPoints();
		InputPointsContainerConstIterator pointItr  = points->Begin();

		this->m_Data = inputMesh->GetGeometryData();
		if ( this->m_Data.IsNull() || this->m_Data->Size()==0)
		{
			std::cout<<"Geometry data is empty!!!"<<std::endl;
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
	}
	
	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
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
	
	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
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

	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
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

		step_sum = shape_force_factor*step_shape + (1.0-shape_force_factor)*step_Target;

		vec_for = vec_normal * step_sum;

		data->externalForce = vec_for * this->GetKappa() + vec_for_gradient * this->GetBeta();
	}
	
	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::Intervene()
	{
		
	}

}/* end namespace itk. */

#endif //__itkDeformableSimplexMesh3DWithShapePriorFilter_hxx
