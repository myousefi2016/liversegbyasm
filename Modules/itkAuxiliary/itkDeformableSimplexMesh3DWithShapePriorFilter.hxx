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
#include "vtkSurface.h"
#include "vtkIsotropicDiscreteRemeshing.h"

#include <set>

namespace itk
{
	/* Constructor. */
	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::DeformableSimplexMesh3DWithShapePriorFilter()

	{
		m_Kappa = 0.1;
		m_NumberOfShapeClusters = 1;
		m_MinShapeDifference = 0.000001;//0.08;
		m_IterationsLv1 = 10;
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

		if (this->GetStatisticalModel() == NULL){
			std::cerr<<"SSM is NULL."<<std::endl;
			return;
		}
		if (this->GetRigidTransform() == NULL){
			std::cerr<<"Rigid transform is NULL."<<std::endl;
			return;
		}
		if (this->GetShapeTransform() == NULL){
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
		m_Pointdata = nonconstInputMesh->GetPointData();
		if (m_Pointdata->Size() <= 0)
		{
			m_Pointdata->Reserve( nonconstInputMesh->GetNumberOfPoints() );
			for (int i=0;i<nonconstInputMesh->GetNumberOfPoints();i++)
			{
				m_Pointdata->SetElement( i, 0.0 );	
			}
		}

		//Initialize shape mesh
		m_ReferenceShapeMesh = m_Model->GetRepresenter()->GetReference();
		m_UpdatedShapeMesh = km::cloneMesh<InputMeshType, InputMeshType>( inputMesh );
		km::copyMeshToMeshGeometry<InputMeshType, InputMeshType>(inputMesh, m_UpdatedShapeMesh);
		km::ComputeGeometry<TOutputMesh>(m_UpdatedShapeMesh, true);

		km::g_liverCentroid.CastFrom(km::getMeshCentroid<InputMeshType>(m_UpdatedShapeMesh));

		m_NumberOfShapeClusters = km::g_number_clusters;

		//Initialize SSM utils.
		this->m_SSMUtils.SetSSM(this->GetStatisticalModel());
		this->m_SSMUtils.SetRigidTransform(this->GetRigidTransform());
		this->m_SSMUtils.SetShapeTransform(this->GetShapeTransform());
		this->m_SSMUtils.SetNumberOfClusters(m_NumberOfShapeClusters);
		this->m_SSMUtils.Initialize();

		//Initialize transforms
		this->GetRigidTransform()->SetCenter(km::g_liverCentroid);
		this->m_CompositeTransform = CompositeTransformType::New();
		this->m_CompositeTransform->AddTransform( this->GetRigidTransform() );
		this->m_CompositeTransform->AddTransform( this->GetShapeTransform() );

		this->m_ClassifierUtils.SetBoundaryClassifier(m_BoundaryClassifier);
		this->m_ClassifierUtils.SetRegionClassifier(m_LiverClassifier);
		this->m_ClassifierUtils.SetProfileExtractor(&m_ProfileExtractor);

		this->m_Phase = Lv1;

		for (int i=0;i<this->GetInput(0)->GetNumberOfPoints();i++)
		{
			m_Forces[i] = 0.0;
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

			std::cout<<".........."<<progress<<"........."<<std::endl;

			this->UpdateProgress(progress);

			this->IntervenePre();

			this->ComputeGeometry();

			if (m_Phase == Lv2)
			{
				if ( m_Step % 10 == 0 && m_Step > 0 )
				{
					this->UpdateReferenceMetrics();
				}
			}

			this->ComputeDisplacement();

			this->IntervenePost();

			m_Step++;
		}
	}
	
	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::ComputeDisplacement()
	{
		if (m_Phase == Lv1)
		{
			VectorType displacement;
			SimplexMeshGeometry *data;
			unsigned int idx;
			InputPointsContainer *nonConstPoints = const_cast< InputPointsContainer * >( this->GetInput(0)->GetPoints() );
			GeometryMapType::Iterator dataIt = this->m_Data->Begin();
			while ( dataIt != this->m_Data->End() )
			{
				idx = dataIt.Index();
				data = dataIt.Value();

				PointType tmpPt;
				this->m_ClassifierUtils.FindNextRegionPoint(tmpPt, data->pos, data, idx, 1.5);
				m_Forces[idx] = dot_product((tmpPt-data->pos).GetVnlVector(), data->normal.GetVnlVector());

				//m_Forces[idx] *= (1.0 - std::abs(m_Errors[idx])); 

				dataIt++;
			}

			double smoothFactor = 0.2; //Indicate how much the external force for current landmark weights.
			dataIt = this->m_Data->Begin();
			while ( dataIt != this->m_Data->End() )
			{
				idx = dataIt.Index();
				data = dataIt.Value();

				double force = 0;
				NeighborSetIterator neighborSetIt = data->neighborSet->begin();
				while(neighborSetIt != data->neighborSet->end())
				{
					int neighborIdx = *neighborSetIt;
					force += m_Forces[neighborIdx];
					neighborSetIt++;
				}
				force = smoothFactor*m_Forces[idx]+(1.0-smoothFactor)*(force/data->neighborSet->size());
				m_Forces[idx] = force;

				displacement.SetVnlVector(data->normal.Get_vnl_vector()*m_Forces[idx]);
				
				data->pos += displacement;
				nonConstPoints->InsertElement(idx, data->pos);

				dataIt++;
			}

		}
		else if (m_Phase == Lv2)
		{
			// Filters should not modify their input...
			// There is a design flaw here.
			const InputMeshType *inputMesh = this->GetInput(0);
			InputPointsContainer *nonConstPoints = const_cast< InputPointsContainer * >( inputMesh->GetPoints() );
			InputPointDataContainer *pointdata = const_cast< InputPointDataContainer * >( inputMesh->GetPointData() );

			VectorType displacement;
			SimplexMeshGeometry *data;
			unsigned int idx;
			typename GeometryMapType::Iterator dataIt = this->m_Data->Begin();
			while ( dataIt != this->m_Data->End() )
			{
				idx = dataIt.Index();
				data = dataIt.Value();

				std::vector<double> features;
				this->m_ProfileExtractor.extractFeatureSet(features, this->m_LiverClassifier->profileCategory, data, data->pos);
				double force = 2*this->m_LiverClassifier->classify(features, idx)-1.0;

				switch( m_SSMUtils.GetLandmarkStatus(idx) )
				{
				case Leaking:
					force = 0;
					break;
				case Abnormal:
					force = 0;
					break;
				}

				m_Forces[idx] = force;
				dataIt++;
			}

			dataIt = this->m_Data->Begin();
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
	}

	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::ComputeExternalForce(SimplexMeshGeometry *data, unsigned int idx)
	{
		VectorType vec_for_classify;
		double force = m_Forces[idx];
		NeighborSetIterator neighborSetIt = data->neighborSet->begin();
		while(neighborSetIt != data->neighborSet->end())
		{
			int neighborIdx = *neighborSetIt;
			force += m_Forces[neighborIdx];
			neighborSetIt++;
		}
		force /= (data->neighborSet->size()+1);
		m_Forces[idx] = force;
		vec_for_classify.SetVnlVector( data->normal.GetVnlVector() * force );

		VectorType vec_for_curvature;
		vec_for_curvature.Fill(0);
		SimplexMeshGeometry *data_shape = this->m_UpdatedShapeMesh->GetGeometryData()->GetElement(idx);
		double phi_shape = data_shape->phi;
		double phi_cur = data->phi;
		double curvatureDiff = (180.0/(20.0*3.14))*(phi_cur - phi_shape); //20 degree as a threshold
		if (curvatureDiff > 1.0){
			curvatureDiff = 1.0;
		}else if (curvatureDiff < -1.0){
			curvatureDiff = -1.0;
		}
		vec_for_curvature.SetVnlVector( data->normal.GetVnlVector() * curvatureDiff );

		data->externalForce = vec_for_classify * this->GetKappa() + vec_for_curvature * this->GetBeta();
	}
	
	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::IntervenePre()
	{
		//Clear cluster force.
		for (int i=0;i<this->GetInput(0)->GetNumberOfPoints();i++)
		{
			m_Forces[i] = 0.0;
		}
	}

	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::IntervenePost()
	{
		if (m_Phase == Lv1)
		{
			//Update shape
			this->UpdateShape();

			this->UpdateLandmarkStatus();

			for (int i=0;i<m_Pointdata->size();i++)
			{
				m_Pointdata->InsertElement(i, (int)m_SSMUtils.GetLandmarkStatus(i));
			}
			km::writeMesh<InputMeshType>(km::g_output_dir, "internalDeformedMesh", m_Step, ".vtk", this->GetInput(0));
			km::writeMesh<InputMeshType>(km::g_output_dir, "internalShapeMesh", m_Step, ".vtk", this->m_UpdatedShapeMesh);

			//Copy updated shape points.
			const InputMeshType *inputMesh = this->GetInput(0);
			InputPointsContainer *nonConstPoints = const_cast< InputPointsContainer * >( inputMesh->GetPoints() );
			typename GeometryMapType::Iterator dataIt = this->m_Data->Begin();
			SimplexMeshGeometry *data;
			unsigned int idx;
			while ( dataIt != this->m_Data->End() )
			{
				idx = dataIt.Index();
				data = dataIt.Value();
				data->pos = this->m_UpdatedShapeMesh->GetPoint(idx);
				nonConstPoints->InsertElement(idx, data->pos);
				dataIt++;
			}

			//Check shape updated difference.
			//double shapeParamDiff = this->m_SSMUtils.CalShapeParaDiff();
			//std::cout<<"Iteration "<<m_Step<<", shape difference: "<<shapeParamDiff<<std::endl;
			if(m_Step >= m_IterationsLv1)
			{
				this->m_Phase = Lv2;
				//m_Step = m_Iterations;
			}
		}
		else if (m_Phase == Lv2)
		{
			//Update shape
			if (this->GetStep()%10 == 9)
			{
				this->UpdateShape();

				for (int i=0;i<m_Pointdata->size();i++)
				{
					m_Pointdata->InsertElement(i, (int)m_SSMUtils.GetLandmarkStatus(i));
				}
				km::writeMesh<InputMeshType>(km::g_output_dir, "internalDeformedMesh", m_Step, ".vtk", this->GetInput(0));
				km::writeMesh<InputMeshType>(km::g_output_dir, "internalShapeMesh", m_Step, ".vtk", this->m_UpdatedShapeMesh);
			}
		}
	}

	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::UpdateShape()
	{
		const InputMeshType *inputMesh = this->GetInput(0);

		//Rigid & shape fitting
		this->m_SSMUtils.Update(inputMesh, m_UpdatedShapeMesh);
		km::ComputeGeometry<TOutputMesh>(m_UpdatedShapeMesh, true);
	}

	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::UpdateLandmarkStatus()
	{
		if (!km::g_landmark_status_evalution)
		{
			return;
		}

		std::vector<double> fittedErrors;
		//fittedErrors.reserve(this->GetInput(0)->GetNumberOfPoints());
		fittedErrors.resize(this->GetInput(0)->GetNumberOfPoints());

		int idx;
		SimplexMeshGeometry *data;
		GeometryMapType::Iterator dataIt = this->m_Data->Begin();
		while ( dataIt != this->m_Data->End() )
		{
			idx = dataIt.Index();
			data = dataIt.Value();

			VectorType diffVector = this->GetInput(0)->GetPoint(idx) - m_UpdatedShapeMesh->GetPoint(idx);

			fittedErrors[idx] = dot_product(diffVector.GetVnlVector(), data->normal.Get_vnl_vector());

			dataIt++;
		}

		double meanFittedError = 0.0;
		double sigmaFittedError = 0.0;
		km::Math::calculateMeanAndSigma(fittedErrors, meanFittedError, sigmaFittedError);

		dataIt = this->m_Data->Begin();
		while ( dataIt != this->m_Data->End() )
		{
			idx = dataIt.Index();
			if (fittedErrors[idx] > meanFittedError+2.5*sigmaFittedError){
				m_SSMUtils.SetLandmarkStatus(idx, Leaking);
			}else if(fittedErrors[idx] < meanFittedError-2.5*sigmaFittedError)
			{
				m_SSMUtils.SetLandmarkStatus(idx, Abnormal);
			}else{
				m_SSMUtils.SetLandmarkStatus(idx, Normal);
			}

			dataIt++;
		}
	}
}/* end namespace itk. */

#endif //__itkDeformableSimplexMesh3DWithShapePriorFilter_hxx
