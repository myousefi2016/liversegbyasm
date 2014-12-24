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
		InputPointDataContainer* pointdata = nonconstInputMesh->GetPointData();
		if (pointdata->Size() <= 0)
		{
			pointdata->Reserve( nonconstInputMesh->GetNumberOfPoints() );
			for (int i=0;i<nonconstInputMesh->GetNumberOfPoints();i++)
			{
				nonconstInputMesh->SetPointData( i, 0.0 );	
			}
		}

		//Initialize shape mesh
		m_ReferenceShapeMesh = this->GetStatisticalModel()->GetRepresenter()->GetReference();
		m_UpdatedShapeMesh = km::cloneMesh<InputMeshType, InputMeshType>( inputMesh );
		km::copyMeshToMeshGeometry<InputMeshType, InputMeshType>(inputMesh, m_UpdatedShapeMesh);
		km::ComputeGeometry<TOutputMesh>(m_UpdatedShapeMesh, true);

		km::g_liverCentroid.CastFrom(km::getMeshCentroid<InputMeshType>(m_UpdatedShapeMesh));

		m_MeanShapeMesh = this->GetStatisticalModel()->DrawMean();
		m_UpdatedMeanShapeMesh = km::cloneMesh<InputMeshType, InputMeshType>( m_MeanShapeMesh );
		km::copyMeshToMeshGeometry<InputMeshType, InputMeshType>(inputMesh, m_UpdatedMeanShapeMesh);
		km::ComputeGeometry<TOutputMesh>(m_UpdatedMeanShapeMesh, false);

		//Initialize SSM utils.
		this->m_SSMUtils.SetSSM(this->GetStatisticalModel());
		this->m_SSMUtils.SetRigidTransform(this->GetRigidTransform());
		this->m_SSMUtils.SetShapeTransform(this->GetShapeTransform());

		//Initialize transforms
		this->GetRigidTransform()->SetCenter(km::g_liverCentroid);
		this->m_CompositeTransform = CompositeTransformType::New();
		this->m_CompositeTransform->AddTransform( this->GetRigidTransform() );
		this->m_CompositeTransform->AddTransform( this->GetShapeTransform() );

		this->m_ClassifierUtils.SetBoundaryClassifier(m_BoundaryClassifier);
		this->m_ClassifierUtils.SetRegionClassifier(m_LiverClassifier);
		this->m_ClassifierUtils.SetProfileExtractor(&m_ProfileExtractor);

		this->m_Phase = Lv1;
	}
	
	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::GenerateData()
	{
		this->Initialize();

		this->Cluster( 7 );

		m_Step = 0;

		while ( m_Step < m_Iterations )
		{
			const float progress = static_cast< float >( m_Step )
				/ static_cast< float >( m_Iterations );

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
			ClusterItem* clusterItem;
			GeometryMapType::Iterator dataIt = this->m_Data->Begin();
			while ( dataIt != this->m_Data->End() )
			{
				idx = dataIt.Index();
				data = dataIt.Value();

				clusterItem = m_ClusterPool.GetClusterItemByPointId(idx);
				PointType tmpPt;
				this->m_ClassifierUtils.FindNextRegionPoint(tmpPt, data->pos, data, idx, 1.5);
				displacement = tmpPt - data->pos;
				clusterItem->force += dot_product(displacement.GetVnlVector(), data->normal.GetVnlVector());

				dataIt++;
			}

			this->m_ClusterPool.Update();

			InputPointsContainer *nonConstPoints = const_cast< InputPointsContainer * >( this->GetInput(0)->GetPoints() );
			dataIt = this->m_Data->Begin();
			while ( dataIt != this->m_Data->End() )
			{
				idx = dataIt.Index();
				data = dataIt.Value();

				clusterItem = m_ClusterPool.GetClusterItemByPointId(idx);
				displacement.SetVnlVector( data->normal.GetVnlVector() );
				displacement *= clusterItem->force;
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
			ClusterItem* clusterItem;
			typename GeometryMapType::Iterator dataIt = this->m_Data->Begin();
			while ( dataIt != this->m_Data->End() )
			{
				idx = dataIt.Index();
				data = dataIt.Value();

				clusterItem = m_ClusterPool.GetClusterItemByPointId(idx);
				std::vector<double> features;
				this->m_ProfileExtractor.extractFeatureSet(features, this->m_LiverClassifier->profileCategory, data, data->pos);
				double force = 2*this->m_LiverClassifier->classify(features, idx)-0.9;
				clusterItem->force += force;

				dataIt++;
			}

			this->m_ClusterPool.Update();

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
		ClusterItem* clusterItem = m_ClusterPool.GetClusterItemByPointId(idx);

		VectorType vec_for_curvature;
		//vec_for_curvature.Fill(0);
		//vec_for_curvature.SetVnlVector(data->normal.GetVnlVector());

		VectorType vec_for_classify;
		vec_for_classify.SetVnlVector(data->normal.GetVnlVector() * clusterItem->force);

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
		this->m_ClusterPool.CleanCache();
	}

	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::IntervenePost()
	{
		if (m_Phase == Lv1)
		{
			this->UpdateShape();

			//km::writeMesh<InputMeshType>(km::g_output_dir, "internalDeformedMesh", m_Step, ".vtk", this->GetInput(0));
			km::writeMesh<InputMeshType>(km::g_output_dir, "internalMeanMesh", m_Step, ".vtk", this->m_UpdatedMeanShapeMesh);
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
			double shapeParamDiff = this->m_SSMUtils.calShapeParaDiff();
			if(shapeParamDiff < 0.02)
			{
				this->m_Phase = Lv2;
			}

			//this->m_ClassifierUtils.updateShapeNormals(this->m_ShapeMeshBeforeFitting, this->m_ShapeMesh);
		}
		else if (m_Phase == Lv2)
		{
			//Update shape
			if (this->GetStep()%10 == 0)
			{
				this->UpdateShape();
			}
		}	
	}

	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::Cluster(int blockSize)
	{
		this->m_ClusterPool.Reset();

		std::cout<<"Start to cluster.."<<std::endl;
		vtkSmartPointer<vtkPolyData> polydata = km::mesh2PolyData<InputMeshType>(this->m_UpdatedShapeMesh);
		std::cout<<"Convert mesh to vtkpolydata done.."<<std::endl;
		vtkSurface* surface=vtkSurface::New();
		surface->CreateFromPolyData(polydata);
		std::cout<<"Vtk surface create successfully."<<std::endl;

		vtkIsotropicDiscreteRemeshing *remesh=vtkIsotropicDiscreteRemeshing::New();
		remesh->SetInput(surface);
		remesh->SetFileLoadSaveOption(0);
		remesh->SetNumberOfClusters( polydata->GetNumberOfPoints()/blockSize );
		remesh->SetConsoleOutput(2);
		remesh->SetSubsamplingThreshold(10);
		remesh->GetMetric()->SetGradation(2);
		remesh->SetDisplay(0);
		//remesh->Remesh();
		std::cout<<"Now call ProcessClustering().."<<std::endl;
		remesh->ProcessClustering();
		std::cout<<"Number of clusters: "<<remesh->GetNumberOfClusters()<<std::endl;
		this->m_UpdatedShapeMesh->GetPointData()->Reserve(m_UpdatedShapeMesh->GetNumberOfPoints());
		for (int pointid = 0; pointid < remesh->GetNumberOfItems(); pointid++)
		{
			int clusterLabel = remesh->GetClustering()->GetValue(pointid);
			this->m_ClusterPool.AddClusterMapping(pointid, clusterLabel);
			this->m_UpdatedShapeMesh->SetPointData(pointid, clusterLabel);
		}
		km::writeMesh<InputMeshType>(km::g_output_dir, "clusteredMesh", blockSize, ".vtk", this->m_UpdatedShapeMesh);
		std::cout<<"Clustering done."<<std::endl;
	}

	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::UpdateShape()
	{
		//Update shape
		const InputMeshType *inputMesh = this->GetInput(0);
		//this->m_SSMUtils.compositeTransformFitting(inputMesh);
		this->m_SSMUtils.rigidTransformFitting(inputMesh);
		this->m_SSMUtils.shapeTransformFitting(inputMesh);

		km::transformMesh<InputMeshType, CompositeTransformType>(m_ReferenceShapeMesh, m_UpdatedShapeMesh, this->m_CompositeTransform);
		km::ComputeGeometry<TOutputMesh>(m_UpdatedShapeMesh, true);

		km::transformMesh<InputMeshType, RigidTransformType>(m_MeanShapeMesh, m_UpdatedMeanShapeMesh, this->m_RigidTransform);
		km::ComputeGeometry<TOutputMesh>(m_UpdatedMeanShapeMesh, false);

		km::g_liverCentroid.CastFrom(km::getMeshCentroid<InputMeshType>(m_UpdatedShapeMesh));
		this->GetRigidTransform()->SetCenter(km::g_liverCentroid);
	}

}/* end namespace itk. */

#endif //__itkDeformableSimplexMesh3DWithShapePriorFilter_hxx
