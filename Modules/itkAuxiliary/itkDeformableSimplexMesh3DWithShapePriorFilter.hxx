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
		m_ShapeMesh = km::cloneMesh<InputMeshType, InputMeshType>( inputMesh );
		km::copyMeshToMeshGeometry<InputMeshType, InputMeshType>(inputMesh, m_ShapeMesh);
		km::ComputeGeometry<TOutputMesh>(m_ShapeMesh);

		m_DeformedMesh = km::cloneMesh<InputMeshType, InputMeshType>( inputMesh );

		//Initialize SSM utils.
		this->m_SSMUtils.SetSSM(this->GetStatisticalModel());
		this->m_SSMUtils.SetRigidTransform(this->GetRigidTransform());
		this->m_SSMUtils.SetShapeTransform(this->GetShapeTransform());

		//Initialize transforms
		km::g_liverCentroid.CastFrom(km::getMeshCentroid<InputMeshType>(m_ShapeMesh));
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

		this->Cluster();

		m_Step = 0;

		while ( m_Step < m_Iterations )
		{
			const float progress = static_cast< float >( m_Step )
				/ static_cast< float >( m_Iterations );

			this->UpdateProgress(progress);

			this->IntervenePre();

			if (m_Phase == Lv2)
			{
				this->ComputeGeometry();

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
			this->m_ClassifierUtils.deformByLiverClassification(m_DeformedMesh, this->m_ShapeMesh, 1.5, 20);
		}
		else if (m_Phase == Lv2)
		{
			this->ComputeClusteredForce();

			// Filters should not modify their input...
			// There is a design flaw here.
			const InputMeshType *inputMesh = this->GetInput(0);
			InputPointsContainer *nonConstPoints = const_cast< InputPointsContainer * >( inputMesh->GetPoints() );
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
	}

	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::ComputeExternalForce(SimplexMeshGeometry *data, unsigned int idx)
	{
		VectorType vec_for;
		VectorType vec_for_gradient;
		VectorType vec_normal;
		vec_normal.SetVnlVector( data->normal.GetVnlVector() );
		vec_for_gradient.Fill( 0 );
		ClusterItem* clusterItem = m_ClusterPool.GetClusterItemByPointId(idx);
		vec_for = vec_normal*(clusterItem->clusterForce/clusterItem->clusterCount);
		data->externalForce = vec_for * this->GetKappa() + vec_for_gradient * this->GetBeta();
	}

	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::ComputeClusteredForce()
	{
		InputPointDataContainer *pointdata = const_cast< InputPointDataContainer * >( this->GetInput(0)->GetPointData() );
		SimplexMeshGeometry *data;
		unsigned int idx;
		typename GeometryMapType::Iterator dataIt = this->m_Data->Begin();
		while ( dataIt != this->m_Data->End() )
		{
			idx = dataIt.Index();
			data = dataIt.Value();

			ClusterItem* clusterItem = m_ClusterPool.GetClusterItemByPointId(idx);
			std::vector<double> features;
			this->m_ProfileExtractor.extractFeatureSet(features, this->m_LiverClassifier->profileCategory, data, data->pos);
			double force = 2*this->m_LiverClassifier->classify(features, idx)-1.0;
			clusterItem->clusterForce += force;
			//pointdata->InsertElement(idx, clusterItem->clusterForce);
			dataIt++;
		}
	}
	
	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::IntervenePre()
	{
		if(m_Phase == Lv2)
		{
			//Clear cluster force.
			this->m_ClusterPool.ClearForce();
		}
	}

	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::IntervenePost()
	{
		if (m_Phase == Lv1)
		{
			this->UpdateShape();

			InputMeshType *inputMesh = const_cast<InputMeshType*>(this->GetInput(0));
			km::copyMeshToMeshPoints<InputMeshType, InputMeshType>(this->m_ShapeMesh, inputMesh);

			double shapeParamDiff = this->m_SSMUtils.calShapeParaDiff();
			if(shapeParamDiff < 0.02)
			{
				this->m_Phase = Lv2;
			}
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
		::Cluster()
	{
		std::cout<<"Start to cluster.."<<std::endl;
		vtkSmartPointer<vtkPolyData> polydata = km::mesh2PolyData<InputMeshType>(this->m_ShapeMesh);
		std::cout<<"Convert mesh to vtkpolydata done.."<<std::endl;
		vtkSurface* surface=vtkSurface::New();
		surface->CreateFromPolyData(polydata);
		std::cout<<"Vtk surface create successfully."<<std::endl;

		vtkIsotropicDiscreteRemeshing *remesh=vtkIsotropicDiscreteRemeshing::New();
		remesh->SetInput(surface);
		remesh->SetFileLoadSaveOption(0);
		remesh->SetNumberOfClusters( polydata->GetNumberOfPoints()/3 );
		remesh->SetConsoleOutput(2);
		remesh->SetSubsamplingThreshold(10);
		remesh->GetMetric()->SetGradation(1);
		remesh->SetDisplay(0);
		//remesh->Remesh();
		std::cout<<"Now call ProcessClustering().."<<std::endl;
		remesh->ProcessClustering();
		std::cout<<"Number of clusters: "<<remesh->GetNumberOfClusters()<<std::endl;
		//this->m_ShapeMesh->GetPointData()->Reserve(m_ShapeMesh->GetNumberOfPoints());
		for (int pointid = 0; pointid < remesh->GetNumberOfItems(); pointid++)
		{
			int clusterLabel = remesh->GetClustering()->GetValue(pointid);
			this->m_ClusterPool.AddClusterMapping(pointid, clusterLabel);
			//this->m_ShapeMesh->SetPointData(pointid, clusterLabel);
		}
		//km::writeMesh<InputMeshType>("D:\\clusteredMesh.vtk", this->m_ShapeMesh);
		std::cout<<"Clustering done."<<std::endl;
	}

	template< class TInputMesh, class TOutputMesh, class TInputImage, class TStatisticalModel, class TRigidTransform, class TShapeTransform>
	void
		DeformableSimplexMesh3DWithShapePriorFilter< TInputMesh, TOutputMesh, TInputImage, TStatisticalModel, TRigidTransform, TShapeTransform>
		::UpdateShape()
	{
		//Update shape
		const InputMeshType *inputMesh = this->GetInput(0);
		this->m_SSMUtils.compositeTransformFitting(inputMesh);

		km::transformMesh<InputMeshType, CompositeTransformType>(m_ReferenceShapeMesh, m_ShapeMesh, this->m_CompositeTransform);
		km::ComputeGeometry<TOutputMesh>(m_ShapeMesh);

		km::g_liverCentroid.CastFrom(km::getMeshCentroid<InputMeshType>(m_ShapeMesh));
		this->GetRigidTransform()->SetCenter(km::g_liverCentroid);
	}

}/* end namespace itk. */

#endif //__itkDeformableSimplexMesh3DWithShapePriorFilter_hxx
