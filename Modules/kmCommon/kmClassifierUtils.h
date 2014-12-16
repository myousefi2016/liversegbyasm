#ifndef __kmClassifierUtils_h
#define __kmClassifierUtils_h

//#include "itkListSample.h"
//#include "itkKdTreeGenerator.h"
//#include "itkWeightedCentroidKdTreeGenerator.h"
#include <map>

#include "itkSimplexMeshGeometry.h"

#include "kmProfileExtractor.h"
#include "kmProfileClassifier.h"

namespace km
{
	template< class TMesh, class TProfileExtractor>
	class ClassifierUtils
	{
	public:
		typedef typename km::ProfileClassifier ClassifierType;
		typedef typename TMesh MeshType;
		typedef typename MeshType::PointType PointType;
		typedef typename PointType::VectorType VectorType;

		//typedef typename Statistics::ListSample<PointType>                       SampleType;
		//typedef typename Statistics::WeightedCentroidKdTreeGenerator<SampleType> TreeGeneratorType;
		//typedef typename TreeGeneratorType::KdTreeType                           KdTreeType;
		//typedef typename KdTreeType::InstanceIdentifierVectorType                NeighborhoodIdentifierType;
		//typedef std::map<int, SampleType::Pointer>                               SampleMapType;
		//typedef std::map<int, TreeGeneratorType::Pointer>                        TreeGeneratorMapType;

		typedef typename TProfileExtractor ProfileExtractorType;

		typedef itk::SimplexMeshGeometry SimplexMeshGeometryType;

		typedef typename MeshType::GeometryMapType GeometryMapType;
		typedef typename GeometryMapType::Iterator GeometryMapIterator;
		
		ClassifierUtils()
		{
		}
		~ClassifierUtils()
		{
		}

		void SetProfileExtractor(ProfileExtractorType* extractor)
		{
			m_ProfileExtractor = extractor;
		}
		
		void SetRegionClassifier(ClassifierType* classifier)
		{
			m_RegionClassifier = classifier;
		}

		void SetBoundaryClassifier(ClassifierType* classifier)
		{
			m_BoundaryClassifier = classifier;
		}

		void deformByLiverClassification( typename MeshType* outputMesh, const typename MeshType* liverMesh, double searchStep = 1.5, unsigned int maxSearchPoints = 20)
		{
			typedef MeshType::PointsContainer PointsContainer;
			typedef MeshType::PointsContainerConstPointer PointsContainerConstPointer;
			typedef MeshType::PointsContainerPointer PointsContainerPointer;
			typedef MeshType::PointsContainerIterator PointsContainerIterator;
			typedef MeshType::PointDataContainer PointDataContainer;
			typedef MeshType::GeometryMapType GeometryMapType;
			typedef GeometryMapType::Pointer GeometryMapPointer;
			typedef GeometryMapType::Iterator GeometryMapIterator;

			km::assigneMesh<MeshType>(outputMesh, 0.0);

			GeometryMapIterator geoIt = liverMesh->GetGeometryData()->Begin();
			GeometryMapIterator geoItEnd = liverMesh->GetGeometryData()->End();

			itk::SimplexMeshGeometry *geodata;
			while (geoIt!=geoItEnd)
			{
				MeshType::PointIdentifier idx = geoIt.Index();
				geodata = geoIt.Value();

				PointType mpoint = liverMesh->GetPoint(idx);
				VectorType normal;
				normal.Set_vnl_vector(geodata->normal.Get_vnl_vector());

				double best_offset = 0.0;
				PointType closestNextRegionPoint = mpoint;
				bool foundNextRegion = true;
				this->FindNextRegionPoint(closestNextRegionPoint, mpoint, geodata, idx, 1.5);

				double distToNextRegion = itk::NumericTraits<double>::max();
				if (foundNextRegion){
					VectorType vec = closestNextRegionPoint - mpoint;
					distToNextRegion = dot_product(vec.GetVnlVector(), normal.GetVnlVector());
				}

				if (/*distToNextRegion<0*/true){
					best_offset = distToNextRegion;
				}else{
					PointType closestBoundaryPoint = mpoint;
					bool foundBoundary = this->FindFirstBoundaryPoint(closestBoundaryPoint, mpoint, geodata, idx, -1.5);

					double distToBoundary = itk::NumericTraits<double>::max();
					if (foundBoundary){
						VectorType vec = closestBoundaryPoint - mpoint;
						distToBoundary = dot_product(vec.GetVnlVector(), normal.GetVnlVector());
					}
					best_offset = std::abs(distToNextRegion)<std::abs(distToBoundary)?distToNextRegion:distToBoundary;
				}
				outputMesh->SetPointData(idx, best_offset);
				geoIt++;
			}

			//Remove noise point.
			km::smoothMeshData<MeshType>(outputMesh, 0);

			geoIt = liverMesh->GetGeometryData()->Begin();
			geoItEnd = liverMesh->GetGeometryData()->End();

			while (geoIt!=geoItEnd)
			{
				unsigned int idx = geoIt.Index();
				geodata = geoIt.Value();

				VectorType normal;
				normal.Set_vnl_vector(geodata->normal.Get_vnl_vector());

				double offsetVal = 0.0;
				outputMesh->GetPointData(idx, &offsetVal);

				PointType oldPos = liverMesh->GetPoint(idx);
				outputMesh->SetPoint(idx, oldPos + normal*offsetVal);

				geoIt++;
			}
		}

		//void UpdateRegionTree(const MeshType * mesh )
		//{
		//	if (m_RegionClassifier == NULL)
		//	{
		//		std::cout<<"Boundary classifier has not been set yet."<<std::endl;
		//		return;
		//	}

		//	if (m_ProfileExtractor == NULL)
		//	{
		//		std::cout<<"Profile extractor has not been set yet."<<std::endl;
		//		return;
		//	}

		//	//Clean samples
		//	for (SampleMapType::iterator it=m_RegionSampleMap.begin(); it!=m_RegionSampleMap.end(); ++it)
		//	{
		//		SampleType::Pointer samples = it->second;
		//		samples->Clear();
		//	}

		//	GeometryMapIterator geoIt = mesh->GetGeometryData()->Begin();
		//	GeometryMapIterator geoItEnd = mesh->GetGeometryData()->End();

		//	SimplexMeshGeometryType *geodata;
		//	while (geoIt!=geoItEnd)
		//	{
		//		MeshType::PointIdentifier idx = geoIt.Index();
		//		PointType curPos = mesh->GetPoint(idx);
		//		geodata = geoIt.Value();

		//		//int clusterlabel = this->m_RegionClassifier->getClusterLabel(idx);
		//		int clusterlabel = 0;
		//		SampleType::Pointer samples = m_RegionSampleMap[clusterlabel];

		//		PointType nextRegionPoint = curPos;
		//		this->FindNextRegionPoint(nextRegionPoint, curPos, geodata, idx, 1.5);

		//		samples->PushBack(nextRegionPoint);

		//		geoIt++;
		//	}

		//	//Update tree.
		//	for (TreeGeneratorMapType::iterator it=m_RegionTreeGeneratorMap.begin(); it!=m_RegionTreeGeneratorMap.end(); ++it)
		//	{
		//		int clusterlabel = it->first;

		//		this->m_RegionTreeGeneratorMap[clusterlabel]->SetSample( this->m_RegionSampleMap[clusterlabel] );
		//		this->m_RegionTreeGeneratorMap[clusterlabel]->SetBucketSize( 4 );
		//		this->m_RegionTreeGeneratorMap[clusterlabel]->Update();
		//	}
		//}

		//void UpdateBoundaryTree(const MeshType * mesh )
		//{
		//	if (m_BoundaryClassifier == NULL)
		//	{
		//		std::cout<<"Boundary classifier has not been set yet."<<std::endl;
		//		return;
		//	}

		//	if (m_ProfileExtractor == NULL)
		//	{
		//		std::cout<<"Profile extractor has not been set yet."<<std::endl;
		//		return;
		//	}

		//	//Clean samples
		//	for (SampleMapType::iterator it=m_BoundarySampleMap.begin(); it!=m_BoundarySampleMap.end(); ++it)
		//	{
		//		SampleType::Pointer samples = it->second;
		//		samples->Clear();
		//	}

		//	GeometryMapIterator geoIt = mesh->GetGeometryData()->Begin();
		//	GeometryMapIterator geoItEnd = mesh->GetGeometryData()->End();

		//	SimplexMeshGeometryType *geodata;
		//	while (geoIt!=geoItEnd)
		//	{
		//		MeshType::PointIdentifier idx = geoIt.Index();
		//		PointType curPos = mesh->GetPoint(idx);
		//		geodata = geoIt.Value();

		//		int clusterlabel = this->m_BoundaryClassifier->getClusterLabel(idx);
		//		SampleType::Pointer samples = m_BoundarySampleMap[clusterlabel];

		//		PointType boudaryPoint = curPos;
		//		bool found;

		//		//Towards outside.
		//		found = this->FindFirstBoundaryPoint(boudaryPoint, curPos, geodata, idx, 1.5);
		//		if (found)
		//		{
		//			samples->PushBack(boudaryPoint);
		//		}

		//		//Towards inside;
		//		found = this->FindFirstBoundaryPoint(boudaryPoint, curPos, geodata, idx, -1.5);
		//		if (found)
		//		{
		//			samples->PushBack(boudaryPoint);
		//		}

		//		geoIt++;
		//	}

		//	//Update tree.
		//	for (TreeGeneratorMapType::iterator it=m_BoundaryTreeGeneratorMap.begin(); it!=m_BoundaryTreeGeneratorMap.end(); ++it)
		//	{
		//		int clusterlabel = it->first;

		//		this->m_BoundaryTreeGeneratorMap[clusterlabel]->SetSample( this->m_BoundarySampleMap[clusterlabel] );
		//		this->m_BoundaryTreeGeneratorMap[clusterlabel]->SetBucketSize( 4 );
		//		this->m_BoundaryTreeGeneratorMap[clusterlabel]->Update();
		//	}
		//}

		//bool FindClosestNextRegionPoint(PointType & pt_next_region, PointType & pt_cur, int pointId = 0)
		//{
		//	//int clusterlabel = this->m_RegionClassifier->getClusterLabel(pointId);
		//	int clusterlabel = 0;

		//	TreeGeneratorMapType::iterator itr = m_RegionTreeGeneratorMap.find(clusterlabel);
		//	if (itr == m_RegionTreeGeneratorMap.end())
		//	{
		//		itr = m_RegionTreeGeneratorMap.find(0);
		//		if (itr == m_RegionTreeGeneratorMap.end())
		//		{
		//			std::cout<<"No tree found."<<std::endl;
		//			return false;
		//		}
		//	}
		//	TreeGeneratorType::Pointer treeGenerator = itr->second;
		//	NeighborhoodIdentifierType neighbors;

		//	try
		//	{
		//		//treeGenerator->GetOutput()->Search( pt_cur, 1u, neighbors );
		//		neighbors.push_back(pointId);
		//	}
		//	catch (...)
		//	{
		//		//std::cout<<"Cannot find any closest boundary point for cluster: "<<clusterlabel<<std::endl;
		//		return false;
		//	}

		//	if (neighbors.size() > 0)
		//	{
		//		pt_next_region = treeGenerator->GetOutput()->GetMeasurementVector( neighbors[0] );

		//		return true;
		//	}
		//	else
		//	{
		//		return false;
		//	}
		//}

		//bool FindClosestBoundaryPoint(PointType & pt_boundary, PointType & pt_cur, int pointId = 0)
		//{
		//	int clusterlabel = this->m_BoundaryClassifier->getClusterLabel(pointId);

		//	TreeGeneratorMapType::iterator itr = m_BoundaryTreeGeneratorMap.find(clusterlabel);
		//	if (itr == m_BoundaryTreeGeneratorMap.end())
		//	{
		//		itr = m_BoundaryTreeGeneratorMap.find(0);
		//		if (itr == m_BoundaryTreeGeneratorMap.end())
		//		{
		//			std::cout<<"No tree found."<<std::endl;
		//			return false;
		//		}
		//	}
		//	TreeGeneratorType::Pointer treeGenerator = itr->second;
		//	NeighborhoodIdentifierType neighbors;

		//	try
		//	{
		//		treeGenerator->GetOutput()->Search( pt_cur, 1u, neighbors );
		//	}
		//	catch (...)
		//	{
		//		//std::cout<<"Cannot find any closest boundary point for cluster: "<<clusterlabel<<std::endl;
		//		return false;
		//	}
		//	
		//	if (neighbors.size() > 0)
		//	{
		//		pt_boundary = treeGenerator->GetOutput()->GetMeasurementVector( neighbors[0] );

		//		return true;
		//	}
		//	else
		//	{
		//		return false;
		//	}
		//}

		void FindNextRegionPoint(PointType & pt, PointType & curPos, SimplexMeshGeometryType* geodata, int pointId, double searchStep)
		{
			bool found = false;

			double searchedDist = 0;
			unsigned int maxSearchPoints = 20;

			VectorType normal;
			normal.Set_vnl_vector(geodata->normal.Get_vnl_vector());

			if (m_ShapeNormalMap.size()>0)
			{
				VectorType shapeNormal = m_ShapeNormalMap[pointId];
				double dotproduct = dot_product(shapeNormal.GetVnlVector(), normal.GetVnlVector());
				if (dotproduct>0){
					shapeNormal.Normalize();
					normal += shapeNormal;
					normal.Normalize();
				}else if (dotproduct<0){
					shapeNormal.Normalize();
					normal -= shapeNormal;
					normal.Normalize();
				}
			}

			double firstDirection = 1.0;
			unsigned int searchedPoints = 1;
			double next_offset = 0;
			while(searchedPoints<=maxSearchPoints)
			{
				PointType pttest = curPos + normal*next_offset;
				double tmpDirection = -1.0;
				if (m_ProfileExtractor->isInsideBuffer(pttest))
				{
					std::vector<FLOATTYPE> feature;
					m_ProfileExtractor->extractFeatureSet(feature, m_RegionClassifier->profileCategory, geodata, pttest);

					tmpDirection = 2.0*m_RegionClassifier->classify(feature,pointId) - 1.0;
				}
				next_offset += tmpDirection * searchStep;
				if (searchedPoints == 1){
					firstDirection = tmpDirection>0?1.0:-1.0;
				}else if (tmpDirection*firstDirection<0){
					break;//Change direction. Break from here.
				}
				searchedPoints++;
			}
			pt = curPos + normal*next_offset;
		}

		bool FindFirstBoundaryPoint(PointType & pt, PointType & curPos, SimplexMeshGeometryType* geodata, int pointId, double searchStep)
		{
			bool found = false;

			double searchedDist = 0;
			unsigned int maxSearchPoints = 20;

			VectorType normal;
			normal.Set_vnl_vector(geodata->normal.Get_vnl_vector());

			double boundaryProbabilityMax = 0.0;
			unsigned int searchedPoints = 1;
			while(searchedPoints<maxSearchPoints)
			{
				PointType pttest = curPos + normal*searchedDist;
				if (this->m_ProfileExtractor->isInsideBuffer(pttest))
				{
					std::vector<FLOATTYPE> feature;
					m_ProfileExtractor->extractFeatureSet(feature, this->m_BoundaryClassifier->profileCategory, geodata, pttest);

					double boundaryProbability = this->m_BoundaryClassifier->classify(feature,pointId);
					if (boundaryProbability>0.8 && boundaryProbability>boundaryProbabilityMax){
						found = true;
						pt = pttest;
						boundaryProbabilityMax = boundaryProbability;
						//break;
					}
				}
				else
				{
					break;
				}
				searchedDist += searchStep;
				searchedPoints++;
			}
			return found;
		}

		void
			updateShapeNormals(const MeshType* meshBeforeShapeFitting, const MeshType* meshAfterShapeFitting)
		{
			if (meshBeforeShapeFitting->GetNumberOfPoints() != meshAfterShapeFitting->GetNumberOfPoints()){
				std::cout<<"Error. Numbers of points are not equal."<<std::endl;
			}

			for (int i=0;i<meshBeforeShapeFitting->GetNumberOfPoints();i++){
				m_ShapeNormalMap[i] = meshAfterShapeFitting->GetPoint(i) - meshBeforeShapeFitting->GetPoint(i);
			}
		}
		
	private:
		typename ProfileExtractorType* m_ProfileExtractor;

		typename ClassifierType* m_BoundaryClassifier;
		//typename SampleMapType m_BoundarySampleMap;
		//typename TreeGeneratorMapType m_BoundaryTreeGeneratorMap;

		typename ClassifierType* m_RegionClassifier;
		//typename SampleMapType m_RegionSampleMap;
		//typename TreeGeneratorMapType m_RegionTreeGeneratorMap;

		std::map<int, VectorType> m_ShapeNormalMap;
	};
}
#endif