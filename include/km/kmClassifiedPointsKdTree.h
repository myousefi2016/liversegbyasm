#ifndef __kmClassifiedPointsKdTree_h
#define __kmClassifiedPointsKdTree_h

#include "itkListSample.h"
#include "itkKdTreeGenerator.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include <map>

#include "kmGlobal.h"
#include "kmProfileClassifier.h"
#include "kmProfileExtractor.h"

namespace km
{
	template< typename TMesh, typename TProfileExtractor>
	class ClassifiedPointsKdTree
	{
	public:
		typedef typename km::ProfileClassifier ClassifierType;
		typedef typename TMesh MeshTpye;
		typedef typename MeshType::PointType PointType;

		typedef typename Statistics::ListSample<PointType>                       SampleType;
		typedef typename Statistics::WeightedCentroidKdTreeGenerator<SampleType> TreeGeneratorType;
		typedef typename TreeGeneratorType::KdTreeType                           KdTreeType;
		typedef typename KdTreeType::InstanceIdentifierVectorType                NeighborhoodIdentifierType;
		typedef std::map<int, SampleType::Pointer>                               SampleMapType;
		typedef std::map<int, TreeGeneratorType::Pointer>                        TreeGeneratorMapType;

		typedef typename TProfileExtractor ProfileExtractorType;

		typedef itk::SimplexMeshGeometry SimplexMeshGeometryType;
		typedef itk::SimplexMeshGeometry::VectorType VectorType;

		typedef MeshType::GeometryMapType GeometryMapType;
		typedef GeometryMapType::Iterator GeometryMapIterator;
		
		ClassifiedPointsKdTree()
		{
		}
		~ClassifiedPointsKdTree()
		{
		}

		void SetProfileExtractor(ProfileExtractorType* extractor)
		{
			m_ProfileExtractor = extractor;
		}
		
		void SetBoundaryClassifier(ClassifierType* classifier)
		{
			m_BoundaryClassifier = classifier;
			for(int i=0;i<classifier->getNumberOfClusters();i++)
			{
				SampleType::Pointer samples = SampleType::New();
				samples->SetMeasurementVectorSize( PointType::PointDimension );
				
				TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
				treeGenerator->SetSample( samples );
				treeGenerator->SetBucketSize( 4 );
				
				m_BoundarySampleMap[i] = samples;
				m_BoundaryTreeGeneratorMap[i] = treeGenerator;
			}
		}
		
		void UpdateBoundaryTree(const MeshType * mesh )
		{
			if (m_BoundaryClassifier == NULL)
			{
				std::cout<<"Boundary classifier has not been set yet."<<std::endl;
				return;
			}

			if (m_ProfileExtractor == NULL)
			{
				std::cout<<"Profile extractor has not been set yet."<<std::endl;
				return;
			}

			//Clean samples
			for (SampleMapType::iterator it=m_BoundarySampleMap.begin(); it!=m_BoundarySampleMap.end(); ++it)
			{
				SampleType::Pointer samples = it->second;
				samples->Clear();
			}

			GeometryMapIterator geoIt = mesh->GetGeometryData()->Begin();
			GeometryMapIterator geoItEnd = mesh->GetGeometryData()->End();

			SimplexMeshGeometryType *geodata;
			while (geoIt!=geoItEnd)
			{
				MeshType::PointIdentifier idx = geoIt.Index();
				PointType curPos = mesh->GetPoint(idx);
				geodata = geoIt.Value();

				int clusterlabel = this->m_BoundaryClassifier->getClusterLabel(idx);
				SampleType::Pointer samples = m_BoundarySampleMap[clusterlabel];

				PointType boudaryPoint = curPos;
				bool found;
				
				//Towards outside.
				found = this->FindFirstBoundaryPoint(boudaryPoint, curPos, geodata, idx, 1.0);
				if (found)
				{
					samples->PushBack(boudaryPoint);
				}

				//Towards inside;
				found = this->FindFirstBoundaryPoint(boudaryPoint, curPos, geodata, idx, -1.0);
				if (found)
				{
					samples->PushBack(boudaryPoint);
				}

				geoIt++;
			}
			
			//Update tree.
			for (TreeGeneratorMapType::iterator it=m_BoundaryTreeGeneratorMap.begin(); it!=m_BoundaryTreeGeneratorMap.end(); ++it)
			{
				int clusterlabel = it->first;
				
				this->m_BoundaryTreeGeneratorMap[clusterlabel]->SetSample( this->m_BoundarySampleMap[clusterlabel] );
				this->m_BoundaryTreeGeneratorMap[clusterlabel]->SetBucketSize( 4 );
				this->m_BoundaryTreeGeneratorMap[clusterlabel]->Update();
			}
		}
		
		bool FindClosetBoundaryPoint(PointType & pt_boundary, PointType & pt_cur, int pointId = 0)
		{
			int clusterlabel = this->m_BoundaryClassifier->getClusterLabel(pointId);

			TreeGeneratorMapType::iterator itr = m_BoundaryTreeGeneratorMap.find(clusterlabel);
			if (itr == m_BoundaryTreeGeneratorMap.end())
			{
				itr = m_BoundaryTreeGeneratorMap.find(0);
				if (itr == m_BoundaryTreeGeneratorMap.end())
				{
					std::cout<<"No tree found."<<std::endl;
					return false;
				}
			}
			TreeGeneratorType::Pointer treeGenerator = itr->second;
			NeighborhoodIdentifierType neighbors;

			try
			{
				treeGenerator->GetOutput()->Search( pt_cur, 1u, neighbors );
			}
			catch (...)
			{
				//std::cout<<"Cannot find any closest boundary point for cluster: "<<clusterlabel<<std::endl;
				return false;
			}
			
			if (neighbors.size() > 0)
			{
				pt_boundary = treeGenerator->GetOutput()->GetMeasurementVector( neighbors[0] );

				return true;
			}
			else
			{
				return false;
			}
		}
		
	private:
		typename ProfileExtractorType* m_ProfileExtractor;
		typename ClassifierType* m_BoundaryClassifier;
		typename SampleMapType m_BoundarySampleMap;
		typename TreeGeneratorMapType m_BoundaryTreeGeneratorMap;

		bool FindFirstBoundaryPoint(PointType & pt, PointType & curPos, SimplexMeshGeometryType* geodata, int pointId, double searchStep)
		{
			bool found = false;

			double maxDist = g_varianceMap[pointId];
			double searchedDist = 0;
			unsigned int maxSearchPoints = 20;

			VectorType normal;
			normal.Set_vnl_vector(geodata->normal.Get_vnl_vector());

			unsigned int searchedPoints = 1;
			while(searchedPoints<maxSearchPoints && std::abs(searchedDist)<maxDist)
			{
				PointType pttest = curPos + normal*searchedDist;
				if (this->m_ProfileExtractor->isInsideBuffer(pttest))
				{
					std::vector<FLOATTYPE> feature;
					m_ProfileExtractor->extractFeatureSet(feature, this->m_BoundaryClassifier->profile_category, geodata, pttest);

					double boundaryProbability = this->m_BoundaryClassifier->classify(feature,pointId);
					if (boundaryProbability>0.5)
					{
						found = true;
						pt = pttest;
						break;
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
		
	};
}
#endif