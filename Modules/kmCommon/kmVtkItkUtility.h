#ifndef __kmVtkItkUtility_h
#define __kmVtkItkUtility_h

#include <iostream>
#include <fstream>
#include <string>
#include <set>

#include "itkSimplexMesh.h"
#include "vtkPolyData.h""
#include "itkLineCell.h"
#include "itkTriangleCell.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"

#include "itkIdentityTransform.h"
#include "itkTransformMeshFilter.h"

#include "itkAdaptSimplexMesh3DFilter.h"

#include "kmUtility.h"

using namespace std;

namespace km
{
	typedef itk::CovariantVector<unsigned int, 3> GeometryVectorType;
	typedef itk::Image<GeometryVectorType, 1>     GeometryImageType;

	template<class MeshType>
	void
		copySimplexMesh(
		const typename MeshType* meshFrom, 
		typename MeshType::Pointer & meshTo, 
		bool copyGeoFlag = true, 
		bool copyPointsFlag = false,
		bool copyCellFlag = false,
		bool computeGeometry = false)
	{
		const unsigned int numberOfPoints = meshFrom->GetNumberOfPoints();

		if (copyGeoFlag)
		{
			typedef typename MeshType::GeometryMapType           GeometryMapType;
			typedef typename MeshType::GeometryMapPointer        GeometryMapPointer;
			typedef typename MeshType::GeometryMapConstIterator  GeometryMapConstIterator;

			GeometryMapPointer inputGeometryData = meshFrom->GetGeometryData();
			GeometryMapPointer outputGeometryData = meshTo->GetGeometryData();

			outputGeometryData->Reserve( numberOfPoints );

			GeometryMapConstIterator inputGeometryItr = inputGeometryData->Begin();

			for( unsigned int pointId = 0; pointId < numberOfPoints; pointId++ )
			{
				SimplexMeshGeometry * outputGeometryDataItem = new SimplexMeshGeometry;
				SimplexMeshGeometry * inputGeometryDataItem = inputGeometryItr.Value();
				outputGeometryDataItem->CopyFrom( *inputGeometryDataItem );
				outputGeometryData->InsertElement( pointId, outputGeometryDataItem );
				++inputGeometryItr;
			}
		}

		meshTo->SetLastCellId( meshFrom->GetLastCellId() );
		meshTo->BuildCellLinks();

		if(computeGeometry)
		{
			ComputeGeometry<MeshType>( meshTo );
		}
	}

	template< typename TInputMesh, typename TOutputMesh >
	void
		CopyInputMeshToOutputMeshPoints(const typename TInputMesh *inputMesh, typename TOutputMesh *outputMesh)
	{
		typedef typename TOutputMesh::PointsContainer OutputPointsContainer;
		typedef typename TInputMesh::PointsContainer  InputPointsContainer;

		typename OutputPointsContainer::Pointer outputPoints = OutputPointsContainer::New();
		const InputPointsContainer *inputPoints = inputMesh->GetPoints();

		if ( inputPoints )
		{
			outputPoints->Reserve( inputPoints->Size() );

			typename InputPointsContainer::ConstIterator inputItr = inputPoints->Begin();
			typename InputPointsContainer::ConstIterator inputEnd = inputPoints->End();

			typename OutputPointsContainer::Iterator outputItr = outputPoints->Begin();

			while ( inputItr != inputEnd )
			{
				outputItr.Value() = inputItr.Value();
				++inputItr;
				++outputItr;
			}

			outputMesh->SetPoints(outputPoints);
		}
	}

	template< typename TInputMesh, typename TOutputMesh >
	void
		CopyInputMeshToOutputMeshCells(const typename TInputMesh *inputMesh, typename TOutputMesh *outputMesh)
	{
		typedef typename TOutputMesh::CellsContainer  OutputCellsContainer;
		typedef typename TInputMesh::CellsContainer   InputCellsContainer;
		typedef typename TOutputMesh::CellAutoPointer CellAutoPointer;

		outputMesh->SetCellsAllocationMethod(TOutputMesh::CellsAllocatedDynamicallyCellByCell);

		typename OutputCellsContainer::Pointer outputCells = OutputCellsContainer::New();
		const InputCellsContainer *inputCells = inputMesh->GetCells();

		if ( inputCells )
		{
			outputCells->Reserve( inputCells->Size() );

			typename InputCellsContainer::ConstIterator inputItr = inputCells->Begin();
			typename InputCellsContainer::ConstIterator inputEnd = inputCells->End();

			typename OutputCellsContainer::Iterator outputItr = outputCells->Begin();

			CellAutoPointer clone;

			while ( inputItr != inputEnd )
			{
				//      outputItr.Value() = inputItr.Value();
				// BUG: FIXME: Here we are copying a pointer, which is a mistake. What we
				// should do is to clone the cell.
				inputItr.Value()->MakeCopy(clone);
				outputItr.Value() = clone.ReleaseOwnership();

				++inputItr;
				++outputItr;
			}

			outputMesh->SetCells(outputCells);
		}
	}

	template< typename TInputMesh, typename TOutputMesh >
	void
		CopyInputMeshToOutputMeshCellLinks(const typename TInputMesh *inputMesh, typename TOutputMesh *outputMesh)
	{
		typedef typename TOutputMesh::CellLinksContainer OutputCellLinksContainer;
		typedef typename TInputMesh::CellLinksContainer  InputCellLinksContainer;

		typename OutputCellLinksContainer::Pointer outputCellLinks = OutputCellLinksContainer::New();
		const InputCellLinksContainer *inputCellLinks = inputMesh->GetCellLinks();

		if ( inputCellLinks )
		{
			outputCellLinks->Reserve( inputCellLinks->Size() );

			typename InputCellLinksContainer::ConstIterator inputItr = inputCellLinks->Begin();
			typename InputCellLinksContainer::ConstIterator inputEnd = inputCellLinks->End();

			typename OutputCellLinksContainer::Iterator outputItr = outputCellLinks->Begin();

			while ( inputItr != inputEnd )
			{
				outputItr.Value() = inputItr.Value();
				++inputItr;
				++outputItr;
			}

			outputMesh->SetCellLinks(outputCellLinks);
		}
	}

	template<class InputMeshType, class OutputMeshType>
	typename OutputMeshType::Pointer
		castMesh( const InputMeshType * inputmesh )
	{
		typedef itk::IdentityTransform<double, InputMeshType::PointDimension> TransformType;
		TransformType::Pointer transform = TransformType::New();

		typedef itk::TransformMeshFilter<InputMeshType, OutputMeshType, TransformType> TransformMeshFilterType;
		TransformMeshFilterType::Pointer transformFilter = TransformMeshFilterType::New();

		transformFilter->SetInput( inputmesh );
		transformFilter->SetTransform( transform );
		transformFilter->Update();

		return transformFilter->GetOutput();
	}

	template<class MeshType>
	GeometryImageType::Pointer
		generateGeoImage( const typename MeshType * mesh  )
	{
#define DEBUG_LOCAL(X) std::cout<<"[DEBUG LOCAL]"<<X<<std::endl;
		unsigned int numberOfPoints = mesh->GetNumberOfPoints();

		GeometryImageType::Pointer geoImage = GeometryImageType::New();
		GeometryImageType::IndexType start;
		GeometryImageType::SizeType  size;
		start.Fill( 0 );
		size.Fill( numberOfPoints );

		GeometryImageType::RegionType region;
		region.SetSize( size );
		region.SetIndex( start );
		geoImage->SetRegions( region );
		geoImage->Allocate();
		GeometryImageType::PixelType zeroVector;
		zeroVector.Fill( 0 );
		geoImage->FillBuffer( zeroVector );

		if( mesh->GetGeometryData().IsNull() )
		{
			KM_DEBUG_ERROR( "Geometry data is NULL in mesh!" );
			return geoImage;
		}

		typedef itk::ImageRegionIterator<GeometryImageType> IteratorType;
		IteratorType it( geoImage, geoImage->GetLargestPossibleRegion() );
		it.GoToBegin();
		unsigned int ppp = 0;
		while (!it.IsAtEnd())
		{
			MeshType::IndexArray neibours = mesh->GetNeighbors( ppp );

			if (neibours.Size()==3)
			{
				GeometryImageType::PixelType pix;
				for (int k=0;k<3;k++)
				{
					pix[k] = neibours[k];
				}
				it.Set( pix );
			}
			else
			{
				std::cout<<"No neibours!"<<std::endl;
			}

			++it;
			++ppp;
		}

		return geoImage;
#undef DEBUG_LOCAL(X)
	}

	template<class MeshType>
	void
		writeSimplexMeshGeometryData(const char* filename, typename MeshType::Pointer mesh)
	{		
		GeometryImageType::Pointer geoImage = km::generateGeoImage<MeshType>( mesh );		

		km::writeImage<GeometryImageType>( filename, geoImage );
	}

	template<class MeshType>
	void
		readSimplexMeshGeometryData(const char* filename, typename MeshType::Pointer mesh)
	{
		unsigned int numberOfPoints = mesh->GetNumberOfPoints();

		mesh->GetGeometryData()->Reserve( numberOfPoints );

		GeometryImageType::Pointer geoImage = km::readImage<GeometryImageType>( filename );
		typedef itk::ImageRegionConstIterator<GeometryImageType> IteratorType;
		IteratorType it( geoImage, geoImage->GetLargestPossibleRegion() );
		it.GoToBegin();
		unsigned int p = 0;
		while (!it.IsAtEnd())
		{
			GeometryImageType::PixelType pix = it.Get();
			SimplexMeshGeometry *data = new itk::SimplexMeshGeometry();
			data->pos = mesh->GetPoint( p );
			for (int k=0;k<3;k++)
			{	
				data->neighborIndices[k] = pix[k];
			}
			mesh->SetGeometryData( p, data );

			++it;
			++p;
		}

		mesh->BuildCellLinks();

		ComputeGeometry<MeshType>( mesh );
	}

	template<class MeshType>
	void
		loadSimplexMeshGeometryData(const GeometryImageType* geoImage, typename MeshType::Pointer mesh)
	{
		unsigned int numberOfPoints = mesh->GetNumberOfPoints();

		mesh->GetGeometryData()->Reserve( numberOfPoints );

		typedef itk::ImageRegionConstIterator<GeometryImageType> IteratorType;
		IteratorType it( geoImage, geoImage->GetLargestPossibleRegion() );
		it.GoToBegin();
		unsigned int p = 0;
		while (!it.IsAtEnd())
		{
			GeometryImageType::PixelType pix = it.Get();
			SimplexMeshGeometry *data = new itk::SimplexMeshGeometry();
			data->pos = mesh->GetPoint( p );
			for (int k=0;k<3;k++)
			{	
				data->neighborIndices[k] = pix[k];
			}
			mesh->SetGeometryData( p, data );

			++it;
			++p;
		}

		mesh->BuildCellLinks();

		ComputeGeometry<MeshType>( mesh );
	}

	/************************************************************************/
	// Conver ITK Simplex Mesh to VTK PolyData
	/************************************************************************/
	template<typename ITKMeshType>
	vtkSmartPointer<vtkPolyData>
		mesh2PolyData( const typename ITKMeshType* mesh )
	{
		typedef typename ITKMeshType::CellType CellType;
		typedef itk::Point<vtkFloatingPointType, 3> ItkPoint;

		typedef typename CellType::PointIdIterator PointIdIterator;

		typedef typename ITKMeshType::CellsContainer::ConstIterator CellIterator;
		typedef typename ITKMeshType::PointsContainer::ConstIterator PointIterator;
		//typedef typename ITKMeshType::LinesContainer

		vtkSmartPointer<vtkPolyData> newPolyData = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

		PointIterator pntIterator = mesh->GetPoints()->Begin();
		PointIterator pntItEnd = mesh->GetPoints()->End();

		for (int i = 0; pntIterator != pntItEnd; ++i, ++pntIterator)
		{
			ItkPoint pnt = pntIterator.Value();
			points->InsertPoint(i, pnt[0], pnt[1], pnt[2]);
		}
		newPolyData->SetPoints(points);
		//points->Delete();

		vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
		vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

		CellIterator cellIt = mesh->GetCells()->Begin();
		CellIterator cellItEnd = mesh->GetCells()->End();

		unsigned long numberoflines = 0;
		unsigned long numberoftriangles = 0;
		unsigned long numberofpolygons = 0;
		unsigned long numberofunknown = 0;

		for (int it = 0; cellIt != cellItEnd; ++it, ++cellIt)
		{
			CellType * cellptr = cellIt.Value();

			PointIdIterator pntIdIter = cellptr->PointIdsBegin();
			PointIdIterator pntIdEnd = cellptr->PointIdsEnd();
			vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();

			for (; pntIdIter != pntIdEnd; ++pntIdIter)
			{
				pts->InsertNextId( *pntIdIter );
			}

			switch( cellptr->GetType() )
			{
			case CellType::LINE_CELL:
				{
					//lines->InsertNextCell(pts);
					//numberoflines++;
					break;
				}
			case CellType::TRIANGLE_CELL:
				{
					cells->InsertNextCell(pts);
					numberoftriangles++;
					break;
				}
			case CellType::POLYGON_CELL:
				{
					cells->InsertNextCell(pts);
					numberofpolygons++;
					break;
				}
			default:
				{
					cells->InsertNextCell(pts);
					numberofunknown++;
					break;
				}
			}

			//cells->InsertNextCell(pts);
		}

		KM_DEBUG_PRINT( "numberoflines: ", numberoflines );
		KM_DEBUG_PRINT( "numberofpolygons: ", numberofpolygons );
		KM_DEBUG_PRINT( "numberoftriangles: ", numberoftriangles );
		KM_DEBUG_PRINT( "numberofunknown: ", numberofunknown );

		if (numberofpolygons>0 || numberoftriangles>0)
		{
			newPolyData->SetPolys(cells);
		}

		if(numberoflines>0)
		{
			newPolyData->SetLines(lines);
		}

		//triangle->Delete();

		return newPolyData;
	}

	template<typename InputMeshType, typename OutputMeshType>
	void
		copyPointDataFromMeshToMesh( typename InputMeshType* inputMesh, typename OutputMeshType* outputMesh )
	{
		unsigned int numberOfPointsInSource = inputMesh->GetNumberOfPoints();
		unsigned int numberOfPointsInTarget = outputMesh->GetNumberOfPoints();

		if (numberOfPointsInTarget!=numberOfPointsInSource)
		{
			std::cerr<<"Number of points is different between source mesh and target mesh !"<<std::endl;
			return;
		}

		typedef typename InputMeshType::PointDataContainerPointer         InputPointDataContainerPointer;
		typedef typename InputMeshType::PointDataContainer::ConstIterator InputPointDataContainerConstIterator;

		typedef typename OutputMeshType::PointDataContainerPointer         OutputPointDataContainerPointer;

		InputPointDataContainerPointer allInputPointData = inputMesh->GetPointData();
		InputPointDataContainerConstIterator inputPointDataIt = allInputPointData->Begin();
		InputPointDataContainerConstIterator inputPointDataItEnd = allInputPointData->End();

		OutputPointDataContainerPointer allOutputPointData = outputMesh->GetPointData();

		while( inputPointDataIt != inputPointDataItEnd )
		{
			//allOutputPoints->InsertElement( inputPointIt.Index(), inputPointIt.Value() );

			allOutputPointData->InsertElement( inputPointDataIt.Index(), inputPointDataIt.Value() );

			inputPointDataIt++;
		}

		// 			outputMesh->SetLastCellId( inputMesh->GetLastCellId() );
		// 			outputMesh->BuildCellLinks();
	}

	template<typename InputMeshType, typename OutputMeshType>
	void
		copyPointsFromMeshToMesh( typename InputMeshType* inputMesh, typename OutputMeshType* outputMesh )
	{
		unsigned int numberOfPointsInSource = inputMesh->GetNumberOfPoints();
		unsigned int numberOfPointsInTarget = outputMesh->GetNumberOfPoints();

		if (numberOfPointsInTarget!=numberOfPointsInSource)
		{
			std::cerr<<"Number of points is different between source mesh and target mesh !"<<std::endl;
			return;
		}

		typedef typename InputMeshType::PointsContainerPointer         InputPointsContainerPointer;
		typedef typename InputMeshType::PointsContainer::ConstIterator InputPointsContainerConstIterator;

		typedef typename OutputMeshType::PointsContainerPointer         OutputPointsContainerPointer;
		typedef typename OutputMeshType::PointsContainer::Iterator      OutputPointsContainerIterator;

		InputPointsContainerPointer allInputPoints = inputMesh->GetPoints();
		InputPointsContainerConstIterator inputPointIt = allInputPoints->Begin();
		InputPointsContainerConstIterator inputPointItEnd = allInputPoints->End();

		//OutputPointsContainerPointer allOutputPoints = outputMesh->GetPoints();

		while( inputPointIt != inputPointItEnd )
		{
			//allOutputPoints->InsertElement( inputPointIt.Index(), inputPointIt.Value() );

			outputMesh->SetPoint( inputPointIt.Index(), inputPointIt.Value() );

			inputPointIt++;
		}

		outputMesh->SetLastCellId( inputMesh->GetLastCellId() );
		outputMesh->BuildCellLinks();
	}

	template<typename MeshType>
	void
		copyPointsFromPolydataToMesh(vtkPolyData* m_PolyData, typename MeshType::Pointer mesh)
	{
		const unsigned int numberOfPointsInPolydata = m_PolyData->GetNumberOfPoints();
		const unsigned int numberOfPointsInMesh = mesh->GetNumberOfPoints();

		if (numberOfPointsInPolydata!=numberOfPointsInMesh)
		{
			std::cout<<"number of points in mesh is not equal to polydata!"<<std::endl;
			return;
		}

		vtkPoints * vtkpoints =  m_PolyData->GetPoints();

		for(unsigned int p =0; p < numberOfPointsInPolydata; p++)
		{
			double* apoint = vtkpoints->GetPoint( p );
			MeshType::PointType pt;
			for(unsigned int i=0;i<3; i++)
			{
				pt[i] = static_cast<vtkFloatingPointType>(apoint[i]);
			}
			mesh->SetPoint( p, pt);
		}
	}

	/************************************************************************/
	// Conver VTK PolyData to ITK Mesh
	/************************************************************************/
	template<typename MeshType>
	typename MeshType::Pointer
		polyData2Mesh( vtkPolyData* m_PolyData, typename MeshType::Pointer refmesh = NULL)
	{
		typedef typename MeshType::PixelType           PixelType;
		typedef typename MeshType::CellType            CellType;
		typedef itk::VertexCell< CellType >            VertexCellType;
		typedef itk::LineCell< CellType >              LineCellType;
		typedef itk::TriangleCell< CellType >          TriangleCellType;
		typedef itk::PolygonCell< CellType >           PolygonCellType;

		typedef typename MeshType::CellAutoPointer     CellAutoPointer;
		typedef typename MeshType::PointIdentifier     PointIdentifier;

		// Create a new mesh
		MeshType::Pointer m_itkMesh = MeshType::New();

		//
		//Points and point data
		//
		const unsigned int numberOfPointsInPolydata = m_PolyData->GetNumberOfPoints();
		unsigned int numberOfPointsInMesh = numberOfPointsInPolydata;
		vtkPoints * vtkpoints =  m_PolyData->GetPoints();
		KM_DEBUG_PRINT( "numberOfPoints: ", numberOfPointsInMesh );
		m_itkMesh->GetPoints()->Reserve( numberOfPointsInMesh );
		m_itkMesh->GetPointData()->Reserve( numberOfPointsInMesh );
		m_itkMesh->GetGeometryData()->Reserve( numberOfPointsInMesh );
		for(unsigned int p =0; p < numberOfPointsInPolydata; p++)
		{
			double* apoint = vtkpoints->GetPoint( p );
			MeshType::PointType pt;
			for(unsigned int i=0;i<3; i++)
			{
				pt[i] = static_cast<vtkFloatingPointType>(apoint[i]);
			}
			m_itkMesh->SetPoint( p, pt);
			m_itkMesh->SetPointData( p, static_cast<PixelType>( 0 ) );
			//m_itkMesh->SetGeometryData( p, new itk::SimplexMeshGeometry() );
			if (NULL!=refmesh)
			{
				SimplexMeshGeometry *refdata = refmesh->GetGeometryData()->GetElement(p);
				SimplexMeshGeometry *data = new itk::SimplexMeshGeometry();
				data->CopyFrom( *refdata );
				m_itkMesh->SetGeometryData( p, data );
			}
		}

		//KM_DEBUG_INFO( "Insert points done!" );


		//
		//Polygons and cell data
		//
		const unsigned int numberOfPolysInPolydata = m_PolyData->GetNumberOfPolys();
		const unsigned int numberOfLinesInPolydata = m_PolyData->GetNumberOfLines();

		//Count the number of all cells which will be inserted into mesh.
		unsigned int numberOfCellsInMesh = numberOfPolysInPolydata;
		vtkCellArray * lines = m_PolyData->GetLines();
		lines->InitTraversal();
		vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
		while( lines->GetNextCell( cellPoints ) )
		{
			numberOfCellsInMesh += ( cellPoints->GetNumberOfIds() - 1);
		}
		KM_DEBUG_PRINT( "numberOfCells: ", numberOfCellsInMesh );
		m_itkMesh->GetCells()->Reserve( numberOfCellsInMesh );
		m_itkMesh->GetCellData()->Reserve( numberOfCellsInMesh );

		//std::cout<<numberOfCellsInMesh<<std::endl;

		//Start to insert cells into mesh.
		vtkCellArray * polygons = m_PolyData->GetPolys();
		polygons->InitTraversal();

		vtkIdType cellId = 0;
		while( polygons->GetNextCell( cellPoints ) )
		{
			MeshType::CellAutoPointer c;

			PolygonCellType *newCell = new PolygonCellType;
			for ( unsigned int jj = 0; jj < cellPoints->GetNumberOfIds(); jj++ )
			{
				newCell->SetPointId( jj, cellPoints->GetId(jj) );
			}
			c.TakeOwnership(newCell);
			//m_itkMesh->SetCell( cellId, c );
			m_itkMesh->AddFace( c );
			m_itkMesh->SetCellData( cellId, static_cast<PixelType>( 0 ) );

			cellId++;
		}

		lines->InitTraversal();

		while( lines->GetNextCell( cellPoints ) )
		{
			MeshType::CellAutoPointer c;

			unsigned int numberOfIdsInLine = cellPoints->GetNumberOfIds();

			for (int i=0;i<numberOfIdsInLine-1;i++)
			{
				//LineCellType *  newCell = new LineCellType;
				vtkIdType pt1 = cellPoints->GetId( i );
				vtkIdType pt2 = cellPoints->GetId( i+1 );
				//newCell->SetPointId(0, pt1);
				//newCell->SetPointId(1, pt2);
				//c.TakeOwnership( newCell );

				//m_itkMesh->SetCell( cellId, c );
				m_itkMesh->AddEdge( pt1, pt2 );
				m_itkMesh->SetCellData( cellId, static_cast<PixelType>( 0 ) );

				cellId++;
			}
		}

		//std::cout<<cellId<<std::endl;

		m_itkMesh->BuildCellLinks();

		return m_itkMesh;
	}

	template<class ItkImageType>
	vtkSmartPointer<vtkImageData>
		castItkImageToVtkImage( const typename ItkImageType* itkImage )
	{
		typedef itk::ImageToVTKImageFilter<ItkImageType> FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput( itkImage );
		filter->Update();

		return filter->GetOutput();
	}

	template<class MeshType>
	typename MeshType::Pointer
		decimateMesh( typename MeshType* inputMesh , unsigned int numberOfPoints = 3000)
	{
		vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();

		polydata = km::mesh2PolyData<MeshType>(inputMesh);
		polydata = km::smoothPolyData( polydata, 30 );
		polydata = km::decimatePolydata(polydata, static_cast<double>(numberOfPoints)/inputMesh->GetNumberOfPoints());

		MeshType::Pointer outputMesh = polyData2Mesh<MeshType>( polydata );

		return outputMesh;
	}

	//If computeAll=false, then only normal will be computed. Otherwise all geometry data like mean curve, sphere radius etc. will be computed.
	template<class MeshType>
	void
		ComputeGeometry(typename MeshType* inputMesh, bool computeAll = false)
	{
		typedef typename MeshType::GeometryMapType      GeometryMapType;
		typedef typename GeometryMapType::Pointer       GeometryMapPointer;
		typedef typename GeometryMapType::Iterator      GeometryMapIterator;

		typename MeshType::PointsContainerPointer allpoints = inputMesh->GetPoints();
		typename MeshType::GeometryMapPointer     allgeodata   = inputMesh->GetGeometryData();
		typename MeshType::PointDataContainerPointer allpointdata = inputMesh->GetPointData();

		if (allpointdata->Size() <= 0)
		{
			allpointdata->Reserve( allpoints->Size() );
		}

		typedef typename SimplexMeshGeometry::PointType PointType;
		typedef typename PointType::VectorType          VectorType;
		typedef CovariantVector< typename VectorType::ValueType, 3 > CovariantVectorType;

		PointType           Foot;
		CovariantVectorType normal;
		CovariantVectorType z;
		VectorType          tmp;

		GeometryMapType::Iterator dataIt = allgeodata->Begin();

		SimplexMeshGeometry *data;

		while ( dataIt != allgeodata->End() )
		{
			//      idx = dataIt.Index();
			data = dataIt.Value();

			data->neighbors[0] = allpoints->GetElement(data->neighborIndices[0]);
			data->neighbors[1] = allpoints->GetElement(data->neighborIndices[1]);
			data->neighbors[2] = allpoints->GetElement(data->neighborIndices[2]);

			// compute normal
			normal.Fill(0.0);

			z.Set_vnl_vector( itk_cross_3d( ( data->neighbors[1] - data->neighbors[0] ).Get_vnl_vector(),
				( data->neighbors[2] - data->neighbors[0] ).Get_vnl_vector() ) );
			z.Normalize();
			normal += z;

			// copy normal
			data->normal = normal;
			data->pos = allpoints->GetElement( dataIt.Index() );

			// compute the simplex angle
			data->ComputeGeometry();

			if ( computeAll )
			{
				tmp = data->neighbors[0] - data->pos;

				double D = 1.0 / ( 2 * data->sphereRadius ); /* */

				double tmpNormalProd = dot_product( tmp.GetVnlVector(), data->normal.GetVnlVector() );

				double sinphi =  2 *data->circleRadius *D *vnl_math_sgn(tmpNormalProd);
				double phi = vcl_asin(sinphi);

				data->phi = phi;
				data->meanCurvature = vcl_abs(sinphi / data->circleRadius);
				tmp = data->pos - data->neighbors[0];

				//compute the foot of p projection of p onto the triangle spanned by its
				// neighbors
				double distance = -tmpNormalProd;
				tmp.SetVnlVector( ( data->pos ).GetVnlVector() - distance * normal.GetVnlVector() );
				Foot.Fill(0.0);
				Foot += tmp;

				data->distance = ( ( data->circleCenter ) - Foot ).GetNorm();

				{
					PointType a, b, c;

					a = data->neighbors[0];
					b = data->neighbors[1];
					c = data->neighbors[2];

					VectorType n, na, nb, nc;
					n.SetVnlVector( itk_cross_3d( ( b - a ).GetVnlVector(), ( c - a ).GetVnlVector() ) );
					na.SetVnlVector( itk_cross_3d( ( c - b ).GetVnlVector(), ( Foot - b ).GetVnlVector() ) );
					nb.SetVnlVector( itk_cross_3d( ( a - c ).GetVnlVector(), ( Foot - c ).GetVnlVector() ) );
					nc.SetVnlVector( itk_cross_3d( ( b - a ).GetVnlVector(), ( Foot - a ).GetVnlVector() ) );

					PointType eps;
					eps[0] = dot_product( n.GetVnlVector(), na.GetVnlVector() ) / n.GetSquaredNorm();
					eps[1] = dot_product( n.GetVnlVector(), nb.GetVnlVector() ) / n.GetSquaredNorm();
					eps[2] = dot_product( n.GetVnlVector(), nc.GetVnlVector() ) / n.GetSquaredNorm();

					data->eps = eps;
				}
			}

			allpointdata->InsertElement( dataIt->Index(), data->phi );

			dataIt.Value() = data;

			dataIt++;
		}
	}

	template<class MeshType>
	void
		adaptMesh(typename MeshType* inputMesh, double m_Gamma, unsigned int m_Iterations)
	{
		typedef itk::AdaptSimplexMesh3DFilter<MeshType, MeshType> AdaptSimplexMesh3DFilterType;
		AdaptSimplexMesh3DFilterType::Pointer adaptorFilter = AdaptSimplexMesh3DFilterType::New();
		adaptorFilter->SetInput( inputMesh );
		adaptorFilter->SetGamma( m_Gamma );
		adaptorFilter->SetTangentFactor( 0.8 );
		adaptorFilter->SetNormalFactor( 0.2 );
		adaptorFilter->SetRigidity( 1 );
		adaptorFilter->SetIterations( m_Iterations );
		adaptorFilter->Update(); //This is just for initialization.
	}

	template<class MeshType>
	void
		smoothMesh(typename MeshType* inputMesh, double factor, unsigned int m_Iterations)
	{
		typedef itk::AdaptSimplexMesh3DFilter<MeshType, MeshType> AdaptSimplexMesh3DFilterType;
		AdaptSimplexMesh3DFilterType::Pointer adaptorFilter = AdaptSimplexMesh3DFilterType::New();
		adaptorFilter->SetInput( inputMesh );
		adaptorFilter->SetGamma( 0.0 );
		adaptorFilter->SetTangentFactor( 0.0 );
		adaptorFilter->SetNormalFactor( factor );
		adaptorFilter->SetRigidity( 1 );
		adaptorFilter->SetIterations( m_Iterations );
		adaptorFilter->Update(); //This is just for initialization.
	}

	template<class MeshType>
	void
		assigneMesh(typename MeshType::Pointer mesh, typename MeshType::PixelType value)
	{
		typedef MeshType::PointsContainer::Iterator      PointsIterator;

		mesh->GetPointData()->Reserve( mesh->GetNumberOfPoints() );

		PointsIterator pointItr = mesh->GetPoints()->Begin();
		PointsIterator pointEnd = mesh->GetPoints()->End();

		while ( pointItr != pointEnd )
		{
			mesh->SetPointData( pointItr.Index(), value );	

			++pointItr;
		}
	}

	template<class MeshType>
	void
		smoothMeshData(typename MeshType::Pointer mesh, unsigned int radius)
	{
		MeshType::PointDataContainerPointer allpointdata = mesh->GetPointData();
		if (allpointdata.IsNull() || allpointdata->Size() <= 0)
		{
			return;
		}

		typedef MeshType::PointDataContainer::ElementIdentifier ElementIdentifier;
		typedef MeshType::PointsContainer::Iterator PointsIterator;
		PointsIterator pointItr = mesh->GetPoints()->Begin();
		PointsIterator pointEnd = mesh->GetPoints()->End();
		while ( pointItr != pointEnd )
		{
			MeshType::PointIdentifier idx = pointItr.Index();
			double datasum = allpointdata->GetElement( idx );

			MeshType::NeighborListType * neighborlist = mesh->GetNeighbors(idx, radius);
			for (int i=0;i<neighborlist->size();i++)
			{
				datasum += allpointdata->GetElement( (*neighborlist)[i] );
			}
			datasum /= (1+neighborlist->size());

			allpointdata->InsertElement( idx, datasum );

			++pointItr;
		}
	}

	//Copy mesh data from sourceMesh to targetMesh.
	template<class TInputMesh, class TOutputMesh>
	void
		copyMeshData(const TInputMesh* sourceMesh, typename TOutputMesh* targetMesh)
	{
		typedef TOutputMesh::PointDataContainer OutputPointDataContainer;
		typedef TInputMesh::PointDataContainer  InputPointDataContainer;

		OutputPointDataContainer::Pointer outputPointData = targetMesh->GetPointData();
		const InputPointDataContainer *inputPointData = sourceMesh->GetPointData();

		if ( inputPointData )
		{
			outputPointData->Reserve( inputPointData->Size() );

			InputPointDataContainer::ConstIterator inputItr = inputPointData->Begin();
			InputPointDataContainer::ConstIterator inputEnd = inputPointData->End();

			OutputPointDataContainer::Iterator outputItr = outputPointData->Begin();

			while ( inputItr != inputEnd )
			{
				outputItr.Value() = inputItr.Value();
				++inputItr;
				++outputItr;
			}
		}
	}
}

#endif