#ifndef __KM_UTILITY_H
#define __KM_UTILITY_H

/*********************************************************/
/*********************************************************/
/*********************************************************/
// 一些常用API
//
// API列表：
//          读入图像，写入图像
//          图像重采样
//          图像平滑，图像切割
//          ... ...
/*********************************************************/
/*********************************************************/
/*********************************************************/

#define kmStaticImageMacro(T) typedef typename T::PixelType PixelType;         \
	typedef typename T::IndexType IndexType;         \
	typedef typename T::PointType PointType;         \
	typedef typename T::SizeType  SizeType;          \
	typedef typename T::SpacingType SpacingType;     \
	typedef typename T::RegionType RegionType;       \
	typedef typename T::DirectionType DirectionType; \

#define KM_DEBUG_ERROR(X) std::cout<<"[DEBUG ERROR]"<<(X)<<std::endl;
#define KM_DEBUG_INFO(X) std::cout<<"[DEBUG INFO] "<<(X)<<std::endl;
#define KM_DEBUG_PRINT(X,Y) std::cout<<"[DEBUG PRINT] "<<X<<":"<<(Y)<<std::endl;
#define KM_DEBUG_PRINT_VALUE(X) std::cout<<"[DEBUG PRINT] "<<"X:"<<(X)<<std::endl;

#define KM_PRINT_EXCEPTION(E) std::cout<<"[EXCEPTION] "<<(E)<<std::endl;

#define KM_ASSERT(X) if(!(X)) std::cout<<"[ASSERT FAIL] "<<std::endl; 

//#include <itkDirectory.h>
#include <string>
#include <fstream>
#include <itkImage.h>
#include <itkSmartPointer.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIterator.h>

#include <itkMesh.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include "itkLineCell.h"
#include "itkTriangleCell.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"

#include <itkVTKPolyDataWriter.h>
#include <itkVTKPolyDataReader.h>

#include <itkMeshFileReader.h>
#include <itkMeshFileWriter.h>

#include <itkResampleImageFilter.h>
#include <itkPadImageFilter.h>
#include <itkConstantBoundaryCondition.h>
#include <itkCastImageFilter.h>
#include <itkThresholdImageFilter.h>
#include "itkTransformFileWriter.h"

#include <itkTranslationTransform.h>

#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryMorphologicalClosingImageFilter.h>
#include <itkBinaryMorphologicalOpeningImageFilter.h>
#include <itkBinaryMedianImageFilter.h>

#include "itkBinaryMask3DMeshSource.h"
#include "itkTriangleMeshToBinaryImageFilter.h"
#include "itkTriangleMeshToSimplexMeshFilter.h"
#include "itkSimplexMeshToTriangleMeshFilter2.h"
#include "itkSimplexMesh.h"

#include <itkWarpMeshFilter.h>
#include <itkTransformMeshFilter.h>

#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkMultiplyImageFilter.h>

#include <vtkSmoothPolyDataFilter.h>

#include <vtkPolyDataNormals.h>
#include <vtkDecimatePro.h>

#include "itkTransformFileReader.h"
#include "itkTransformFactoryBase.h"

#include <itkMinMaxCurvatureFlowImageFilter.h>
#include <itkCurvatureAnisotropicDiffusionImageFilter.h>
#include <itkGradientAnisotropicDiffusionImageFilter.h>
#include <itkMedianImageFilter.h>
#include "itkDiscreteGaussianImageFilter.h"
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <itkGradientImageFilter.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkGradientRecursiveGaussianImageFilter.h>
#include <itkGradientVectorFlowImageFilter.h>
#include <itkSigmoidImageFilter.h>
#include <itkChangeInformationImageFilter.h>
#include "itkMinimumMaximumImageCalculator.h"
#include "itkRescaleIntensityImageFilter.h"
#include <itkHistogramMatchingImageFilter.h>
#include <itkIntensityWindowingImageFilter.h>

#include "itkImageSliceIteratorWithIndex.h"

#include <itkShiftScaleImageFilter.h>

#include "itkAddImageFilter.h"
#include "itkAbsImageFilter.h"
#include <itkPowImageFilter.h>

#include "itkBinaryContourImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include <itkImageDuplicator.h>
#include <itkWatershedImageFilter.h>

#include <vtkImageMarchingCubes.h>

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include "itkRegionOfInterestImageFilter.h"

#include "itkNormalizedCorrelationImageFilter.h"
#include "itkSimplexMeshAdaptTopologyFilter.h"

#include <sstream>

#include "kmAlgorithm.h"

using namespace itk;

namespace km
{
	struct HistogramNode 
	{
		HistogramNode()
		{
			lower = 0;
			upper = 0;
			percentage = 0;
			count = 0;
		};
		double lower;
		double upper;
		double percentage;
		unsigned count;
	};

	struct Histogram
	{
		double statisticMin;
		double statisticMax;
		unsigned int binNum;
		double binWidth;
		HistogramNode * nodes;

		//Histogram(double _staMin, double _staMax, unsigned int _binNum)
		//{
		//	statisticMin = _staMin;
		//	statisticMax = _staMax;
		//	binNum = _binNum;

		//	KM_ASSERT(statisticMax>=statisticMin);
		//	KM_ASSERT(binNum>0);

		//	binWidth = (statisticMax-statisticMin)/binNum;

		//	nodes = new HistogramNode[binNum];
		//}

		Histogram(double _staMin, double _staMax, double _binWidth)
		{
			statisticMin = _staMin;
			statisticMax = _staMax;
			binWidth = _binWidth;

			KM_ASSERT(statisticMax>=statisticMin);
			KM_ASSERT(binWidth>0);
			
			binNum = static_cast<unsigned int>(std::ceil((statisticMax-statisticMin)/binWidth)) + 1;

			nodes = new HistogramNode[binNum];
		}

		~Histogram()
		{
			delete[] nodes;
		}

		void clear()
		{
			for (int i=0;i<binNum;i++)
			{
				this->nodes[i].count = 0;
				this->nodes[i].percentage = 0.0;
			}
		}

		void calcThresholdOfRatio( double ratio, double & thresholdLower, double & thresholdUpper )
		{
			if (ratio<=0)
			{
				return;
			}

			double totalPercentage = 0;
			double max=this->nodes[0].percentage;

			int binNumM;

			for(int i=0;i<this->binNum;i++)
			{
				if (this->nodes[i].percentage>max)
				{
					max=this->nodes[i].percentage;
					binNumM = i;
				}
			}

			double sum=max;
			int i=binNumM;
			int j=binNumM;
			int flag=0;

			bool foundflag = false;

			while(i>=0&&j<this->binNum)
			{
				if (flag==0)
				{
					i--;
					flag=1;
					sum+=this->nodes[i].percentage;
				}
				else
				{
					j++; 
					flag=0; 
					sum+=this->nodes[j].percentage;
				}

				if (sum>ratio)
				{
					thresholdLower = this->nodes[i].lower;
					thresholdUpper = this->nodes[j].upper; 

					foundflag = true;

					break;
				}
			}

			if (!foundflag)
			{
				if(i<0) i+=1;
				if(j>binNumM) j-=1;

				thresholdLower = this->nodes[i].lower;
				thresholdUpper = this->nodes[j].upper;
			}
		}

		double calcRatioOfThreshold( double thresholdLower, double thresholdUpper )
		{
			thresholdLower = std::max(thresholdLower, this->statisticMin);
			thresholdUpper = std::min(thresholdUpper, this->statisticMax);

			if(thresholdUpper < thresholdLower)
			{
				return 0.0;
			}

			unsigned int binIdxLower, binIdxUpper;
			binIdxLower = static_cast<unsigned int>(std::floor( (thresholdLower - this->statisticMin)/binWidth ));
			binIdxUpper = static_cast<unsigned int>(std::ceil(  (thresholdUpper - this->statisticMin)/binWidth ));

			double ratio = 0;
			for (unsigned int i=binIdxLower;i<binIdxUpper;i++)
			{
				ratio += this->nodes[i].percentage;
			}

			return ratio;
		}

		void print()
		{

		}
	};

	/************************************************************************/
	/*                                                                      */
	/************************************************************************/
	template<typename ImageType>
	void 
		calculateHistogram( const ImageType* image, Histogram* histogram )
	{
		/*double binSize = (staUpper-staLower+1)/binNum;

		if (!histogram)
		{
			histogram = new HistogramNode[binNum];
		}*/

		histogram->clear();

		unsigned int totalCount = 0; 
		itk::ImageRegionConstIterator<ImageType> it( image, image->GetLargestPossibleRegion() );
		it.GoToBegin();
		while(!it.IsAtEnd())
		{
			double val = it.Get();
			if (val<histogram->statisticMin || val>histogram->statisticMax)
			{
				++it;
				continue;
			}
			else
			{
				unsigned int binIndex = static_cast<unsigned int>(std::floor( (val - histogram->statisticMin)/histogram->binWidth ));

				KM_ASSERT( binIndex >= 0 && binIndex < histogram->binNum );

				histogram->nodes[ binIndex ].count++;

				totalCount++;
				++it;
			}
		}

		double totalPercentage = 0;
		for(int i=0;i<histogram->binNum;i++)
		{
			histogram->nodes[i].lower = histogram->statisticMin+i*histogram->binWidth;
			histogram->nodes[i].upper = histogram->nodes[i].lower + histogram->binWidth - 1;
			histogram->nodes[i].percentage = static_cast<double>(histogram->nodes[i].count)/totalCount;

			totalPercentage += histogram->nodes[i].percentage;
		}
		KM_DEBUG_PRINT( "Total percentage: ", totalPercentage );
		KM_DEBUG_INFO( "Calculate histogram DONE!" );
	}

	/************************************************************************/
	/* Decorate file name                                                   */
	/************************************************************************/
	const std::string decorateFilename( const char* oldFilename, const char* extension, const char* addname  )
	{
		std::stringstream ss;

		std::string oldnamestr(oldFilename);
		ss << oldnamestr.substr(0, oldnamestr.find_last_of( '.' )) << addname << extension;

		return ss.str();
	}

	/************************************************************************/
	/* Get Data List                                                        */
	/************************************************************************/
	int getDataList( std::string dataListFileName, std::vector<std::string> &files )
	{
		std::ifstream  fin(dataListFileName.c_str(), std::ios::in);  
		char  filename[1024]={0};
		int num = 0;
		while(fin.getline(filename, sizeof(filename)))  
		{
			std::string str_filename( filename );

			//std::cout<<">>>>"<<filename<<std::endl;

			if(str_filename.size()==0)
			{
				continue;
			}
			
			if( str_filename.find('#') != std::string::npos )
			{
				continue;
			}
			
			if ( str_filename.find( "END" ) != std::string::npos )
			{
				break;
			}

			files.push_back(str_filename);
			num ++;
		}  
		fin.clear();  
		fin.close();

		return num;
	}

	template<class ImageType1, class ImageType2, class OutputImageType>
	typename OutputImageType::Pointer
		addImage(const typename ImageType1* image1,
		typename ImageType2* image2)
	{
		typedef itk::AddImageFilter<ImageType1, ImageType2, OutputImageType> AddImageFilterType;
		AddImageFilterType::Pointer adder = AddImageFilterType::New();
		adder->SetInput1( image1 );
		adder->SetInput2( image2 );
		adder->Update();

		return adder->GetOutput();
	}


	/************************************************************************/
	/* Read ITK Image                                                        */
	/************************************************************************/
	template<typename TImageType>
	typename TImageType::Pointer
		readImage(const char* filename)
	{
		typedef itk::ImageFileReader<TImageType> ImageFileReaderType;
		typedef typename ImageFileReaderType::Pointer ImageFileReaderPointer;
		ImageFileReaderPointer reader = ImageFileReaderType::New();
		reader->SetFileName( filename );
		try
		{
			reader->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return NULL; // Since the goal of the example is to catch the exception, we declare this a success.
		}
		return reader->GetOutput();
	}

	template<typename TImageType>
	typename TImageType::Pointer
		readImage(const std::string filename)
	{
		return readImage<TImageType>(filename.c_str());
	}

	/************************************************************************/
	/* Write ITK Image                                                      */
	/************************************************************************/
	template<typename TImageType>
	void
		writeImage(const char* filename, const typename TImageType* image)
	{
		typedef itk::ImageFileWriter<TImageType> ImageFileWriterType;
		typedef typename ImageFileWriterType::Pointer ImageFileWriterPointer;
		ImageFileWriterPointer writer = ImageFileWriterType::New();
		writer->SetFileName( filename );
		writer->SetInput( image );
		try
		{
			writer->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
		}
	}

	template<typename TImageType>
	void
		writeImage(const std::string filename, const typename TImageType* image)
	{
		writeImage<TImageType>(filename.c_str(), image);
	}

	template<typename TImageType>
	void
		writeImage(const std::string filenameprefix, 
								int index, 
								const char* extension, 
								const typename TImageType* image)
	{
		std::stringstream ss;
		ss  <<  filenameprefix   <<  "."  <<  index+1  <<  extension;

		km::writeImage<TImageType>( ss.str().c_str(), image );
	}

	/************************************************************************/
	/* Read Mesh                                                            */
	/************************************************************************/
	template<typename TMeshType>
	typename TMeshType::Pointer
		readMesh( const char* filename )
	{
		//typedef itk::VTKPolyDataReader<TMeshType>        VTKPolyDataReaderType;
		//VTKPolyDataReaderType::Pointer reader = VTKPolyDataReaderType::New();
		typedef itk::MeshFileReader<TMeshType> MeshReaderType;
		MeshReaderType::Pointer reader = MeshReaderType::New();
		reader->SetFileName( filename );
		try
		{
			reader->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return NULL;
		}

		return reader->GetOutput();
	}

	template<typename TMeshType>
	typename TMeshType::Pointer
		readMesh( const std::string filename )
	{
		return readMesh<TMeshType>(filename.c_str());
	}

	/************************************************************************/
	/* Write Mesh                                                           */
	/************************************************************************/
	template<typename TMeshType>
	void
		writeMesh( const char* filename, const typename TMeshType* mesh )
	{
		//std::cout<<"here!!!!!"<<std::endl;
		//typedef itk::VTKPolyDataWriter<TMeshType>       VTKPolyDataWriterType;
		//VTKPolyDataWriterType::Pointer writer = VTKPolyDataWriterType::New();
		typedef itk::MeshFileWriter<TMeshType> MeshWriterType;
		MeshWriterType::Pointer writer = MeshWriterType::New();
		writer->SetInput( mesh );
		writer->SetFileName( filename );
		try
		{
			writer->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
		}
		catch(...)
		{
			std::cerr << "Shit happens!"<<std::endl;
		}
	}

	template<typename TMeshType>
	void
		writeMesh( const std::string filename, const TMeshType* mesh )
	{
		writeMesh<TMeshType>(filename.c_str(), mesh);
	}

	/************************************************************************/
	/* Read VTKPolyData                                                     */
	/************************************************************************/
	vtkSmartPointer<vtkPolyData>
		readPolyData( const char* filename )
	{
		vtkSmartPointer<vtkPolyDataReader> inputPolyDataReader = vtkSmartPointer<vtkPolyDataReader>::New();
		inputPolyDataReader->SetFileName( filename );
		try
		{
			inputPolyDataReader->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return NULL;
		}

		return inputPolyDataReader->GetOutput();
	}

	vtkSmartPointer<vtkPolyData>
		readPolyData( const std::string filename )
	{
		return readPolyData(filename.c_str());
	}

	/************************************************************************/
	/* Write VTKPolyData                                                    */
	/************************************************************************/
	void
		writePolyData( const char* filename, vtkPolyData* polydata )
	{
		vtkSmartPointer<vtkPolyDataWriter> outputPolyDataWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
		outputPolyDataWriter->SetFileName( filename );
		outputPolyDataWriter->SetInput( polydata );
		try
		{
			outputPolyDataWriter->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
		}
	}

	void
		writePolyData( const std::string filename, vtkPolyData* polydata )
	{
		writePolyData(filename.c_str(), polydata);
	}

	template<typename TImageType>
	void
		writeImage(const char* dirname, 
		const char* filename, 
		const typename TImageType* image)
	{
		std::stringstream ss;
		ss  <<  dirname  <<  "/"  <<  filename;

		km::writeImage<TImageType>( ss.str().c_str(), image );
	}

	template<typename TImageType>
	void
		writeImage(const char* dirname, 
		const std::string filenameprefix, 
		unsigned int index, 
		const char* extension, 
		const typename TImageType* image)
	{
		std::stringstream ss;
		ss  <<  dirname  <<  "/"  <<  filenameprefix   <<  "."  <<  index+1  <<  extension;

		km::writeImage<TImageType>( ss.str().c_str(), image );
	}

	template<typename TTransformType>
	void
		writeTransform(const char* dirname, 
		const std::string filenameprefix, 
		unsigned int index, 
		const char* extension, 
		const typename TTransformType* transform)
	{
		std::stringstream ss;
		ss  <<  dirname  <<  "/"  <<  filenameprefix   <<  "."  <<  index+1  <<  extension;

		km::writeTransform< TTransformType >( ss.str().c_str(), transform );
	}

	void
		writePolyData(const char* dirname, 
		const std::string filenameprefix, 
		unsigned int index, 
		const char* extension, 
		vtkPolyData* polydata)
	{
		std::stringstream ss;
		ss  <<  dirname  <<  "/"  <<  filenameprefix   <<  "."  <<  index+1  <<  extension;

		km::writePolyData( ss.str().c_str(), polydata );
	}

	template<typename TMeshType>
	void
		writeMesh(const char* dirname, 
		const std::string filename, 
		const TMeshType* mesh)
	{
		std::stringstream ss;
		ss  <<  dirname  <<  "/"  <<  filename;

		km::writeMesh< TMeshType >( ss.str().c_str(), mesh );
	}

	template<typename TMeshType>
	void
		writeMesh(const char* dirname, 
		const std::string filenameprefix, 
		unsigned int index, 
		const char* extension, 
		const TMeshType* mesh)
	{
		std::stringstream ss;
		ss  <<  dirname  <<  "/"  <<  filenameprefix   <<  "."  <<  index+1  <<  extension;

		km::writeMesh< TMeshType >( ss.str().c_str(), mesh );
	}

	/************************************************************************/
	/* Decimate a vtkpolydata                                               */
	/************************************************************************/
	vtkSmartPointer<vtkPolyData>
		decimatePolydata( vtkSmartPointer<vtkPolyData> inputPolyData, double newRatio )
	{
		//std::cout<<"decimation ratio: "<<newRatio<<std::endl;
		vtkSmartPointer<vtkDecimatePro> decimate =
			vtkSmartPointer<vtkDecimatePro>::New();
#if VTK_MAJOR_VERSION <= 5
		decimate->SetInputConnection( inputPolyData->GetProducerPort() );
#else
		decimate->SetInputData( inputPolyData );
#endif
		decimate->SetFeatureAngle( 60 );
		decimate->SplittingOff();
		decimate->PreSplitMeshOff();
		decimate->BoundaryVertexDeletionOff();
		decimate->SetTargetReduction( static_cast<double> ( 1.0-newRatio )); //10% reduction (if there was 100 triangles, now there will be 90)
		decimate->PreserveTopologyOn();
		decimate->Update();

		return decimate->GetOutput();
	}

	/*
	//////////////////////////////////////////////////////////////////////////
	// Conver ITK Mesh to VTK PolyData
	//////////////////////////////////////////////////////////////////////////
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

			//if( cellptr->GetType() == CellType::LINE_CELL )
			//{
			//	LineType * line = static_cast<LineType *>( cell );
			//	std::cout << "dimension = " << line->GetDimension();
			//	std::cout << " # points = " << line->GetNumberOfPoints();
			//	std::cout << std::endl;
			//}

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
					lines->InsertNextCell(pts);
					numberoflines++;
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

	//////////////////////////////////////////////////////////////////////////
	// Conver VTK PolyData to ITK Mesh
	//////////////////////////////////////////////////////////////////////////
	template<typename MeshType>
	typename MeshType::Pointer
		polyData2Mesh( vtkPolyData* m_PolyData, int meshType = 0 ) //0:triangle 1:polygon 2:quadedge
	{
		typedef typename MeshType::CellType            CellType;
		typedef itk::VertexCell< CellType >            VertexCellType;
		typedef itk::LineCell< CellType >              LineCellType;
		typedef itk::TriangleCell< CellType >          TriangleCellType;
		typedef itk::PolygonCell< CellType >           PolygonCellType;
		typedef itk::TetrahedronCell< CellType >       TetrahedronCellType;
		typedef itk::HexahedronCell< CellType >        HexahedronCellType;
		typedef itk::QuadrilateralCell< CellType >     QuadrilateralCellType;
		typedef itk::QuadraticEdgeCell< CellType >     QuadraticEdgeCellType;
		typedef itk::QuadraticTriangleCell< CellType > QuadraticTriangleCellType;

		typedef typename MeshType::CellAutoPointer     CellAutoPointer;
		typedef typename MeshType::PointIdentifier     PointIdentifier;

		// Create a new mesh
		MeshType::Pointer m_itkMesh = MeshType::New();

		const unsigned int numberOfPointsInPolydata = m_PolyData->GetNumberOfPoints();
		const unsigned int numberOfCellsInPolydata = m_PolyData->GetNumberOfPolys();
		const unsigned int numberOfLinesInPolydata = m_PolyData->GetNumberOfLines();
		const unsigned int numberOfStipsInPolydata = m_PolyData->GetNumberOfStrips();

		KM_DEBUG_INFO( "Information in polydata: " );
		KM_DEBUG_PRINT( "numberOfPoints: ", numberOfPointsInPolydata );
		KM_DEBUG_PRINT( "numberOfCells: ", numberOfCellsInPolydata );
		KM_DEBUG_PRINT( "numberOfLines: ", numberOfLinesInPolydata );
		KM_DEBUG_PRINT( "numberOfStips: ", numberOfStipsInPolydata );

		unsigned int numberOfPointsInMesh = 0;
		unsigned int numberOfCellsInMesh = 0;
		unsigned int numberOfLinesInMesh = 0;
		//
		//Points
		//
		vtkPoints * vtkpoints =  m_PolyData->GetPoints();
		m_itkMesh->GetPoints()->Reserve(numberOfPointsInPolydata );
		for(unsigned int p =0; p < numberOfPointsInPolydata; p++)
		{
			double* apoint = vtkpoints->GetPoint( p );
			MeshType::PointType pt;
			for(unsigned int i=0;i<3; i++)
			{
				pt[i] = static_cast<vtkFloatingPointType>(apoint[i]);
			}
			m_itkMesh->SetPoint( p, pt);
			numberOfPointsInMesh++;
		}
		//KM_DEBUG_INFO( "Insert points done!" );

		vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
		m_itkMesh->GetCells()->Reserve( numberOfCellsInPolydata + numberOfLinesInPolydata );

		//
		//Lines
		//
		vtkCellArray * lines = m_PolyData->GetLines();

		int cellId = 0;
		lines->InitTraversal();

		vtkIdType numberOfCellPoints = 0;
		while( lines->GetNextCell( cellPoints ) )
		{
			MeshType::CellAutoPointer c;

			unsigned int numberOfIdsInLine = cellPoints->GetNumberOfIds();

			for (int i=0;i<numberOfIdsInLine-1;i++)
			{
				LineCellType *  newCell = new LineCellType;
				newCell->SetPointId(0, cellPoints->GetId( i ));
				newCell->SetPointId(1, cellPoints->GetId( i+1));
				c.TakeOwnership( newCell );
				m_itkMesh->SetCell( cellId++, c );
				numberOfLinesInMesh++;
			}
		} 

		//KM_DEBUG_INFO( "Insert lines done!" );

		//
		//Cells
		//
		vtkCellArray * polygons = m_PolyData->GetPolys();
		polygons->InitTraversal();

		numberOfCellPoints = 0;
		while( polygons->GetNextCell( cellPoints ) )
		{
			MeshType::CellAutoPointer c;

			switch(meshType)
			{
			case 0:
				{
					TriangleCellType *newCell = new TriangleCellType;
					for ( unsigned int jj = 0; jj < cellPoints->GetNumberOfIds(); jj++ )
					{
						newCell->SetPointId( jj, cellPoints->GetId(jj) );
					}
					c.TakeOwnership(newCell);
					break;
				}
			case 1:
				{
					PolygonCellType *newCell = new PolygonCellType;
					for ( unsigned int jj = 0; jj < cellPoints->GetNumberOfIds(); jj++ )
					{
						newCell->SetPointId( jj, cellPoints->GetId(jj) );
					}
					c.TakeOwnership(newCell);
					break;
				}
			case 2:
				{
					PolygonCellType *newCell = new PolygonCellType;
					for ( unsigned int jj = 0; jj < cellPoints->GetNumberOfIds(); jj++ )
					{
						newCell->SetPointId( jj, cellPoints->GetId(jj) );
					}
					c.TakeOwnership(newCell);
					break;
				}
			default:
				{
					std::cout<<"Unkown mesh type!"<<std::endl;
					TriangleCellType *newCell = new TriangleCellType;
					for ( unsigned int jj = 0; jj < cellPoints->GetNumberOfIds(); jj++ )
					{
						newCell->SetPointId( jj, cellPoints->GetId(jj) );
					}
					c.TakeOwnership(newCell);
					break;
				}
			}

			m_itkMesh->SetCell( cellId++, c );

			numberOfCellsInMesh++;
		} 

		//KM_DEBUG_INFO( "Insert polygons done!" );

		KM_DEBUG_INFO( "Information in mesh: " );
		KM_DEBUG_PRINT( "numberofpoints: ", numberOfPointsInMesh );
		KM_DEBUG_PRINT( "numberofcells: ", numberOfCellsInMesh );
		KM_DEBUG_PRINT( "numberoflines: ", numberOfLinesInMesh );

		return m_itkMesh;
	}

	//template<class ItkImageType>
	//vtkSmartPointer<vtkImageData>
	//	castItkImageToVtkImage( const typename ItkImageType* itkImage )
	//{
	//	typedef itk::ImageToVTKImageFilter<ItkImageType> FilterType;
	//	FilterType::Pointer filter = FilterType::New();
	//	filter->SetInput( itkImage );
	//	filter->Update();

	//	return filter->GetOutput();
	//} */

	//vtkSmartPointer<vtkPolyData>
	//	generatePolyDataFromVtkImage( vtkImageData* image )
	//{
	//	vtkSmartPointer<vtkImageMarchingCubes>  vmarchingcubes = vtkSmartPointer<vtkImageMarchingCubes>::New();
	//	vmarchingcubes->SetInput( image );
	//	vmarchingcubes->SetValue(0, 0.5);
	//	vmarchingcubes->ComputeScalarsOff();
	//	vmarchingcubes->ComputeNormalsOff();
	//	vmarchingcubes->ComputeGradientsOff();
	//	vmarchingcubes->SetInputMemoryLimit(1000);
	//	vmarchingcubes->Update();

	//	return vmarchingcubes->GetOutput();
	//}

	/************************************************************************/
	/* Image Casting                                                        */
	/************************************************************************/
	template<typename TImageType1, typename TImageType2>
	typename TImageType2::Pointer
		castImage( const typename TImageType1* inputImage )
	{
		typedef itk::CastImageFilter<TImageType1, TImageType2> CastImageFilterType;
		typename CastImageFilterType::Pointer caster = CastImageFilterType::New();
		caster->SetInput( inputImage );
		caster->Update();

		return caster->GetOutput();
	}

	template<typename TImageType>
	typename TImageType::PointType
		getCenter(const typename TImageType* inputImage)
	{
		TImageType::RegionType region = inputImage->GetLargestPossibleRegion();
		TImageType::IndexType  origin = region.GetIndex();
		TImageType::SizeType   size   = region.GetSize();
		TImageType::IndexType  centerIndex;
		for(int i=0;i<TImageType::ImageDimension;i++)
		{
			centerIndex[i] = origin[i] + size[i]/2;
		}
		TImageType::PointType  centerPoint;
		inputImage->TransformIndexToPhysicalPoint(centerIndex, centerPoint);

		return centerPoint;
	}

	/************************************************************************/
	/* Get Centroid                                                         */
	/************************************************************************/
	template<typename TImageType>
	void
		getCentroid(const typename TImageType* inputImage, 
		typename TImageType::PointType &center_pt,
		typename TImageType::IndexType &center_idx)
	{
		kmStaticImageMacro(TImageType);
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);

		center_idx.Fill(0.0);
		int n = 0;
		RegionType requestedRegion = inputImage->GetLargestPossibleRegion();
		itk::ImageRegionConstIterator<TImageType> it(inputImage, requestedRegion);
		int count = 0;
		while( !it.IsAtEnd() )
		{
			PixelType val = it.Get();
			IndexType ind = it.GetIndex();
			if( val != 0 )
			{
				n++;
				for(int d = 0; d < Dimension; d++)
				{
					center_idx[d] += ind[d];
				}
			}
			++it;
			++count;
		}

		for(int d = 0; d < Dimension; d++)
		{
			center_idx[d] /= n;
		}

		inputImage->TransformIndexToPhysicalPoint(center_idx, center_pt); 
	}

	template<typename TImageType>
	typename TImageType::PointType
		getCentroid(const typename TImageType* inputImage)
	{
		TImageType::PointType center_pt;
		TImageType::IndexType center_idx;
		getCentroid<TImageType>( inputImage, center_pt, center_idx );

		return center_pt;
	}

	/************************************************************************/
	template<typename TImageType>
	typename TImageType::PointType
		getCentroid(const typename TImageType* inputImage, double lowThreshold, double highThreshold)
	{
		return getCentroid(inputImage, inputImage->GetLargestPossibleRegion(), lowThreshold, highThreshold);
	}

	/************************************************************************/
	template<typename TImageType>
	typename TImageType::PointType
		getCentroid(const typename TImageType* inputImage, const typename TImageType::RegionType& searchRegion, double lowThreshold, double highThreshold)
	{
		kmStaticImageMacro(TImageType);
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);

		IndexType center_idx;
		center_idx.Fill(0.0);
		int n = 0;
		itk::ImageRegionConstIterator<TImageType> it(inputImage, searchRegion);
		int count = 0;
		while( !it.IsAtEnd() )
		{
			PixelType val = it.Get();
			IndexType ind = it.GetIndex();
			if( val>=lowThreshold && val<=highThreshold )
			{
				n++;
				for(int d = 0; d < Dimension; d++)
				{
					center_idx[d] += ind[d];
				}
			}
			++it;
			++count;
		}

		for(int d = 0; d < Dimension; d++)
		{
			center_idx[d] /= n;
		}

		PointType center_pt;
		inputImage->TransformIndexToPhysicalPoint(center_idx, center_pt); 
		return center_pt;
	}

	template<typename TImageType>
	int
		countPixels(const typename TImageType* inputImage, double lowThreshold, double highThreshold)
	{
		return countPixels<TImageType>( inputImage, inputImage->GetLargestPossibleRegion(), lowThreshold, highThreshold );
	}

	template<typename TImageType>
	int
		countPixels(const typename TImageType* inputImage, const typename TImageType::RegionType & searchRegion, double lowThreshold, double highThreshold)
	{
		kmStaticImageMacro(TImageType);
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);

		TImageType::RegionType croppedRegion(searchRegion);
		croppedRegion.Crop( inputImage->GetLargestPossibleRegion() );

		itk::ImageRegionConstIterator<TImageType> it(inputImage, croppedRegion);
		int count = 0;
		while( !it.IsAtEnd() )
		{
			PixelType val = it.Get();
			
			if (val>=lowThreshold && val<=highThreshold)
			{
				count++;
			}

			++it;
		}

		return count;
	}

	void
		calculatePolydataCentroid(vtkPolyData* inputPolydata, itk::FixedArray< double, 3 > & centroid)
	{
		centroid.Fill(0);
		const unsigned int numberOfPoints = inputPolydata->GetNumberOfPoints();
		vtkPoints * vtkpoints =  inputPolydata->GetPoints();

		for(unsigned int p =0; p < numberOfPoints; p++)
		{
			double* apoint = vtkpoints->GetPoint( p );
			for(unsigned int i=0;i<3; i++)
			{
				centroid[i] += apoint[i];
			}
		}
		for(unsigned int i=0;i<3; i++)
		{
			centroid[i] /= numberOfPoints;
		}
		return;
	}

	template<class MeshType>
	typename MeshType::PointType
		getMeshCentroid(const MeshType* pointset)
	{
		typedef MeshType::PointType PointType;
		PointType centroid;
		centroid.Fill(0);

		unsigned int Dimension = MeshType::PointDimension;

		unsigned long numberOfPoints = pointset->GetNumberOfPoints();
		typedef MeshType::PointsContainer::ConstIterator PointsIterator;

		PointsIterator pointItr = pointset->GetPoints()->Begin();
		PointsIterator pointEnd = pointset->GetPoints()->End();

		while(pointItr!=pointEnd)
		{
			PointType pt = pointItr.Value();

			for(unsigned int i=0;i<Dimension; i++)
			{
				centroid[i] += pt[i];
			}

			pointItr++;
		}

		for(unsigned int i=0;i<Dimension; i++)
		{
			centroid[i] /= numberOfPoints;
		}
		return centroid;
	}

	/************************************************************************/
	/* Resample Image                                                       */
	/************************************************************************/
	template<typename TImageType>
	typename TImageType::Pointer
		resampleImage( 
		const typename TImageType* inputImage, 
		const typename TImageType::SpacingType spacing, 
		double defaultValue = 0,
		int interpolatorType = 0)
	{
		kmStaticImageMacro(TImageType);
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);

		SizeType inputSize = inputImage->GetLargestPossibleRegion().GetSize();
		SpacingType inputSpacing = inputImage->GetSpacing();

		itk::Vector<double> inputLenth;
		for(int i = 0;i<Dimension; i++)
		{
			inputLenth[i] = ((float)inputSize[i])*inputSpacing[i];
		}

		SizeType outputSize;
		SpacingType outputSpacing;
		for(int i = 0;i<Dimension; i++)
		{
			outputSize[i] = inputLenth[i]/spacing[i];
			outputSpacing[i] = spacing[i];
		}

		typedef itk::ResampleImageFilter<TImageType, TImageType> ResampleImageFilterType;
		typedef typename ResampleImageFilterType::Pointer        ResampleImageFilterPointer;
		ResampleImageFilterPointer resampler = ResampleImageFilterType::New();
		resampler->SetInput( inputImage );
		resampler->SetOutputParametersFromImage( inputImage );
		resampler->SetOutputOrigin( inputImage->GetOrigin() );
		resampler->SetSize( outputSize );
		resampler->SetOutputSpacing( outputSpacing );
		resampler->SetDefaultPixelValue( defaultValue );

		typedef itk::LinearInterpolateImageFunction<TImageType> LinearInterpolateType;
		typedef itk::NearestNeighborInterpolateImageFunction<TImageType> NearestNeighborInterpolateType;
		switch (interpolatorType)
		{
		case 0:
			{
				LinearInterpolateType::Pointer interpolator = LinearInterpolateType::New();
				resampler->SetInterpolator( interpolator );
				break;
			}
		case 1:
			{
				NearestNeighborInterpolateType::Pointer interpolator2 = NearestNeighborInterpolateType::New();
				resampler->SetInterpolator( interpolator2 );
				break;
			}

		default:
			{
				std::cout<<"Unknown interpolator type!"<<std::endl;
				break;
			}
		}

		resampler->Update();

		return resampler->GetOutput();
	}

	/************************************************************************/
	/* Resample Image                                                       */
	/************************************************************************/
	template<typename TImageType>
	typename TImageType::Pointer
		resampleImage( const typename TImageType* inputImage, const typename TImageType::SpacingType spacing, const typename TImageType::SizeType size, double defaultValue )
	{
		typedef itk::ResampleImageFilter<TImageType, TImageType> ResampleImageFilterType;
		typedef typename ResampleImageFilterType::Pointer        ResampleImageFilterPointer;
		ResampleImageFilterPointer resampler = ResampleImageFilterType::New();
		resampler->SetInput( inputImage );
		resampler->SetOutputParametersFromImage( inputImage );
		resampler->SetOutputOrigin( inputImage->GetOrigin() );
		resampler->SetSize( size );
		resampler->SetOutputSpacing( spacing );
		resampler->SetDefaultPixelValue( defaultValue );
		resampler->Update();

		return resampler->GetOutput();
	}

	/************************************************************************/
	/* Pading Image                                                         */
	/************************************************************************/

	template<typename TImageType>
	typename TImageType::Pointer
		padImage( const typename TImageType* inputImage, typename TImageType::SizeType lowerBound, typename TImageType::SizeType upperBound, double defaultValue=0)
	{
		kmStaticImageMacro(TImageType);
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);

		typedef itk::PadImageFilter<TImageType, TImageType> PadImageFilterType;
		typedef typename PadImageFilterType::Pointer        PadImageFilterPointer;

		PadImageFilterPointer filter = PadImageFilterType::New();
		filter->SetInput( inputImage );
		filter->SetPadLowerBound( lowerBound );
		filter->SetPadUpperBound( upperBound );
		itk::ConstantBoundaryCondition< TImageType > bc;
		bc.SetConstant( defaultValue );
		filter->SetBoundaryCondition( &bc );
		filter->Update();

		return filter->GetOutput();
	}

	/************************************************************************/
	/* Translate Image                                                      */
	/************************************************************************/
	template<typename TImageType>
	typename TImageType::Pointer
		translateImage( const typename TImageType* inputImage, itk::Vector<double> & offsetVector, const typename TImageType* referenceImage )
	{
		//kmStaticImageMacro(TImageType);
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);

		typedef itk::TranslationTransform<double, Dimension> TranslationTransformType;
		typedef typename TranslationTransformType::Pointer TranslationTransformPointer;

		TranslationTransformPointer translationTransform = TranslationTransformType::New();
		translationTransform->Translate(offsetVector);

		typedef itk::ResampleImageFilter<TImageType, TImageType> ResampleImageFilterType;
		typedef typename ResampleImageFilterType::Pointer                 ResampleImageFilterPointer;
		ResampleImageFilterPointer resampleFilter = ResampleImageFilterType::New();
		resampleFilter->SetTransform(translationTransform.GetPointer());
		resampleFilter->SetInput(inputImage);
		resampleFilter->SetOutputParametersFromImage( referenceImage );
		resampleFilter->Update();

		return resampleFilter->GetOutput();
	}

	/************************************************************************/
	/* Binary Close Image                                                   */
	/************************************************************************/
	template<typename TImageType>
	typename TImageType::Pointer
		binaryClose( const typename TImageType* inputImage, typename TImageType::SizeType kernelRadius )
	{
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);

		typedef itk::BinaryBallStructuringElement<unsigned char, 3> BinaryBallStructuringElementType;
		typedef itk::BinaryMorphologicalClosingImageFilter<TImageType, TImageType, BinaryBallStructuringElementType> BinaryMorphologicalClosingImageFilterType;
		typedef typename BinaryMorphologicalClosingImageFilterType::Pointer BinaryMorphologicalClosingImageFilterPointer;
		BinaryMorphologicalClosingImageFilterPointer filter = BinaryMorphologicalClosingImageFilterType::New();
		BinaryBallStructuringElementType structuringElement;
		BinaryBallStructuringElementType::RadiusType radius;
		for(int i=0;i<Dimension;i++)
		{
			radius[i] = kernelRadius[i];
		}
		structuringElement.SetRadius(radius);
		structuringElement.CreateStructuringElement();
		filter->SetInput( inputImage );
		filter->SetKernel( structuringElement );
		filter->SetForegroundValue( 1 );
		filter->Update();

		return filter->GetOutput();
	}

	/************************************************************************/
	/* Binary Open Image                                                   */
	/************************************************************************/
	template<typename TImageType>
	typename TImageType::Pointer
		binaryOpen( const typename TImageType* inputImage, typename TImageType::SizeType kernelRadius )
	{
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);

		typedef itk::BinaryBallStructuringElement<unsigned char, 3> BinaryBallStructuringElementType;
		typedef itk::BinaryMorphologicalOpeningImageFilter<TImageType, TImageType, BinaryBallStructuringElementType> BinaryMorphologicalClosingImageFilterType;
		typedef typename BinaryMorphologicalClosingImageFilterType::Pointer BinaryMorphologicalClosingImageFilterPointer;
		BinaryMorphologicalClosingImageFilterPointer filter = BinaryMorphologicalClosingImageFilterType::New();
		BinaryBallStructuringElementType structuringElement;
		BinaryBallStructuringElementType::RadiusType radius;
		for(int i=0;i<Dimension;i++)
		{
			radius[i] = kernelRadius[i];
		}
		structuringElement.SetRadius(radius);
		structuringElement.CreateStructuringElement();
		filter->SetInput( inputImage );
		filter->SetKernel( structuringElement );
		filter->SetForegroundValue( 1 );
		filter->Update();

		return filter->GetOutput();
	}

	/************************************************************************/
	/* Binary Dilate Image                                                   */
	/************************************************************************/
	template<typename TImageType>
	typename TImageType::Pointer
		binaryDilate( const typename TImageType* inputImage, typename TImageType::SizeType kernelRadius )
	{
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);

		typedef itk::BinaryBallStructuringElement<unsigned char, 3> BinaryBallStructuringElementType;
		typedef itk::BinaryDilateImageFilter<TImageType, TImageType, BinaryBallStructuringElementType> FilterType;
		FilterType::Pointer filter = FilterType::New();
		BinaryBallStructuringElementType structuringElement;
		BinaryBallStructuringElementType::RadiusType radius;
		for(int i=0;i<Dimension;i++)
		{
			radius[i] = kernelRadius[i];
		}
		structuringElement.SetRadius(radius);
		structuringElement.CreateStructuringElement();
		filter->SetInput( inputImage );
		filter->SetKernel( structuringElement );
		filter->SetForegroundValue( 1 );
		filter->Update();

		return filter->GetOutput();
	}

	/************************************************************************/
	/* Binary Median Image                                                   */
	/************************************************************************/
	template<typename TImageType>
	typename TImageType::Pointer
		binaryMedianSmooth( const typename TImageType* inputImage, typename TImageType::SizeType kernelRadius )
	{
		typedef itk::BinaryMedianImageFilter<TImageType, TImageType> BinaryMedianImageFilterType;
		BinaryMedianImageFilterType::Pointer filter = BinaryMedianImageFilterType::New();
		filter->SetInput( inputImage );
		filter->SetForegroundValue( 1 );
		filter->SetBackgroundValue( 0 );
		filter->SetRadius( kernelRadius );
		filter->Update();

		return filter->GetOutput();
	}

	/************************************************************************/
	/* Generate VTKPolyData from Binary Image                               */
	/************************************************************************/
	template<typename TImageType, typename TMeshType>
	typename TMeshType::Pointer
		generateMeshFromBinary( const typename TImageType* inputImage )
	{
		//kmStaticMeshMacro(TMeshType);

		typedef itk::BinaryMask3DMeshSource<TImageType, TMeshType> BinaryMask3DMeshSourceType;
		typedef typename BinaryMask3DMeshSourceType::Pointer       BinaryMask3DMeshSourcePointer;
		BinaryMask3DMeshSourcePointer meshSource = BinaryMask3DMeshSourceType::New();
		meshSource->SetInput( inputImage );
		meshSource->SetObjectValue( 1 );
		meshSource->SetRegionOfInterest( inputImage->GetLargestPossibleRegion() );
		meshSource->Update();

		return meshSource->GetOutput();
	}


	/************************************************************************/
	/*                                                                      */
	/************************************************************************/
	template<typename TMeshType, typename TImageType, typename TReferenceImageType>
	typename TImageType::Pointer
		generateBinaryFromMesh( 
		const typename TMeshType* inputMesh, 
		const typename TReferenceImageType* referenceImage )
	{
		kmStaticImageMacro(TImageType);
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);

		typedef itk::TriangleMeshToBinaryImageFilter<TMeshType, TImageType> TriangleMeshToBinaryImageFilterType;
		typedef typename TriangleMeshToBinaryImageFilterType::Pointer       TriangleMeshToBinaryImageFilterPointer;
		TriangleMeshToBinaryImageFilterPointer filter = TriangleMeshToBinaryImageFilterType::New();
		filter->SetInput( const_cast<TMeshType*>(inputMesh) );

		filter->SetSize( referenceImage->GetLargestPossibleRegion().GetSize() );
		filter->SetSpacing( referenceImage->GetSpacing() );
		filter->SetIndex( referenceImage->GetLargestPossibleRegion().GetIndex() );
		filter->SetOrigin( referenceImage->GetOrigin() );
		filter->SetInsideValue( 1 );
		filter->SetOutsideValue( 0 );
		filter->Update();

		return filter->GetOutput();
	}

	template<class ImageType, class PointSetType>
	void
		binaryImageToPointSet(const typename ImageType* image, typename PointSetType::Pointer & pointSet)
	{
		itk::ImageRegionConstIterator< ImageType >  it( image, image->GetRequestedRegion() );

		typedef PointSetType::PointType PointType;
		PointType point;

		unsigned long pointId = 0;

		it.GoToBegin();
		while(!it.IsAtEnd())
		{
			if (it.Get() != 0)
			{
				image->TransformIndexToPhysicalPoint( it.GetIndex() , point );
				pointSet->SetPoint( pointId, point );

				// Transfer the pixel data to the value associated with the point.
				pointSet->SetPointData( pointId, it.Get() );

				++pointId;

			}

			++it;
		}
	}

	template<typename TMeshType, typename TImageType>
	typename TImageType::Pointer
		generateBinaryFromMesh( 
														const typename TMeshType* inputMesh, 
														typename TImageType::SizeType outputSize, 
														typename TImageType::SpacingType outputSpacing )
	{
		kmStaticImageMacro(TImageType);
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);

		typedef itk::TriangleMeshToBinaryImageFilter<TMeshType, TImageType> TriangleMeshToBinaryImageFilterType;
		typedef typename TriangleMeshToBinaryImageFilterType::Pointer       TriangleMeshToBinaryImageFilterPointer;
		TriangleMeshToBinaryImageFilterPointer filter = TriangleMeshToBinaryImageFilterType::New();
		filter->SetInput( const_cast<TMeshType*>(inputMesh) );

		PointType origin;
		IndexType index;
		for( int i=0;i<Dimension;i++ )
		{
			origin[i] = 0.0;
			index[i] = 0;
		}

		filter->SetSize( outputSize );
		filter->SetSpacing( outputSpacing );
		filter->SetIndex( index );
		filter->SetOrigin( origin );
		filter->SetInsideValue( 1 );
		filter->SetOutsideValue( 0 );
		filter->Update();

		return filter->GetOutput();
	}

	/************************************************************************/
	/* Calculate Mauer Distance Map of Binary Image                         */
	/************************************************************************/
	template<typename TImageType, typename TDistanceMapType>
	typename TDistanceMapType::Pointer
		calculateDistanceMap( const typename TImageType* inputImage )
	{
		typedef itk::SignedMaurerDistanceMapImageFilter<TImageType, TDistanceMapType> SignedMaurerDistanceMapImageFilterType;
		typedef typename SignedMaurerDistanceMapImageFilterType::Pointer              SignedMaurerDistanceMapImageFilterPointer;
		SignedMaurerDistanceMapImageFilterPointer distanceFilter = SignedMaurerDistanceMapImageFilterType::New();
		distanceFilter->SetInput( inputImage );
		distanceFilter->Update();

		return distanceFilter->GetOutput();
	}

	template<typename TImageType, typename TGradientImageType>
	typename TGradientImageType::Pointer
		calculateRecursiveGradientImage( const typename TImageType* inputImage, double sigma )
	{
		typedef itk::GradientRecursiveGaussianImageFilter<TImageType, TGradientImageType> GradientRecursiveGaussianImageFilterType;
		GradientRecursiveGaussianImageFilterType::Pointer filter = GradientRecursiveGaussianImageFilterType::New();
		filter->SetInput( inputImage );
		filter->SetSigma( sigma );
		filter->SetNormalizeAcrossScale( true );
		filter->Update();

		return filter->GetOutput();
	}

	template<typename TImageType, typename TGradientImageType>
	typename TGradientImageType::Pointer
		calculateGradientMagnitudeImage( const typename TImageType* inputImage, double sigma )
	{
		typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<TImageType, TGradientImageType> GradientMagnitudeRecursiveGaussianImageFilterType;
		GradientMagnitudeRecursiveGaussianImageFilterType::Pointer filter = GradientMagnitudeRecursiveGaussianImageFilterType::New();
		filter->SetInput( inputImage );
		filter->SetSigma( sigma );
		filter->Update();

		return filter->GetOutput();
	}

	template<typename TImageType, typename TSigmoidImageType>
	typename TSigmoidImageType::Pointer
		calculateSigmoidImage( const typename TImageType* inputImage, double outputMinimum, double outputMaximum, double alpha, double beta )
	{
		typedef itk::SigmoidImageFilter<TImageType, TSigmoidImageType> SigmoidImageFilterType;
		SigmoidImageFilterType::Pointer filter = SigmoidImageFilterType::New();
		filter->SetInput( inputImage );
		filter->SetOutputMinimum( outputMinimum );
		filter->SetOutputMaximum( outputMaximum );
		filter->SetAlpha( alpha );
		filter->SetBeta( beta );
		filter->Update();

		return filter->GetOutput();
	} 

	template<typename TGradientImageType>
	typename TGradientImageType::Pointer
		calcuateGVF(const typename TGradientImageType* inputGradient, int iterations, double noiseLevel)
	{
		typedef itk::GradientVectorFlowImageFilter< TGradientImageType, TGradientImageType >  FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput( inputGradient );
		filter->SetIterationNum( iterations );
		filter->SetNoiseLevel( noiseLevel );
		filter->Update();

		return filter->GetOutput();
	}

	/************************************************************************/
	/* Warp/Transform a Mesh                                                */
	/************************************************************************/
	template<typename TMeshType, typename TDisplacementType>
	typename TMeshType::Pointer
		warpMesh( const typename TMeshType* inputMesh, const typename TDisplacementType* displacement  )
	{
		typedef itk::WarpMeshFilter<TMeshType, TMeshType, TDisplacementType> WarpMeshFilterType;
		typedef typename WarpMeshFilterType::Pointer                         WarpMeshFilterPointer;
		WarpMeshFilterPointer meshWarper = WarpMeshFilterType::New();
		meshWarper->SetInput( inputMesh );
		meshWarper->SetDisplacementField( displacement );
		meshWarper->Update();

		return meshWarper->GetOutput();
	}

	template<typename TMeshType, typename TTransformType>
	typename TMeshType::Pointer
		transformMesh( const typename TMeshType* inputMesh, typename TTransformType* transform  )
	{
		typedef itk::TransformMeshFilter<TMeshType, TMeshType, TTransformType> TransformMeshFilterType;
		TransformMeshFilterType::Pointer filter = TransformMeshFilterType::New();
		filter->SetInput( inputMesh );
		filter->SetTransform( transform );
		filter->Update();

		return filter->GetOutput();
	}

	template<typename TMeshType, typename TTransformType>
	void
		transformMesh( const typename TMeshType* inputMesh, typename TMeshType* outputMesh, typename TTransformType* transform  )
	{
		typedef TMeshType::PointType PointType;
		typedef TMeshType::PointsContainer::ConstPointer  PointsContainerConstPointer;
		typedef TMeshType::PointsContainer::Pointer       PointsContainerPointer;
		typedef TMeshType::PointsContainer::ConstIterator PointsContainerConstIterator;

		PointsContainerConstPointer inputPoints = inputMesh->GetPoints();
		PointsContainerPointer      outputPoints = outputMesh->GetPoints();

		if (outputPoints->Size() != inputPoints->Size())
		{
			outputPoints->Reserve( inputPoints->Size() );
		}

		PointsContainerConstIterator inputPointIt = inputPoints->Begin();
		PointsContainerConstIterator outputPointIt = outputPoints->Begin();

		while( (inputPointIt!=inputPoints->End()) && (outputPointIt!=outputPoints->End()) )
		{
			if ( inputPointIt->Index() != outputPointIt->Index() )
			{
				KM_DEBUG_ERROR( "Index of input mesh != index of output mesh" );
			}

			PointType p_old = inputPointIt->Value();
			PointType p_new = transform->TransformPoint( p_old );

			outputPoints->InsertElement( inputPointIt->Index(), p_new );

			inputPointIt++;
			outputPointIt++;
		}
	}

	/************************************************************************/
	/* Warp/Transform a Image                                               */
	/************************************************************************/
	template<typename TImageType, typename TDisplacementType>
	typename TImageType::Pointer
		warpImage( const typename TImageType* inputImage, const typename TDisplacementType* displacement  )
	{
		typedef itk::WarpImageFilter<TImageType, TImageType, TDisplacementType> FilterType;
		FilterType::Pointer meshWarper = FilterType::New();
		meshWarper->SetInput( inputImage );
		meshWarper->SetDisplacementField( displacement );
		meshWarper->SetOutputSpacing( inputImage->GetSpacing() );
		meshWarper->SetOutputOrigin(  inputImage->GetOrigin() );
		meshWarper->Update();

		return meshWarper->GetOutput();
	}

	template<typename TImageType, typename TTransformType>
	typename TImageType::Pointer
		transformImage( const typename TImageType* inputImage, const typename TTransformType* transform  )
	{
		typedef itk::ResampleImageFilter<TImageType, TImageType> FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput( inputImage );
		filter->SetTransform( transform );
		filter->SetSize(    inputImage->GetLargestPossibleRegion().GetSize() );
		filter->SetOutputOrigin(  inputImage->GetOrigin() );
		filter->SetOutputSpacing( inputImage->GetSpacing() );
		filter->SetOutputDirection( inputImage->GetDirection() );
		filter->Update();

		return filter->GetOutput();
	}


	/************************************************************************/
	/* Smooth a PolyData                                                    */
	/************************************************************************/
	vtkSmartPointer<vtkPolyData>
		smoothPolyData( vtkPolyData * inputPolyData, int iterations, double relaxfactor = 0.025 )
	{
		vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
		smoothFilter->SetInput( inputPolyData );
		smoothFilter->SetNumberOfIterations( iterations );
		smoothFilter->BoundarySmoothingOn();
		smoothFilter->SetFeatureAngle( 45 );
		smoothFilter->SetEdgeAngle( 15 );
		smoothFilter->SetRelaxationFactor( 0.05 );

		//smoothFilter->Print(std::cout);

		smoothFilter->Update();

		return smoothFilter->GetOutput();
	}

	/************************************************************************/
	/* Calculate the normals of a PolyData                                  */
	/************************************************************************/
	vtkSmartPointer<vtkPolyData>
		calculateNormals( vtkPolyData* inputPolyData )
	{
		vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
		normalGenerator->SetInput(inputPolyData);
		normalGenerator->ComputePointNormalsOn();
		normalGenerator->ComputeCellNormalsOff();
		normalGenerator->Update();

		return normalGenerator->GetOutput();
	}

	/************************************************************************/
	/* Theshold Image                                                       */
	/************************************************************************/
	template<typename TImageType>
	typename TImageType::Pointer
		thresholdImage( const TImageType* inputImage, double lowValue, double highValue, double outsideValue = 0 )
	{
		typedef itk::ThresholdImageFilter<TImageType> ThresholdImageFilterType;
		typename ThresholdImageFilterType::Pointer thresholdFilter = ThresholdImageFilterType::New();
		thresholdFilter->SetInput( inputImage );
		thresholdFilter->ThresholdOutside( static_cast<typename TImageType::PixelType>(lowValue), static_cast<typename TImageType::PixelType>(highValue) );
		thresholdFilter->SetOutsideValue( outsideValue );
		thresholdFilter->Update();

		return thresholdFilter->GetOutput();
	}

	/************************************************************************/
	/* Theshold Image                                                       */
	/************************************************************************/
	template<typename TImageType>
	typename TImageType::Pointer
		thresholdAboveImage( const TImageType* inputImage, double highValue, double outsideValue = highValue )
	{
		typedef itk::ThresholdImageFilter<TImageType> ThresholdImageFilterType;
		typename ThresholdImageFilterType::Pointer thresholdFilter = ThresholdImageFilterType::New();
		thresholdFilter->SetInput( inputImage );
		thresholdFilter->ThresholdAbove ( static_cast<typename TImageType::PixelType>(highValue) );
		thresholdFilter->SetOutsideValue( outsideValue );
		thresholdFilter->Update();

		return thresholdFilter->GetOutput();
	}

	template<typename TImageType, typename TBinaryImageType>
	typename TBinaryImageType::Pointer
		binaryThresholdImage( const TImageType* inputImage, double lowValue, double highValue, double insideValue, double outsideValue )
	{
		typedef itk::BinaryThresholdImageFilter<TImageType, TBinaryImageType> ThresholdImageFilterType;
		typename ThresholdImageFilterType::Pointer thresholdFilter = ThresholdImageFilterType::New();
		thresholdFilter->SetInput( inputImage );
		thresholdFilter->SetLowerThreshold(lowValue);
		thresholdFilter->SetUpperThreshold(highValue);
		thresholdFilter->SetInsideValue(insideValue);
		thresholdFilter->SetOutsideValue(outsideValue);
		thresholdFilter->Update();

		return thresholdFilter->GetOutput();
	}

	template<typename TImageType, typename TBinaryImageType>
	typename TBinaryImageType::Pointer
		binaryContourImage( const TImageType* inputImage, double foregroundValue, double backgourndValue )
	{
		typedef itk::BinaryContourImageFilter <TImageType, TBinaryImageType> BinaryContourImageFilterType;
		typename BinaryContourImageFilterType::Pointer filter = BinaryContourImageFilterType::New();
		filter->SetInput( inputImage );
		filter->SetForegroundValue(foregroundValue);
		filter->SetBackgroundValue(backgourndValue);
		filter->FullyConnectedOn();
		filter->Update();

		return filter->GetOutput();
	}

	template<typename TriangleMeshType, typename SimplexMeshType>
	typename SimplexMeshType::Pointer
		triangleMeshToSimplexMesh(const typename TriangleMeshType* inputMesh)
	{
		typedef itk::TriangleMeshToSimplexMeshFilter<TriangleMeshType, SimplexMeshType> FilterType;
		typename FilterType::Pointer filter = FilterType::New();
		try
		{
			filter->SetInput( inputMesh );
			filter->Update();
		}
		catch ( itk::ExceptionObject & e  )
		{
			KM_PRINT_EXCEPTION(e);
			//std::cout<<"Exception occurs when converting triangle mesh to simplex mesh!"<<std::endl;
			//std::cout<<
			//std::cout<<e.GetFile()<<std::endl;
			//std::cout<<e.GetLine()<<std::endl;
			//std::cout<<e.GetLocation()<<std::endl;
		}
		catch (...)
		{
			std::cout<<"Exception occurs when converting triangle mesh to simplex mesh!"<<std::endl;
		}

		SimplexMeshType::Pointer outputMesh = filter->GetOutput();
		//outputMesh->DisconnectPipeline();

		// Do not delete the filter because it is a smart pointer
		//filter->Delete();

		return outputMesh;
	}

	template<typename SimplexMeshType, typename TriangleMeshType>
	typename TriangleMeshType::Pointer
		simplexMeshToTriangleMesh(const typename SimplexMeshType* inputMesh)
	{
		typedef itk::SimplexMeshToTriangleMeshFilter<SimplexMeshType, TriangleMeshType> FilterType;
		typename FilterType::Pointer filter = FilterType::New();
		try
		{
			filter->SetInput( inputMesh );
			filter->Update();
		}
		catch ( itk::ExceptionObject & e  )
		{
			KM_PRINT_EXCEPTION(e);
			//std::cout<<"Exception occurs when converting triangle mesh to simplex mesh!"<<std::endl;
			//std::cout<<
			//std::cout<<e.GetFile()<<std::endl;
			//std::cout<<e.GetLine()<<std::endl;
			//std::cout<<e.GetLocation()<<std::endl;
		}
		catch (...)
		{
			std::cout<<"Exception occurs when converting triangle mesh to simplex mesh!"<<std::endl;
		}

		TriangleMeshType::Pointer outputMesh = filter->GetOutput();
		//outputMesh->DisconnectPipeline();

		// Do not delete the filter because it is a smart pointer
		//filter->Delete();

		return outputMesh;
	}

	template<typename TTransform>
	void
		writeTransform(const char* filename, const typename TTransform* transform)
	{
		typedef itk::TransformFileWriter WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetInput(transform);
		writer->SetFileName(filename);
		writer->Update();
	}

	template<typename TTransform>
	typename TTransform::Pointer
		readTransform(const char* filename)
	{
		TTransform::Pointer transform = TTransform::New();
		const char* transformname = transform->GetNameOfClass();
		std::cout<<transformname<<std::endl;

		itk::TransformFactoryBase::RegisterDefaultTransforms();
		itk::TransformFileReader::Pointer reader = itk::TransformFileReader::New();
		reader->SetFileName(filename);
		try
		{
			reader->Update();
		}
		catch( itk::ExceptionObject & excp )
		{
			std::cerr << "Error while reading the transform file" << std::endl;
			std::cerr << excp << std::endl;
			std::cerr << "[FAILED]" << std::endl;
			return transform;
		}

		typedef itk::TransformFileReader::TransformListType * TransformListType;
		TransformListType transforms = reader->GetTransformList();
		std::cout << "Number of transforms = " << transforms->size() << std::endl;

		itk::TransformFileReader::TransformListType::const_iterator it = transforms->begin();
		itk::TransformFileReader::TransformListType::const_iterator itEnd = transforms->end();

		while(it!=itEnd)
		{
			if(!strcmp((*it)->GetNameOfClass(),transformname))
			{
				transform = static_cast<TTransform*>((*it).GetPointer());
			}
			it++;
		}

		return transform;
	}

	template<typename TImageType>
	typename TImageType::Pointer
		gaussSmooth(const typename TImageType* inputImage, double sigma)
	{
		typedef itk::DiscreteGaussianImageFilter<TImageType, TImageType> FilterType;
		FilterType::Pointer filter=FilterType::New();
		filter->SetInput(inputImage);
		filter->SetVariance(sigma);
		filter->Update();

		return filter->GetOutput();
	}


	template<typename TImageType>
	typename TImageType::Pointer
		minMaxSmooth(const typename TImageType* inputImage, int numberOfIterations, double timeStep, double radius)
	{
		typedef itk::Image<float, TImageType::ImageDimension> FloatImageType;
		FloatImageType::Pointer image = km::castImage<TImageType, FloatImageType>( inputImage );

		typedef itk::MinMaxCurvatureFlowImageFilter<FloatImageType, FloatImageType> MinMaxCurvatureFlowImageFilterType;
		typedef MinMaxCurvatureFlowImageFilterType::Pointer MinMaxCurvatureFlowImageFilterPointer;
		MinMaxCurvatureFlowImageFilterPointer minMaxFilter = MinMaxCurvatureFlowImageFilterType::New();
		minMaxFilter->SetInput( image );
		minMaxFilter->SetTimeStep( timeStep );
		minMaxFilter->SetNumberOfIterations( numberOfIterations );
		minMaxFilter->SetStencilRadius( radius );
		minMaxFilter->Update();

		FloatImageType::Pointer smoothed = minMaxFilter->GetOutput();

		TImageType::Pointer smoothedImage = km::castImage<FloatImageType, TImageType>( smoothed );

		return smoothedImage;
	}

	template<typename TImageType>
	typename TImageType::Pointer
		anisotropicSmooth(const typename TImageType* inputImage, int numberOfIterations, double timeStep, double conductance)
	{
		typedef itk::GradientAnisotropicDiffusionImageFilter<  TImageType,  TImageType >  AnisotropicSmoothFilterType;
		typedef AnisotropicSmoothFilterType::Pointer AnisotropicSmoothFilterPointer;
		AnisotropicSmoothFilterPointer anisoSmoothfilter = AnisotropicSmoothFilterType::New();
		anisoSmoothfilter->SetNumberOfIterations( numberOfIterations );
		anisoSmoothfilter->SetTimeStep( timeStep );
		anisoSmoothfilter->SetConductanceParameter( conductance );
		anisoSmoothfilter->SetInput( inputImage );
		anisoSmoothfilter->UseImageSpacingOn();
		anisoSmoothfilter->Update();

		return anisoSmoothfilter->GetOutput();
	}

	template<typename TImageType>
	typename TImageType::Pointer
		medianSmooth(const typename TImageType* inputImage, int radiusValue)
	{
		typedef itk::MedianImageFilter<TImageType, TImageType > FilterType;
		FilterType::Pointer medianFilter = FilterType::New();
		FilterType::InputSizeType radius;
		radius.Fill(radiusValue);
		medianFilter->SetRadius(radius);
		medianFilter->SetInput( inputImage );
		medianFilter->Update();

		return medianFilter->GetOutput();
	}

	template<typename TImageType>
	typename TImageType::Pointer
		medianSmooth(const typename TImageType* inputImage, typename TImageType::SizeType radius)
	{
		typedef itk::MedianImageFilter<TImageType, TImageType > FilterType;
		FilterType::Pointer medianFilter = FilterType::New();
		medianFilter->SetRadius(radius);
		medianFilter->SetInput( inputImage );
		medianFilter->Update();

		return medianFilter->GetOutput();
	}

	template<typename TImageType>
	typename TImageType::Pointer
		discreteGaussSmooth(const typename TImageType* inputImage, double variance)
	{
		typedef itk::DiscreteGaussianImageFilter<TImageType,TImageType> DiscreteGaussianImageFilterType;
		DiscreteGaussianImageFilterType::Pointer filter = DiscreteGaussianImageFilterType::New();
		filter->SetInput( inputImage );
		filter->SetVariance( variance );
		filter->Update();

		return filter->GetOutput();
	}

	template<class TImageType, class ReferenceImageType, class TTransformType>
	typename TImageType::Pointer
		transformImageByReference(
		const typename TImageType* inputimage, 
		const typename ReferenceImageType* referenceimage, 
		const typename TTransformType* transform,
		double defaultValue = 0,
		int   interpolatorType = 0) //0:liner interpolator 1:nearest...
	{
		typedef itk::ResampleImageFilter< TImageType, TImageType,double >  ResampleFilterType;
		ResampleFilterType::Pointer resample = ResampleFilterType::New();
		resample->SetOutputOrigin( referenceimage->GetOrigin() );
		resample->SetOutputSpacing( referenceimage->GetSpacing() );
		resample->SetOutputDirection ( referenceimage->GetDirection() );
		resample->SetOutputStartIndex ( referenceimage->GetLargestPossibleRegion().GetIndex() );
		resample->SetSize ( referenceimage->GetLargestPossibleRegion().GetSize() );
		resample->SetDefaultPixelValue( defaultValue );
		resample->SetTransform( transform );
		resample->SetInput( inputimage );

		typedef itk::LinearInterpolateImageFunction<TImageType> LinearInterpolateType;
		typedef itk::NearestNeighborInterpolateImageFunction<TImageType> NearestNeighborInterpolateType;
		switch (interpolatorType)
		{
		case 0:
			{
				LinearInterpolateType::Pointer interpolator = LinearInterpolateType::New();
				resample->SetInterpolator( interpolator );
				break;
			}
		case 1:
			{
				NearestNeighborInterpolateType::Pointer interpolator2 = NearestNeighborInterpolateType::New();
				resample->SetInterpolator( interpolator2 );
				break;
			}
			
		default:
			{
				std::cout<<"Unknown interpolator type!"<<std::endl;
				break;
			}
		}

		resample->Update();

		return resample->GetOutput();
	}

	template<class TImageType, class ReferenceImageType>
	typename TImageType::Pointer
		resampleImageByReference(
		const typename TImageType* inputimage, 
		const typename ReferenceImageType* referenceimage, 
		double defaultValue = 0,
		int interpolatorType = 0)
	{
		typedef itk::ResampleImageFilter< TImageType, TImageType,double >  ResampleFilterType;
		ResampleFilterType::Pointer resample = ResampleFilterType::New();
		//resample->SetOutputParametersFromImage( referenceimage );
		resample->SetOutputOrigin( referenceimage->GetOrigin() );
		resample->SetOutputSpacing( referenceimage->GetSpacing() );
		resample->SetOutputDirection ( referenceimage->GetDirection() );
		resample->SetOutputStartIndex ( referenceimage->GetLargestPossibleRegion().GetIndex() );
		resample->SetSize ( referenceimage->GetLargestPossibleRegion().GetSize() );
		resample->SetDefaultPixelValue( defaultValue );
		//resample->SetTransform( transform );
		resample->SetInput( inputimage );

		typedef itk::LinearInterpolateImageFunction<TImageType> LinearInterpolateType;
		typedef itk::NearestNeighborInterpolateImageFunction<TImageType> NearestNeighborInterpolateType;
		switch (interpolatorType)
		{
		case 0:
			{
				LinearInterpolateType::Pointer interpolator = LinearInterpolateType::New();
				resample->SetInterpolator( interpolator );
				break;
			}
		case 1:
			{
				NearestNeighborInterpolateType::Pointer interpolator2 = NearestNeighborInterpolateType::New();
				resample->SetInterpolator( interpolator2 );
				break;
			}

		default:
			{
				std::cout<<"Unknown interpolator type!"<<std::endl;
				break;
			}
		}

		resample->Update();

		return resample->GetOutput();
	}

	template<class TImageType, class TTransformType>
	typename TImageType::Pointer
		transformImage(const typename TImageType* inputimage, const typename TImageType* referenceimage, const typename TTransformType* transform, double defaultValue = 0)
	{
		typedef itk::ResampleImageFilter< TImageType, TImageType,double >  ResampleFilterType;
		ResampleFilterType::Pointer resample = ResampleFilterType::New();
		//resample->SetOutputParametersFromImage( referenceimage );
		resample->SetOutputOrigin( referenceimage->GetOrigin() );
		resample->SetOutputSpacing( referenceimage->GetSpacing() );
		resample->SetOutputDirection ( referenceimage->GetDirection() );
		resample->SetOutputStartIndex ( referenceimage->GetLargestPossibleRegion().GetIndex() );
		resample->SetSize ( referenceimage->GetLargestPossibleRegion().GetSize() );
		resample->SetDefaultPixelValue( defaultValue );
		resample->SetTransform( transform );
		resample->SetInput( inputimage );
		resample->Update();

		return resample->GetOutput();
	}

	template<typename TImageType, typename TBinaryImageType>
	typename TImageType::Pointer
		maskImage(const typename TImageType* inputimage, const typename TBinaryImageType* maskimage)
	{
		typedef itk::MultiplyImageFilter< TImageType, TBinaryImageType >  FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput1( inputimage );
		filter->SetInput2( maskimage );
		try
		{
			filter->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return NULL; 
		}

		return filter->GetOutput();
	}

	template<typename ImageType1, typename ImageType2, typename OutputImageType>
	typename OutputImageType::Pointer
		multiplyImage( const ImageType1* image1, const ImageType2* image2 )
	{
		typedef itk::MultiplyImageFilter<ImageType1, ImageType2, OutputImageType> FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput1( image1 );
		filter->SetInput2( image2 );
		filter->Update();

		return filter->GetOutput();
	}

	template<typename TImageType>
	typename TImageType::Pointer
		shiftScale(const typename TImageType* inputImage, double shiftValue, double scaleValue)
	{
		typedef itk::ShiftScaleImageFilter<TImageType, TImageType> FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput( inputImage );
		filter->SetScale( scaleValue );
		filter->SetShift( shiftValue );
		filter->Update();

		return filter->GetOutput();
	}

	template<typename InputImageType, typename OutputImageType>
	typename OutputImageType::Pointer
		rescaleIntensity(const typename InputImageType* inputImage, double outputMinimum, double outputMaximum)
	{
		typedef itk::RescaleIntensityImageFilter<InputImageType, OutputImageType> FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput(inputImage);
		filter->SetOutputMinimum(outputMinimum);
		filter->SetOutputMaximum(outputMaximum);
		filter->Update();

		return filter->GetOutput();
	}

	template<typename ImageType>
	typename ImageType::Pointer
		intensityWindow( typename ImageType* inputImage, double windowsmin, double windowmax, double outputmin, double outputmax )
	{
		typedef itk::IntensityWindowingImageFilter<ImageType> FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput( inputImage );
		filter->SetWindowMinimum(windowsmin);
		filter->SetWindowMaximum(windowmax);
		filter->SetOutputMinimum(outputmin);
		filter->SetOutputMaximum(outputmax);
		filter->Update();

		return filter->GetOutput();
	}

	template<typename TImageType>
	void
		calculateMinAndMax(const typename TImageType* inputImage, double& minimumValue, double& maximumValue)
	{
		typedef itk::MinimumMaximumImageCalculator<TImageType> MinimumMaximumImageCalculatorType;
		MinimumMaximumImageCalculatorType::Pointer calculator = MinimumMaximumImageCalculatorType::New();
		calculator->SetImage( inputImage );
		calculator->Compute();

		minimumValue = calculator->GetMinimum();
		maximumValue = calculator->GetMaximum();
	}

	template<typename TImageType>
	void
		calculateMinAndMax(const typename TImageType* inputImage, double& minimumValue, double& maximumValue, typename TImageType::IndexType & minIndex, typename TImageType::IndexType & maxIndex)
	{
		typedef itk::MinimumMaximumImageCalculator<TImageType> MinimumMaximumImageCalculatorType;
		MinimumMaximumImageCalculatorType::Pointer calculator = MinimumMaximumImageCalculatorType::New();
		calculator->SetImage( inputImage );
		calculator->Compute();

		minimumValue = calculator->GetMinimum();
		maximumValue = calculator->GetMaximum();

		minIndex = calculator->GetIndexOfMinimum();
		maxIndex = calculator->GetIndexOfMaximum();
	}

	template<typename InputImageType, typename ReferenceImageType>
	typename InputImageType::Pointer
		changeByCopying(const typename InputImageType* inputImage, const typename ReferenceImageType* referenceImage)
	{
		typedef itk::ChangeInformationImageFilter<InputImageType> ChangeInformationImageFilterType;
		ChangeInformationImageFilterType::Pointer changer = ChangeInformationImageFilterType::New();
		changer->SetInput( inputImage );
		changer->SetOutputOrigin( referenceImage->GetOrigin() );
		changer->ChangeOriginOn();
		//changer->SetOutputSpacing( referenceImage->GetSpacing() );
		//changer->ChangeSpacingOn();
		changer->Update();

		return changer->GetOutput();
	}

	template<typename InputImageType, typename ReferenceImageType>
	typename InputImageType::Pointer
		resampleByCopying(const typename InputImageType* inputImage, const typename ReferenceImageType* referenceImage)
	{
		typedef itk::ResampleImageFilter<InputImageType, InputImageType> ResampleImageFilterType;
		ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
		resampler->SetInput( inputImage );
		resampler->SetOutputOrigin ( referenceImage->GetOrigin() );
		resampler->SetOutputSpacing ( referenceImage->GetSpacing() );
		resampler->SetOutputDirection ( referenceImage->GetDirection() );
		resampler->SetOutputStartIndex ( referenceImage->GetLargestPossibleRegion().GetIndex() );
		resampler->SetSize ( referenceImage->GetLargestPossibleRegion().GetSize() );
		resampler->Update();

		return resampler->GetOutput();
	}

	/************************************************************************/
	/* 得到最小包围盒的上下index                                            */
	/************************************************************************/
	template<typename TImageType>
	void
		getBoundSliceIndex(const typename TImageType* image, 
											 int & lowestSliceIndex,
											 int & highestSliceIndex,
											 double backgroundValue = 0)
	{
		typedef itk::ImageSliceConstIteratorWithIndex< TImageType > ConstIteratorType;
		ConstIteratorType it( image, image->GetLargestPossibleRegion() );
		it.GoToBegin();
		it.SetFirstDirection( 0 );  // 0=x, 1=y, 2=z
		it.SetSecondDirection( 1 ); // 0=x, 1=y, 2=z

		lowestSliceIndex = 9999;
		highestSliceIndex = -9999;

		while( !it.IsAtEnd() )
		{
			while( !it.IsAtEndOfSlice() )
			{
				while( !it.IsAtEndOfLine() )
				{
					TImageType::IndexType ind = it.GetIndex();
					TImageType::PixelType val = it.Get();

					if ( val != backgroundValue )
					{
						if( ind[2]>highestSliceIndex )
						{
							highestSliceIndex = ind[2];
						}
						if( ind[2]<lowestSliceIndex )
						{
							lowestSliceIndex = ind[2];
						}
						it.NextSlice();
					}
					else
					{
						++it;
					}
				}
				it.NextLine();
			}
			it.NextSlice();
		}
	}


	/************************************************************************/
	/* 得到最小包围盒                                                       */
	/************************************************************************/
	template<typename TImageType>
	void
		getBoundRegion(const typename TImageType* image, 
		typename TImageType::RegionType& region, 
		typename TImageType::PixelType backgroundValue = 0,
		double padding = 0)
	{
		kmStaticImageMacro(TImageType);
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);
		typedef itk::ImageRegionConstIteratorWithIndex<TImageType> IteratorType;
		IteratorType it( image, image->GetLargestPossibleRegion() );
		it.GoToBegin();
		IndexType startIndex;
		startIndex.Fill(9999);
		IndexType endIndex;
		endIndex.Fill(0);
		while(!it.IsAtEnd())
		{
			PixelType val = it.Get();

			if(val !=backgroundValue )
			{
				IndexType idx = it.GetIndex();
				for(unsigned int i=0;i<Dimension;i++)
				{
					if(idx[i]<startIndex[i])
					{
						startIndex[i] = idx[i];
					}
					
					if(idx[i]>endIndex[i])
					{
						endIndex[i] = idx[i];
					}
				}
			}

			it++;
		}

// 		std::cout<<startIndex<<std::endl;
// 		std::cout<<endIndex<<std::endl;

		SpacingType spacing = image->GetSpacing();
		SizeType size;
		for(unsigned int i=0;i<Dimension;i++)
		{
			startIndex[i] -= static_cast<int>(padding/spacing[i]);
			endIndex[i] += static_cast<int>(padding/spacing[i]);
		}

		region.SetIndex( startIndex );
		region.SetUpperIndex( endIndex );

		region.Crop( image->GetLargestPossibleRegion() );

		return ;
	}

	template<typename TMeshType, typename TImageType>
	void
		getBoundRegion(
		const typename TMeshType* mesh,
		const typename TImageType* referenceimage, 
		typename TImageType::RegionType& region,
		double padding = 0.0)
	{
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);
		typedef TMeshType::PointsContainerConstPointer     PointsContainer;
		typedef TMeshType::PointsContainerConstIterator    PointsIterator;

		TMeshType::PointType minPoint;
		minPoint.Fill(9999);
		TMeshType::PointType maxPoint;
		maxPoint.Fill(-9999);

		PointsIterator pointItr = mesh->GetPoints()->Begin();
		PointsIterator pointEnd = mesh->GetPoints()->End();

		while ( pointItr != pointEnd )
		{
			TMeshType::PointType pt = pointItr.Value();

			for(int i=0;i<Dimension;i++)
			{
				if(pt[i]<minPoint[i])
				{
					minPoint[i] = pt[i];
				}

				if(pt[i]>maxPoint[i])
				{
					maxPoint[i] = pt[i];
				}
			}

			++pointItr;
		}

		TImageType::PointType startPoint;
		TImageType::PointType endPoint;
		TImageType::IndexType startIndex;
		TImageType::IndexType endIndex;
		TImageType::SizeType  boundSize;
		for (int i=0;i<Dimension;i++)
		{
			startPoint[i] = (minPoint[i]-padding);
			endPoint[i]   = (maxPoint[i]+padding);
		}

		referenceimage->TransformPhysicalPointToIndex( startPoint, startIndex );
		referenceimage->TransformPhysicalPointToIndex( endPoint, endIndex );

		region.SetIndex( startIndex );
		region.SetUpperIndex( endIndex );
	}

	template<typename TImageType>
	void
		boundRegions( 
		typename TImageType* image, 
		typename TImageType::RegionType region1,
		typename TImageType::RegionType region2,
		typename TImageType::RegionType regionResult)
	{
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);
		TImageType::IndexType startIndex1 = region1.GetIndex();
		TImageType::IndexType endIndex1 = region1.GetUpperIndex();
		TImageType::IndexType startIndex2 = region2.GetIndex();
		TImageType::IndexType endIndex2 = region2.GetUpperIndex();

		for(int i=0;i<Dimension;i++)
		{
			if(startIndex1[i]>startIndex2[i])
			{
				startIndex1[i] = startIndex2[i];
			}

			if(endIndex1[i]<endIndex2[i])
			{
				endIndex1[i] = endIndex2[i];
			}
		}

		regionResult.SetIndex( startIndex1 );
		regionResult.SetUpperIndex( endIndex1 );

		//regionResult.Crop( image->GetLargestPossibleRegion() );
	}


	template<typename TImageType>
	unsigned long
		calculateNumberOfPixels(const typename TImageType* image, 
		typename TImageType::PixelType lowerThreshold,
		typename TImageType::PixelType upperThreshold)
	{
		typedef itk::ImageRegionConstIteratorWithIndex<TImageType> IteratorType;
		IteratorType it( image, image->GetLargestPossibleRegion() );
		it.GoToBegin();

		unsigned long counts = 0;
		while(!it.IsAtEnd())
		{
			double val = it.Get();
			if(val>=lowerThreshold && val<=upperThreshold)
			{
				counts++;
			}
			++it;
		}
		
		return counts;
	}

	template<typename ImageType>
	typename ImageType::Pointer
		extractByBoundRegion(const typename ImageType* inputImage, typename ImageType::RegionType roi)
	{
		typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > ExtractFilterType;
		ExtractFilterType::Pointer extractfilter = ExtractFilterType::New();
		extractfilter->SetRegionOfInterest(roi);
		extractfilter->SetInput( inputImage );
		//extractfilter->SetDirectionCollapseToIdentity(); // This is required.
		extractfilter->Update();

		//extractfilter->UpdateLargestPossibleRegion();

		return extractfilter->GetOutput();
	}

	template<class TImageType>
	void
		hightlight(typename TImageType* image, typename TImageType::PointType pt, double hvalue, int radius = 0 )
	{
		TImageType::IndexType index;
		image->TransformPhysicalPointToIndex( pt, index );
		TImageType::SizeType rad;
		rad.Fill( 2*radius + 1 );

		TImageType::IndexType idx;
		for(int i=0;i<TImageType::ImageDimension;i++)
		{
			idx[i] = index[i] - radius;
		}

		TImageType::RegionType hRegion( idx, rad );
		hRegion.Crop( image->GetLargestPossibleRegion() );

		itk::ImageRegionIteratorWithIndex<TImageType> it( image, hRegion );
		it.GoToBegin();
		while(!it.IsAtEnd())
		{
			it.Set( hvalue );
			it++;
		}
	}

	template< typename ImageType, typename TransformType,typename DisplacementFieldType >
	void
		generateDisplacementFieldByTransform( const typename ImageType * fixedImage, typename TransformType* transform, typename DisplacementFieldType::Pointer & field )
	{
			//DisplacementFieldType::Pointer field = DisplacementFieldType::New();
			field->SetRegions( fixedImage->GetLargestPossibleRegion() );
			field->SetOrigin( fixedImage->GetOrigin() );
			field->SetSpacing( fixedImage->GetSpacing() );
			field->SetDirection( fixedImage->GetDirection() );
			field->Allocate();

			typedef itk::ImageRegionIterator< DisplacementFieldType > FieldIterator;
			FieldIterator fi( field, fixedImage->GetLargestPossibleRegion() );

			fi.GoToBegin();

			TransformType::InputPointType  fixedPoint;
			TransformType::OutputPointType movingPoint;
			DisplacementFieldType::IndexType index;

			DisplacementFieldType::PixelType displacement;

			while( ! fi.IsAtEnd() )
			{
				index = fi.GetIndex();
				field->TransformIndexToPhysicalPoint( index, fixedPoint );
				movingPoint = transform->TransformPoint( fixedPoint );
				displacement = movingPoint - fixedPoint;
				fi.Set( displacement );
				++fi;
			}

	}

	template<typename ImageType>
	typename ImageType::Pointer
		cloneImage (const typename ImageType* image)
	{
		typedef itk::ImageDuplicator<ImageType> ImageDuplicatorType;
		ImageDuplicatorType::Pointer filter = ImageDuplicatorType::New();
		filter->SetInputImage( image );
		filter->Update();

		return filter->GetOutput();
	}

	template<typename FilterType>
	typename FilterType::OutputImageType::Pointer
		watershedSeg( const typename FilterType::InputImageType * inputImage, double threshold, double level )
	{

		typedef itk::WatershedImageFilter<FilterType::InputImageType> WatershedImageFilterType;
		typedef WatershedImageFilterType::OutputImageType LabelImageType;
		WatershedImageFilterType::Pointer segFilter = WatershedImageFilterType::New();
		segFilter->SetInput( inputImage );
		segFilter->SetThreshold( 0.005 );
		segFilter->SetLevel( 0.5 );
		segFilter->Update();

		return segFilter->GetOutput();
	}

	template<typename ImageType>
	void
		extendBinaryByDist( typename ImageType::Pointer & inputImage, double dist )
	{

		typedef itk::Image<float, ImageType::ImageDimension> DistMapType;
		DistMapType::Pointer distmap = km::calculateDistanceMap<ImageType, DistMapType>( inputImage );

		itk::ImageRegionIteratorWithIndex<DistMapType> it( distmap, distmap->GetLargestPossibleRegion() );
		it.GoToBegin();

		while(!it.IsAtEnd())
		{
			if( it.Get() > dist )
			{
				inputImage->SetPixel( it.GetIndex(), 0 );
			}
			else
			{
				inputImage->SetPixel( it.GetIndex(), 1 );
			}

			++it;
		}
	}


	template<typename ImageType, typename LabelType>
	typename LabelType::Pointer
		generateLabelImage( const ImageType* inputImage, unsigned int & numberOfComponents)
	{
		typedef itk::ConnectedComponentImageFilter<ImageType, LabelType> LabelFilterType;
		LabelFilterType::Pointer labelFilter = LabelFilterType::New();
		labelFilter->SetInput( inputImage );
		//labelFilter->SetBackgroundValue( 0 );
		//labelFilter->SetFullyConnected(true);
		labelFilter->Update();

		numberOfComponents = labelFilter->GetObjectCount();

		return labelFilter->GetOutput();
	}

	template<typename ImageType>
	void
		thresholdOutImage(typename ImageType::Pointer & image, double thresholdLower, double thresholdUpper, double replaceValue)
	{
		itk::ImageRegionIterator<ImageType> it( image, image->GetLargestPossibleRegion() );
		it.GoToBegin();
		while(!it.IsAtEnd())
		{
			ImageType::PixelType val = it.Get();
			if (val>=thresholdLower && val<=thresholdUpper)
			{
				it.Set( replaceValue );
			}
			it++;
		}
	}

	template<typename ImageType>
	void
		thresholdLabels(typename ImageType::Pointer & image, 
										const std::vector<int> labels , 
										typename int insideValue, 
										typename int outsideValue )
	{
		if (labels.size() == 0)
		{
			image->FillBuffer(0);
			return;
		}

		itk::ImageRegionIterator<ImageType> it( image, image->GetLargestPossibleRegion() );
		it.GoToBegin();
		while(!it.IsAtEnd())
		{
			ImageType::PixelType val = it.Get();
		
			for (int i=0;i<labels.size();i++)
			{
				if (val == labels[i])
				{
					it.Set( insideValue );
				}
				else
				{
					it.Set( outsideValue );
				}
			}

			it++;
		}
	}

	template<class ImageType>
	typename ImageType::Pointer
		absImage(const typename ImageType* inputImage)
	{
		typedef itk::AbsImageFilter<ImageType, ImageType> AbsImageFilterType;
		AbsImageFilterType::Pointer filter = AbsImageFilterType::New();
		filter->SetInput( inputImage );
		filter->Update();

		return filter->GetOutput();
	}

	template<class ImageType>
	typename ImageType::Pointer
		powImage(const typename ImageType * inputImage, const double alpha)
	{
		typedef itk::PowImageFilter<ImageType> PowImageFilterType;
		PowImageFilterType::Pointer filter = PowImageFilterType::New();
		filter->SetInput1(inputImage);
		filter->SetConstant2(alpha);
		filter->Update();

		return filter->GetOutput();
	}

	template<class ImageType>
	typename ImageType::Pointer
		histogramMatch(const ImageType* inputImage, const ImageType* refImage)
	{
		typedef itk::HistogramMatchingImageFilter<ImageType,ImageType> HEFilterType;
		HEFilterType::Pointer IntensityEqualizeFilter = HEFilterType::New();
		IntensityEqualizeFilter->SetReferenceImage( refImage );
		IntensityEqualizeFilter->SetInput( inputImage );
		IntensityEqualizeFilter->SetNumberOfHistogramLevels( 300);
		IntensityEqualizeFilter->SetNumberOfMatchPoints( 15);
		IntensityEqualizeFilter->ThresholdAtMeanIntensityOn();
		IntensityEqualizeFilter->Update();

		return IntensityEqualizeFilter->GetOutput();
	}
}

#endif
