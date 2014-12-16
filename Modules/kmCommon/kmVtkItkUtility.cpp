#ifndef __kmVtkItkUtility_cpp
#define __kmVtkItkUtility_cpp

#include "kmVtkItkUtility.h"
namespace km
{
	//-------------------------------------------------------//
	//                      VTK related                      //
	//-------------------------------------------------------//
	//Calculate the normals of a PolyData
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
	
	//Smooth a PolyData 
	vtkSmartPointer<vtkPolyData>
		smoothPolyData( vtkPolyData * inputPolyData, int iterations, double relaxfactor )
	{
		vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
		smoothFilter->SetInput( inputPolyData );
		smoothFilter->SetNumberOfIterations( iterations );
		smoothFilter->BoundarySmoothingOn();
		smoothFilter->SetFeatureAngle( 45 );
		smoothFilter->SetEdgeAngle( 15 );
		smoothFilter->SetRelaxationFactor( 0.05 );
		smoothFilter->Update();
		return smoothFilter->GetOutput();
	}
	
	//Decimate a vtkpolydata
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
	
	//Read VTKPolyData
	vtkSmartPointer<vtkPolyData>
		readPolyData( const char* filename )
	{
		vtkSmartPointer<vtkPolyDataReader> inputPolyDataReader = vtkSmartPointer<vtkPolyDataReader>::New();
		inputPolyDataReader->SetFileName( filename );
		try{
			inputPolyDataReader->Update();
		}catch( itk::ExceptionObject & err ){
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return NULL;
		}
		return inputPolyDataReader->GetOutput();
	}
	
	//Read VTKPolyData
	vtkSmartPointer<vtkPolyData>
		readPolyData( const std::string filename )
	{
		return readPolyData(filename.c_str());
	}
	
	//Write VTKPolyData
	void writePolyData( const char* filename, vtkPolyData* polydata )
	{
		vtkSmartPointer<vtkPolyDataWriter> outputPolyDataWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
		outputPolyDataWriter->SetFileName( filename );
		outputPolyDataWriter->SetInput( polydata );
		try{
			outputPolyDataWriter->Update();
		}catch( itk::ExceptionObject & err ){
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
		}
	}
	
	//Write VTKPolyData
	void writePolyData( const std::string filename, vtkPolyData* polydata )
	{
		writePolyData(filename.c_str(), polydata);
	}
	
	//Write a vtk polydata
	void writePolyData(const char* dirname, const std::string filenameprefix, unsigned int index, const char* extension, vtkPolyData* polydata)
	{
		std::stringstream ss;
		ss  <<  dirname  <<  "/"  <<  filenameprefix   <<  "."  <<  index+1  <<  extension;
		km::writePolyData( ss.str().c_str(), polydata );
	}
}

#endif