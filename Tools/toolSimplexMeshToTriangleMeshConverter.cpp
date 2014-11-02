#include <itkTriangleMeshToSimplexMeshFilter.h>
#include <itkSimplexMesh.h>
#include "itkSimplexMeshAdaptTopologyFilter.h"

#include "kmUtility.h"
#include "kmVtkItkUtility.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

#include <vtkPLYWriter.h>
#include <vtkSTLWriter.h>
#include <vtkCellData.h>
#include "vtkSurface.h"
#include "vtkIsotropicDiscreteRemeshing.h"

using namespace std;
using namespace km;

const unsigned Dimensions = 3;
typedef itk::SimplexMesh<float, Dimensions> MeshType;

int main(int argc, char* argv[]) {

	if (argc < 4) {
		std::cout << "usage " << argv[0] << " inputMesh geoImage outputMesh " << std::endl;
		exit(-1);
	}

	itk::TimeProbesCollectorBase chronometer;
	itk::MemoryProbesCollectorBase memorymeter;
	memorymeter.Start( "converting" );
	chronometer.Start( "converting" );

	MeshType::Pointer inputMesh = readMesh<MeshType>( argv[1] );
	km::readSimplexMeshGeometryData<MeshType>(argv[2], inputMesh);

	MeshType::Pointer triMesh = km::simplexMeshToTriangleMesh<MeshType, MeshType>(inputMesh);
	km::writeMesh<MeshType>( "triMesh.vtk", triMesh );

	//Check orientation
	{
		vtkSurface *Mesh=vtkSurface::New();
		Mesh->CreateFromFile("triMesh.vtk");

		vtkPolyDataNormals *Normals=vtkPolyDataNormals::New();
		Normals->SetInput(Mesh);
		Normals->SplittingOff();
		Normals->FlipNormalsOn();
		Normals->Update();

		vtkPolyData *Output=Normals->GetOutput();
		Output->GetPointData()->SetScalars(0);
		Output->GetCellData()->SetScalars(0);

		vtkPolyDataWriter *Writer=vtkPolyDataWriter::New();
		Writer->SetInput(Output);
		Writer->SetFileName("good_orientation.vtk");
		Writer->Write();
		Writer->Delete();
		Normals->Delete();
		Mesh->Delete();
	}

	//Remesh
	{
		vtkSurface *Mesh=vtkSurface::New();
		Mesh->CreateFromFile("good_orientation.vtk");
		Mesh->GetCellData()->Initialize();
		Mesh->GetPointData()->Initialize();
		Mesh->DisplayMeshProperties();

		vtkIsotropicDiscreteRemeshing *Remesh=vtkIsotropicDiscreteRemeshing::New();
		Remesh->SetInput(Mesh);
		Remesh->SetFileLoadSaveOption(0);
		Remesh->SetNumberOfClusters( triMesh->GetNumberOfPoints() );
		Remesh->SetConsoleOutput(2);
		Remesh->SetSubsamplingThreshold(10);
		Remesh->GetMetric()->SetGradation(1);
		Remesh->Remesh();

		Remesh->GetOutput()->WriteToFile("FUCK.vtk");
		Remesh->Delete();
		Mesh->Delete();
	}

	triMesh = km::readMesh<MeshType>( "FUCK.vtk" );

	MeshType::Pointer simplexMesh = km::triangleMeshToSimplexMesh<MeshType, MeshType>( triMesh );
	km::writeMesh<MeshType>( "SimplexMesh.vtk", simplexMesh );

	chronometer.Stop( "converting" );
	memorymeter.Stop( "converting" );
	chronometer.Report( std::cout );
	memorymeter.Report( std::cout );

	return EXIT_SUCCESS;
}

