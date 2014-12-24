#include <iostream>
#include "itkSimplexMesh.h"
#include "kmCommon.h"
#include "vtkSmartPointer.h"
#include "vtkSurface.h"
#include "vtkIsotropicDiscreteRemeshing.h"

int main(int argc, char* argv[])
{
	const char* meshfile = "D:\\Workspace\\ASM\\projects\\LiverSegbyASM\\experiments\\modelFitting\\output_20140815\\initializedMesh.vtk";
	typedef itk::SimplexMesh<double,3> SimplexMeshType;
	SimplexMeshType::Pointer inputMesh = km::readMesh<SimplexMeshType>(meshfile);
	
	vtkSmartPointer<vtkPolyData> polydata = km::mesh2PolyData<SimplexMeshType>(inputMesh);
	//km::writePolyData("polydata.vtk", polydata);

	vtkSurface* surface=vtkSurface::New();
	surface->CreateFromPolyData(polydata);
	//surface->WriteToFile("surface.vtk");

	vtkIsotropicDiscreteRemeshing *remesh=vtkIsotropicDiscreteRemeshing::New();
	remesh->SetInput(surface);
	remesh->SetFileLoadSaveOption(0);
	remesh->SetNumberOfClusters( 18 );
	remesh->SetConsoleOutput(2);
	remesh->SetSubsamplingThreshold(10);
	remesh->GetMetric()->SetGradation(2);
	remesh->SetDisplay(1);
	//remesh->Remesh();
	std::cout<<"Now call ProcessClustering().."<<std::endl;
	remesh->ProcessClustering();

	//SimplexMeshType::PointDataContainerPointer pointDatas = inputMesh->GetPointData();
	//pointDatas->Reserve(inputMesh->GetNumberOfPoints());
	//for (int pointid = 0; pointid < remesh->GetNumberOfItems(); pointid++)
	//{
	//	int clusterLabel = remesh->GetClustering()->GetValue(pointid);
	//	pointDatas->InsertElement(pointid, clusterLabel);
	//}

	//km::writeMesh<SimplexMeshType>("labeledMesh.vtk", inputMesh);

	return 0;
}

