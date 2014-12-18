#include <iostream>
#include "itkSimplexMesh.h"
#include "kmCommon.h"
#include "vtkSmartPointer.h"
#include "vtkSurface.h"

int main(int argc, char* argv[])
{
	const char* meshfile = "D:\\Workspace\\ASM\\projects\\LiverSegbyASM\\experiments\\modelFitting\\output_20140815\\initializedMesh.vtk";
	typedef itk::SimplexMesh<double,3> SimplexMeshType;
	SimplexMeshType::Pointer inputMesh = km::readMesh<SimplexMeshType>(meshfile);
	
	vtkSmartPointer<vtkPolyData> polydata = km::mesh2PolyData<SimplexMeshType>(inputMesh);
	km::writePolyData("polydata.vtk", polydata);

	vtkSurface* surface=vtkSurface::New();
	surface->CreateFromPolyData(polydata);
	surface->WriteToFile("surface.vtk");
	return 0;
}

