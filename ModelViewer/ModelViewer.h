#ifndef _GUI_h
#define _GUI_h

#include <iostream>

#include <QMainWindow>
#include "ui_ModelViewer.h"

#include "itkSimplexMesh.h"
#include "Representers/ITK/itkSimplexMeshRepresenter.h"
#include "statismo_ITK/itkStatisticalModel.h"
#include "statismo_ITK/itkStatisticalShapeModelTransform.h"
#include "itkMesh.h"
#include "Representers/ITK/itkMeshRepresenter.h"
#include "vtkPolyData.h"
#include "Representers/VTK/vtkPolyDataRepresenter.h"
#include "vtkSmartPointer.h"
#include <itkMesh.h>
#include <memory>
class vtkRenderer;
class vtkEventQtSlotConnect;
class vtkObject;
class vtkCommand;
class vtkActor;
class vtkPolyDataMapper;

using namespace statismo;

const unsigned int Dimension = 3;
typedef double MeshPixelType;
typedef itk::SimplexMeshRepresenter<MeshPixelType, Dimension> MeshRepresenterType;
typedef StatisticalModel<MeshRepresenterType>                 MeshStatisticalModelType;

typedef MeshRepresenterType::MeshType                         MeshType;
typedef MeshType::Pointer                                     MeshPointer;

typedef vtkPolyDataRepresenter                     PolyDataRepresenterType;
typedef StatisticalModel<PolyDataRepresenterType>  PolyDataStatisticalModelType;

class ModelViewer : public QMainWindow, public Ui::MainWindow
{
  Q_OBJECT
public:
  ModelViewer();
  ~ModelViewer();

public slots:
  //void updateCoords(vtkObject*);
  //void popup(vtkObject * obj, unsigned long, 
  //           void * client_data, void *,
  //           vtkCommand * command);
  void loadITKShapeModel();
  void loadVTKShapeModel();
  void exit();
  void setLamda1(int v);
  void setLamda2(int v);
  void setLamda3(int v);
  void setLamda4(int v);
	void setLamdaX(int v);
	void setLamdaXIndex(int v);
	void save();

public:
  enum ModelType {
    ITKShape = 0,
    VTKShape = 1,
    ITKIntensity = 2,
    VTKIntensity = 3
  };

protected:
  vtkSmartPointer<vtkRenderer>  m_Renderer;
  vtkSmartPointer<vtkActor> m_PolyDataActor;
  vtkSmartPointer<vtkPolyDataMapper > m_PolyDataMapper;

private:
  double m_LamdaValue1;
  double m_LamdaValue2;
  double m_LamdaValue3;
  double m_LamdaValue4;
	double m_LamdaValueX;
	unsigned int lamdaX;

  std::string m_ModelFileName;
  void getModelFileName();
  void updateModelSample();

  MeshStatisticalModelType* meshModel;
  PolyDataStatisticalModelType* polyDataModel;

  vtkSmartPointer<vtkPolyData> sampledPolyData;
	MeshPointer sampledMesh;


  ModelType m_CurrentModelType;
};

#endif // _GUI_h

