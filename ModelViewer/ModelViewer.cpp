#include "ModelViewer.h"

#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkCommand.h"
#include "vtkEventQtSlotConnect.h"
#include "vtkConeSource.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkTDxInteractorStyleCamera.h"
#include "vtkTDxInteractorStyleSettings.h"

#include <QMenu>
#include "QVTKInteractor.h"
#include <QtGui>
#include <QtCore>
#include <QString>

#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkCamera.h>

#include "kmUtility.h"
#include "kmVtkItkUtility.h"


ModelViewer::ModelViewer():
m_LamdaValue1(0), 
m_LamdaValue2(0), 
m_LamdaValue3(0), 
m_LamdaValue4(0), 
m_LamdaValueX(0), 
lamdaX(1),
m_CurrentModelType(ModelType::VTKShape)
{
  this->setupUi(this); 

  const double angleSensitivity=0.02;
  const double translationSensitivity=0.001;

  QVTKInteractor *iren=qvtkWidget->GetInteractor();
  vtkInteractorStyle *s=
    static_cast<vtkInteractorStyle *>(iren->GetInteractorStyle());
  vtkTDxInteractorStyleCamera *t=
    static_cast<vtkTDxInteractorStyleCamera *>(s->GetTDxStyle());

  t->GetSettings()->SetAngleSensitivity(angleSensitivity);
  t->GetSettings()->SetTranslationXSensitivity(translationSensitivity);
  t->GetSettings()->SetTranslationYSensitivity(translationSensitivity);
  t->GetSettings()->SetTranslationZSensitivity(translationSensitivity);

  // add a renderer
  m_Renderer = vtkSmartPointer< vtkRenderer >::New();
  qvtkWidget->GetRenderWindow()->AddRenderer(m_Renderer);

  connect(horizontalSlider1, SIGNAL(valueChanged(int)), this, SLOT(setLamda1(int)));
  connect(horizontalSlider2, SIGNAL(valueChanged(int)), this, SLOT(setLamda2(int)));
  connect(horizontalSlider3, SIGNAL(valueChanged(int)), this, SLOT(setLamda3(int)));
  connect(horizontalSlider4, SIGNAL(valueChanged(int)), this, SLOT(setLamda4(int)));
	connect(horizontalSliderX, SIGNAL(valueChanged(int)), this, SLOT(setLamdaX(int)));
	connect(spinBox,           SIGNAL(valueChanged(int)), this, SLOT(setLamdaXIndex(int)));
	connect(pushButton,        SIGNAL(clicked(void)),     this, SLOT(save()));

  connect(actionOpenITKShapeModel, SIGNAL(triggered()),  this, SLOT(loadITKShapeModel()));
  connect(actionOpenVTKShapeModel, SIGNAL(triggered()),  this, SLOT(loadVTKShapeModel()));
  connect(actionExit,              SIGNAL(triggered()),  this, SLOT(exit()));

	sampledMesh = MeshType::New();
	sampledPolyData = vtkSmartPointer<vtkPolyData>::New();

}

ModelViewer::~ModelViewer()
{
	if (meshModel)
	{
		meshModel->Delete();
	}

	if (polyDataModel)
	{
		polyDataModel->Delete();
	}
}

void ModelViewer::getModelFileName()
{
  std::cout<<"load model file"<<std::endl;
  QString filter;
  filter = "Model file (*.h5)";

  //file path
  QDir dir;
  QString fileName = QFileDialog::getOpenFileName(this,tr("Open Model File"), dir.absolutePath() , filter );
  m_ModelFileName = fileName.toStdString();

  std::cout<<"Model Name: "<<m_ModelFileName<<std::endl;
}

void ModelViewer::loadITKShapeModel()
{
  m_CurrentModelType = ModelType::ITKShape;
  getModelFileName();

  if (m_ModelFileName.empty())
  {
	  return;
  }

  //meshModel = MeshStatisticalModelType::New();
  //meshModel->Load( m_ModelFileName.c_str() );

	meshModel = MeshStatisticalModelType::Load( m_ModelFileName );

	unsigned int numberOfPrincipalComponents = meshModel->GetNumberOfPrincipalComponents();
  std::cout << "loaded model with " << numberOfPrincipalComponents << " Principal Components" << std::endl;

	spinBox->setMaximum ( numberOfPrincipalComponents );
	spinBox->setMinimum ( 1 );
	spinBox->setValue( numberOfPrincipalComponents );

  //typedef MeshStatisticalModelType::VectorType VectorType;
  sampledMesh = meshModel->DrawMean();

  sampledPolyData = km::mesh2PolyData<MeshType>(sampledMesh);
	std::cout<<"number of points: "<<sampledPolyData->GetNumberOfPoints()<<std::endl;
	std::cout<<"number of cells: "<<sampledPolyData->GetNumberOfCells()<<std::endl;
	std::cout<<"number of verts: "<<sampledPolyData->GetNumberOfVerts()<<std::endl;

  m_PolyDataActor = vtkSmartPointer< vtkActor >::New();
  m_PolyDataMapper = vtkSmartPointer< vtkPolyDataMapper >::New();

  m_PolyDataMapper->SetInput(sampledPolyData);
  m_PolyDataActor->SetMapper(m_PolyDataMapper);

  m_Renderer->AddViewProp(m_PolyDataActor);  
  m_Renderer->ResetCamera();
  qvtkWidget->GetRenderWindow()->Render();
  qvtkWidget->update();
}

void ModelViewer::loadVTKShapeModel()
{
  m_CurrentModelType = ModelType::VTKShape;
  getModelFileName();

  if (m_ModelFileName.empty())
  {
	  return;
  }

  polyDataModel = PolyDataStatisticalModelType::Load( m_ModelFileName );

	unsigned int numberOfPrincipalComponents = polyDataModel->GetNumberOfPrincipalComponents();
	std::cout << "loaded model with " << numberOfPrincipalComponents << " Principal Components" << std::endl;

	spinBox->setMaximum ( numberOfPrincipalComponents );
	spinBox->setMinimum ( 1 );
	spinBox->setValue( numberOfPrincipalComponents );

  sampledPolyData = polyDataModel->DrawMean();
	std::cout<<"number of points: "<<sampledPolyData->GetNumberOfPoints()<<std::endl;
	std::cout<<"number of cells: "<<sampledPolyData->GetNumberOfCells()<<std::endl;
	std::cout<<"number of verts: "<<sampledPolyData->GetNumberOfVerts()<<std::endl;

  m_PolyDataActor = vtkSmartPointer< vtkActor >::New();
  m_PolyDataMapper = vtkSmartPointer< vtkPolyDataMapper >::New();

  m_PolyDataMapper->SetInput(sampledPolyData);
  m_PolyDataActor->SetMapper(m_PolyDataMapper);

  m_Renderer->AddViewProp(m_PolyDataActor);  
  m_Renderer->ResetCamera();
  qvtkWidget->GetRenderWindow()->Render();
  qvtkWidget->update();
}

void ModelViewer::exit()
{
  std::cout<<"exit"<<std::endl;
}

void ModelViewer::setLamda1(int v)
{
  m_LamdaValue1 = static_cast<double>(v) / 10.0; 
  updateModelSample();
}

void ModelViewer::setLamda2(int v)
{
  m_LamdaValue2 = static_cast<double>(v) / 10.0; 
  updateModelSample();
}

void ModelViewer::setLamda3(int v)
{
  m_LamdaValue3 = static_cast<double>(v) / 10.0; 
  updateModelSample();
}

void ModelViewer::setLamda4(int v)
{
  m_LamdaValue4 = static_cast<double>(v) / 10.0; 
  updateModelSample();
}

void ModelViewer::setLamdaX(int v)
{
	m_LamdaValue1 = 0;
	m_LamdaValue2 = 0;
	m_LamdaValue3 = 0;
	m_LamdaValue4 = 0;
	m_LamdaValueX = static_cast<double>(v) / 10.0; 
	updateModelSample();
}

void ModelViewer::setLamdaXIndex(int v)
{
	lamdaX = static_cast<unsigned int>(v);
}

void ModelViewer::updateModelSample()
{
  std::cout<<"Lamda1: "<<m_LamdaValue1<<",  ""Lamda2: "<<m_LamdaValue2<<",  ""Lamda3: "<<m_LamdaValue3<<",  ""Lamda4: "<<m_LamdaValue4<<std::endl;
  if( m_CurrentModelType == ModelType::VTKShape )
  {
		VectorType coefficients = VectorType::Zero(polyDataModel->GetNumberOfPrincipalComponents());
		//std::cout<<"here!"<<std::endl;

    coefficients(0) = m_LamdaValue1;
    coefficients(1) = m_LamdaValue2;
    coefficients(2) = m_LamdaValue3;
    coefficients(3) = m_LamdaValue4;
		coefficients(lamdaX-1) = m_LamdaValueX;
		//std::cout<<"here!"<<std::endl;
   /* vtkSmartPointer<vtkPolyData> */sampledPolyData = polyDataModel->DrawSample(coefficients);
		std::cout<<"number of points: "<<sampledPolyData->GetNumberOfPoints()<<std::endl;
		std::cout<<"number of cells: "<<sampledPolyData->GetNumberOfCells()<<std::endl;
		std::cout<<"number of verts: "<<sampledPolyData->GetNumberOfVerts()<<std::endl;
		//std::cout<<"here!"<<std::endl;
    m_PolyDataMapper->SetInput(sampledPolyData);
		//std::cout<<"here!"<<std::endl;
    qvtkWidget->update();
		//std::cout<<"here!"<<std::endl;
  }
  else if( m_CurrentModelType == ModelType::ITKShape )
	{
		//typedef MeshStatisticalModelType::VectorType VectorType;
		//VectorType coefficients;
		VectorType coefficients = VectorType::Zero(meshModel->GetNumberOfPrincipalComponents());
		//std::cout<<"here!"<<std::endl;
		std::cout<<coefficients<<std::endl;
    
    coefficients(0) = m_LamdaValue1;
    coefficients(1) = m_LamdaValue2;
    coefficients(2) = m_LamdaValue3;
    coefficients(3) = m_LamdaValue4;
		coefficients(lamdaX-1) = m_LamdaValueX;
		//std::cout<<"here!"<<std::endl;
    /*MeshType::Pointer */sampledMesh = meshModel->DrawSample(coefficients);
		//std::cout<<"here!"<<std::endl;
    /*vtkSmartPointer<vtkPolyData>*/ sampledPolyData = km::mesh2PolyData<MeshType>(sampledMesh);
		std::cout<<"number of points: "<<sampledPolyData->GetNumberOfPoints()<<std::endl;
		std::cout<<"number of cells: "<<sampledPolyData->GetNumberOfCells()<<std::endl;
		std::cout<<"number of verts: "<<sampledPolyData->GetNumberOfVerts()<<std::endl;
		//std::cout<<"here!"<<std::endl;
    m_PolyDataMapper->SetInput(sampledPolyData);
		//std::cout<<"here!"<<std::endl;
    qvtkWidget->update();
		//std::cout<<"here!"<<std::endl;
  }
}

void ModelViewer::save()
{
	QString filter;
	filter = "Polydata file (*.vtk)";
	QDir dir;
	QString fileName = QFileDialog::getSaveFileName(this,tr("Open Model File"), dir.absolutePath() , filter );
	std::string savefilename = fileName.toStdString();

	if( m_CurrentModelType == ModelType::VTKShape )
	{
		km::writePolyData( savefilename, sampledPolyData );
	}
	else if( m_CurrentModelType == ModelType::ITKShape )
	{
		km::writeMesh<MeshType>( savefilename, sampledMesh );
	}

	std::cout<<"Saved file name: "<<savefilename<<std::endl;
}