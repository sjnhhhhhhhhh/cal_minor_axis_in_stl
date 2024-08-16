#include <vtkSmartPointer.h>
#include <vtkRendererCollection.h>
#include <vtkPointPicker.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkActor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkProperty.h>
#include <vtkSTLReader.h>
#include <vtkCylinderSource.h>
#include <vtkCubeAxesActor.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAutoInit.h>
#include "MovableAxesRepresentation.h"

VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType)

int main(int, char *[])
{
    // 创建读取器
    vtkSmartPointer<vtkPolyData> input1 = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkSTLReader> reader1 = vtkSmartPointer<vtkSTLReader>::New();
    reader1->SetFileName("C:/code/extract3d/stl/2.stl");
    reader1->Update();
    input1->DeepCopy(reader1->GetOutput());

    // 创建映射器和演员
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(reader1->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetOpacity(0.7);

    // 创建渲染器、渲染窗口和交互器
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->SetWindowName("STL Model Display");
    renderWindow->AddRenderer(renderer);
    
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // 将演员添加到场景中
    renderer->AddActor(actor);
    renderer->SetBackground(0.1, 0.2, 0.3);

    // 添加可移动的坐标轴
    vtkMovableAxesRepresentation* axesRepresentation = new vtkMovableAxesRepresentation();
    for (int i = 0; i < 3; ++i) {
        renderer->AddActor(axesRepresentation->GetAxisCircle(i));
        renderer->AddActor(axesRepresentation->GetHandleAxis(i));
    }

    // 调整摄像机视角，使得模型和坐标轴都能正确显示
    renderer->ResetCamera();
    renderWindow->Render();
    renderWindowInteractor->Start();

    delete axesRepresentation; // 清理自定义对象

    return EXIT_SUCCESS;
}
