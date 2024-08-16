#include <vtkSmartPointer.h>
#include <vtkRendererCollection.h>
#include <vtkSTLReader.h>
#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>

#include <nlohmann/json.hpp>
#include <fstream>

int main(int, char *[])
{
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName("C:/code/extract3d/stl/2.stl");
    reader->Update();

    double bounds[6];
    reader->GetOutput()->GetBounds(bounds);

    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
    cutter->SetCutFunction(plane);
    cutter->SetInputConnection(reader->GetOutputPort());

    nlohmann::json root;

    int sliceIndex = 0;
    int counter=0;
    for (double z = bounds[4]+0.5; z <= bounds[5]; z += 0.5)
    {
        counter += 1;
        plane->SetOrigin(0, 0, z);
        cutter->Update();

        vtkSmartPointer<vtkPolyData> cutData = cutter->GetOutput();
        vtkSmartPointer<vtkPoints> points = cutData->GetPoints();

        nlohmann::json sliceJson;
        sliceJson["z"] = z;
        sliceJson["id"] = counter;

        for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
        {
            double point[3];
            points->GetPoint(i, point);

            nlohmann::json pointJson = {
                {"x", point[0]},
                {"y", point[1]},
                {"z", point[2]}
            };

            sliceJson["contour"].push_back(pointJson);
        }

        root["slices"].push_back(sliceJson);
    }

    // 将JSON数据写入文件
    std::ofstream file("C:/code/extract3d/src/output.json");
    file << root.dump(4); // 4个空格缩进
    file.close();

    return EXIT_SUCCESS;
}
