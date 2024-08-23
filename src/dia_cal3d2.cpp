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
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <Eigen/Dense>
#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkLineSource.h>
#include <vtkConvexHull2D.h>
#include <vtkDelaunay3D.h>
#include <vtkGeometryFilter.h>
#include <vtkConvexPointSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCommand.h>
#include <vtkCamera.h>
#include <chrono>

using Eigen::Vector2d;


VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType)

// 三维点结构体
struct Point {
    double x;
    double y;
    double z;
    Point() : x(0), y(0), z(0) {}
    Point(double x1, double y1, double z1) : x(x1), y(y1), z(z1) {}
};

// 马后炮，Point转换成vector3d
Eigen::Vector3d pointToVector3d(const Point& p) {
    return Eigen::Vector3d(p.x, p.y, p.z);
}

// 平面结构体
struct Plane {
    double A, B, C, D; // 平面方程的系数
};

Plane calculate_normal_plane(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2) {
    Eigen::Vector3d major_axis = p2 - p1;
    major_axis.normalize();

    Eigen::Vector3d normal;
    if (std::abs(major_axis.x()) > 1e-6) {
        normal = Eigen::Vector3d(-major_axis.y(), major_axis.x(), 0).normalized();
    } else {
        normal = Eigen::Vector3d(0, -major_axis.z(), major_axis.y()).normalized();
    }

    double A = normal.x();
    double B = normal.y();
    double C = normal.z();
    double D = -(A * p1.x() + B * p1.y() + C * p1.z());

    return Plane{A, B, C, D};
}

// 三维两点距离
double dis_cal(const Point &p1, const Point &p2) {
    double a = std::pow((p1.x - p2.x), 2) + std::pow((p1.y - p2.y), 2) + std::pow((p1.z - p2.z), 2);
    return std::sqrt(a);
}

// 提取stl点
std::vector<Point> extractpoints(vtkSmartPointer<vtkPolyData> data) {
    vtkSmartPointer<vtkPoints> points = data->GetPoints();

    std::vector<Point> pointsvector;
    vtkIdType num = points->GetNumberOfPoints();
    for (vtkIdType i = 0; i < num; i++) {
        double temp[3];
        points->GetPoint(i, temp);
        Point p(temp[0], temp[1], temp[2]);
        pointsvector.push_back(p);
    }
    return pointsvector;
}

// 计算长径
std::tuple<double, Point, Point> cal_majoraxis(const std::vector<Point>& pointsvector) {
    int sizeofvec = pointsvector.size();
    double max = 0;  // 初始化最大距离为0
    Point p1max;
    Point p2max;

    for (int i = 0; i < sizeofvec - 1; i++) {
        for (int j = i + 1; j < sizeofvec; j++) {
            double dis = dis_cal(pointsvector[i], pointsvector[j]);
            if (dis > max) {
                p1max = pointsvector[i];
                p2max = pointsvector[j];
                max = dis;
            }
        }
    }
    return std::tuple<double, Point, Point>(max, p1max, p2max);
}

// 计算短径
std::tuple<double,Point,Point> cal_minor(const vtkSmartPointer<vtkPolyData>& cutData)
{   
    std::vector<Point> Points;
    vtkSmartPointer<vtkPoints> points = cutData->GetPoints();
    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i){
        double point[3];
        points->GetPoint(i,point);
        Point point1(point[0],point[1],point[2]);
        Points.push_back(point1);
    }
    auto [minor,p3,p4] = cal_majoraxis(Points);
    return std::tuple<double,Point,Point> (minor,p3,p4);
}

// 输出径信息
void output1(const double &dis_long, const Point &p1, const Point &p2) {
    std::cout << "point1:(" << p1.x << "," << p1.y << "," << p1.z << ")\n";
    std::cout << "point2:(" << p2.x << "," << p2.y << "," << p2.z << ")\n";
    std::cout << "distance:" << dis_long << std::endl;
}

// 输出径信息
void output2(const double &dis_long, const Point &p1, const Point &p2) {
    std::cout << "point3:(" << p1.x << "," << p1.y << "," << p1.z << ")\n";
    std::cout << "point4:(" << p2.x << "," << p2.y << "," << p2.z << ")\n";
    std::cout << "distance:" << dis_long << std::endl;
}

// 切割
vtkSmartPointer<vtkPolyData> cutWithPlane(const vtkSmartPointer<vtkPolyData>& inputData, const Eigen::Vector3d& pointOnPlane, const Eigen::Vector3d& normal) {
    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    plane->SetOrigin(pointOnPlane.x(), pointOnPlane.y(), pointOnPlane.z());
    plane->SetNormal(normal.x(), normal.y(), normal.z());

    vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
    cutter->SetCutFunction(plane);
    cutter->SetInputData(inputData);
    cutter->Update();

    return cutter->GetOutput();
}

// 优化切割+短径计算
vtkSmartPointer<vtkPolyData> improved_cut(const vtkSmartPointer<vtkPolyData>& inputData, const Eigen::Vector3d& pointOnPlane, const Eigen::Vector3d& normal, double* mbounds,const double& bounds_size) {
    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    plane->SetOrigin(pointOnPlane.x(), pointOnPlane.y(), pointOnPlane.z());
    plane->SetNormal(normal.x(), normal.y(), normal.z());

    vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
    cutter->SetCutFunction(plane);
    cutter->SetInputData(inputData);
    cutter->Update();
    // 提取切片包围盒
    vtkSmartPointer<vtkPolyData> cutData = cutter->GetOutput();

    
    




    return cutter->GetOutput();
}


Eigen::Vector3d cal_normal(Point p1,Point p2,Point p3,Point p4)
{
    Eigen::Vector3d p1_vec = pointToVector3d(p1);
    Eigen::Vector3d p2_vec = pointToVector3d(p2);
    Eigen::Vector3d p3_vec = pointToVector3d(p3);
    Eigen::Vector3d p4_vec = pointToVector3d(p4);
    Eigen::Vector3d orient_major = p2_vec - p1_vec;
    Eigen::Vector3d orient_minor = p4_vec - p3_vec;
    Eigen::Vector3d normal = orient_major.cross(orient_minor);
    normal.normalize();

    return normal;
}


int main(int, char*[]) {
    // 获取开始时间点
    auto start = std::chrono::high_resolution_clock::now();
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName("C:/code/extract3d/stl/3.stl");
    reader->Update();

    vtkSmartPointer<vtkPolyData> data = reader->GetOutput();

    // 提取stl数据集所有点
    auto pointsvector = extractpoints(data);

    // 计算长径
    auto [dis_long, p1, p2] = cal_majoraxis(pointsvector);

    /*                         --以上没有问题                               --*/                  

    // 转换长径端点
    Eigen::Vector3d p1_vec = pointToVector3d(p1);
    Eigen::Vector3d p2_vec = pointToVector3d(p2);

    // 计算长径的法平面
    Plane normal_plane = calculate_normal_plane(p1_vec, p2_vec);

    // 切割
    Eigen::Vector3d major_axis = p2_vec - p1_vec;
    major_axis.normalize();
    double step = 0.1;
    double max_minor = 0;
    Point p3max;
    Point p4max;
    std::vector<double> max_t2_list;
    Eigen::Vector3d vecQ1,vecQ2;
    double bounds_size[3];


    for (double t = 0.0; t <= (p2_vec - p1_vec).norm(); t += step) {

        Eigen::Vector3d pointOnPlane = p1_vec + t * major_axis;

        // 切割平面法向量设置为长径向量的法向量
        Eigen::Vector3d normal = major_axis;  // 这里使用长径的方向作为法向量

        vtkSmartPointer<vtkPolyData> cutData = cutWithPlane(data, pointOnPlane, normal);

        double bounds[6];
        cutData->GetBounds(bounds);
        double size_x = bounds[1] - bounds[0];
        double size_y = bounds[3] - bounds[2];
        double size_z = bounds[5] - bounds[4];
        //如果切片整体大于目前最大切片，就更新最大切片
        //如果切片整体小于最大切片，那就直接跳过
        //其他情况不变
        if (size_x > bounds_size[0] && size_y > bounds_size[1] && size_z > bounds_size[2]){
            bounds_size[0] = size_x;
            bounds_size[1] = size_y;
            bounds_size[2] = size_z;
        }
        else if (size_x < bounds_size[0] && size_y < bounds_size[1] && size_z < bounds_size[2]){
            continue;
        }

//2047529

        auto [minor,p3,p4] = cal_minor(cutData);


        if (minor > max_minor) {
            max_minor = minor;
            p3max = p3;
            p4max = p4;
        }
        std::cout<<"max_minor="<<max_minor<<"\n";
    }

    output1(dis_long, p1, p2);
    output2(max_minor, p3max, p4max);
    auto normal2 = cal_normal(p1,p2,p3max, p4max);
    // 输出文件路径
    std::string output_file_path = "C:/code/extract3d/src/result_3d.txt";
    std::ofstream output_file(output_file_path);
    if (!output_file.is_open()) {
        std::cerr << "Error opening output file: " << output_file_path << std::endl;
        return 1;
    }
    output_file << "normal:" << "(" << normal2.x() << "," << normal2.y() << "," << normal2.z() << ")" << "\n"
                << "major_axis p1:" << "(" << p1.x << "," << p1.y << "," << p1.z << ")" << "\n"
                << "major_axis p2:" << "(" << p2.x << "," << p2.y << "," << p2.z << ")" << "\n"
                << "minor_axis p3:" << "(" << p3max.x << "," << p3max.y << "," << p3max.z << ")" << "\n"
                << "minor_axis p4:" << "(" << p4max.x << "," << p4max.y << "," << p4max.z << ")" << "\n" ;
    // 获取结束时间点
    auto end = std::chrono::high_resolution_clock::now();

    // 计算耗时（以微秒为单位）
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // 打印耗时
    std::cout << "Elapsed time: " << duration.count() << " microseconds" << std::endl;

    // 创建映射器和演员
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(reader->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetOpacity(0.7);

    // 创建渲染器、渲染窗口和交互器
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->SetWindowName("本程序基于一个离谱的Bug运行");
    renderWindow->AddRenderer(renderer);
    
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // 添加STL模型演员到场景中
    renderer->AddActor(actor);
    renderer->SetBackground(0.1, 0.2, 0.3);

    // 创建长径的线段
    vtkSmartPointer<vtkLineSource> majorLine = vtkSmartPointer<vtkLineSource>::New();
    majorLine->SetPoint1(p1.x, p1.y, p1.z);
    majorLine->SetPoint2(p2.x, p2.y, p2.z);
    
    vtkSmartPointer<vtkPolyDataMapper> majorLineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    majorLineMapper->SetInputConnection(majorLine->GetOutputPort());
    
    vtkSmartPointer<vtkActor> majorLineActor = vtkSmartPointer<vtkActor>::New();
    majorLineActor->SetMapper(majorLineMapper);
    majorLineActor->GetProperty()->SetColor(1.0, 0.0, 0.0); // 红色

    // 创建短径的线段
    vtkSmartPointer<vtkLineSource> minorLine = vtkSmartPointer<vtkLineSource>::New();
    minorLine->SetPoint1(p3max.x, p3max.y, p3max.z);
    minorLine->SetPoint2(p4max.x, p4max.y, p4max.z);
    
    vtkSmartPointer<vtkPolyDataMapper> minorLineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    minorLineMapper->SetInputConnection(minorLine->GetOutputPort());
    
    vtkSmartPointer<vtkActor> minorLineActor = vtkSmartPointer<vtkActor>::New();
    minorLineActor->SetMapper(minorLineMapper);
    minorLineActor->GetProperty()->SetColor(0.0, 1.0, 0.0); // 绿色

    // 添加长短径演员到渲染器
    renderer->AddActor(majorLineActor);
    renderer->AddActor(minorLineActor);

    // 调整摄像机视角，使得模型和坐标轴都能正确显示
    //renderer->ResetCamera();
    renderWindow->Render();
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}