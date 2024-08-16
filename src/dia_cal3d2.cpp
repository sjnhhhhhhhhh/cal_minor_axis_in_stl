/*
短径算法仍然有问题，主要是延长倍数t2计算有问题，但是在t2+0.5后很离奇的正确了，并且两个stl文件都适用
短径计算我先后尝试了三种办法：
    （1） 切片后直接找切片内部最大距离，但我很快发现这种办法找到的线段不与长径相交，原理上就有问题
    （2） 第二个办法，我想着切片后对每个切片单独运行2d凸包短径算法，但是遇到了如下问题
        1.cutter生成的polydata里面包含的points全是三维的，这导致我们必须投影再套算法 
        2.长径作为法向量并不平行于任一坐标轴，因此投影后会造成无法映射回三维造成信息丢失，但我们可以建立映射表
        3.这个问题没能解决，2d情况下我们也知道，短径计算时大部分情况无法做到点对点，必须计算交点，但是新产生的交点不在映射表中，我们只能找一个最近的投影后点映射回去，这又导致了误差，短径无法与长径相交
    （3） 直接建立三维凸壳，由于我们的cutter数据虽然是三维的，但是其本质上只是一个位于三维坐标系的二维平面，因此可以计算，只不过向量计算比较复杂，也因此遇到了目前的问题
        1.计算的短径绝对是对的，但是没有正确延伸到另一边，而是与长径相交后就不延伸了，这一定有蹊跷
        2.在t2+0.5后很离奇的正确了，并且两个stl文件都适用，我目前正在debug


*/



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

// 切割并计算最长距离（短径）
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

std::vector<Point> convex_hull(const vtkSmartPointer<vtkPolyData>& cutData) {
    vtkSmartPointer<vtkPoints> points = cutData->GetPoints();
    if (!points || points->GetNumberOfPoints() < 4) {
        std::cerr << "Not enough points to compute a 3D convex hull." << std::endl;
        return std::vector<Point>();  // 返回一个空的点集
    }

    // 创建一个 ConvexPointSet
    vtkSmartPointer<vtkConvexPointSet> convexPointSet = vtkSmartPointer<vtkConvexPointSet>::New();
    
    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
        convexPointSet->GetPointIds()->InsertNextId(i);
    }

    vtkSmartPointer<vtkUnstructuredGrid> convexHull = vtkSmartPointer<vtkUnstructuredGrid>::New();
    convexHull->Allocate(1, 1);
    convexHull->InsertNextCell(convexPointSet->GetCellType(), convexPointSet->GetPointIds());
    convexHull->SetPoints(points);

    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(convexHull->GetPoints());

    // 获取凸包点集
    std::vector<Point> pointsOnHull;
    for (vtkIdType i = 0; i < polyData->GetNumberOfPoints(); ++i) {
        double p[3];
        polyData->GetPoint(i, p);
        pointsOnHull.emplace_back(p[0], p[1], p[2]);
    }

    return pointsOnHull;
}



double compute_shortest_axis_with_pointOnPlane(const std::vector<Point>& pointsOnHull, const Eigen::Vector3d& pointOnPlane, Point& short_p1, Point& short_p2) {
    double max_distance = 0.0;
    double max_t2;
    for (const auto& p : pointsOnHull) {
        Eigen::Vector3d vecP(p.x, p.y, p.z);//起点
        //Eigen::Vector3d vecPointOnPlane(pointOnPlane.x, pointOnPlane.y, pointOnPlane.z);
        Eigen::Vector3d direction = pointOnPlane - vecP;//线段方向（点到Pointonplane）

        // 遍历凸包的所有边，找出与延长线的交点
        for (size_t i = 0; i < pointsOnHull.size(); ++i) {
            const Point& q1 = pointsOnHull[i];
            const Point& q2 = pointsOnHull[(i + 1) % pointsOnHull.size()];  // Next point (circular)

            Eigen::Vector3d vecQ1(q1.x, q1.y, q1.z);
            Eigen::Vector3d vecQ2(q2.x, q2.y, q2.z);

            Eigen::Vector3d segmentDir = vecQ2 - vecQ1; //边向量

            Eigen::Vector3d crossProduct = direction.cross(segmentDir); //计算二者叉积

            // 如果方向与边平行（即叉积接近零），跳过
            if (crossProduct.norm() < 1e-6) continue;

            // 计算延长线与边的交点参数t1和t2

            //初步确认t2计算有无导致延长线没能正确抵达另一边，停在了长径上
            double t1 = (vecQ1 - pointOnPlane).cross(segmentDir).dot(crossProduct) / crossProduct.squaredNorm();//延长线是否与边产生交点
            double t2 = (vecQ1 - vecP).cross(direction).dot(crossProduct) / crossProduct.squaredNorm();//交点是否在延长线上
            //std::cout<<"start point: "<< vecP <<"\n";
            //if(t2>max_t2) max_t2 = t2;
            // 如果t1在[0,1]范围内，表示延长线与边相交
            if (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1) {
                Eigen::Vector3d intersection = vecP + (t2+0.5) * direction;
                double dist = (intersection - vecP).norm();
                if (dist > max_distance) {
                    max_distance = dist;
                    short_p1 = p;
                    short_p2 = Point(intersection.x(), intersection.y(), intersection.z());
                }
            }
        }
        
    }
    std::cout<<"start_point:"<<"("<<short_p1.x<<","<<short_p1.y<<","<<short_p1.z<<")\n";
    std::cout<<"inter_point:"<<"("<<short_p2.x<<","<<short_p2.y<<","<<short_p2.z<<")\n";
    //if (max_t2>=0.1) std::cout<<"max_t2="<<max_t2<<"\n";
    
    return max_distance;
}


int main(int, char*[]) {
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName("C:/code/extract3d/stl/3.stl");
    reader->Update();

    vtkSmartPointer<vtkPolyData> data = reader->GetOutput();

    // 提取stl数据集所有点
    auto pointsvector = extractpoints(data);

    // 计算长径
    auto [dis_long, p1, p2] = cal_majoraxis(pointsvector);

    

    // 转换长径
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

    for (double t = 0.0; t <= (p2_vec - p1_vec).norm(); t += step) {

        Eigen::Vector3d pointOnPlane = p1_vec + t * major_axis;

        // 切割平面法向量设置为长径向量的法向量
        Eigen::Vector3d normal = major_axis;  // 这里使用长径的方向作为法向量

        vtkSmartPointer<vtkPolyData> cutData = cutWithPlane(data, pointOnPlane, normal);

        auto pointsOnHull = convex_hull(cutData);

        Point p3,p4;

        auto minor = compute_shortest_axis_with_pointOnPlane(   pointsOnHull, 
                                                                pointOnPlane, 
                                                                p3, 
                                                                p4);

        if (minor > max_minor) {
            max_minor = minor;
            p3max = p3;
            p4max = p4;
        }
        std::cout<<"max_minor="<<max_minor<<"\n";
    }

    output1(dis_long, p1, p2);
    output2(max_minor, p3max, p4max);

    // 创建映射器和演员
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(reader->GetOutputPort());

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
    renderer->ResetCamera();
    renderWindow->Render();
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}