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

// 计算二维向量的叉积
double cross(const Vector2d& a, const Vector2d& b) {
    return a.x() * b.y() - a.y() * b.x();
}

// 计算两个点之间的距离
double distance(const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
    return (a - b).norm();
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
/*vtkSmartPointer<vtkPolyData> cutWithPlane(const vtkSmartPointer<vtkPolyData>& inputData, const Eigen::Vector3d& pointOnPlane, const Eigen::Vector3d& normal) {
    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    plane->SetOrigin(pointOnPlane.x(), pointOnPlane.y(), pointOnPlane.z());
    plane->SetNormal(normal.x(), normal.y(), normal.z());

    vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
    cutter->SetCutFunction(plane);
    cutter->SetInputData(inputData);
    cutter->Update();

    return cutter->GetOutput();
}*/

// 计算直线与线段的交点
Vector2d line_segment_intersection(const Vector2d& a1, const Vector2d& a2, const Vector2d& b1, const Vector2d& b2, bool& intersects) {
    Vector2d r = a2 - a1; //代表线段a1a2的方向
    Vector2d s = b2 - b1; //代表线段b1b2的方向
    double rxs = cross(r, s);
    
    if (rxs == 0) {
        // 线段平行或共线
        intersects = false;
        return Vector2d(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());
    }

    double t = cross(b1 - a1, s) / rxs;
    double u = cross(b1 - a1, r) / rxs;

    if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
        intersects = true;
        return a1 + t * r; //如果存在交点，就返回交点
    } else {
        intersects = false;
        return Vector2d(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());
    }
}

// 分别求上下凸包并合成
std::vector<Eigen::Vector2d> compute_convex_hull(std::vector<Eigen::Vector2d> points) {
    std::sort(points.begin(), points.end(), [](const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
        return (p1.x() < p2.x()) || (p1.x() == p2.x() && p1.y() < p2.y());
    });

    std::vector<Eigen::Vector2d> hull;
    for (const auto& point : points) {
        while (hull.size() >= 2 && cross(hull[hull.size() - 1] - hull[hull.size() - 2], point - hull[hull.size() - 1]) <= 0) {
            hull.pop_back();
        }
        hull.push_back(point);
    }

    size_t lower_hull_size = hull.size();
    for (auto it = points.rbegin(); it != points.rend(); ++it) {
        while (hull.size() > lower_hull_size && cross(hull[hull.size() - 1] - hull[hull.size() - 2], *it - hull[hull.size() - 1]) <= 0) {
            hull.pop_back();
        }
        hull.push_back(*it);
    }

    hull.pop_back(); // 删除最后一个点，因为它与第一个点相同
    return hull;
}

double cal_minor_axis(const std::vector<Eigen::Vector2d>& hull, const Eigen::Vector2d& pointOnPlane2D, const std::vector<std::pair<Eigen::Vector2d, Eigen::Vector3d>>& projectedPoints, Eigen::Vector3d& short_p1, Eigen::Vector3d& short_p2) {
    double max_length = 0;
    std::cout << "-- start cal_minor --";
    for (const auto& point_pair : projectedPoints) {
        const Eigen::Vector2d& projectedPoint = point_pair.first;
        const Eigen::Vector3d& originalPoint = point_pair.second;

        // 沿着 pointOnPlane2D 和 projectedPoint 方向延伸线段
        for (size_t i = 0; i < hull.size(); ++i) {
            size_t next_i = (i + 1) % hull.size();

            // 计算延伸线段与凸包边界的交点
            bool intersects;
            Eigen::Vector2d intersection = line_segment_intersection(projectedPoint, projectedPoint + 100 * (pointOnPlane2D - projectedPoint), hull[i], hull[next_i], intersects);
            if (intersects && intersection != projectedPoint) {  // 确保交点与pointOnPlane2D不同
                double dist = (projectedPoint - intersection).norm();

                if (dist > max_length) {
                    max_length = dist;
                    short_p1 = point_pair.second; // 短径起点对应的三维点
                    // 找到与交点最近的三维点，由于要找最近的点，导致长短径仍然不相交，这个要改
                    Eigen::Vector3d closest_point = originalPoint;
                    double min_dist = (intersection - projectedPoint).norm();
                    for (const auto& candidate : projectedPoints) {
                        double candidate_dist = (intersection - candidate.first).norm();
                        if (candidate_dist < min_dist) {
                            closest_point = candidate.second;
                            min_dist = candidate_dist;
                        }
                    }
                    short_p2 = closest_point; // 短径终点对应的三维点
                }
            }
        }
    }

    return max_length;
}


Point vector3dToPoint(const Eigen::Vector3d& vec) {
    return Point(vec.x(), vec.y(), vec.z());
}


int main(int, char*[]) {
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName("C:/code/extract3d/stl/2.stl");
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
    Eigen::Vector3d p33d,p43d;

    for (double t = 0.0; t <= (p2_vec - p1_vec).norm(); t += step) {
        Eigen::Vector3d pointOnPlane = p1_vec + t * major_axis;
        //--计算切割平面--

        // 切割平面法向量设置为长径向量的法向量
        Eigen::Vector3d normal = major_axis;  // 这里使用长径的方向作为法向量

        vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
        plane->SetOrigin(pointOnPlane.x(), pointOnPlane.y(), pointOnPlane.z());
        plane->SetNormal(normal.x(), normal.y(), normal.z());
        

        vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
        cutter->SetCutFunction(plane);
        cutter->SetInputData(data);
        cutter->Update();

       vtkSmartPointer<vtkPolyData> cutData = cutter->GetOutput();
        Eigen::Vector3d u, v;

        // 假设 normal 不平行于 z 轴
        if (std::abs(normal.z()) > 1e-6) {
            u = normal.cross(Eigen::Vector3d(1, 0, 0)).normalized(); // normal 与 x 轴的叉积
        } else {
            u = normal.cross(Eigen::Vector3d(0, 1, 0)).normalized(); // normal 与 y 轴的叉积
        }

        v = normal.cross(u).normalized(); // v 与 u 及 normal 垂直

        //投影并保存投影前后点的映射关系

        std::vector<std::pair<Eigen::Vector2d, Eigen::Vector3d>> projectedPoints;
        vtkSmartPointer<vtkPoints> points = cutData->GetPoints();
        std::vector<Vector2d> points2d;

        for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
            double p[3];
            points->GetPoint(i, p);
            Eigen::Vector3d point(p[0], p[1], p[2]);

            // 计算投影
            double x_proj = (point - pointOnPlane).dot(u);
            double y_proj = (point - pointOnPlane).dot(v);
            Eigen::Vector2d projectedPoint(x_proj, y_proj);

            // 保存二维投影点和三维点的映射关系
            projectedPoints.emplace_back(projectedPoint, point);
            points2d.push_back(projectedPoint);
        }

        // 投影 pointOnPlane 到 u-v 平面
        double x_proj_onPlane = pointOnPlane.dot(u);
        double y_proj_onPlane = pointOnPlane.dot(v);

        Vector2d pointonplane_projected(x_proj_onPlane,y_proj_onPlane);

        auto hull = compute_convex_hull(points2d);

        Eigen::Vector3d short_p1,short_p2;
        auto temp_max = cal_minor_axis(hull,pointonplane_projected,projectedPoints,short_p1,short_p2);
        if (temp_max > max_minor)
        {
            max_minor = temp_max;
            std::cout<<"    current_max="<<max_minor<<"\n";
            p33d = short_p1;
            p43d = short_p2;
        }
        else if (temp_max <= max_minor){std::cout<<"     temp_max_too_low="<<temp_max<<"\n";} 




        /*std::vector<Point> pointsvector2 = extractpoints(cutData);

        if (pointsvector2.size() < 2) continue;

        auto [minor, p3, p4] = cal_majoraxis(pointsvector2);

        if (minor > max_minor) {
            max_minor = minor;
            p3max = p3;
            p4max = p4;
        }*/
    }
    //std::cout<<p33d<<"\n"<<p43d<<"\n";

    p3max = vector3dToPoint(p33d);
    p4max = vector3dToPoint(p43d);

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
