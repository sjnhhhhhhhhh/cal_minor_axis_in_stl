/*
短径算法仍然有问题，主要是延长倍数t2计算有问题
短径计算我先后尝试了三种办法：
    （1） 切片后直接找切片内部最大距离，但我很快发现这种办法找到的线段不与长径相交，原理上就有问题
    （2） 第二个办法，我想着切片后对每个切片单独运行2d凸包短径算法，但是遇到了如下问题
        1.cutter生成的polydata里面包含的points全是三维的，这导致我们必须投影再套算法 
        2.长径作为法向量并不平行于任一坐标轴，因此投影后会造成无法映射回三维造成信息丢失，但我们可以建立映射表
        3.这个问题没能解决，2d情况下我们也知道，短径计算时大部分情况无法做到点对点，必须计算交点，但是新产生的交点不在映射表中，我们只能找一个最近的投影后点映射回去，这又导致了误差，短径无法与长径相交
    （3） 直接建立三维凸壳，由于我们的cutter数据虽然是三维的，但是其本质上只是一个位于三维坐标系的二维平面，因此可以计算，只不过向量计算比较复杂，也因此遇到了目前的问题
        1.计算的短径绝对是对的，但是没有正确延伸到另一边，超出去了，但我的t2是解方程求得的表达式，不可能出错
        2.在限定t2范围到[0,1]，然后再t2+0.5后很离奇的正确了，最离奇的是直接设定1.5 x direction并不对，说明t2并不是简单的=1了，我目前正在debug

        问题修正1：凸包合成的其实是错的，现在以及按照顺时针生成了有序的切片边界点集
        问题修正2：t2和t1的逻辑顺序起初有一点点浪费内存，现在改进了
        问题修正3：交点的计算逻辑现在正确了，之前由于存在误差，导致短径并没有和边真正相交，这就导致短径怎么延伸都找不到交点，现在改成了2d找算法，相当于是重新定了一个相对坐标系，然后忽略z轴进行计算
        

575959  17905892
525441  15592840

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
//合成三维凸包
std::vector<Point> convex_hull(const vtkSmartPointer<vtkPolyData>& cutData) {
    vtkSmartPointer<vtkPoints> points = cutData->GetPoints();
    if (!points || points->GetNumberOfPoints() < 4) {
        std::cerr << "Not enough points to compute a 3D convex hull." << std::endl;
        return std::vector<Point>();  // 返回一个空的点集
    }

    // 使用 Delaunay3D 生成三维凸包
    vtkSmartPointer<vtkDelaunay3D> delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
    delaunay->SetInputData(cutData);
    delaunay->Update();

    vtkSmartPointer<vtkUnstructuredGrid> convexHull = delaunay->GetOutput();

    // 提取凸包点集
    std::vector<Point> pointsOnHull;
    vtkSmartPointer<vtkPoints> hullPoints = convexHull->GetPoints();
    for (vtkIdType i = 0; i < hullPoints->GetNumberOfPoints(); ++i) {
        double p[3];
        hullPoints->GetPoint(i, p);
        pointsOnHull.emplace_back(p[0], p[1], p[2]);
    }

    return pointsOnHull;
}

std::vector<Point> ordered_points(const vtkSmartPointer<vtkPolyData>& cutData)//返回顺时针排列的切片边界点集
{
    vtkSmartPointer<vtkPoints> points = cutData->GetPoints();
    vtkSmartPointer<vtkCellArray> lines = cutData->GetLines();

    // 使用 unordered_map 存储连接关系
    std::unordered_map<vtkIdType, std::vector<vtkIdType>> pointAdjacency;

    lines->InitTraversal();
    vtkIdType npts;
    vtkIdType const* pts;
    while (lines->GetNextCell(npts, pts)) {
        for (vtkIdType i = 0; i < npts - 1; i++) {
            pointAdjacency[pts[i]].push_back(pts[i+1]);
            pointAdjacency[pts[i+1]].push_back(pts[i]);
        }
    }

    // 从起点开始排序
    vtkIdType startPoint = pointAdjacency.begin()->first; // 从任意一个点开始
    std::vector<vtkIdType> orderedPoints;
    orderedPoints.push_back(startPoint);

    vtkIdType currentPoint = startPoint;
    vtkIdType previousPoint = -1;

    while (true) {
        const auto& neighbors = pointAdjacency[currentPoint];
        vtkIdType nextPoint = -1;

        for (auto neighbor : neighbors) {
            if (neighbor != previousPoint) {
                nextPoint = neighbor;
                break;
            }
        }

        if (nextPoint == -1 || nextPoint == startPoint) break;

        orderedPoints.push_back(nextPoint);
        previousPoint = currentPoint;
        currentPoint = nextPoint;
    }

    // 计算中心点
    double centroid[3] = {0.0, 0.0, 0.0};
    for (auto ptId : orderedPoints) {
        double p[3];
        points->GetPoint(ptId, p);
        centroid[0] += p[0];
        centroid[1] += p[1];
        centroid[2] += p[2];
    }
    centroid[0] /= orderedPoints.size();
    centroid[1] /= orderedPoints.size();
    centroid[2] /= orderedPoints.size();

    // 确定角度并排序
    std::sort(orderedPoints.begin(), orderedPoints.end(),
        [&points, &centroid](vtkIdType a, vtkIdType b) {
            double pa[3], pb[3];
            points->GetPoint(a, pa);
            points->GetPoint(b, pb);
            double angleA = atan2(pa[1] - centroid[1], pa[0] - centroid[0]);
            double angleB = atan2(pb[1] - centroid[1], pb[0] - centroid[0]);
            return angleA < angleB;  // 顺时针方向
        }
    );

    // 转换为 std::vector<Point> 返回
    std::vector<Point> result;
    for (auto ptId : orderedPoints) {
        double p[3];
        points->GetPoint(ptId, p);
        result.emplace_back(p[0], p[1], p[2]);
    }

    return result;
}




bool point_equal(Point p1, Eigen::Vector3d p2){
    Eigen::Vector3d p1_vec(p1.x,p1.y,p1.z);
    if (p1_vec.x() == p2.x() && p1_vec.y() == p2.y() && p1_vec.z() == p2.z()) return true;
    else return false;
}

std::tuple<bool, Eigen::Vector3d> line_segment_intersection_3D(const Eigen::Vector3d& P1, const Eigen::Vector3d& P2,
                                                               const Eigen::Vector3d& Q1, const Eigen::Vector3d& Q2) {
    Eigen::Vector3d r = P2 - P1;  // 线段P1P2的方向向量
    Eigen::Vector3d s = Q2 - Q1;  // 线段Q1Q2的方向向量
    Eigen::Vector3d PQ = Q1 - P1;

    // 构造矩阵A和向量B
    Eigen::Matrix2d A;
    A << r.x(), -s.x(),
         r.y(), -s.y();

    Eigen::Vector2d B(PQ.x(), PQ.y());

    // 解方程组，判断是否有交点
    if (A.determinant() != 0) {
        Eigen::Vector2d t_u = A.inverse() * B;
        double t = t_u.x();
        double u = t_u.y();

        if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
            Eigen::Vector3d intersection = P1 + t * r;
            return std::make_tuple(true, intersection);
        }
    }

    // 如果没有交点，返回false和一个无效的交点
    return std::make_tuple(false, Eigen::Vector3d(std::numeric_limits<double>::lowest(), 
                                                  std::numeric_limits<double>::lowest(), 
                                                  std::numeric_limits<double>::lowest()));
}



std::tuple<double,Eigen::Vector3d,Eigen::Vector3d> compute_shortest_axis_with_pointOnPlane(const std::vector<Point>& pointsOnSlice, const Eigen::Vector3d& pointOnPlane, Point& short_p1, Point& short_p2) {
    double max_distance = 0.0;
    Eigen::Vector3d mvecQ1,mvecQ2;
    //double max_t2 = 0.0;
    for (const auto& p : pointsOnSlice) {
        Eigen::Vector3d vecP(p.x, p.y, p.z);//起点
        //Eigen::Vector3d vecPointOnPlane(pointOnPlane.x, pointOnPlane.y, pointOnPlane.z);
        Eigen::Vector3d direction = pointOnPlane - vecP;//线段方向（点到Pointonplane）

        // 遍历凸包的所有边，找出与延长线的交点
        for (size_t i = 0; i < pointsOnSlice.size(); ++i) {
            const Point& q1 = pointsOnSlice[i];
            if (point_equal(q1 , vecP)) continue; //保证当前遍历边和起点无关
            const Point& q2 = pointsOnSlice[(i + 1) % pointsOnSlice.size()];  // Next point (circular)
            if (point_equal(q2 , vecP)) continue;

            Eigen::Vector3d vecQ1(q1.x, q1.y, q1.z);
            Eigen::Vector3d vecQ2(q2.x, q2.y, q2.z);

            Eigen::Vector3d segmentDir = vecQ2 - vecQ1; //边向量

            Eigen::Vector3d crossProduct = direction.cross(segmentDir); //计算二者叉积

            // 如果方向与边平行（即叉积接近零），跳过
            if (crossProduct.norm() < 0.1) continue;
            auto [intersects,vecinter] = line_segment_intersection_3D(vecP - 100*direction,pointOnPlane + 100*direction,vecQ1,vecQ2);
            if (intersects){
                double dist = (vecinter - vecP).norm();
                if (dist > max_distance) {
                    //std::cout << "t1 = " << t1 << "\n";
                    //std::cout << "crossproduct = " << crossProduct << "\n";
                    max_distance = dist;
                    short_p1 = p;
                    short_p2 = Point(vecinter.x(), vecinter.y(), vecinter.z());
                    mvecQ1 = vecQ1;
                    mvecQ2 = vecQ2;
                    
                
                }
            }

                
            
        }
        
    }

    //std::cout<<"start_point:"<<"("<<short_p1.x<<","<<short_p1.y<<","<<short_p1.z<<")\n";
    //std::cout<<"inter_point:"<<"("<<short_p2.x<<","<<short_p2.y<<","<<short_p2.z<<")\n";

    
    return std::make_tuple(max_distance, mvecQ1, mvecQ2);
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
    double bounds_size[3];
    Eigen::Vector3d vecQ1,vecQ2;


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

        auto points = ordered_points(cutData); //切片的顺时针有序点集
        //auto points_hull = convex_hull(cutData); 

        Point p3,p4;
        

        auto [minor,mvecQ1,mvecQ2] = compute_shortest_axis_with_pointOnPlane(   points, 
                                                                                pointOnPlane, 
                                                                                p3, 
                                                                                p4);

        

        if (minor > max_minor) {
            max_minor = minor;
            p3max = p3;
            p4max = p4;
            vecQ1 = mvecQ1;
            vecQ2 = mvecQ2;
        }
        std::cout<<"max_minor="<<max_minor<<"\n";
    }

    output1(dis_long, p1, p2);
    output2(max_minor, p3max, p4max);

   
    std::cout<<"vecQ1= "<<"("<<vecQ1.x()<<","<<vecQ1.y()<<","<<vecQ1.z()<<")\n";
    std::cout<<"vecQ2= "<<"("<<vecQ2.x()<<","<<vecQ2.y()<<","<<vecQ2.z()<<")\n";
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