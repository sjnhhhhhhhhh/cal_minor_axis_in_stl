# CMakeLists.txt: test1 的 CMake 项目，在此处包括源代码并定义项目特定的逻辑。

# 将源代码添加到此项目的可执行文件。
add_executable(extract "dia_cal3d2.cpp")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET extract PROPERTY CXX_STANDARD 20)
endif()

# 查找并链接VTK库
find_package(VTK REQUIRED)
include(${VTK_USE_FILE}) # 确保 VTK 配置正确
target_link_libraries(extract PRIVATE ${VTK_LIBRARIES})

# 在目标创建后使用 vtk_module_autoinit
if (VTK_VERSION VERSION_GREATER "9.00.0")
    vtk_module_autoinit(
        TARGETS extract
        MODULES ${VTK_LIBRARIES}
    )
endif ()

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${EIGEN3_INCLUDE_DIR})