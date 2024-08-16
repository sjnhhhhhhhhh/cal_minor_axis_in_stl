# CMakeList.txt: 顶层 CMake 项目文件，在此处执行全局配置
# 并包含子项目。
#
cmake_minimum_required (VERSION 3.8)

# 如果支持，请为 MSVC 编译器启用热重载。
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()


set(CMAKE_CXX_STANDARD 17)  
if(CMAKE_BUILD_TYPE MATCHES Debug)  
    message(STATUS "Building in Debug mode.")  
    # 在这里添加Debug模式下的特定设置或命令  
    set(VTK_DIR C:/source/VTK/lib/cmake/vtk-9.3)
    #set(ITK_DIR D:/openlib/ITK/lib/cmake/ITK-5.3)

elseif(CMAKE_BUILD_TYPE MATCHES Release)  
    message(STATUS "Building in Release mode.")  
    # 在这里添加Release模式下的特定设置或命令  
    set(VTK_DIR C:/source/VTK/lib/cmake/vtk-9.3)
    #set(ITK_DIR D:/openlib/ITK-R/lib/cmake/ITK-5.3)
else()  
    message(STATUS "Building in an unknown mode: ${CMAKE_BUILD_TYPE}")  
    # 处理其他可能的构建类型或未知构建类型的情况  
endif()
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
include_directories(${CMAKE_SOURCE_DIR})

include_directories(${PROJECT_SOURCE_DIR}/include)

# 设置 Eigen 库的路径（假设 Eigen 库已下载并解压到指定目录）
set(EIGEN3_INCLUDE_DIR "C:/source/eigen-3.4.0")
include_directories(${EIGEN3_INCLUDE_DIR})

project ("extract3d")

# 包含子项目。
add_subdirectory ("src")
