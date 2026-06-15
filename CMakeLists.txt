cmake_minimum_required(VERSION 3.24)
project(kraken_gpu LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_ARCHITECTURES 86) # Targeting RTX A5000 (Ampere)

# CUDA setup
find_package(CUDAToolkit REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

# CUDA flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -lineinfo")

include_directories(include ${ZLIB_INCLUDE_DIRS})

file(GLOB_RECURSE SOURCES src/*.cpp src/*.cu)

add_executable(kraken_gpu main.cu ${SOURCES})

target_link_libraries(kraken_gpu PRIVATE CUDA::cudart ${ZLIB_LIBRARIES} Threads::Threads)
