cmake_minimum_required(VERSION 3.5)

project(FRET_Embedlab)
enable_language(Fortran)

set(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/cmake/math
    /usr/local/share/cmake/Modules
    )

#set(CMAKE_MODULE_PATH "/opt/cmake-3.14.0-rc1-Linux-x86_64/share/cmake-3.14/Modules/;${CMAKE_MODULE_PATH}")
option(ENABLE_AUTO_BLAS      "Enable CMake to autodetect BLAS"                  ON)
option(ENABLE_AUTO_LAPACK    "Enable CMake to autodetect LAPACK"                ON)
option(ENABLE_THREADED_MKL   "Enable OpenMP parallelization in MKL"             ON)
option(ENABLE_OMP            "Enable OpenMP parallelization"                    OFF)


if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -O3 -funroll-loops -Wall")
    set(bounds "-fbounds-check")
endif()
#if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
#    set(dialect "-free ")
#    set(bounds "-check bounds")
#endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
# 
#   Set standard flags 
# 
    set(CMAKE_Fortran_FLAGS "-fpp -O3 -W1 -xHost -no-wrap-margin -heap-arrays")
    #set(CMAKE_Fortran_FLAGS "-fpp -O3 -W1 -xHost -no-wrap-margin -i4 -g -traceback -check all -fp-stack-check")
# 
#   Enable openmp if requested (default) 
# 
    if(ENABLE_OMP)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -qopenmp -parallel"
            )
    endif()
#
endif()


if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
endif()


if(ENABLE_OMP)
   find_package(OpenMP COMPONENTS Fortran)
endif()
#find_package(BLAS REQUIRED)
#find_package(LAPACK REQUIRED)

## FOR THE TESTS
#include(ConfigTesting)

# Math detection
set(BLAS_LANG "Fortran")
set(LAPACK_LANG "Fortran")
set(MKL_COMPILER_BINDINGS ${CMAKE_Fortran_COMPILER_ID})
include(ConfigMath)


if(LAPACK_FOUND AND BLAS_FOUND)
   set(lapackblas_libraries ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

#
# Compile.
#
add_executable(FRET_Embedlab "") 
add_subdirectory(src) 
if(ENABLE_OMP)
   target_link_libraries(FRET_Embedlab 
     PRIVATE 
     ${LAPACK_LIBRARIES} 
     ${BLAS_LIBRARIES} 
     PUBLIC
     OpenMP::OpenMP_Fortran
   )
else()
   target_link_libraries(FRET_Embedlab 
     PRIVATE 
     ${LAPACK_LIBRARIES} 
     ${BLAS_LIBRARIES} 
   )
endif()
  
## THIS IS FOR THE TESTS
#file(COPY ${CMAKE_SOURCE_DIR}/tests DESTINATION ${CMAKE_BINARY_DIR})
#
#if(EXISTS ${CMAKE_SOURCE_DIR}/submodules/runtest/runtest)
#   file(COPY ${CMAKE_SOURCE_DIR}/submodules/runtest/runtest DESTINATION ${CMAKE_BINARY_DIR}/tests)
#endif()
