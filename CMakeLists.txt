cmake_minimum_required(VERSION 2.6)

# Enable C99
if(CMAKE_VERSION VERSION_LESS "3.1")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
else()
  set(CMAKE_C_STANDARD 99)
endif()

# Find OPENMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")

# include source code directory
set(SRC src)
include_directories(${SRC})

# include external code directory
set(EXTERNAL ${SRC}/external)
include_directories(${EXTERNAL})

# include NFFT
add_subdirectory(${EXTERNAL}/nfft-3.5.0)
include_directories(${NFFT3_INCPATH})
link_directories(   ${NFFT3_LIBPATH})

# set libraries
set(LIBS m jpeg png tiff fftw3 fftw3_threads)

# geometric transformation using NFFT
add_executable(homographic_transform ${SRC}/main_transform.c ${SRC}/homographic_transform.c ${SRC}/fft_core.c ${SRC}/homography_core.c ${EXTERNAL}/iio.c)
add_dependencies(homographic_transform nfft-3.5.0)
target_link_libraries(homographic_transform ${LIBS} nfft3_threads)

# shift interpolation
add_executable(shift_transform ${SRC}/main_shift.c ${SRC}/shift_transform.c ${SRC}/fft_core.c ${EXTERNAL}/iio.c)
target_link_libraries(shift_transform ${LIBS})

# upsampling
add_executable(upsampling ${SRC}/main_upsampling.c ${SRC}/resampling.c ${SRC}/fft_core.c ${EXTERNAL}/iio.c)
target_link_libraries(upsampling ${LIBS})

# downsampling
add_executable(downsampling ${SRC}/main_downsampling.c ${SRC}/resampling.c ${SRC}/fft_core.c ${EXTERNAL}/iio.c)
target_link_libraries(downsampling ${LIBS})
