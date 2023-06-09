# Copyright (c) 2011-2017, NVIDIA CORPORATION. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of NVIDIA CORPORATION nor the names of its
#    contributors may be used to endorse or promote products derived
#    from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
# OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

cmake_minimum_required (VERSION 3.18)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../include" "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")

GET_FILENAME_COMPONENT(CMAKE_C_COMPILER_NAME "${CMAKE_C_COMPILER}" NAME)
IF(CMAKE_C_COMPILER_NAME MATCHES cl AND NOT CMAKE_C_COMPILER_NAME MATCHES clang)
  set(libs_all CUDA::cusparse CUDA::cusolver)
  set(dyn_libs amgxsh CUDA::cudart CUDA::cublas)
ELSE(CMAKE_C_COMPILER_NAME MATCHES cl AND NOT CMAKE_C_COMPILER_NAME MATCHES clang)
  set(libs_all rt dl CUDA::cusparse CUDA::cusolver)
  set(dyn_libs amgxsh rt dl CUDA::cudart CUDA::cublas)
ENDIF(CMAKE_C_COMPILER_NAME MATCHES cl AND NOT CMAKE_C_COMPILER_NAME MATCHES clang)

ADD_EXECUTABLE(eigensolver eigensolver.c)

target_link_libraries(eigensolver ${dyn_libs})

IF(WIN32)
    IF(MSVC)
      set_source_files_properties( eigensolver.c PROPERTIES LANGUAGE CXX)
      set_target_properties ( eigensolver PROPERTIES LINKER_LANGUAGE CXX )
    ENDIF(MSVC)
ELSE(WIN32)
  SET_SOURCE_FILES_PROPERTIES( eigensolver.c PROPERTIES COMPILE_FLAGS -std=c99 )
ENDIF(WIN32)

if(MPI_FOUND)
  ADD_EXECUTABLE(eigensolver_mpi eigensolver_mpi.c)

  IF(WIN32)
      IF(MSVC)
        set_source_files_properties( eigensolver_mpi.c PROPERTIES LANGUAGE CXX)
        set_target_properties ( eigensolver_mpi PROPERTIES COMPILE_FLAGS "${OpenMP_C_FLAGS}" LINKER_LANGUAGE CXX )
        set_target_properties ( eigensolver_mpi PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
        target_link_libraries(eigensolver_mpi ${MPI_C_LIBRARIES} ${dyn_libs})
      ENDIF(MSVC)
    ELSE(WIN32)
      SET_SOURCE_FILES_PROPERTIES( eigensolver_mpi.c PROPERTIES COMPILE_FLAGS -std=c99 )
      set_target_properties ( eigensolver_mpi PROPERTIES COMPILE_FLAGS "${OpenMP_C_FLAGS} -pthread" )
      set_target_properties ( eigensolver_mpi PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
      target_link_libraries(eigensolver_mpi amgxsh ${omp_lib} ${MPI_C_LIBRARIES} ${dyn_libs})
    ENDIF(WIN32) 
    
    if(MPI_COMPILE_FLAGS)
      set_target_properties(eigensolver_mpi PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    endif(MPI_COMPILE_FLAGS)

endif(MPI_FOUND)
