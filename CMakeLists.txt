project( elemesholve )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  
endif()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
find_package(Eigen3 REQUIRED)
message(STATUS "Found Eigen3 include dir '${EIGEN3_INCLUDE_DIR}'.")
 
# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()

find_package( TBB )
if( TBB_FOUND )
    include(${TBB_USE_FILE})
    list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DCGAL_CONCURRENT_MESH_3")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DCGAL_CONCURRENT_MESH_3")
endif()

# include for local directory

# include for local package
include_directories(BEFORE ../include $ENV{MPMUTILS}/Matrix/ $ENV{MPMUTILS}/Visualization/ $ENV{MPMUTILS}/GeneralUtils/ ${EIGEN3_INCLUDE_DIR} /usr/include/suitesparse/ )
link_directories($ENV{MPMUTILS}/Visualization/ $ENV{MPMUTILS}/GeneralUtils/ )

#################
# Documentation with Doxygen
#################
find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)


# Creating entries for target: elemesholve
# ############################

add_executable( elemesholve elemesholve.cc CalculationProcess.cc CellMatrix.cc EigenSparse.cc FEMesh2.cc FEMesh3.cc FEMesh3Slice.cc GeomDef.cc GeomSetup.cc MeshSlice.cc MeshVis.cc UmfSparse.cc )
target_compile_features(elemesholve PRIVATE cxx_auto_type)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS elemesholve )

# Link the executable to CGAL and third-party libraries
target_link_libraries(elemesholve ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} MPMGeneralUtils MPMVis GL glut blas pthread umfpack cholmod )

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_BUILD_TYPE debug)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3 -Wextra -Wpedantic -Woverloaded-virtual")
