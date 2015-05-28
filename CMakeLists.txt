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

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()

# include for local directory

# include for local package
include_directories( BEFORE ../include $ENV{MPMUTILS}/Matrix/ $ENV{MPMUTILS}/Visualization/ )
link_directories($ENV{MPMUTILS}/Visualization/)

# Creating entries for target: elemesholve
# ############################

add_executable( elemesholve elemesholve.cc UmfSparse.cc FEMesh.cc )
target_compile_features(elemesholve PRIVATE cxx_auto_type)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS elemesholve )

# Link the executable to CGAL and third-party libraries
target_link_libraries(elemesholve ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} MPMVis GL glut blas umfpack )
