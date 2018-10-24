INCLUDE(FindPkgConfig)
PKG_CHECK_MODULES(PC_DIGITALHF digitalhf)

FIND_PATH(
    DIGITALHF_INCLUDE_DIRS
    NAMES digitalhf/api.h
    HINTS $ENV{DIGITALHF_DIR}/include
        ${PC_DIGITALHF_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    DIGITALHF_LIBRARIES
    NAMES gnuradio-digitalhf
    HINTS $ENV{DIGITALHF_DIR}/lib
        ${PC_DIGITALHF_LIBDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(DIGITALHF DEFAULT_MSG DIGITALHF_LIBRARIES DIGITALHF_INCLUDE_DIRS)
MARK_AS_ADVANCED(DIGITALHF_LIBRARIES DIGITALHF_INCLUDE_DIRS)

