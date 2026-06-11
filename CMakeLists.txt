cmake_minimum_required(VERSION 3.16)

# ── Automatic version from git tags ────────────────────────────────────────────
# Create the first tag with: git tag v0.1.0
execute_process(
  COMMAND git describe --tags --always --dirty
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GAPSTATS_GIT_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
if(NOT GAPSTATS_GIT_VERSION OR GAPSTATS_GIT_VERSION STREQUAL "")
  set(GAPSTATS_GIT_VERSION "0.0.0-unversioned")
endif()
# Extract a clean x.y.z for project(VERSION)
string(REGEX REPLACE "^v([0-9]+\\.[0-9]+\\.[0-9]+).*" "\\1"
       _SEMVER "${GAPSTATS_GIT_VERSION}")
if(NOT _SEMVER MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+$")
  set(_SEMVER "0.0.0")
endif()

project(gapstats VERSION "${_SEMVER}" LANGUAGES CXX C)
message(STATUS "gapstats ${GAPSTATS_GIT_VERSION}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Build type ──────────────────────────────────────────────────────────────────
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# ── Coverage ──────────────────────────────────────────────────────────────────
option(COVERAGE "Enable gcov coverage instrumentation (-O0 --coverage)" OFF)
if(COVERAGE)
  add_compile_options(-O0 -g --coverage)
  add_link_options(--coverage)
endif()

# ── Optimisation (Release only, not during coverage profiling) ───────────────
if(CMAKE_BUILD_TYPE STREQUAL "Release" AND NOT COVERAGE)
  add_compile_options(-O3 -march=native)
endif()
add_compile_options(-Wall -Wextra)

# ── htslib: prefer system install (pixi/conda/apt), fall back to ExternalProject
find_package(PkgConfig QUIET)
set(_HTSLIB_SYSTEM FALSE)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(HTSLIB QUIET htslib>=1.20)
  if(HTSLIB_FOUND)
    set(_HTSLIB_SYSTEM TRUE)
    message(STATUS "Using system htslib ${HTSLIB_VERSION}")
    add_library(htslib_dep INTERFACE)
    target_include_directories(htslib_dep INTERFACE ${HTSLIB_INCLUDE_DIRS})
    target_link_libraries(htslib_dep INTERFACE ${HTSLIB_LINK_LIBRARIES})
    target_link_directories(htslib_dep INTERFACE ${HTSLIB_LIBRARY_DIRS})
  endif()
endif()

if(NOT _HTSLIB_SYSTEM)
  message(STATUS "System htslib not found — building 1.21 from source")
  include(ExternalProject)
  find_package(ZLIB REQUIRED)
  find_package(Threads REQUIRED)
  find_package(OpenSSL REQUIRED)
  find_library(BZ2_LIB  NAMES bz2  REQUIRED)
  find_library(LZMA_LIB NAMES lzma REQUIRED)
  find_library(CURL_LIB NAMES curl REQUIRED)

  set(HTSLIB_PREFIX "${CMAKE_BINARY_DIR}/htslib-prefix")
  file(MAKE_DIRECTORY "${HTSLIB_PREFIX}/include")
  file(MAKE_DIRECTORY "${HTSLIB_PREFIX}/lib")

  ExternalProject_Add(htslib_ext
    GIT_REPOSITORY  https://github.com/samtools/htslib.git
    GIT_TAG         1.21
    GIT_SHALLOW     TRUE
    PREFIX          "${HTSLIB_PREFIX}"
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND
      autoreconf -i
      COMMAND ./configure
        --prefix=${HTSLIB_PREFIX}
        --enable-libcurl
        --with-libdeflate
    BUILD_COMMAND   make -j${CMAKE_BUILD_PARALLEL_LEVEL}
    INSTALL_COMMAND make install
    UPDATE_COMMAND  ""
  )

  add_library(htslib_dep INTERFACE)
  add_dependencies(htslib_dep htslib_ext)
  target_include_directories(htslib_dep INTERFACE "${HTSLIB_PREFIX}/include")
  target_link_libraries(htslib_dep INTERFACE
    "${HTSLIB_PREFIX}/lib/libhts.a"
    ZLIB::ZLIB ${BZ2_LIB} ${LZMA_LIB} ${CURL_LIB} deflate
    Threads::Threads OpenSSL::Crypto
  )
endif()

# ── CLI11 via FetchContent (header-only) ──────────────────────────────────────
include(FetchContent)
FetchContent_Declare(
  cli11
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
  GIT_TAG        v2.4.1
  GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(cli11)

# ── Main executable ───────────────────────────────────────────────────────────
add_executable(gapstats src/main.cpp)
target_include_directories(gapstats PRIVATE src)
target_compile_definitions(gapstats PRIVATE
  GAPSTATS_VERSION="${GAPSTATS_GIT_VERSION}"
)
target_link_libraries(gapstats PRIVATE htslib_dep CLI11::CLI11)
install(TARGETS gapstats RUNTIME DESTINATION bin)

# ── Tests ─────────────────────────────────────────────────────────────────────
include(CTest)
if(BUILD_TESTING)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.6.0
    GIT_SHALLOW    TRUE
  )
  FetchContent_MakeAvailable(Catch2)
  add_subdirectory(tests)
endif()
