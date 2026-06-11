add_executable(test_gapstats test_core.cpp)
target_include_directories(test_gapstats PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_gapstats PRIVATE Catch2::Catch2WithMain)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(test_gapstats)

# ── Coverage report target (only active when -DCOVERAGE=ON) ──────────────────
if(COVERAGE)
  find_program(LCOV_CMD    lcov    REQUIRED)
  find_program(GENHTML_CMD genhtml REQUIRED)

  add_custom_target(coverage
    # 1. Capture all .gcda data
    COMMAND ${LCOV_CMD}
            --directory "${CMAKE_BINARY_DIR}"
            --capture
            --output-file coverage.info
    # 2. Remove third-party and system paths
    COMMAND ${LCOV_CMD}
            --remove coverage.info
            "${CMAKE_BINARY_DIR}/_deps/*"
            "${CMAKE_BINARY_DIR}/htslib-prefix/*"
            "/usr/*"
            --output-file coverage.info
    # 3. HTML report
    COMMAND ${GENHTML_CMD} coverage.info
            --output-directory coverage_html
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    COMMENT "Generating lcov coverage report → build-cov/coverage_html/index.html"
    VERBATIM
  )
endif()
