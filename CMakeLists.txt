cmake_minimum_required(VERSION 3.16)

project(LowMC VERSION 0.1.0 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(LOWMC_ENABLE_HEPMC3 "Build the optional HepMC3 adapter target" OFF)
option(LOWMC_ENABLE_PYTHIA8 "Build the optional Pythia8 adapter target" OFF)
option(LOWMC_ENABLE_PHOTOS "Build the optional PHOTOS HepMC3 benchmark" OFF)

set(PHOTOS_ROOT "" CACHE PATH "PHOTOS installation prefix")

enable_testing()

find_package(Python3 COMPONENTS Interpreter)

if(LOWMC_ENABLE_HEPMC3 OR LOWMC_ENABLE_PHOTOS)
    find_package(HepMC3 CONFIG QUIET)
    if(NOT HepMC3_FOUND)
        message(FATAL_ERROR
            "The enabled HepMC3-dependent target requires HepMC3, but CMake could not find "
            "a HepMC3 package configuration. Install HepMC3 or pass its prefix "
            "with -DCMAKE_PREFIX_PATH=/path/to/hepmc3."
        )
    endif()
endif()

if(LOWMC_ENABLE_PHOTOS)
    if(PHOTOS_ROOT)
        find_path(
            PHOTOS_INCLUDE_DIR
            NAMES Photos/Photos.h
            PATHS "${PHOTOS_ROOT}/include"
            NO_DEFAULT_PATH
        )
    else()
        find_path(
            PHOTOS_INCLUDE_DIR
            NAMES Photos/Photos.h
            PATH_SUFFIXES include
        )
    endif()
    find_library(
        PHOTOS_CORE_LIBRARY
        NAMES Photospp libPhotospp
        HINTS "${PHOTOS_ROOT}"
        PATH_SUFFIXES lib
    )
    find_library(
        PHOTOS_HEPMC3_LIBRARY
        NAMES PhotosppHepMC3 libPhotosppHepMC3
        HINTS "${PHOTOS_ROOT}"
        PATH_SUFFIXES lib
    )
    if(NOT PHOTOS_INCLUDE_DIR OR NOT PHOTOS_CORE_LIBRARY OR NOT PHOTOS_HEPMC3_LIBRARY)
        message(FATAL_ERROR
            "LOWMC_ENABLE_PHOTOS=ON requires PHOTOS with its HepMC3 interface. "
            "Pass the installation prefix with -DPHOTOS_ROOT=/path/to/photospp."
        )
    endif()
endif()

if(LOWMC_ENABLE_PYTHIA8)
    find_program(
        PYTHIA8_CONFIG_EXECUTABLE
        NAMES pythia8-config
        HINTS "$ENV{HOME}/pythia8317/bin"
    )
    if(NOT PYTHIA8_CONFIG_EXECUTABLE)
        message(FATAL_ERROR
            "LOWMC_ENABLE_PYTHIA8=ON requires Pythia8, but CMake could not "
            "find pythia8-config. Install Pythia8, add pythia8-config to PATH, "
            "or use a standard local prefix such as $HOME/pythia8317."
        )
    endif()
endif()

include(cmake/LowMCTargets.cmake)

install(
    TARGETS lowmc_core
    EXPORT LowMCTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
    EXPORT LowMCTargets
    FILE LowMCTargets.cmake
    NAMESPACE LowMC::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LowMC
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/LowMCConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/LowMCConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LowMC
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/LowMCConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(
    FILES
        ${CMAKE_CURRENT_BINARY_DIR}/LowMCConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/LowMCConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LowMC
)

install(
    FILES ${CMAKE_CURRENT_SOURCE_DIR}/examples/README.md
    DESTINATION ${CMAKE_INSTALL_DATADIR}/lowmc/examples
)

install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples/basic/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/lowmc/examples/basic
    FILES_MATCHING PATTERN "*.cpp"
)

install(
    FILES ${CMAKE_CURRENT_SOURCE_DIR}/examples/advanced/shape_unweighting.cpp
    DESTINATION ${CMAKE_INSTALL_DATADIR}/lowmc/examples/advanced
)

if(Python3_Interpreter_FOUND)
    add_test(
        NAME lowmc_plot_eta_phi_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/plot_eta_phi.py --help
    )
    add_test(
        NAME lowmc_compare_fixed_grid_weighted_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_fixed_grid_weighted.py --help
    )
    add_test(
        NAME lowmc_plot_sampler_study_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/plot_sampler_study.py --help
    )
    add_test(
        NAME lowmc_run_grid_weighted_stress_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_grid_weighted_stress.py --help
    )
    add_test(
        NAME lowmc_scan_adaptive_tuning_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/scan_adaptive_tuning.py --help
    )
    add_test(
        NAME lowmc_plot_adaptive_tuning_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/plot_adaptive_tuning.py --help
    )
    add_test(
        NAME lowmc_run_p100_focus_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_p100_focus.py --help
    )
    add_test(
        NAME lowmc_scan_grid_tuning_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/scan_grid_tuning.py --help
    )
    add_test(
        NAME lowmc_run_p100_reference_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_p100_reference.py --help
    )
    add_test(
        NAME lowmc_retabulate_grid_tuning_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/retabulate_grid_tuning.py --help
    )
    add_test(
        NAME lowmc_compare_p100_proposals_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_p100_proposals.py --help
    )
    add_test(
        NAME lowmc_scan_smooth_tuning_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/scan_smooth_tuning.py --help
    )
    add_test(
        NAME lowmc_compare_variance_grid_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_variance_grid.py --help
    )
    add_test(
        NAME lowmc_multiseed_variance_grid_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/multiseed_variance_grid.py --help
    )
    add_test(
        NAME lowmc_compare_full_weighted_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_full_weighted.py --help
    )
    add_test(
        NAME lowmc_multirange_full_weighted_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/multirange_full_weighted.py --help
    )
    add_test(
        NAME lowmc_plot_sampler_validation_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/plot_sampler_validation.py --help
    )
    add_test(
        NAME lowmc_plot_alice_dielectron_comparison_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/plot_alice_dielectron_comparison.py --help
    )
    add_test(
        NAME lowmc_plot_alice_dielectron_publication_comparison_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/plot_alice_dielectron_publication_comparison.py --help
    )
    add_test(
        NAME lowmc_scan_unweighting_envelope_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/scan_unweighting_envelope.py --help
    )
    add_test(
        NAME lowmc_stress_unweighting_bound_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/stress_unweighting_bound.py --help
    )
    add_test(
        NAME lowmc_compare_shape_unweighting_diagnostics_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_shape_unweighting_diagnostics.py --help
    )
    add_test(
        NAME lowmc_validate_event_yield_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_event_yield.py --help
    )
    add_test(
        NAME lowmc_validate_count_shape_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_count_shape.py --help
    )
    add_test(
        NAME lowmc_validate_release_candidate_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_release_candidate.py --help
    )
    add_test(
        NAME lowmc_validate_release_candidate_tests
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_release_validation.py
    )
    add_test(
        NAME lowmc_generate_pythia_hepmc3_sample_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/generate_pythia_hepmc3_sample.py --help
    )
    add_test(
        NAME lowmc_validate_pythia_afterburner_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_pythia_afterburner.py --help
    )
    add_test(
        NAME lowmc_summarize_pythia_afterburner_timing_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/summarize_pythia_afterburner_timing.py --help
    )
    add_test(
        NAME lowmc_tune_event_yield_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/tune_event_yield.py --help
    )
    add_test(
        NAME lowmc_validate_shape_unweighting_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_shape_unweighting.py --help
    )
    add_test(
        NAME lowmc_validate_shape_unweighting_candidates_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_shape_unweighting_candidates.py --help
    )
    add_test(
        NAME lowmc_tune_shape_unweighting_bound_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/tune_shape_unweighting_bound.py --help
    )
    add_test(
        NAME lowmc_plot_delphi_muon_benchmark_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/plot_delphi_muon_benchmark.py --help
    )
    add_test(
        NAME lowmc_plot_photos_zmumu_comparison_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/plot_photos_zmumu_comparison.py --help
    )
    add_test(
        NAME lowmc_run_photos_zmumu_soft_scan_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_photos_zmumu_soft_scan.py --help
    )
    add_test(
        NAME lowmc_compare_leg_centered_proposal_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare_leg_centered_proposal.py --help
    )
    add_test(
        NAME lowmc_summarize_pythia_leg_tail_study_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/summarize_pythia_leg_tail_study.py --help
    )
    add_test(
        NAME lowmc_summarize_pythia_leg_envelope_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/summarize_pythia_leg_envelope.py --help
    )
    add_test(
        NAME lowmc_summarize_delphi_hadronic_benchmark_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/summarize_delphi_hadronic_benchmark.py --help
    )
    add_test(
        NAME lowmc_summarize_wa102_benchmark_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/summarize_wa102_benchmark.py --help
    )
    add_test(
        NAME lowmc_summarize_delphi_jet_accounting_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/summarize_delphi_jet_accounting.py --help
    )
    add_test(
        NAME lowmc_summarize_analytic_bound_audit_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/summarize_analytic_bound_audit.py --help
    )
    add_test(
        NAME lowmc_summarize_rapidity_center_audit_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/summarize_rapidity_center_audit.py --help
    )
endif()
