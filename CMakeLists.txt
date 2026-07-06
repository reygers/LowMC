cmake_minimum_required(VERSION 3.16)

project(LowMC VERSION 0.1.0 LANGUAGES CXX)

include(GNUInstallDirs)

option(LOWMC_ENABLE_HEPMC3 "Build the optional HepMC3 adapter target" OFF)
option(LOWMC_ENABLE_PYTHIA8 "Build the optional Pythia8 adapter target" OFF)

enable_testing()

find_package(Python3 COMPONENTS Interpreter)

if(LOWMC_ENABLE_HEPMC3)
    find_package(HepMC3 CONFIG QUIET)
    if(NOT HepMC3_FOUND)
        message(FATAL_ERROR
            "LOWMC_ENABLE_HEPMC3=ON requires HepMC3, but CMake could not find "
            "a HepMC3 package configuration. Install HepMC3 or pass its prefix "
            "with -DCMAKE_PREFIX_PATH=/path/to/hepmc3."
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
        NAME lowmc_generate_pythia_hepmc3_sample_help
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/generate_pythia_hepmc3_sample.py --help
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
endif()
