execute_process(
    COMMAND ${PYTHIA8_CONFIG_EXECUTABLE} --includedir
    OUTPUT_VARIABLE LOWMC_PYTHIA8_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE LOWMC_PYTHIA8_INCLUDE_RESULT
)
if(NOT LOWMC_PYTHIA8_INCLUDE_RESULT EQUAL 0)
    message(FATAL_ERROR "pythia8-config --includedir failed")
endif()

execute_process(
    COMMAND ${PYTHIA8_CONFIG_EXECUTABLE} --libdir
    OUTPUT_VARIABLE LOWMC_PYTHIA8_LIBRARY_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE LOWMC_PYTHIA8_LIBRARY_RESULT
)
if(NOT LOWMC_PYTHIA8_LIBRARY_RESULT EQUAL 0)
    message(FATAL_ERROR "pythia8-config --libdir failed")
endif()

execute_process(
    COMMAND ${PYTHIA8_CONFIG_EXECUTABLE} --xmldoc
    OUTPUT_VARIABLE LOWMC_PYTHIA8_XMLDOC_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE LOWMC_PYTHIA8_XMLDOC_RESULT
)
if(NOT LOWMC_PYTHIA8_XMLDOC_RESULT EQUAL 0)
    message(FATAL_ERROR "pythia8-config --xmldoc failed")
endif()

find_library(
    LOWMC_PYTHIA8_LIBRARY
    NAMES pythia8
    HINTS ${LOWMC_PYTHIA8_LIBRARY_DIR}
    REQUIRED
)

add_library(lowmc_pythia8_adapter STATIC
    src/Pythia8Adapter.cpp
)

target_compile_features(lowmc_pythia8_adapter PUBLIC cxx_std_17)

target_include_directories(lowmc_pythia8_adapter
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${LOWMC_PYTHIA8_INCLUDE_DIR}
)

target_link_libraries(lowmc_pythia8_adapter
    PUBLIC
        lowmc_core
        ${LOWMC_PYTHIA8_LIBRARY}
        ${CMAKE_DL_LIBS}
)

if(APPLE)
    set_target_properties(lowmc_pythia8_adapter PROPERTIES
        INSTALL_RPATH ${LOWMC_PYTHIA8_LIBRARY_DIR}
        BUILD_RPATH ${LOWMC_PYTHIA8_LIBRARY_DIR}
    )
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_pythia8_adapter PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_pythia8_adapter PRIVATE /W4)
endif()

add_executable(lowmc_pythia8_adapter_tests
    tests/test_pythia8_adapter.cpp
)

target_link_libraries(lowmc_pythia8_adapter_tests
    PRIVATE
        lowmc_pythia8_adapter
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_pythia8_adapter_tests PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_pythia8_adapter_tests PRIVATE /W4)
endif()

add_test(NAME lowmc_pythia8_adapter_tests COMMAND lowmc_pythia8_adapter_tests)

add_executable(lowmc_pythia8_adapter_diagnostic
    diagnostics/pythia8_adapter_diagnostic.cpp
)

target_link_libraries(lowmc_pythia8_adapter_diagnostic
    PRIVATE
        lowmc_pythia8_adapter
)

target_compile_definitions(lowmc_pythia8_adapter_diagnostic
    PRIVATE
        LOWMC_PYTHIA8_XMLDOC="${LOWMC_PYTHIA8_XMLDOC_DIR}"
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_pythia8_adapter_diagnostic PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_pythia8_adapter_diagnostic PRIVATE /W4)
endif()

add_test(
    NAME lowmc_pythia8_adapter_diagnostic_help
    COMMAND lowmc_pythia8_adapter_diagnostic --help
)

add_test(
    NAME lowmc_pythia8_adapter_diagnostic_toy
    COMMAND lowmc_pythia8_adapter_diagnostic --mode toy
)
set_tests_properties(
    lowmc_pythia8_adapter_diagnostic_toy
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "event_index=0 mode=toy selected_incoming_charged=2 selected_outgoing_charged=2 ignored_neutral_candidates=1 rejected_charged_candidates=0 signed_charge_sum=0 ordinary_charge_sum=4"
)

add_test(
    NAME lowmc_pythia8_adapter_diagnostic_softqcd
    COMMAND lowmc_pythia8_adapter_diagnostic --mode softqcd --events 1 --seed 12345
)
set_tests_properties(
    lowmc_pythia8_adapter_diagnostic_softqcd
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "event_index=0 mode=softqcd selected_incoming_charged=2 .* signed_charge_sum=0"
)

add_executable(lowmc_pythia8_eta_phi_map
    diagnostics/pythia8_eta_phi_map.cpp
)

target_link_libraries(lowmc_pythia8_eta_phi_map
    PRIVATE
        lowmc_pythia8_adapter
)

target_compile_definitions(lowmc_pythia8_eta_phi_map
    PRIVATE
        LOWMC_PYTHIA8_XMLDOC="${LOWMC_PYTHIA8_XMLDOC_DIR}"
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_pythia8_eta_phi_map PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_pythia8_eta_phi_map PRIVATE /W4)
endif()

add_executable(lowmc_pythia8_weighted_photon_diagnostic
    diagnostics/pythia8_weighted_photon_diagnostic.cpp
)

target_link_libraries(lowmc_pythia8_weighted_photon_diagnostic
    PRIVATE
        lowmc_pythia8_adapter
)

target_compile_definitions(lowmc_pythia8_weighted_photon_diagnostic
    PRIVATE
        LOWMC_PYTHIA8_XMLDOC="${LOWMC_PYTHIA8_XMLDOC_DIR}"
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_pythia8_weighted_photon_diagnostic PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_pythia8_weighted_photon_diagnostic PRIVATE /W4)
endif()

add_test(
    NAME lowmc_pythia8_weighted_photon_diagnostic_help
    COMMAND lowmc_pythia8_weighted_photon_diagnostic --help
)

add_test(
    NAME lowmc_pythia8_weighted_photon_diagnostic_toy
    COMMAND lowmc_pythia8_weighted_photon_diagnostic --mode toy --samples 8 --seed 12345
)
set_tests_properties(
    lowmc_pythia8_weighted_photon_diagnostic_toy
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "event_index=0 mode=toy selected_incoming_charged=2 selected_outgoing_charged=2 ignored_neutral_candidates=1 rejected_charged_candidates=0 signed_charge_sum=0 ordinary_charge_sum=4 photon_status=generated samples=8"
)

add_test(
    NAME lowmc_pythia8_weighted_photon_diagnostic_nonconserved_toy
    COMMAND lowmc_pythia8_weighted_photon_diagnostic --mode toy-nonconserved --samples 8
)
set_tests_properties(
    lowmc_pythia8_weighted_photon_diagnostic_nonconserved_toy
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "event_index=0 mode=toy-nonconserved .* photon_status=skipped_nonconserved_charge"
)

add_test(
    NAME lowmc_pythia8_weighted_photon_diagnostic_softqcd
    COMMAND lowmc_pythia8_weighted_photon_diagnostic --mode softqcd --events 1 --samples 4 --seed 12345
)

add_executable(lowmc_pythia8_alice_dielectron_benchmark
    diagnostics/pythia8_alice_dielectron_benchmark.cpp
)

target_link_libraries(lowmc_pythia8_alice_dielectron_benchmark
    PRIVATE
        lowmc_pythia8_adapter
)

target_compile_definitions(lowmc_pythia8_alice_dielectron_benchmark
    PRIVATE
        LOWMC_PYTHIA8_XMLDOC="${LOWMC_PYTHIA8_XMLDOC_DIR}"
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(
        lowmc_pythia8_alice_dielectron_benchmark PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_pythia8_alice_dielectron_benchmark PRIVATE /W4)
endif()

add_test(
    NAME lowmc_pythia8_alice_dielectron_benchmark_help
    COMMAND lowmc_pythia8_alice_dielectron_benchmark --help
)

add_executable(lowmc_pythia8_yield_convergence
    diagnostics/pythia8_yield_convergence.cpp
)

target_link_libraries(lowmc_pythia8_yield_convergence
    PRIVATE
        lowmc_pythia8_adapter
)

target_compile_definitions(lowmc_pythia8_yield_convergence
    PRIVATE
        LOWMC_PYTHIA8_XMLDOC="${LOWMC_PYTHIA8_XMLDOC_DIR}"
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(
        lowmc_pythia8_yield_convergence PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_pythia8_yield_convergence PRIVATE /W4)
endif()

add_test(
    NAME lowmc_pythia8_yield_convergence_help
    COMMAND lowmc_pythia8_yield_convergence --help
)

add_test(
    NAME lowmc_pythia8_alice_dielectron_benchmark_smoke
    COMMAND lowmc_pythia8_alice_dielectron_benchmark
            ${CMAKE_CURRENT_BINARY_DIR}/alice_dielectron_smoke.csv
            --events 1
            --samples 16
            --direct-samples 16
            --seed 12345
)
set_tests_properties(
    lowmc_pythia8_alice_dielectron_benchmark_smoke
    PROPERTIES PASS_REGULAR_EXPRESSION "aggregate_status=ok.*accepted_events=1"
)

add_test(
    NAME lowmc_pythia8_alice_dielectron_benchmark_publication_smoke
    COMMAND lowmc_pythia8_alice_dielectron_benchmark
            ${CMAKE_CURRENT_BINARY_DIR}/alice_dielectron_publication_smoke.csv
            --events 1
            --samples 16
            --direct-samples 16
            --seed 12345
            --comparison publication
)
set_tests_properties(
    lowmc_pythia8_alice_dielectron_benchmark_publication_smoke
    PROPERTIES PASS_REGULAR_EXPRESSION
        "aggregate_status=ok.*comparison_mode=published_cross_section"
)

add_executable(
    lowmc_example_pythia8_weighted_photons
    ${CMAKE_CURRENT_SOURCE_DIR}/../../examples/advanced/pythia8_weighted_photons.cpp
)

target_link_libraries(
    lowmc_example_pythia8_weighted_photons
    PRIVATE
        lowmc_pythia8_adapter
)

target_compile_definitions(
    lowmc_example_pythia8_weighted_photons
    PRIVATE
        LOWMC_PYTHIA8_XMLDOC="${LOWMC_PYTHIA8_XMLDOC_DIR}"
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(
        lowmc_example_pythia8_weighted_photons PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_example_pythia8_weighted_photons PRIVATE /W4)
endif()

add_test(
    NAME lowmc_example_pythia8_weighted_photons
    COMMAND lowmc_example_pythia8_weighted_photons
)
set_tests_properties(
    lowmc_example_pythia8_weighted_photons
    PROPERTIES PASS_REGULAR_EXPRESSION "charge_conserved=1.*weighted_photons=100"
)

add_executable(lowmc_pythia8_leg_proposal_tail_diagnostic
    diagnostics/pythia8_leg_proposal_tail_diagnostic.cpp
)

target_link_libraries(lowmc_pythia8_leg_proposal_tail_diagnostic
    PRIVATE
        lowmc_pythia8_adapter
)

target_compile_definitions(lowmc_pythia8_leg_proposal_tail_diagnostic
    PRIVATE
        LOWMC_PYTHIA8_XMLDOC="${LOWMC_PYTHIA8_XMLDOC_DIR}"
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(
        lowmc_pythia8_leg_proposal_tail_diagnostic PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_pythia8_leg_proposal_tail_diagnostic PRIVATE /W4)
endif()

add_test(
    NAME lowmc_pythia8_leg_proposal_tail_diagnostic_help
    COMMAND lowmc_pythia8_leg_proposal_tail_diagnostic --help
)

add_test(
    NAME lowmc_pythia8_leg_proposal_tail_diagnostic_ee_smoke
    COMMAND lowmc_pythia8_leg_proposal_tail_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/pythia8_leg_proposal_tail_ee.csv
            --mode ee-hadrons
            --events 1
            --samples 8
            --seed 12345
)

add_test(
    NAME lowmc_pythia8_leg_proposal_tail_diagnostic_softqcd_smoke
    COMMAND lowmc_pythia8_leg_proposal_tail_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/pythia8_leg_proposal_tail_softqcd.csv
            --mode softqcd
            --events 1
            --samples 8
            --seed 12345
)

add_executable(lowmc_pythia8_leg_envelope_diagnostic
    diagnostics/pythia8_leg_envelope_diagnostic.cpp
)

target_link_libraries(lowmc_pythia8_leg_envelope_diagnostic
    PRIVATE
        lowmc_pythia8_adapter
)

target_compile_definitions(lowmc_pythia8_leg_envelope_diagnostic
    PRIVATE
        LOWMC_PYTHIA8_XMLDOC="${LOWMC_PYTHIA8_XMLDOC_DIR}"
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(
        lowmc_pythia8_leg_envelope_diagnostic PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_pythia8_leg_envelope_diagnostic PRIVATE /W4)
endif()

add_test(
    NAME lowmc_pythia8_leg_envelope_diagnostic_help
    COMMAND lowmc_pythia8_leg_envelope_diagnostic --help
)

add_test(
    NAME lowmc_pythia8_leg_envelope_diagnostic_ee_smoke
    COMMAND lowmc_pythia8_leg_envelope_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/pythia8_leg_envelope_ee.csv
            --mode ee-hadrons
            --events 1
            --calibration-samples 8
            --trial-samples 16
            --seed 12345
)

add_test(
    NAME lowmc_pythia8_leg_envelope_diagnostic_softqcd_smoke
    COMMAND lowmc_pythia8_leg_envelope_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/pythia8_leg_envelope_softqcd.csv
            --mode softqcd
            --events 1
            --calibration-samples 8
            --trial-samples 16
            --seed 12345
)

add_executable(lowmc_pythia8_delphi_hadronic_benchmark
    diagnostics/pythia8_delphi_hadronic_benchmark.cpp
)

target_link_libraries(lowmc_pythia8_delphi_hadronic_benchmark
    PRIVATE
        lowmc_pythia8_adapter
)

target_compile_definitions(lowmc_pythia8_delphi_hadronic_benchmark
    PRIVATE
        LOWMC_PYTHIA8_XMLDOC="${LOWMC_PYTHIA8_XMLDOC_DIR}"
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(
        lowmc_pythia8_delphi_hadronic_benchmark PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_pythia8_delphi_hadronic_benchmark PRIVATE /W4)
endif()

add_test(
    NAME lowmc_pythia8_delphi_hadronic_benchmark_help
    COMMAND lowmc_pythia8_delphi_hadronic_benchmark --help
)

add_test(
    NAME lowmc_pythia8_delphi_hadronic_benchmark_smoke
    COMMAND lowmc_pythia8_delphi_hadronic_benchmark
            ${CMAKE_CURRENT_BINARY_DIR}/pythia8_delphi_hadronic_smoke.csv
            --events 1
            --mc-samples 64
            --direct-loge-bins 2
            --direct-eta-bins 4
            --direct-phi-bins 4
            --direct-refinement-levels 0
            --direct-subbins 1
            --seed 12345
)
set_tests_properties(
    lowmc_pythia8_delphi_hadronic_benchmark_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "aggregate_status=ok accepted_events=1"
)

add_executable(lowmc_pythia8_wa102_benchmark
    diagnostics/pythia8_wa102_benchmark.cpp
)

target_link_libraries(lowmc_pythia8_wa102_benchmark
    PRIVATE
        lowmc_pythia8_adapter
)

target_compile_definitions(lowmc_pythia8_wa102_benchmark
    PRIVATE
        LOWMC_PYTHIA8_XMLDOC="${LOWMC_PYTHIA8_XMLDOC_DIR}"
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(
        lowmc_pythia8_wa102_benchmark PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_pythia8_wa102_benchmark PRIVATE /W4)
endif()

add_test(
    NAME lowmc_pythia8_wa102_benchmark_help
    COMMAND lowmc_pythia8_wa102_benchmark --help
)

add_test(
    NAME lowmc_pythia8_wa102_benchmark_smoke
    COMMAND lowmc_pythia8_wa102_benchmark
            ${CMAKE_CURRENT_BINARY_DIR}/pythia8_wa102_smoke.csv
            --events 1
            --mc-samples 64
            --direct-loge-bins 2
            --direct-eta-bins 4
            --direct-phi-bins 4
            --direct-refinement-levels 0
            --direct-subbins 1
            --seed 12345
)
set_tests_properties(
    lowmc_pythia8_wa102_benchmark_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "aggregate_status=ok accepted_events=1"
)

add_test(
    NAME lowmc_pythia8_delphi_hadronic_jet_accounting_smoke
    COMMAND lowmc_pythia8_delphi_hadronic_benchmark
            ${CMAKE_CURRENT_BINARY_DIR}/pythia8_delphi_hadronic_jet_accounting_smoke.csv
            --events 1
            --jet-accounting-only
            --seed 12345
)
set_tests_properties(
    lowmc_pythia8_delphi_hadronic_jet_accounting_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "jet_accounting_status=ok seed=12345 events=1"
)

foreach(LOWMC_DELPHI_JET_ALGORITHM IN ITEMS Durham Jade)
    string(
        TOLOWER
        ${LOWMC_DELPHI_JET_ALGORITHM}
        LOWMC_DELPHI_JET_ALGORITHM_LOWER
    )
    add_test(
        NAME
            lowmc_pythia8_delphi_hadronic_benchmark_${LOWMC_DELPHI_JET_ALGORITHM_LOWER}_smoke
        COMMAND
            lowmc_pythia8_delphi_hadronic_benchmark
            ${CMAKE_CURRENT_BINARY_DIR}/pythia8_delphi_hadronic_${LOWMC_DELPHI_JET_ALGORITHM_LOWER}_smoke.csv
            --events 1
            --mc-samples 64
            --direct-loge-bins 2
            --direct-eta-bins 4
            --direct-phi-bins 4
            --direct-refinement-levels 0
            --direct-subbins 1
            --jet-algorithm ${LOWMC_DELPHI_JET_ALGORITHM}
            --seed 12345
    )
    set_tests_properties(
        lowmc_pythia8_delphi_hadronic_benchmark_${LOWMC_DELPHI_JET_ALGORITHM_LOWER}_smoke
        PROPERTIES
            PASS_REGULAR_EXPRESSION
                "aggregate_status=ok accepted_events=1 accepted_jets=.*jet_algorithm=${LOWMC_DELPHI_JET_ALGORITHM}"
    )
endforeach()

add_executable(lowmc_pythia8_analytic_bound_audit
    diagnostics/pythia8_analytic_bound_audit.cpp
)

target_link_libraries(lowmc_pythia8_analytic_bound_audit
    PRIVATE
        lowmc_pythia8_adapter
)

target_compile_definitions(lowmc_pythia8_analytic_bound_audit
    PRIVATE
        LOWMC_PYTHIA8_XMLDOC="${LOWMC_PYTHIA8_XMLDOC_DIR}"
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(
        lowmc_pythia8_analytic_bound_audit
        PRIVATE -Wall -Wextra -pedantic
    )
elseif(MSVC)
    target_compile_options(
        lowmc_pythia8_analytic_bound_audit PRIVATE /W4
    )
endif()

add_test(
    NAME lowmc_pythia8_analytic_bound_audit_help
    COMMAND lowmc_pythia8_analytic_bound_audit --help
)

add_test(
    NAME lowmc_pythia8_analytic_bound_audit_smoke
    COMMAND lowmc_pythia8_analytic_bound_audit
            ${CMAKE_CURRENT_BINARY_DIR}/pythia8_analytic_bound_audit_smoke.csv
            --events 1
            --trials 64
            --calibration-trials 16
            --seed 12345
)
set_tests_properties(
    lowmc_pythia8_analytic_bound_audit_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "analytic_bound_audit_status=ok mode=ee-hadrons seed=12345 eta_center=pseudorapidity events=1"
)

add_test(
    NAME lowmc_pythia8_true_rapidity_center_smoke
    COMMAND lowmc_pythia8_analytic_bound_audit
            ${CMAKE_CURRENT_BINARY_DIR}/pythia8_true_rapidity_center_smoke.csv
            --events 1
            --trials 64
            --calibration-trials 16
            --seed 12345
            --true-rapidity-centers
)
set_tests_properties(
    lowmc_pythia8_true_rapidity_center_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "analytic_bound_audit_status=ok mode=ee-hadrons seed=12345 eta_center=rapidity events=1"
)

if(Python3_Interpreter_FOUND)
    add_test(
        NAME lowmc_pythia8_delphi_hadronic_summary_smoke
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/scripts/summarize_delphi_hadronic_benchmark.py
                ${CMAKE_CURRENT_BINARY_DIR}/pythia8_delphi_hadronic_smoke.csv
                --summary
                ${CMAKE_CURRENT_BINARY_DIR}/pythia8_delphi_hadronic_smoke_summary.csv
                --joint-inputs
                ${CMAKE_CURRENT_BINARY_DIR}/pythia8_delphi_hadronic_smoke_joint.csv
                --theta-summary
                ${CMAKE_CURRENT_BINARY_DIR}/pythia8_delphi_hadronic_smoke_theta_summary.csv
    )
    set_tests_properties(
        lowmc_pythia8_delphi_hadronic_summary_smoke
        PROPERTIES
            DEPENDS lowmc_pythia8_delphi_hadronic_benchmark_smoke
            PASS_REGULAR_EXPRESSION
                "theta_events=1 theta_jets=.*theta_pt400_direct_per_jet=.*theta_pt80_diagnostic_direct_per_jet="
    )

    add_test(
        NAME lowmc_pythia8_wa102_summary_smoke
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/scripts/summarize_wa102_benchmark.py
                ${CMAKE_CURRENT_BINARY_DIR}/pythia8_wa102_smoke.csv
                --summary
                ${CMAKE_CURRENT_BINARY_DIR}/pythia8_wa102_smoke_summary.csv
    )
    set_tests_properties(
        lowmc_pythia8_wa102_summary_smoke
        PROPERTIES
            DEPENDS lowmc_pythia8_wa102_benchmark_smoke
            PASS_REGULAR_EXPRESSION
                "events=1 .*lowmc_over_direct=.*comparison_lowmc_over_direct="
    )
endif()

set_tests_properties(
    lowmc_pythia8_weighted_photon_diagnostic_softqcd
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "event_index=0 mode=softqcd selected_incoming_charged=2 .* photon_status=generated samples=4"
)

if(LOWMC_ENABLE_HEPMC3)
    add_executable(lowmc_pythia8_hepmc3_crosscheck
        diagnostics/pythia8_hepmc3_crosscheck.cpp
    )

    target_link_libraries(lowmc_pythia8_hepmc3_crosscheck
        PRIVATE
            lowmc_pythia8_adapter
            lowmc_hepmc3_adapter
    )

    target_compile_definitions(lowmc_pythia8_hepmc3_crosscheck
        PRIVATE
            LOWMC_PYTHIA8_XMLDOC="${LOWMC_PYTHIA8_XMLDOC_DIR}"
    )

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        target_compile_options(lowmc_pythia8_hepmc3_crosscheck PRIVATE -Wall -Wextra -pedantic)
    elseif(MSVC)
        target_compile_options(lowmc_pythia8_hepmc3_crosscheck PRIVATE /W4)
    endif()

    add_test(
        NAME lowmc_pythia8_hepmc3_crosscheck_help
        COMMAND lowmc_pythia8_hepmc3_crosscheck --help
    )

    add_test(
        NAME lowmc_pythia8_hepmc3_crosscheck_softqcd
        COMMAND lowmc_pythia8_hepmc3_crosscheck --events 3 --samples 4 --seed 12345
    )
    set_tests_properties(
        lowmc_pythia8_hepmc3_crosscheck_softqcd
        PROPERTIES
            PASS_REGULAR_EXPRESSION
                "aggregate_status=match events_compared=3 .* mismatch_count=0 .* max_relative_weighted_summary_difference=[0-9.e+-]+"
    )
endif()
