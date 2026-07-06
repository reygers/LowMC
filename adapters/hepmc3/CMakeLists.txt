add_library(lowmc_hepmc3_adapter STATIC
    src/HepMC3Adapter.cpp
)

target_compile_features(lowmc_hepmc3_adapter PUBLIC cxx_std_17)

target_include_directories(lowmc_hepmc3_adapter
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if(TARGET HepMC3::HepMC3)
    set(LOWMC_HEPMC3_TARGET HepMC3::HepMC3)
elseif(TARGET HepMC3::hepmc3)
    set(LOWMC_HEPMC3_TARGET HepMC3::hepmc3)
elseif(TARGET HepMC3::HepMC3_static)
    set(LOWMC_HEPMC3_TARGET HepMC3::HepMC3_static)
else()
    message(FATAL_ERROR
        "HepMC3 was found, but LowMC does not recognize the imported CMake "
        "target name. Expected one of HepMC3::HepMC3, HepMC3::hepmc3, or "
        "HepMC3::HepMC3_static."
    )
endif()

target_link_libraries(lowmc_hepmc3_adapter
    PUBLIC
        lowmc_core
        ${LOWMC_HEPMC3_TARGET}
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_hepmc3_adapter PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_hepmc3_adapter PRIVATE /W4)
endif()

add_executable(lowmc_hepmc3_adapter_tests
    tests/test_hepmc3_adapter.cpp
)

target_link_libraries(lowmc_hepmc3_adapter_tests
    PRIVATE
        lowmc_hepmc3_adapter
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_hepmc3_adapter_tests PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_hepmc3_adapter_tests PRIVATE /W4)
endif()

add_test(NAME lowmc_hepmc3_adapter_tests COMMAND lowmc_hepmc3_adapter_tests)

add_executable(lowmc_hepmc3_adapter_diagnostic
    diagnostics/hepmc3_adapter_diagnostic.cpp
)

target_link_libraries(lowmc_hepmc3_adapter_diagnostic
    PRIVATE
        lowmc_hepmc3_adapter
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_hepmc3_adapter_diagnostic PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_hepmc3_adapter_diagnostic PRIVATE /W4)
endif()

add_test(
    NAME lowmc_hepmc3_adapter_diagnostic_help
    COMMAND lowmc_hepmc3_adapter_diagnostic --help
)

add_test(
    NAME lowmc_hepmc3_adapter_diagnostic_tiny_fixture
    COMMAND lowmc_hepmc3_adapter_diagnostic
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/tiny_event.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_adapter_diagnostic_tiny_fixture
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "event_index=0 momentum_unit=MEV selected_incoming_charged=2 selected_outgoing_charged=1 ignored_neutral_candidates=2 rejected_charged_candidates=0 signed_charge_sum=-1 ordinary_charge_sum=3"
)

add_test(
    NAME lowmc_hepmc3_adapter_diagnostic_charge_conserving_fixture
    COMMAND lowmc_hepmc3_adapter_diagnostic
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/charge_conserving_event.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_adapter_diagnostic_charge_conserving_fixture
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "event_index=0 momentum_unit=MEV selected_incoming_charged=2 selected_outgoing_charged=2 ignored_neutral_candidates=1 rejected_charged_candidates=0 signed_charge_sum=0 ordinary_charge_sum=4"
)

add_test(
    NAME lowmc_hepmc3_adapter_diagnostic_validation_fixture
    COMMAND lowmc_hepmc3_adapter_diagnostic
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/validation_event.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_adapter_diagnostic_validation_fixture
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "event_index=0 momentum_unit=MEV selected_incoming_charged=2 selected_outgoing_charged=2 ignored_neutral_candidates=1 rejected_charged_candidates=0 signed_charge_sum=0 ordinary_charge_sum=4"
)

add_executable(lowmc_hepmc3_selected_particle_text
    diagnostics/hepmc3_selected_particle_text.cpp
)

target_link_libraries(lowmc_hepmc3_selected_particle_text
    PRIVATE
        lowmc_hepmc3_adapter
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_hepmc3_selected_particle_text PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_hepmc3_selected_particle_text PRIVATE /W4)
endif()

add_test(
    NAME lowmc_hepmc3_selected_particle_text_help
    COMMAND lowmc_hepmc3_selected_particle_text --help
)

add_test(
    NAME lowmc_hepmc3_selected_particle_text_validation_fixture
    COMMAND lowmc_hepmc3_selected_particle_text
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/validation_event.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/validation_event_particles.lowmc
)
set_tests_properties(
    lowmc_hepmc3_selected_particle_text_validation_fixture
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "event_index=0 momentum_unit=MEV selected_incoming_charged=2 selected_outgoing_charged=2 ignored_neutral_candidates=1 rejected_charged_candidates=0 signed_charge_sum=0 ordinary_charge_sum=4 output=4"
)

add_test(
    NAME lowmc_hepmc3_selected_particle_text_matches_checked_in_fixture
    COMMAND ${CMAKE_COMMAND} -E compare_files
            ${CMAKE_CURRENT_BINARY_DIR}/validation_event_particles.lowmc
            ${CMAKE_CURRENT_SOURCE_DIR}/../../fixtures/hepmc3_validation_event_particles.lowmc
)
set_tests_properties(
    lowmc_hepmc3_selected_particle_text_matches_checked_in_fixture
    PROPERTIES
        DEPENDS lowmc_hepmc3_selected_particle_text_validation_fixture
)

add_executable(lowmc_hepmc3_weighted_photon_diagnostic
    diagnostics/hepmc3_weighted_photon_diagnostic.cpp
)

target_link_libraries(lowmc_hepmc3_weighted_photon_diagnostic
    PRIVATE
        lowmc_hepmc3_adapter
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_hepmc3_weighted_photon_diagnostic PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_hepmc3_weighted_photon_diagnostic PRIVATE /W4)
endif()

add_test(
    NAME lowmc_hepmc3_weighted_photon_diagnostic_help
    COMMAND lowmc_hepmc3_weighted_photon_diagnostic --help
)

add_test(
    NAME lowmc_hepmc3_weighted_photon_diagnostic_charge_conserving_fixture
    COMMAND lowmc_hepmc3_weighted_photon_diagnostic
            --samples 8
            --seed 12345
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/charge_conserving_event.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_weighted_photon_diagnostic_charge_conserving_fixture
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "event_index=0 momentum_unit=MEV selected_incoming_charged=2 selected_outgoing_charged=2 ignored_neutral_candidates=1 rejected_charged_candidates=0 signed_charge_sum=0 ordinary_charge_sum=4 photon_status=generated samples=8"
)

add_executable(lowmc_hepmc3_afterburner_diagnostic
    diagnostics/hepmc3_afterburner_diagnostic.cpp
)

target_link_libraries(lowmc_hepmc3_afterburner_diagnostic
    PRIVATE
        lowmc_hepmc3_adapter
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_hepmc3_afterburner_diagnostic PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_hepmc3_afterburner_diagnostic PRIVATE /W4)
endif()

add_executable(lowmc_hepmc3_output_summary_diagnostic
    diagnostics/hepmc3_output_summary_diagnostic.cpp
)

target_link_libraries(lowmc_hepmc3_output_summary_diagnostic
    PRIVATE
        lowmc_hepmc3_adapter
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_hepmc3_output_summary_diagnostic PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_hepmc3_output_summary_diagnostic PRIVATE /W4)
endif()

add_executable(lowmc_hepmc3_output_file_tests
    tests/check_hepmc3_low_photon_output.cpp
)

target_link_libraries(lowmc_hepmc3_output_file_tests
    PRIVATE
        lowmc_hepmc3_adapter
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_hepmc3_output_file_tests PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_hepmc3_output_file_tests PRIVATE /W4)
endif()

add_executable(lowmc_hepmc3_multi_event_output_tests
    tests/check_hepmc3_multi_event_output.cpp
)

target_link_libraries(lowmc_hepmc3_multi_event_output_tests
    PRIVATE
        lowmc_hepmc3_adapter
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_hepmc3_multi_event_output_tests PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_hepmc3_multi_event_output_tests PRIVATE /W4)
endif()

add_executable(lowmc_hepmc3_realistic_multi_event_output_tests
    tests/check_hepmc3_realistic_multi_event_output.cpp
)

target_link_libraries(lowmc_hepmc3_realistic_multi_event_output_tests
    PRIVATE
        lowmc_hepmc3_adapter
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(
        lowmc_hepmc3_realistic_multi_event_output_tests PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_hepmc3_realistic_multi_event_output_tests PRIVATE /W4)
endif()

add_executable(lowmc_hepmc3_generate_realistic_sample
    tests/generate_realistic_sample.cpp
)

target_link_libraries(lowmc_hepmc3_generate_realistic_sample
    PRIVATE
        lowmc_hepmc3_adapter
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_hepmc3_generate_realistic_sample PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_hepmc3_generate_realistic_sample PRIVATE /W4)
endif()

add_executable(lowmc_hepmc3_generated_sample_output_tests
    tests/check_hepmc3_generated_sample_output.cpp
)

target_link_libraries(lowmc_hepmc3_generated_sample_output_tests
    PRIVATE
        lowmc_hepmc3_adapter
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(
        lowmc_hepmc3_generated_sample_output_tests PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_hepmc3_generated_sample_output_tests PRIVATE /W4)
endif()

add_executable(lowmc_hepmc3_no_lowmc_output_tests
    tests/check_hepmc3_no_lowmc_output.cpp
)

target_link_libraries(lowmc_hepmc3_no_lowmc_output_tests
    PRIVATE
        lowmc_hepmc3_adapter
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_hepmc3_no_lowmc_output_tests PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_hepmc3_no_lowmc_output_tests PRIVATE /W4)
endif()

add_executable(lowmc_hepmc3_charge_conserving_no_lowmc_output_tests
    tests/check_hepmc3_charge_conserving_no_lowmc_output.cpp
)

target_link_libraries(lowmc_hepmc3_charge_conserving_no_lowmc_output_tests
    PRIVATE
        lowmc_hepmc3_adapter
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(
        lowmc_hepmc3_charge_conserving_no_lowmc_output_tests PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_hepmc3_charge_conserving_no_lowmc_output_tests PRIVATE /W4)
endif()

install(
    TARGETS
        lowmc_hepmc3_adapter_diagnostic
        lowmc_hepmc3_selected_particle_text
        lowmc_hepmc3_weighted_photon_diagnostic
        lowmc_hepmc3_afterburner_diagnostic
        lowmc_hepmc3_output_summary_diagnostic
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(
    PROGRAMS
        ${CMAKE_CURRENT_SOURCE_DIR}/../../examples/hepmc3/hepmc3_quickstart.py
    DESTINATION ${CMAKE_INSTALL_DATADIR}/lowmc/examples/hepmc3
)

install(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/../../examples/hepmc3/README.md
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/charge_conserving_event.hepmc3
    DESTINATION ${CMAKE_INSTALL_DATADIR}/lowmc/examples/hepmc3
)

set_target_properties(
    lowmc_hepmc3_adapter_diagnostic
    lowmc_hepmc3_selected_particle_text
    lowmc_hepmc3_weighted_photon_diagnostic
    lowmc_hepmc3_afterburner_diagnostic
    lowmc_hepmc3_output_summary_diagnostic
    PROPERTIES
        INSTALL_RPATH_USE_LINK_PATH TRUE
)

add_test(
    NAME lowmc_hepmc3_afterburner_diagnostic_help
    COMMAND lowmc_hepmc3_afterburner_diagnostic --help
)

add_test(
    NAME lowmc_hepmc3_output_summary_diagnostic_help
    COMMAND lowmc_hepmc3_output_summary_diagnostic --help
)

add_test(
    NAME lowmc_hepmc3_afterburner_rejects_invalid_enhancement
    COMMAND lowmc_hepmc3_afterburner_diagnostic
            --enhancement-factor 0
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/charge_conserving_event.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/afterburner_rejects_invalid_enhancement.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_rejects_invalid_enhancement
    PROPERTIES
        WILL_FAIL TRUE
)

add_test(
    NAME lowmc_hepmc3_afterburner_neutral_only_skip
    COMMAND lowmc_hepmc3_afterburner_diagnostic
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/neutral_only_event.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/neutral_only_skip_output.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_neutral_only_skip
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "event_index=0 momentum_unit=MEV selected_incoming_charged=0 selected_outgoing_charged=0 ignored_neutral_candidates=1 rejected_charged_candidates=0 signed_charge_sum=0 ordinary_charge_sum=0 photon_status=skipped_invalid_selection.*photons_added=0"
)

add_test(
    NAME lowmc_hepmc3_afterburner_neutral_only_output_unchanged
    COMMAND lowmc_hepmc3_no_lowmc_output_tests
            ${CMAKE_CURRENT_BINARY_DIR}/neutral_only_skip_output.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_neutral_only_output_unchanged
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_neutral_only_skip
)

add_test(
    NAME lowmc_hepmc3_output_summary_neutral_only_skip
    COMMAND lowmc_hepmc3_output_summary_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/neutral_only_skip_output.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_output_summary_neutral_only_skip
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_neutral_only_skip
        PASS_REGULAR_EXPRESSION
            "event_index=0 lowmc_photons=0 photon_mode=none enhancement_factor=0 photon_weight_sum=0 photon_weight_attributes=0 event_metadata=0"
)

add_test(
    NAME lowmc_hepmc3_afterburner_neutral_only_fail_policy
    COMMAND lowmc_hepmc3_afterburner_diagnostic
            --failure-policy fail
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/neutral_only_event.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/neutral_only_fail_output.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_neutral_only_fail_policy
    PROPERTIES
        WILL_FAIL TRUE
)

add_test(
    NAME lowmc_hepmc3_afterburner_insufficient_candidates_skip
    COMMAND lowmc_hepmc3_afterburner_diagnostic
            --enhancement-factor 100
            --multiplicity-seed 2
            --shape-seed 67890
            --min-shape-trials 1
            --shape-trials-per-photon 1
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/charge_conserving_event.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/insufficient_candidates_skip_output.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_insufficient_candidates_skip
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "event_index=0 momentum_unit=MEV selected_incoming_charged=2 selected_outgoing_charged=2 ignored_neutral_candidates=1 rejected_charged_candidates=0 signed_charge_sum=0 ordinary_charge_sum=4 photon_status=skipped_insufficient_shape_candidates.*requested_photons=12 retained_photons=5 photons_added=0 shape_status=ok"
)

add_test(
    NAME lowmc_hepmc3_afterburner_insufficient_candidates_output_unchanged
    COMMAND lowmc_hepmc3_charge_conserving_no_lowmc_output_tests
            ${CMAKE_CURRENT_BINARY_DIR}/insufficient_candidates_skip_output.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_insufficient_candidates_output_unchanged
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_insufficient_candidates_skip
)

add_test(
    NAME lowmc_hepmc3_output_summary_insufficient_candidates_skip
    COMMAND lowmc_hepmc3_output_summary_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/insufficient_candidates_skip_output.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_output_summary_insufficient_candidates_skip
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_insufficient_candidates_skip
        PASS_REGULAR_EXPRESSION
            "event_index=0 lowmc_photons=0 photon_mode=none enhancement_factor=0 photon_weight_sum=0 photon_weight_attributes=0 event_metadata=0"
)

add_test(
    NAME lowmc_hepmc3_afterburner_insufficient_candidates_fail_policy
    COMMAND lowmc_hepmc3_afterburner_diagnostic
            --failure-policy fail
            --enhancement-factor 100
            --multiplicity-seed 2
            --shape-seed 67890
            --min-shape-trials 1
            --shape-trials-per-photon 1
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/charge_conserving_event.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/insufficient_candidates_fail_output.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_insufficient_candidates_fail_policy
    PROPERTIES
        WILL_FAIL TRUE
)

add_test(
    NAME lowmc_hepmc3_afterburner_charge_conserving_fixture
    COMMAND lowmc_hepmc3_afterburner_diagnostic
            --max-events 1
            --multiplicity-seed 2
            --shape-seed 67890
            --min-shape-trials 1000
            --shape-trials-per-photon 200
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/charge_conserving_event.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/charge_conserving_with_low_photon.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_charge_conserving_fixture
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "event_index=0 momentum_unit=MEV selected_incoming_charged=2 selected_outgoing_charged=2 ignored_neutral_candidates=1 rejected_charged_candidates=0 signed_charge_sum=0 ordinary_charge_sum=4 photon_status=generated photon_mode=unweighted enhancement_factor=1 photon_weight=1.*requested_photons=1 retained_photons=1 photons_added=1"
)

add_test(
    NAME lowmc_hepmc3_afterburner_output_file_contains_low_photon
    COMMAND lowmc_hepmc3_output_file_tests
            ${CMAKE_CURRENT_BINARY_DIR}/charge_conserving_with_low_photon.hepmc3
            unweighted
            1
            1
            1
)
set_tests_properties(
    lowmc_hepmc3_afterburner_output_file_contains_low_photon
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_charge_conserving_fixture
)

add_test(
    NAME lowmc_hepmc3_afterburner_enhanced_charge_conserving_fixture
    COMMAND lowmc_hepmc3_afterburner_diagnostic
            --max-events 1
            --multiplicity-seed 2
            --shape-seed 67890
            --min-shape-trials 1000
            --shape-trials-per-photon 200
            --enhancement-factor 10
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/charge_conserving_event.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/charge_conserving_with_enhanced_low_photons.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_enhanced_charge_conserving_fixture
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "event_index=0 momentum_unit=MEV selected_incoming_charged=2 selected_outgoing_charged=2 ignored_neutral_candidates=1 rejected_charged_candidates=0 signed_charge_sum=0 ordinary_charge_sum=4 photon_status=generated photon_mode=enhanced enhancement_factor=10 photon_weight=0.1.*requested_photons=4 retained_photons=4 photons_added=4"
)

add_test(
    NAME lowmc_hepmc3_afterburner_enhanced_output_file_contains_weighted_low_photons
    COMMAND lowmc_hepmc3_output_file_tests
            ${CMAKE_CURRENT_BINARY_DIR}/charge_conserving_with_enhanced_low_photons.hepmc3
            enhanced
            10
            0.1
            4
)
set_tests_properties(
    lowmc_hepmc3_afterburner_enhanced_output_file_contains_weighted_low_photons
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_enhanced_charge_conserving_fixture
)

add_test(
    NAME lowmc_hepmc3_afterburner_atomic_success_fixture
    COMMAND lowmc_hepmc3_afterburner_diagnostic
            --atomic-output
            --max-events 1
            --multiplicity-seed 2
            --shape-seed 67890
            --min-shape-trials 1000
            --shape-trials-per-photon 200
            --enhancement-factor 10
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/charge_conserving_event.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/charge_conserving_atomic_low_photons.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_atomic_success_fixture
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "event_index=0 momentum_unit=MEV selected_incoming_charged=2 selected_outgoing_charged=2 ignored_neutral_candidates=1 rejected_charged_candidates=0 signed_charge_sum=0 ordinary_charge_sum=4 photon_status=generated photon_mode=enhanced enhancement_factor=10 photon_weight=0.1 lambda_event=0.144148 requested_photons=4 retained_photons=4 photons_added=4 shape_status=ok"
)

add_test(
    NAME lowmc_hepmc3_afterburner_atomic_success_output_file
    COMMAND lowmc_hepmc3_output_file_tests
            ${CMAKE_CURRENT_BINARY_DIR}/charge_conserving_atomic_low_photons.hepmc3
            enhanced
            10
            0.1
            4
)
set_tests_properties(
    lowmc_hepmc3_afterburner_atomic_success_output_file
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_atomic_success_fixture
)

add_test(
    NAME lowmc_hepmc3_afterburner_boosted_p100_fixture
    COMMAND lowmc_hepmc3_afterburner_diagnostic
            --max-events 1
            --multiplicity-seed 2
            --shape-seed 67890
            --min-shape-trials 1000
            --shape-trials-per-photon 200
            --enhancement-factor 10
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/boosted_transverse_p100_event.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/boosted_transverse_p100_low_photons.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_boosted_p100_fixture
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "event_index=0 momentum_unit=MEV selected_incoming_charged=2 selected_outgoing_charged=2 ignored_neutral_candidates=1 rejected_charged_candidates=0 signed_charge_sum=0 ordinary_charge_sum=4 photon_status=generated photon_mode=enhanced enhancement_factor=10 photon_weight=0.1.*lambda_event=0.362493 requested_photons=5 retained_photons=5 photons_added=5 shape_status=ok"
)

add_test(
    NAME lowmc_hepmc3_afterburner_boosted_p100_output_file
    COMMAND lowmc_hepmc3_output_file_tests
            ${CMAKE_CURRENT_BINARY_DIR}/boosted_transverse_p100_low_photons.hepmc3
            enhanced
            10
            0.1
            5
            0.362493
)
set_tests_properties(
    lowmc_hepmc3_afterburner_boosted_p100_output_file
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_boosted_p100_fixture
)

add_test(
    NAME lowmc_hepmc3_output_summary_boosted_p100
    COMMAND lowmc_hepmc3_output_summary_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/boosted_transverse_p100_low_photons.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_output_summary_boosted_p100
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_boosted_p100_fixture
        PASS_REGULAR_EXPRESSION
            "event_index=0 lowmc_photons=5 photon_mode=enhanced enhancement_factor=10 photon_weight_sum=0.5 photon_weight_attributes=5 event_metadata=1"
)

add_test(
    NAME lowmc_hepmc3_afterburner_shape_envelope_failure_skip
    COMMAND lowmc_hepmc3_afterburner_diagnostic
            --max-events 1
            --multiplicity-seed 2
            --shape-seed 67890
            --min-shape-trials 10000
            --shape-trials-per-photon 200
            --enhancement-factor 10
            --shape-safety-factor 0.1
            --boosted-shape-safety-factor 0.1
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/boosted_transverse_p100_event.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/shape_envelope_failure_skip_output.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_shape_envelope_failure_skip
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "event_index=0 momentum_unit=MEV selected_incoming_charged=2 selected_outgoing_charged=2 ignored_neutral_candidates=1 rejected_charged_candidates=0 signed_charge_sum=0 ordinary_charge_sum=4 photon_status=skipped_failed_shape_envelope.*requested_photons=5 retained_photons=5 photons_added=0 shape_status=envelope_failed max_trial_bound_ratio=1[0-9]"
)

add_test(
    NAME lowmc_hepmc3_afterburner_shape_envelope_failure_output_unchanged
    COMMAND lowmc_hepmc3_charge_conserving_no_lowmc_output_tests
            ${CMAKE_CURRENT_BINARY_DIR}/shape_envelope_failure_skip_output.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_shape_envelope_failure_output_unchanged
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_shape_envelope_failure_skip
)

add_test(
    NAME lowmc_hepmc3_output_summary_shape_envelope_failure
    COMMAND lowmc_hepmc3_output_summary_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/shape_envelope_failure_skip_output.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_output_summary_shape_envelope_failure
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_shape_envelope_failure_skip
        PASS_REGULAR_EXPRESSION
            "event_index=0 lowmc_photons=0 photon_mode=none enhancement_factor=0 photon_weight_sum=0 photon_weight_attributes=0 event_metadata=0"
)

add_test(
    NAME lowmc_hepmc3_afterburner_shape_envelope_failure_fail_policy
    COMMAND lowmc_hepmc3_afterburner_diagnostic
            --failure-policy fail
            --max-events 1
            --multiplicity-seed 2
            --shape-seed 67890
            --min-shape-trials 10000
            --shape-trials-per-photon 200
            --enhancement-factor 10
            --shape-safety-factor 0.1
            --boosted-shape-safety-factor 0.1
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/boosted_transverse_p100_event.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/shape_envelope_failure_fail_output.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_shape_envelope_failure_fail_policy
    PROPERTIES
        WILL_FAIL TRUE
)

add_test(
    NAME lowmc_hepmc3_afterburner_forward_eta5_fixture
    COMMAND lowmc_hepmc3_afterburner_diagnostic
            --max-events 1
            --multiplicity-seed 2
            --shape-seed 67890
            --min-shape-trials 1000
            --shape-trials-per-photon 200
            --enhancement-factor 10
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/forward_pion_pt1_eta5_event.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/forward_pion_pt1_eta5_low_photons.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_forward_eta5_fixture
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "event_index=0 momentum_unit=MEV selected_incoming_charged=2 selected_outgoing_charged=2 ignored_neutral_candidates=1 rejected_charged_candidates=0 signed_charge_sum=0 ordinary_charge_sum=0 photon_status=generated photon_mode=enhanced enhancement_factor=10 photon_weight=0.1.*lambda_event=0.0351479 requested_photons=2 retained_photons=2 photons_added=2 shape_status=ok"
)

add_test(
    NAME lowmc_hepmc3_afterburner_forward_eta5_output_file
    COMMAND lowmc_hepmc3_output_file_tests
            ${CMAKE_CURRENT_BINARY_DIR}/forward_pion_pt1_eta5_low_photons.hepmc3
            enhanced
            10
            0.1
            2
            0.0351479
)
set_tests_properties(
    lowmc_hepmc3_afterburner_forward_eta5_output_file
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_forward_eta5_fixture
)

add_test(
    NAME lowmc_hepmc3_output_summary_forward_eta5
    COMMAND lowmc_hepmc3_output_summary_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/forward_pion_pt1_eta5_low_photons.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_output_summary_forward_eta5
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_forward_eta5_fixture
        PASS_REGULAR_EXPRESSION
            "event_index=0 lowmc_photons=2 photon_mode=enhanced enhancement_factor=10 photon_weight_sum=0.2 photon_weight_attributes=2 event_metadata=1"
)

add_test(
    NAME lowmc_hepmc3_afterburner_extreme_p1000_fixture
    COMMAND lowmc_hepmc3_afterburner_diagnostic
            --max-events 1
            --multiplicity-seed 2
            --shape-seed 67890
            --min-shape-trials 1000
            --shape-trials-per-photon 200
            --enhancement-factor 10
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/extreme_pion_p1000_event.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/extreme_pion_p1000_low_photons.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_extreme_p1000_fixture
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "event_index=0 momentum_unit=MEV selected_incoming_charged=2 selected_outgoing_charged=2 ignored_neutral_candidates=1 rejected_charged_candidates=0 signed_charge_sum=0 ordinary_charge_sum=0 photon_status=generated photon_mode=enhanced enhancement_factor=10 photon_weight=0.1.*lambda_event=0.579429 requested_photons=7 retained_photons=7 photons_added=7 shape_status=ok"
)

add_test(
    NAME lowmc_hepmc3_afterburner_extreme_p1000_output_file
    COMMAND lowmc_hepmc3_output_file_tests
            ${CMAKE_CURRENT_BINARY_DIR}/extreme_pion_p1000_low_photons.hepmc3
            enhanced
            10
            0.1
            7
            0.579429
)
set_tests_properties(
    lowmc_hepmc3_afterburner_extreme_p1000_output_file
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_extreme_p1000_fixture
)

add_test(
    NAME lowmc_hepmc3_output_summary_extreme_p1000
    COMMAND lowmc_hepmc3_output_summary_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/extreme_pion_p1000_low_photons.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_output_summary_extreme_p1000
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_extreme_p1000_fixture
        PASS_REGULAR_EXPRESSION
            "event_index=0 lowmc_photons=7 photon_mode=enhanced enhancement_factor=10 photon_weight_sum=0.7 photon_weight_attributes=7 event_metadata=1"
)

add_test(
    NAME lowmc_hepmc3_output_summary_single_enhanced
    COMMAND lowmc_hepmc3_output_summary_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/charge_conserving_with_enhanced_low_photons.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_output_summary_single_enhanced
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_enhanced_charge_conserving_fixture
        PASS_REGULAR_EXPRESSION
            "event_index=0 lowmc_photons=4 photon_mode=enhanced enhancement_factor=10 photon_weight_sum=0.4 photon_weight_attributes=4 event_metadata=1"
)

add_test(
    NAME lowmc_hepmc3_afterburner_multi_event_physical_fixture
    COMMAND lowmc_hepmc3_afterburner_diagnostic
            --multiplicity-seed 2
            --shape-seed 67890
            --min-shape-trials 1000
            --shape-trials-per-photon 200
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/multi_event_mixed.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/multi_event_physical_low_photons.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_multi_event_physical_fixture
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "events_read=2 events_written=2 photons_added=1"
)

add_test(
    NAME lowmc_hepmc3_afterburner_multi_event_physical_output_file
    COMMAND lowmc_hepmc3_multi_event_output_tests
            ${CMAKE_CURRENT_BINARY_DIR}/multi_event_physical_low_photons.hepmc3
            unweighted
            1
            1
            1
)
set_tests_properties(
    lowmc_hepmc3_afterburner_multi_event_physical_output_file
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_multi_event_physical_fixture
)

add_test(
    NAME lowmc_hepmc3_output_summary_multi_event_physical
    COMMAND lowmc_hepmc3_output_summary_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/multi_event_physical_low_photons.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_output_summary_multi_event_physical
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_multi_event_physical_fixture
        PASS_REGULAR_EXPRESSION
            "event_index=0 lowmc_photons=1 photon_mode=unweighted enhancement_factor=1 photon_weight_sum=1 photon_weight_attributes=0 event_metadata=1"
)

add_test(
    NAME lowmc_hepmc3_afterburner_multi_event_enhanced_fixture
    COMMAND lowmc_hepmc3_afterburner_diagnostic
            --multiplicity-seed 2
            --shape-seed 67890
            --min-shape-trials 1000
            --shape-trials-per-photon 200
            --enhancement-factor 10
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/multi_event_mixed.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/multi_event_enhanced_low_photons.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_multi_event_enhanced_fixture
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "events_read=2 events_written=2 photons_added=4"
)

add_test(
    NAME lowmc_hepmc3_afterburner_multi_event_enhanced_output_file
    COMMAND lowmc_hepmc3_multi_event_output_tests
            ${CMAKE_CURRENT_BINARY_DIR}/multi_event_enhanced_low_photons.hepmc3
            enhanced
            10
            0.1
            4
)
set_tests_properties(
    lowmc_hepmc3_afterburner_multi_event_enhanced_output_file
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_multi_event_enhanced_fixture
)

add_test(
    NAME lowmc_hepmc3_output_summary_multi_event_enhanced
    COMMAND lowmc_hepmc3_output_summary_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/multi_event_enhanced_low_photons.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_output_summary_multi_event_enhanced
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_multi_event_enhanced_fixture
        PASS_REGULAR_EXPRESSION
            "event_index=1 lowmc_photons=0 photon_mode=none enhancement_factor=0 photon_weight_sum=0 photon_weight_attributes=0 event_metadata=0"
)

add_test(
    NAME lowmc_hepmc3_afterburner_realistic_multi_event_fixture
    COMMAND lowmc_hepmc3_afterburner_diagnostic
            --multiplicity-seed 2
            --shape-seed 67890
            --min-shape-trials 1000
            --shape-trials-per-photon 200
            --enhancement-factor 10
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/realistic_mixed_three_event.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/realistic_mixed_three_event_low_photons.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_realistic_multi_event_fixture
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "events_read=3 events_written=3 photons_added=7"
)

add_test(
    NAME lowmc_hepmc3_afterburner_realistic_multi_event_output_file
    COMMAND lowmc_hepmc3_realistic_multi_event_output_tests
            ${CMAKE_CURRENT_BINARY_DIR}/realistic_mixed_three_event_low_photons.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_realistic_multi_event_output_file
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_realistic_multi_event_fixture
)

add_test(
    NAME lowmc_hepmc3_output_summary_realistic_multi_event
    COMMAND lowmc_hepmc3_output_summary_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/realistic_mixed_three_event_low_photons.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_output_summary_realistic_multi_event
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_realistic_multi_event_fixture
        PASS_REGULAR_EXPRESSION
            "event_index=2 lowmc_photons=0 photon_mode=none enhancement_factor=0 photon_weight_sum=0 photon_weight_attributes=0 event_metadata=0"
)

add_test(
    NAME lowmc_hepmc3_generated_realistic_sample_fixture
    COMMAND lowmc_hepmc3_generate_realistic_sample
            ${CMAKE_CURRENT_BINARY_DIR}/generated_realistic_sample.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_generated_realistic_sample_fixture
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "generated_hepmc3_realistic_sample events=10"
)

add_test(
    NAME lowmc_hepmc3_afterburner_generated_realistic_sample
    COMMAND lowmc_hepmc3_afterburner_diagnostic
            --multiplicity-seed 2
            --shape-seed 67890
            --min-shape-trials 800
            --shape-trials-per-photon 120
            --enhancement-factor 10
            ${CMAKE_CURRENT_BINARY_DIR}/generated_realistic_sample.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/generated_realistic_sample_low_photons.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_generated_realistic_sample
    PROPERTIES
        DEPENDS lowmc_hepmc3_generated_realistic_sample_fixture
        PASS_REGULAR_EXPRESSION
            "events_read=10 events_written=10 photons_added=36"
)

add_test(
    NAME lowmc_hepmc3_generated_realistic_sample_output_file
    COMMAND lowmc_hepmc3_generated_sample_output_tests
            ${CMAKE_CURRENT_BINARY_DIR}/generated_realistic_sample_low_photons.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_generated_realistic_sample_output_file
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_generated_realistic_sample
)

add_test(
    NAME lowmc_hepmc3_output_summary_generated_realistic_sample
    COMMAND lowmc_hepmc3_output_summary_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/generated_realistic_sample_low_photons.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_output_summary_generated_realistic_sample
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_generated_realistic_sample
        PASS_REGULAR_EXPRESSION
            "event_index=9"
)

add_test(
    NAME lowmc_hepmc3_afterburner_realistic_multi_event_fail_policy
    COMMAND lowmc_hepmc3_afterburner_diagnostic
            --failure-policy fail
            --multiplicity-seed 2
            --shape-seed 67890
            --min-shape-trials 1000
            --shape-trials-per-photon 200
            --enhancement-factor 10
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/realistic_mixed_three_event.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/realistic_mixed_three_event_fail_policy_output.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_realistic_multi_event_fail_policy
    PROPERTIES
        WILL_FAIL TRUE
)

add_test(
    NAME lowmc_hepmc3_output_summary_realistic_multi_event_fail_policy_partial_output
    COMMAND lowmc_hepmc3_output_summary_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/realistic_mixed_three_event_fail_policy_output.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_output_summary_realistic_multi_event_fail_policy_partial_output
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_realistic_multi_event_fail_policy
        PASS_REGULAR_EXPRESSION
            "event_index=1 lowmc_photons=2 photon_mode=enhanced enhancement_factor=10 photon_weight_sum=0.2 photon_weight_attributes=2 event_metadata=1"
)

add_test(
    NAME lowmc_hepmc3_afterburner_atomic_fail_policy_clean
    COMMAND ${CMAKE_COMMAND} -E rm -f
            ${CMAKE_CURRENT_BINARY_DIR}/realistic_mixed_three_event_atomic_fail_output.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/realistic_mixed_three_event_atomic_fail_output.hepmc3.tmp
)

add_test(
    NAME lowmc_hepmc3_afterburner_atomic_fail_policy
    COMMAND lowmc_hepmc3_afterburner_diagnostic
            --atomic-output
            --failure-policy fail
            --multiplicity-seed 2
            --shape-seed 67890
            --min-shape-trials 1000
            --shape-trials-per-photon 200
            --enhancement-factor 10
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/realistic_mixed_three_event.hepmc3
            ${CMAKE_CURRENT_BINARY_DIR}/realistic_mixed_three_event_atomic_fail_output.hepmc3
)
set_tests_properties(
    lowmc_hepmc3_afterburner_atomic_fail_policy
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_atomic_fail_policy_clean
        WILL_FAIL TRUE
)

add_test(
    NAME lowmc_hepmc3_afterburner_atomic_fail_policy_no_output
    COMMAND ${CMAKE_COMMAND}
            -DCHECK_PATH=${CMAKE_CURRENT_BINARY_DIR}/realistic_mixed_three_event_atomic_fail_output.hepmc3
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_file_absent.cmake
)
set_tests_properties(
    lowmc_hepmc3_afterburner_atomic_fail_policy_no_output
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_atomic_fail_policy
)

add_test(
    NAME lowmc_hepmc3_afterburner_atomic_fail_policy_no_temporary_output
    COMMAND ${CMAKE_COMMAND}
            -DCHECK_PATH=${CMAKE_CURRENT_BINARY_DIR}/realistic_mixed_three_event_atomic_fail_output.hepmc3.tmp
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_file_absent.cmake
)
set_tests_properties(
    lowmc_hepmc3_afterburner_atomic_fail_policy_no_temporary_output
    PROPERTIES
        DEPENDS lowmc_hepmc3_afterburner_atomic_fail_policy
)

if(Python3_Interpreter_FOUND)
    add_test(
        NAME lowmc_hepmc3_count_yield_smoke_help
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/../../scripts/hepmc3_count_yield_smoke.py
                --help
    )

    add_test(
        NAME lowmc_hepmc3_count_yield_highstat_help
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/../../scripts/hepmc3_count_yield_highstat.py
                --help
    )

    add_test(
        NAME lowmc_hepmc3_count_yield_smoke
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/../../scripts/hepmc3_count_yield_smoke.py
                --tool-dir $<TARGET_FILE_DIR:lowmc_hepmc3_afterburner_diagnostic>
                --input ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/charge_conserving_event.hepmc3
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/realistic_mixed_three_event.hepmc3
                --output-dir ${CMAKE_CURRENT_BINARY_DIR}/count_yield_smoke
                --seeds 1 2 3 4
                --enhancement-factors 1 10
    )
    set_tests_properties(
        lowmc_hepmc3_count_yield_smoke
        PROPERTIES
            PASS_REGULAR_EXPRESSION
                "aggregate cases=16 events=32 .* requested_photons="
    )

    add_test(
        NAME lowmc_hepmc3_count_yield_generated_sample
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/../../scripts/hepmc3_count_yield_smoke.py
                --tool-dir $<TARGET_FILE_DIR:lowmc_hepmc3_afterburner_diagnostic>
                --input ${CMAKE_CURRENT_BINARY_DIR}/generated_realistic_sample.hepmc3
                --output-dir ${CMAKE_CURRENT_BINARY_DIR}/count_yield_generated_sample
                --seeds 1 2
                --enhancement-factors 10 25
                --max-events 2
                --max-poisson-sigma 2
    )
    set_tests_properties(
        lowmc_hepmc3_count_yield_generated_sample
        PROPERTIES
            DEPENDS lowmc_hepmc3_generated_realistic_sample_fixture
            PASS_REGULAR_EXPRESSION
                "aggregate cases=4 events=8 .* requested_photons="
    )

    add_test(
        NAME lowmc_hepmc3_quickstart_help
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/../../scripts/hepmc3_quickstart.py
                --help
    )

    add_test(
        NAME lowmc_hepmc3_quickstart_smoke
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/../../scripts/hepmc3_quickstart.py
                --tool-dir $<TARGET_FILE_DIR:lowmc_hepmc3_afterburner_diagnostic>
                --input ${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/charge_conserving_event.hepmc3
                --output-dir ${CMAKE_CURRENT_BINARY_DIR}/quickstart
    )
    set_tests_properties(
        lowmc_hepmc3_quickstart_smoke
        PROPERTIES
            PASS_REGULAR_EXPRESSION
                "atomic_output=1.*event_index=0 lowmc_photons=4 photon_mode=enhanced enhancement_factor=10 photon_weight_sum=0.4 photon_weight_attributes=4 event_metadata=1"
    )

    add_test(
        NAME lowmc_hepmc3_installed_examples_install
        COMMAND ${CMAKE_COMMAND}
                --install ${CMAKE_BINARY_DIR}
                --prefix ${CMAKE_CURRENT_BINARY_DIR}/installed_examples_prefix
    )

    add_test(
        NAME lowmc_hepmc3_installed_quickstart_smoke
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_BINARY_DIR}/installed_examples_prefix/${CMAKE_INSTALL_DATADIR}/lowmc/examples/hepmc3/hepmc3_quickstart.py
                --output-dir ${CMAKE_CURRENT_BINARY_DIR}/installed_examples_quickstart
    )
    set_tests_properties(
        lowmc_hepmc3_installed_quickstart_smoke
        PROPERTIES
            DEPENDS lowmc_hepmc3_installed_examples_install
            PASS_REGULAR_EXPRESSION
                "atomic_output=1.*event_index=0 lowmc_photons=4 photon_mode=enhanced enhancement_factor=10 photon_weight_sum=0.4 photon_weight_attributes=4 event_metadata=1"
    )
endif()
