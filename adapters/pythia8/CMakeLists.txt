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
