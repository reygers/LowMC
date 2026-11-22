add_library(lowmc_core STATIC
    src/AnalyticLowBound.cpp
    src/Balek.cpp
    src/Dielectron.cpp
    src/EventYield.cpp
    src/FixedKTWeightedSampler.cpp
    src/GridEtaPhiProposal.cpp
    src/LegCenteredEtaPhiProposal.cpp
    src/LowCurrent.cpp
    src/ShapeUnweighting.cpp
    src/WeightedEventGenerator.cpp
    src/WeightedPhotonGenerator.cpp
    src/WeightedPhotonSampler.cpp
)

add_library(LowMC::core ALIAS lowmc_core)

target_compile_features(lowmc_core PUBLIC cxx_std_17)

target_include_directories(lowmc_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set_target_properties(lowmc_core PROPERTIES EXPORT_NAME core)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_core PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_core PRIVATE /W4)
endif()

add_executable(lowmc_core_tests
    tests/test_low_current.cpp
)

target_link_libraries(lowmc_core_tests
    PRIVATE
        lowmc_core
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_core_tests PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_core_tests PRIVATE /W4)
endif()

add_test(NAME lowmc_core_tests COMMAND lowmc_core_tests)

add_executable(lowmc_two_particle_diagnostic
    diagnostics/two_particle_eta_phi.cpp
)

target_link_libraries(lowmc_two_particle_diagnostic
    PRIVATE
        lowmc_core
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_two_particle_diagnostic PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_two_particle_diagnostic PRIVATE /W4)
endif()

add_test(
    NAME lowmc_two_particle_diagnostic_smoke
    COMMAND lowmc_two_particle_diagnostic ${CMAKE_CURRENT_BINARY_DIR}/two_particle_eta_phi.csv
)

add_test(
    NAME lowmc_two_particle_diagnostic_custom_grid_smoke
    COMMAND lowmc_two_particle_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/two_particle_eta_phi_custom_grid.csv
            --eta-bins 5
            --phi-bins 7
)

add_test(
    NAME lowmc_two_particle_diagnostic_asymmetric_smoke
    COMMAND lowmc_two_particle_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/two_particle_eta_phi_asymmetric.csv
            --eta-bins 5
            --phi-bins 7
            --momentum 100
            --particle-eta 0.35
            --particle-phi 1.1
            --second-momentum 35
            --second-particle-eta -0.75
            --second-particle-phi 4.0
)

add_test(
    NAME lowmc_two_particle_diagnostic_event4_smoke
    COMMAND lowmc_two_particle_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/two_particle_eta_phi_event4.csv
            --topology event4
            --eta-bins 5
            --phi-bins 7
)

add_test(
    NAME lowmc_two_particle_diagnostic_particles_file_smoke
    COMMAND lowmc_two_particle_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/two_particle_eta_phi_particles.csv
            --particles ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/event4_particles.lowmc
            --eta-bins 5
            --phi-bins 7
)

add_test(
    NAME lowmc_two_particle_diagnostic_clustered_particles_file_smoke
    COMMAND lowmc_two_particle_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/two_particle_eta_phi_clustered_particles.csv
            --particles ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/clustered8_particles.lowmc
            --eta-bins 5
            --phi-bins 7
)

add_test(
    NAME lowmc_two_particle_diagnostic_incoming_outgoing_particles_file_smoke
    COMMAND lowmc_two_particle_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/two_particle_eta_phi_incoming_outgoing_particles.csv
            --particles ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/incoming_outgoing4_particles.lowmc
            --eta-bins 5
            --phi-bins 7
)

add_test(
    NAME lowmc_two_particle_diagnostic_hepmc3_validation_particles_file_smoke
    COMMAND lowmc_two_particle_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/two_particle_eta_phi_hepmc3_validation_particles.csv
            --particles ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/hepmc3_validation_event_particles.lowmc
            --eta-bins 5
            --phi-bins 7
)

add_test(
    NAME lowmc_two_particle_diagnostic_forward4_smoke
    COMMAND lowmc_two_particle_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/two_particle_eta_phi_forward4.csv
            --topology forward4
            --eta-bins 5
            --phi-bins 7
)

add_test(
    NAME lowmc_two_particle_diagnostic_multi8_smoke
    COMMAND lowmc_two_particle_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/two_particle_eta_phi_multi8.csv
            --topology multi8
            --eta-bins 5
            --phi-bins 7
)

add_test(
    NAME lowmc_two_particle_diagnostic_rejects_invalid_phase_space
    COMMAND lowmc_two_particle_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/invalid_two_particle_eta_phi.csv
            --eta-min 1.0
            --eta-max 1.0
)
set_tests_properties(
    lowmc_two_particle_diagnostic_rejects_invalid_phase_space
    PROPERTIES WILL_FAIL TRUE
)

add_test(
    NAME lowmc_two_particle_diagnostic_rejects_invalid_grid
    COMMAND lowmc_two_particle_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/invalid_two_particle_eta_phi_grid.csv
            --eta-bins 1
)
set_tests_properties(
    lowmc_two_particle_diagnostic_rejects_invalid_grid
    PROPERTIES WILL_FAIL TRUE
)

add_executable(lowmc_fixed_kt_weighted_samples
    diagnostics/fixed_kt_weighted_samples.cpp
)

target_link_libraries(lowmc_fixed_kt_weighted_samples
    PRIVATE
        lowmc_core
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_fixed_kt_weighted_samples PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_fixed_kt_weighted_samples PRIVATE /W4)
endif()

add_test(
    NAME lowmc_fixed_kt_weighted_samples_smoke
    COMMAND lowmc_fixed_kt_weighted_samples
            ${CMAKE_CURRENT_BINARY_DIR}/fixed_kt_weighted_samples.csv
            --samples 8
)

add_test(
    NAME lowmc_fixed_kt_adaptive_weighted_samples_smoke
    COMMAND lowmc_fixed_kt_weighted_samples
            ${CMAKE_CURRENT_BINARY_DIR}/fixed_kt_adaptive_weighted_samples.csv
            --samples 8
            --proposal adaptive
            --adaptive-uniform-fraction 0.25
            --adaptive-width-scale 10
            --adaptive-min-half-width 0.02
            --adaptive-max-half-width 0.8
)

add_test(
    NAME lowmc_fixed_kt_smooth_weighted_samples_smoke
    COMMAND lowmc_fixed_kt_weighted_samples
            ${CMAKE_CURRENT_BINARY_DIR}/fixed_kt_smooth_weighted_samples.csv
            --samples 8
            --proposal smooth
            --smooth-uniform-fraction 0.25
            --smooth-width-scale 10
            --smooth-min-half-width 0.02
            --smooth-max-half-width 0.8
)

add_test(
    NAME lowmc_fixed_kt_leg_weighted_samples_smoke
    COMMAND lowmc_fixed_kt_weighted_samples
            ${CMAKE_CURRENT_BINARY_DIR}/fixed_kt_leg_weighted_samples.csv
            --samples 8
            --proposal leg
            --leg-uniform-fraction 0.25
)

add_test(
    NAME lowmc_fixed_kt_leg_particles_file_smoke
    COMMAND lowmc_fixed_kt_weighted_samples
            ${CMAKE_CURRENT_BINARY_DIR}/fixed_kt_leg_particles_file.csv
            --samples 8
            --proposal leg
            --particles ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/event4_particles.lowmc
)

add_test(
    NAME lowmc_fixed_kt_grid_weighted_samples_smoke
    COMMAND lowmc_fixed_kt_weighted_samples
            ${CMAKE_CURRENT_BINARY_DIR}/fixed_kt_grid_weighted_samples.csv
            --samples 8
            --proposal grid
            --grid-eta-bins 4
            --grid-phi-bins 4
            --grid-refinement-levels 1
            --grid-refine-fraction 0.25
            --grid-uniform-fraction 0.2
)

add_executable(lowmc_weighted_sampler_study
    diagnostics/weighted_sampler_study.cpp
)

target_link_libraries(lowmc_weighted_sampler_study
    PRIVATE
        lowmc_core
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_weighted_sampler_study PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_weighted_sampler_study PRIVATE /W4)
endif()

add_test(
    NAME lowmc_weighted_sampler_study_smoke
    COMMAND lowmc_weighted_sampler_study
            ${CMAKE_CURRENT_BINARY_DIR}/weighted_sampler_study.csv
            --samples 8
)

add_test(
    NAME lowmc_weighted_sampler_study_adaptive_smoke
    COMMAND lowmc_weighted_sampler_study
            ${CMAKE_CURRENT_BINARY_DIR}/weighted_sampler_study_adaptive.csv
            --samples 8
            --proposal adaptive
            --adaptive-uniform-fraction 0.25
            --adaptive-width-scale 10
            --adaptive-min-half-width 0.02
            --adaptive-max-half-width 0.8
)

add_test(
    NAME lowmc_weighted_sampler_study_leg_smoke
    COMMAND lowmc_weighted_sampler_study
            ${CMAKE_CURRENT_BINARY_DIR}/weighted_sampler_study_leg.csv
            --samples 8
            --momentum 10
            --proposal leg
            --leg-uniform-fraction 0.25
)

add_test(
    NAME lowmc_weighted_sampler_study_smooth_smoke
    COMMAND lowmc_weighted_sampler_study
            ${CMAKE_CURRENT_BINARY_DIR}/weighted_sampler_study_smooth.csv
            --samples 8
            --proposal smooth
            --smooth-uniform-fraction 0.25
            --smooth-width-scale 10
            --smooth-min-half-width 0.02
            --smooth-max-half-width 0.8
)

add_test(
    NAME lowmc_weighted_sampler_study_grid_smoke
    COMMAND lowmc_weighted_sampler_study
            ${CMAKE_CURRENT_BINARY_DIR}/weighted_sampler_study_grid.csv
            --samples 8
            --proposal grid
            --grid-eta-bins 4
            --grid-phi-bins 4
            --grid-refinement-levels 1
            --grid-refine-fraction 0.25
            --grid-uniform-fraction 0.2
)

add_test(
    NAME lowmc_weighted_sampler_study_grid_variance_smoke
    COMMAND lowmc_weighted_sampler_study
            ${CMAKE_CURRENT_BINARY_DIR}/weighted_sampler_study_grid_variance.csv
            --samples 8
            --proposal grid
            --grid-eta-bins 4
            --grid-phi-bins 4
            --grid-refinement-levels 1
            --grid-refine-fraction 0.25
            --grid-uniform-fraction 0.2
            --grid-refinement-score variance
            --grid-variance-probe-bins 2
)

add_executable(lowmc_full_weighted_sampler_study
    diagnostics/full_weighted_sampler_study.cpp
)

target_link_libraries(lowmc_full_weighted_sampler_study
    PRIVATE
        lowmc_core
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_full_weighted_sampler_study PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_full_weighted_sampler_study PRIVATE /W4)
endif()

add_test(
    NAME lowmc_full_weighted_sampler_study_smoke
    COMMAND lowmc_full_weighted_sampler_study
            ${CMAKE_CURRENT_BINARY_DIR}/full_weighted_sampler_study.csv
            --samples 8
)

add_test(
    NAME lowmc_full_weighted_sampler_study_grid_smoke
    COMMAND lowmc_full_weighted_sampler_study
            ${CMAKE_CURRENT_BINARY_DIR}/full_weighted_sampler_study_grid.csv
            --samples 8
            --proposal grid
            --grid-eta-bins 4
            --grid-phi-bins 4
            --grid-refinement-levels 1
            --grid-refine-fraction 0.25
            --grid-uniform-fraction 0.2
)

add_test(
    NAME lowmc_full_weighted_sampler_study_particles_file_smoke
    COMMAND lowmc_full_weighted_sampler_study
            ${CMAKE_CURRENT_BINARY_DIR}/full_weighted_sampler_study_particles.csv
            --particles ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/event4_particles.lowmc
            --samples 8
            --proposal grid
            --grid-eta-bins 4
            --grid-phi-bins 4
            --grid-refinement-levels 1
            --grid-refine-fraction 0.25
            --grid-uniform-fraction 0.2
)

add_executable(lowmc_weighted_event_example
    diagnostics/weighted_event_example.cpp
)

target_link_libraries(lowmc_weighted_event_example
    PRIVATE
        lowmc_core
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_weighted_event_example PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_weighted_event_example PRIVATE /W4)
endif()

add_test(
    NAME lowmc_weighted_event_example_smoke
    COMMAND lowmc_weighted_event_example
            ${CMAKE_CURRENT_BINARY_DIR}/weighted_event_example.csv
            --samples 8
)

add_executable(lowmc_event_yield_diagnostic
    diagnostics/event_yield_diagnostic.cpp
)

target_link_libraries(lowmc_event_yield_diagnostic
    PRIVATE
        lowmc_core
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_event_yield_diagnostic PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_event_yield_diagnostic PRIVATE /W4)
endif()

add_test(
    NAME lowmc_event_yield_diagnostic_smoke
    COMMAND lowmc_event_yield_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/event_yield_diagnostic.csv
            --log-kt-bins 4
            --eta-bins 8
            --phi-bins 8
            --seed 12345
)
set_tests_properties(
    lowmc_event_yield_diagnostic_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "method=log-kt-eta-phi-midpoint status=ok .* lambda_event="
)

add_test(
    NAME lowmc_event_yield_diagnostic_particles_file_smoke
    COMMAND lowmc_event_yield_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/event_yield_particles.csv
            --particles ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/event4_particles.lowmc
            --log-kt-bins 4
            --eta-bins 8
            --phi-bins 8
            --oversampling-factor 3
            --seed 12345
)
set_tests_properties(
    lowmc_event_yield_diagnostic_particles_file_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "method=log-kt-eta-phi-midpoint status=ok particle_count=4 .* poisson_mean="
)

add_test(
    NAME lowmc_event_yield_diagnostic_adaptive_smoke
    COMMAND lowmc_event_yield_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/event_yield_adaptive.csv
            --method adaptive
            --log-kt-bins 2
            --adaptive-grid-eta-bins 4
            --adaptive-grid-phi-bins 4
            --adaptive-refinement-levels 1
            --adaptive-refine-fraction 0.25
            --adaptive-eta-sub-bins 1
            --adaptive-phi-sub-bins 1
            --seed 12345
)
set_tests_properties(
    lowmc_event_yield_diagnostic_adaptive_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "method=adaptive-log-kt-grid-eta-phi-midpoint status=ok .* grid_cells="
)

add_executable(lowmc_count_shape_diagnostic
    diagnostics/count_shape_diagnostic.cpp
)

target_link_libraries(lowmc_count_shape_diagnostic
    PRIVATE
        lowmc_core
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_count_shape_diagnostic PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_count_shape_diagnostic PRIVATE /W4)
endif()

add_test(
    NAME lowmc_count_shape_diagnostic_help
    COMMAND lowmc_count_shape_diagnostic --help
)

add_test(
    NAME lowmc_count_shape_diagnostic_smoke
    COMMAND lowmc_count_shape_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/count_shape_diagnostic.csv
            --momentum 10
            --oversampling-factor 100
            --multiplicity-seed 12345
            --shape-seed 67890
            --min-shape-trials 1000
            --shape-trials-per-photon 500
)
set_tests_properties(
    lowmc_count_shape_diagnostic_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "status=ok .* lambda_event="
)

add_executable(lowmc_unweighting_diagnostic
    diagnostics/unweighting_diagnostic.cpp
)

target_link_libraries(lowmc_unweighting_diagnostic
    PRIVATE
        lowmc_core
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_unweighting_diagnostic PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_unweighting_diagnostic PRIVATE /W4)
endif()

add_test(
    NAME lowmc_unweighting_diagnostic_help
    COMMAND lowmc_unweighting_diagnostic --help
)

add_test(
    NAME lowmc_unweighting_diagnostic_smoke
    COMMAND lowmc_unweighting_diagnostic
            --calibration-samples 16
            --trial-samples 16
            --safety-factor 100
            --seed 12345
)
set_tests_properties(
    lowmc_unweighting_diagnostic_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "calibration_samples=16 trial_samples=16 .* envelope_status=held"
)

add_test(
    NAME lowmc_unweighting_diagnostic_cell_smoke
    COMMAND lowmc_unweighting_diagnostic
            --calibration-samples 16
            --trial-samples 16
            --safety-factor 100
            --seed 12345
            --envelope-mode cell
)
set_tests_properties(
    lowmc_unweighting_diagnostic_cell_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "calibration_samples=16 trial_samples=16 .* envelope_mode=cell .* envelope_status=held"
)

add_test(
    NAME lowmc_unweighting_diagnostic_probed_cell_smoke
    COMMAND lowmc_unweighting_diagnostic
            --calibration-samples 16
            --trial-samples 16
            --safety-factor 100
            --seed 12345
            --envelope-mode probed-cell
            --probe-bins 2
            --probe-rule edge
)
set_tests_properties(
    lowmc_unweighting_diagnostic_probed_cell_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "calibration_samples=16 trial_samples=16 .* envelope_mode=probed-cell .* envelope_status=held"
)

add_test(
    NAME lowmc_unweighting_diagnostic_forward_smoke
    COMMAND lowmc_unweighting_diagnostic
            --calibration-samples 16
            --trial-samples 16
            --safety-factor 100
            --seed 12345
            --momentum 100
            --particle-eta 4
            --envelope-mode probed-cell
            --probe-bins 2
)
set_tests_properties(
    lowmc_unweighting_diagnostic_forward_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "particle_eta=4 .* envelope_mode=probed-cell .* envelope_status=held"
)

add_executable(lowmc_shape_unweighting_diagnostic
    diagnostics/shape_unweighting_diagnostic.cpp
)

target_link_libraries(lowmc_shape_unweighting_diagnostic
    PRIVATE
        lowmc_core
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_shape_unweighting_diagnostic PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_shape_unweighting_diagnostic PRIVATE /W4)
endif()

add_test(
    NAME lowmc_shape_unweighting_diagnostic_help
    COMMAND lowmc_shape_unweighting_diagnostic --help
)

add_test(
    NAME lowmc_shape_unweighting_diagnostic_smoke
    COMMAND lowmc_shape_unweighting_diagnostic
            --trials 16
            --calibration-trials 16
            --safety-factor 100
            --probe-bins 2
            --seed 12345
            --momentum 100
)
set_tests_properties(
    lowmc_shape_unweighting_diagnostic_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "status=ok .* trials=16 .* accepted_candidates="
)

add_test(
    NAME lowmc_shape_unweighting_true_rapidity_center_smoke
    COMMAND lowmc_shape_unweighting_diagnostic
            --trials 16
            --calibration-trials 16
            --safety-factor 100
            --probe-bins 2
            --probe-policy boosted-cells
            --seed 12345
            --momentum 0.1
            --particle-eta 2
            --true-rapidity-centers
)
set_tests_properties(
    lowmc_shape_unweighting_true_rapidity_center_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "status=ok.*charged_particle_eta_center=rapidity"
)

add_test(
    NAME lowmc_shape_unweighting_analytic_bound_audit_smoke
    COMMAND lowmc_shape_unweighting_diagnostic
            --trials 16
            --calibration-trials 16
            --safety-factor 100
            --probe-bins 2
            --seed 12345
            --momentum 10
            --analytic-bound-audit
)
set_tests_properties(
    lowmc_shape_unweighting_analytic_bound_audit_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "analytic_bound_available=1.*analytic_bound_fallback_cells=0.*analytic_bound_observed_violations=0.*analytic_bound_trial_violations=0"
)

add_test(
    NAME lowmc_shape_unweighting_diagnostic_asymmetric_smoke
    COMMAND lowmc_shape_unweighting_diagnostic
            --trials 16
            --calibration-trials 16
            --safety-factor 100
            --probe-bins 2
            --seed 12345
            --momentum 100
            --particle-eta 0.35
            --particle-phi 1.1
            --second-momentum 35
            --second-particle-eta -0.75
            --second-particle-phi 4.0
)
set_tests_properties(
    lowmc_shape_unweighting_diagnostic_asymmetric_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "status=ok .* signed_charge_sum=0"
)

add_test(
    NAME lowmc_shape_unweighting_diagnostic_event4_smoke
    COMMAND lowmc_shape_unweighting_diagnostic
            --topology event4
            --trials 16
            --calibration-trials 16
            --safety-factor 100
            --probe-bins 2
            --seed 12345
)
set_tests_properties(
    lowmc_shape_unweighting_diagnostic_event4_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "status=ok.*topology=event4.*particle_count=4.*signed_charge_sum=0"
)

add_test(
    NAME lowmc_shape_unweighting_diagnostic_particles_file_smoke
    COMMAND lowmc_shape_unweighting_diagnostic
            --particles ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/event4_particles.lowmc
            --trials 16
            --calibration-trials 16
            --safety-factor 100
            --probe-bins 2
            --seed 12345
)
set_tests_properties(
    lowmc_shape_unweighting_diagnostic_particles_file_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "status=ok.*topology=particles.*particle_count=4.*signed_charge_sum=0"
)

add_test(
    NAME lowmc_shape_unweighting_diagnostic_clustered_particles_file_smoke
    COMMAND lowmc_shape_unweighting_diagnostic
            --particles ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/clustered8_particles.lowmc
            --trials 16
            --calibration-trials 16
            --safety-factor 100
            --probe-bins 2
            --seed 12345
)
set_tests_properties(
    lowmc_shape_unweighting_diagnostic_clustered_particles_file_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "status=ok.*topology=particles.*particle_count=8.*signed_charge_sum=0"
)

add_test(
    NAME lowmc_shape_unweighting_diagnostic_incoming_outgoing_particles_file_smoke
    COMMAND lowmc_shape_unweighting_diagnostic
            --particles ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/incoming_outgoing4_particles.lowmc
            --trials 16
            --calibration-trials 16
            --safety-factor 100
            --probe-bins 2
            --seed 12345
)
set_tests_properties(
    lowmc_shape_unweighting_diagnostic_incoming_outgoing_particles_file_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "status=ok.*topology=particles.*particle_count=4.*signed_charge_sum=0"
)

add_test(
    NAME lowmc_shape_unweighting_diagnostic_hepmc3_validation_particles_file_smoke
    COMMAND lowmc_shape_unweighting_diagnostic
            --particles ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/hepmc3_validation_event_particles.lowmc
            --trials 16
            --calibration-trials 16
            --safety-factor 100
            --probe-bins 2
            --seed 12345
)
set_tests_properties(
    lowmc_shape_unweighting_diagnostic_hepmc3_validation_particles_file_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "status=ok.*topology=particles.*particle_count=4.*signed_charge_sum=0"
)

add_test(
    NAME lowmc_shape_unweighting_diagnostic_forward4_smoke
    COMMAND lowmc_shape_unweighting_diagnostic
            --topology forward4
            --trials 16
            --calibration-trials 16
            --safety-factor 100
            --probe-bins 2
            --seed 12345
)
set_tests_properties(
    lowmc_shape_unweighting_diagnostic_forward4_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "status=ok.*topology=forward4.*particle_count=4.*signed_charge_sum=0"
)

add_test(
    NAME lowmc_shape_unweighting_diagnostic_multi8_smoke
    COMMAND lowmc_shape_unweighting_diagnostic
            --topology multi8
            --trials 16
            --calibration-trials 16
            --safety-factor 100
            --probe-bins 2
            --seed 12345
)
set_tests_properties(
    lowmc_shape_unweighting_diagnostic_multi8_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "status=ok.*topology=multi8.*particle_count=8.*signed_charge_sum=0"
)

add_executable(lowmc_delphi_muon_fsr_diagnostic
    diagnostics/delphi_muon_fsr_diagnostic.cpp
)

target_link_libraries(lowmc_delphi_muon_fsr_diagnostic
    PRIVATE
        lowmc_core
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(lowmc_delphi_muon_fsr_diagnostic PRIVATE -Wall -Wextra -pedantic)
elseif(MSVC)
    target_compile_options(lowmc_delphi_muon_fsr_diagnostic PRIVATE /W4)
endif()

add_test(
    NAME lowmc_delphi_muon_fsr_diagnostic_help
    COMMAND lowmc_delphi_muon_fsr_diagnostic --help
)

add_test(
    NAME lowmc_delphi_muon_fsr_diagnostic_smoke
    COMMAND lowmc_delphi_muon_fsr_diagnostic
            ${CMAKE_CURRENT_BINARY_DIR}/delphi_muon_fsr_samples.csv
            --trials 64
            --calibration-trials 16
            --probe-bins 2
            --safety-factor 100
            --seed 12345
)
set_tests_properties(
    lowmc_delphi_muon_fsr_diagnostic_smoke
    PROPERTIES
        PASS_REGULAR_EXPRESSION
            "status=ok .* gamma=430 .* expected_peak_mrad="
)

if(LOWMC_ENABLE_PHOTOS)
    add_executable(lowmc_photos_zmumu_benchmark
        diagnostics/photos_zmumu_benchmark.cpp
    )

    target_include_directories(lowmc_photos_zmumu_benchmark
        PRIVATE
            ${PHOTOS_INCLUDE_DIR}
    )

    target_link_libraries(lowmc_photos_zmumu_benchmark
        PRIVATE
            lowmc_core
            ${PHOTOS_HEPMC3_LIBRARY}
            ${PHOTOS_CORE_LIBRARY}
            HepMC3::HepMC3
    )

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        target_compile_options(lowmc_photos_zmumu_benchmark PRIVATE -Wall -Wextra -pedantic)
    elseif(MSVC)
        target_compile_options(lowmc_photos_zmumu_benchmark PRIVATE /W4)
    endif()

    add_test(
        NAME lowmc_photos_zmumu_benchmark_help
        COMMAND lowmc_photos_zmumu_benchmark --help
    )

    add_test(
        NAME lowmc_photos_zmumu_benchmark_smoke
        COMMAND lowmc_photos_zmumu_benchmark
                ${CMAKE_CURRENT_BINARY_DIR}/photos_zmumu_photons.csv
                ${CMAKE_CURRENT_BINARY_DIR}/photos_zmumu_summary.csv
                --events 64
                --seed1 1802
                --seed2 9373
                --kt-min 0.05
                --kt-max 0.2
                --eta-max 2
                --infrared-cutoff 0.0003
                --lowmc-samples 128
                --lowmc-refinement-levels 3
                --lowmc-output ${CMAKE_CURRENT_BINARY_DIR}/lowmc_zmumu_photons.csv
    )
    set_tests_properties(
        lowmc_photos_zmumu_benchmark_smoke
        PROPERTIES
            PASS_REGULAR_EXPRESSION "status=ok.*events=64.*multiphoton_events=0"
    )

    add_test(
        NAME lowmc_photos_zmumu_energy_window_smoke
        COMMAND lowmc_photos_zmumu_benchmark
                ${CMAKE_CURRENT_BINARY_DIR}/photos_zmumu_energy_photons.csv
                ${CMAKE_CURRENT_BINARY_DIR}/photos_zmumu_energy_summary.csv
                --events 64
                --energy-max 0.1
                --seed1 1802
                --seed2 9373
                --kt-min 0.05
                --kt-max 0.2
                --eta-max 2
                --infrared-cutoff 0.0003
                --lowmc-samples 64
                --lowmc-refinement-levels 3
                --lowmc-output ${CMAKE_CURRENT_BINARY_DIR}/lowmc_zmumu_energy_photons.csv
    )
    set_tests_properties(
        lowmc_photos_zmumu_energy_window_smoke
        PROPERTIES
            PASS_REGULAR_EXPRESSION "status=ok.*events=64.*multiphoton_events=0"
    )
endif()

function(lowmc_add_core_example target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} PRIVATE lowmc_core)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        target_compile_options(${target} PRIVATE -Wall -Wextra -pedantic)
    elseif(MSVC)
        target_compile_options(${target} PRIVATE /W4)
    endif()
endfunction()

lowmc_add_core_example(
    lowmc_example_low_current
    examples/basic/evaluate_low_current.cpp
)
lowmc_add_core_example(
    lowmc_example_weighted_photons
    examples/basic/generate_weighted_photons.cpp
)
lowmc_add_core_example(
    lowmc_example_weighted_dielectrons
    examples/basic/generate_weighted_dielectrons.cpp
)
lowmc_add_core_example(
    lowmc_example_event_multiplicity
    examples/basic/sample_event_multiplicity.cpp
)
lowmc_add_core_example(
    lowmc_example_shape_unweighting
    examples/advanced/shape_unweighting.cpp
)

add_test(NAME lowmc_example_low_current COMMAND lowmc_example_low_current)
set_tests_properties(
    lowmc_example_low_current
    PROPERTIES PASS_REGULAR_EXPRESSION "charge_conserved=1"
)

add_test(NAME lowmc_example_weighted_photons COMMAND lowmc_example_weighted_photons)
set_tests_properties(
    lowmc_example_weighted_photons
    PROPERTIES PASS_REGULAR_EXPRESSION "weighted_photons=1000"
)

add_test(NAME lowmc_example_weighted_dielectrons COMMAND lowmc_example_weighted_dielectrons)
set_tests_properties(
    lowmc_example_weighted_dielectrons
    PROPERTIES PASS_REGULAR_EXPRESSION "weighted_dielectrons=1000"
)

add_test(NAME lowmc_example_event_multiplicity COMMAND lowmc_example_event_multiplicity)
set_tests_properties(
    lowmc_example_event_multiplicity
    PROPERTIES PASS_REGULAR_EXPRESSION "integration_status=ok"
)

add_test(NAME lowmc_example_shape_unweighting COMMAND lowmc_example_shape_unweighting)
set_tests_properties(
    lowmc_example_shape_unweighting
    PROPERTIES PASS_REGULAR_EXPRESSION "status_ok=1.*envelope_failures=0"
)

if(LOWMC_ENABLE_HEPMC3)
    add_subdirectory(adapters/hepmc3)
endif()

if(LOWMC_ENABLE_PYTHIA8)
    add_subdirectory(adapters/pythia8)
endif()
