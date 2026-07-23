cmake_minimum_required(VERSION 3.16)

project(LowMCInstallConsumer LANGUAGES CXX)

enable_testing()

find_package(LowMC 0.1 CONFIG REQUIRED)

add_executable(lowmc_install_consumer main.cpp)
target_link_libraries(lowmc_install_consumer PRIVATE LowMC::core)

add_test(NAME lowmc_install_consumer_runs COMMAND lowmc_install_consumer)
