target_sources(FretLab PUBLIC 
#${CMAKE_CURRENT_SOURCE_DIR}/general/general_module.f90
${CMAKE_CURRENT_SOURCE_DIR}/tools/parallel_module.f90
${CMAKE_CURRENT_SOURCE_DIR}/algorithm/algorithm_module.f90
${CMAKE_CURRENT_SOURCE_DIR}/density_nanoparticle/integrals.f90
${CMAKE_CURRENT_SOURCE_DIR}/density_nanoparticle/nanoparticle.f90
${CMAKE_CURRENT_SOURCE_DIR}/density_nanoparticle/density.f90
${CMAKE_CURRENT_SOURCE_DIR}/input_output/output_module.f90
${CMAKE_CURRENT_SOURCE_DIR}/input_output/input_module.f90
#${CMAKE_CURRENT_SOURCE_DIR}/tools/array_manipulation_module.f90
${CMAKE_CURRENT_SOURCE_DIR}/tools/string_manipulation_module.f90
${CMAKE_CURRENT_SOURCE_DIR}/tools/time_module.f90
${CMAKE_CURRENT_SOURCE_DIR}/tools/parameters_module.f90
${CMAKE_CURRENT_SOURCE_DIR}/target/target_module.f90
${CMAKE_CURRENT_SOURCE_DIR}/FretLab.f90
)
