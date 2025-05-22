macro(add_FRET_Embedlab_runtest _name _labels)
    add_test(
        ${_name}
        python3 ${PROJECT_BINARY_DIR}/tests/${_name}/test --binary-dir=${PROJECT_BINARY_DIR}  --work-dir=${PROJECT_BINARY_DIR}/tests/${_name} --verbose)
    if(NOT "${_labels}" STREQUAL "")
        set_tests_properties(${_name} PROPERTIES LABELS "${_labels}")
    endif()
endmacro()

# All tests here should contain the label "FRET_Embedlab"

# Add a keyword for the length of the test: 
# 
# 	short < 30 seconds
# 	medium > 30 seconds < 120 seconds
# 	long > 120 seconds < 200 seconds
# 	verylong > 200 seconds
# 
# NEVER comment out tests
add_FRET_Embedlab_runtest(integrate_density                     "FRET_Embedlab;density_integral")
add_FRET_Embedlab_runtest(acceptor_donor_with_overlap_integral  "FRET_Embedlab;aceptor_donor_charges_overlap;")
add_FRET_Embedlab_runtest(acceptor_donor_coulomb                "FRET_Embedlab;aceptor_donor_coulomb;")
add_FRET_Embedlab_runtest(acceptor_np_charges                   "FRET_Embedlab;aceptor_np_charges;")

add_FRET_Embedlab_runtest(acceptor_np_charges_donor_with_overlap_integral     "FRET_Embedlab;aceptor_np_donor_charges_overlap;")

