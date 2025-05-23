macro(add_FretLab_runtest _name _labels)
    add_test(
        ${_name}
        python3 ${PROJECT_BINARY_DIR}/tests/${_name}/test --binary-dir=${PROJECT_BINARY_DIR}  --work-dir=${PROJECT_BINARY_DIR}/tests/${_name} --verbose)
    if(NOT "${_labels}" STREQUAL "")
        set_tests_properties(${_name} PROPERTIES LABELS "${_labels}")
    endif()
endmacro()

# All tests here should contain the label "FretLab"

# Add a keyword for the length of the test: 
# 
add_FretLab_runtest(integrate_density                                "FretLab;density_integral")
add_FretLab_runtest(acceptor_donor_with_overlap_integral             "FretLab;aceptor_donor_charges_overlap;")
add_FretLab_runtest(acceptor_donor_coulomb                           "FretLab;aceptor_donor_coulomb;")
add_FretLab_runtest(acceptor_np_charges                              "FretLab;aceptor_np_charges;")
add_FretLab_runtest(acceptor_np_charges_dipoles_donor_coulomb        "FretLab;acceptor_np_charges_dipoles_donor_coulomb;")
add_FretLab_runtest(acceptor_np_charges_donor_coulomb                "FretLab;aceptor_np_charges_donor_coulomb;")
add_FretLab_runtest(acceptor_np_charges_donor_with_overlap_integral  "FretLab;aceptor_np_donor_charges_overlap;")

