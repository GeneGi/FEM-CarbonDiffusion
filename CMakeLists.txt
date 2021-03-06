cmake_minimum_required(VERSION 3.6)
project(FEM_CarbonDiffusion)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    cpp/include/project1.h
    cpp/lib/googletest/codegear/gtest_all.cc
    cpp/lib/googletest/codegear/gtest_link.cc
    cpp/lib/googletest/include/gtest/internal/custom/gtest-port.h
    cpp/lib/googletest/include/gtest/internal/custom/gtest-printers.h
    cpp/lib/googletest/include/gtest/internal/custom/gtest.h
    cpp/lib/googletest/include/gtest/internal/gtest-death-test-internal.h
    cpp/lib/googletest/include/gtest/internal/gtest-filepath.h
    cpp/lib/googletest/include/gtest/internal/gtest-internal.h
    cpp/lib/googletest/include/gtest/internal/gtest-linked_ptr.h
    cpp/lib/googletest/include/gtest/internal/gtest-param-util-generated.h
    cpp/lib/googletest/include/gtest/internal/gtest-param-util.h
    cpp/lib/googletest/include/gtest/internal/gtest-port-arch.h
    cpp/lib/googletest/include/gtest/internal/gtest-port.h
    cpp/lib/googletest/include/gtest/internal/gtest-string.h
    cpp/lib/googletest/include/gtest/internal/gtest-tuple.h
    cpp/lib/googletest/include/gtest/internal/gtest-type-util.h
    cpp/lib/googletest/include/gtest/gtest-death-test.h
    cpp/lib/googletest/include/gtest/gtest-message.h
    cpp/lib/googletest/include/gtest/gtest-param-test.h
    cpp/lib/googletest/include/gtest/gtest-printers.h
    cpp/lib/googletest/include/gtest/gtest-spi.h
    cpp/lib/googletest/include/gtest/gtest-test-part.h
    cpp/lib/googletest/include/gtest/gtest-typed-test.h
    cpp/lib/googletest/include/gtest/gtest.h
    cpp/lib/googletest/include/gtest/gtest_pred_impl.h
    cpp/lib/googletest/include/gtest/gtest_prod.h
    cpp/lib/googletest/samples/prime_tables.h
    cpp/lib/googletest/samples/sample1.cc
    cpp/lib/googletest/samples/sample1.h
    cpp/lib/googletest/samples/sample10_unittest.cc
    cpp/lib/googletest/samples/sample1_unittest.cc
    cpp/lib/googletest/samples/sample2.cc
    cpp/lib/googletest/samples/sample2.h
    cpp/lib/googletest/samples/sample2_unittest.cc
    cpp/lib/googletest/samples/sample3-inl.h
    cpp/lib/googletest/samples/sample3_unittest.cc
    cpp/lib/googletest/samples/sample4.cc
    cpp/lib/googletest/samples/sample4.h
    cpp/lib/googletest/samples/sample4_unittest.cc
    cpp/lib/googletest/samples/sample5_unittest.cc
    cpp/lib/googletest/samples/sample6_unittest.cc
    cpp/lib/googletest/samples/sample7_unittest.cc
    cpp/lib/googletest/samples/sample8_unittest.cc
    cpp/lib/googletest/samples/sample9_unittest.cc
    cpp/lib/googletest/src/gtest-all.cc
    cpp/lib/googletest/src/gtest-death-test.cc
    cpp/lib/googletest/src/gtest-filepath.cc
    cpp/lib/googletest/src/gtest-internal-inl.h
    cpp/lib/googletest/src/gtest-port.cc
    cpp/lib/googletest/src/gtest-printers.cc
    cpp/lib/googletest/src/gtest-test-part.cc
    cpp/lib/googletest/src/gtest-typed-test.cc
    cpp/lib/googletest/src/gtest.cc
    cpp/lib/googletest/src/gtest_main.cc
    cpp/lib/googletest/test/gtest-death-test_ex_test.cc
    cpp/lib/googletest/test/gtest-death-test_test.cc
    cpp/lib/googletest/test/gtest-filepath_test.cc
    cpp/lib/googletest/test/gtest-linked_ptr_test.cc
    cpp/lib/googletest/test/gtest-listener_test.cc
    cpp/lib/googletest/test/gtest-message_test.cc
    cpp/lib/googletest/test/gtest-options_test.cc
    cpp/lib/googletest/test/gtest-param-test2_test.cc
    cpp/lib/googletest/test/gtest-param-test_test.cc
    cpp/lib/googletest/test/gtest-param-test_test.h
    cpp/lib/googletest/test/gtest-port_test.cc
    cpp/lib/googletest/test/gtest-printers_test.cc
    cpp/lib/googletest/test/gtest-test-part_test.cc
    cpp/lib/googletest/test/gtest-tuple_test.cc
    cpp/lib/googletest/test/gtest-typed-test2_test.cc
    cpp/lib/googletest/test/gtest-typed-test_test.cc
    cpp/lib/googletest/test/gtest-typed-test_test.h
    cpp/lib/googletest/test/gtest-unittest-api_test.cc
    cpp/lib/googletest/test/gtest_all_test.cc
    cpp/lib/googletest/test/gtest_break_on_failure_unittest_.cc
    cpp/lib/googletest/test/gtest_catch_exceptions_test_.cc
    cpp/lib/googletest/test/gtest_color_test_.cc
    cpp/lib/googletest/test/gtest_env_var_test_.cc
    cpp/lib/googletest/test/gtest_environment_test.cc
    cpp/lib/googletest/test/gtest_filter_unittest_.cc
    cpp/lib/googletest/test/gtest_help_test_.cc
    cpp/lib/googletest/test/gtest_list_tests_unittest_.cc
    cpp/lib/googletest/test/gtest_main_unittest.cc
    cpp/lib/googletest/test/gtest_no_test_unittest.cc
    cpp/lib/googletest/test/gtest_output_test_.cc
    cpp/lib/googletest/test/gtest_pred_impl_unittest.cc
    cpp/lib/googletest/test/gtest_premature_exit_test.cc
    cpp/lib/googletest/test/gtest_prod_test.cc
    cpp/lib/googletest/test/gtest_repeat_test.cc
    cpp/lib/googletest/test/gtest_shuffle_test_.cc
    cpp/lib/googletest/test/gtest_sole_header_test.cc
    cpp/lib/googletest/test/gtest_stress_test.cc
    cpp/lib/googletest/test/gtest_throw_on_failure_ex_test.cc
    cpp/lib/googletest/test/gtest_throw_on_failure_test_.cc
    cpp/lib/googletest/test/gtest_uninitialized_test_.cc
    cpp/lib/googletest/test/gtest_unittest.cc
    cpp/lib/googletest/test/gtest_xml_outfile1_test_.cc
    cpp/lib/googletest/test/gtest_xml_outfile2_test_.cc
    cpp/lib/googletest/test/gtest_xml_output_unittest_.cc
    cpp/lib/googletest/test/production.cc
    cpp/lib/googletest/test/production.h
    cpp/lib/googletest/xcode/Samples/FrameworkSample/widget.cc
    cpp/lib/googletest/xcode/Samples/FrameworkSample/widget.h
    cpp/lib/googletest/xcode/Samples/FrameworkSample/widget_test.cc
    cpp/lib/CMakeLists.txt
    cpp/src/main.cpp
    cpp/src/project1.cpp
        cpp/test/test_project1.cpp
    cpp/CMakeLists.txt
    cpp/README.md)

add_executable(FEM_CarbonDiffusion ${SOURCE_FILES})