#------------------------------------------------------------------------------#
# libdetran/src/angle/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran angle tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_GaussLegendre           test_GaussLegendre.cc)
TARGET_LINK_LIBRARIES(test_GaussLegendre    angle utilities)

ADD_EXECUTABLE(test_QuadrupleRange          test_QuadrupleRange.cc)
TARGET_LINK_LIBRARIES(test_QuadrupleRange   angle utilities)

ADD_EXECUTABLE(test_Collocated              test_Collocated.cc)
TARGET_LINK_LIBRARIES(test_Collocated       angle utilities)

ADD_EXECUTABLE(test_Uniform                 test_Uniform.cc)
TARGET_LINK_LIBRARIES(test_Uniform          angle utilities)

ADD_EXECUTABLE(test_TabuchiYamamoto         test_TabuchiYamamoto.cc)
TARGET_LINK_LIBRARIES(test_TabuchiYamamoto  angle utilities)

ADD_EXECUTABLE(test_MomentToDiscrete        test_MomentToDiscrete.cc)
TARGET_LINK_LIBRARIES(test_MomentToDiscrete angle utilities)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_GaussLegendre_basic   test_GaussLegendre      0)
ADD_TEST(test_QuadrupleRange_basic  test_QuadrupleRange     0)
ADD_TEST(test_Collocated            test_Collocated         0)
ADD_TEST(test_Uniform               test_Uniform            0)
ADD_TEST(test_TabuchiYamamoto       test_TabuchiYamamoto    0)
ADD_TEST(test_MomentToDiscrete_1D   test_MomentToDiscrete   0)
ADD_TEST(test_MomentToDiscrete_2D   test_MomentToDiscrete   0)

