#------------------------------------------------------------------------------#
# libdetran/src/transport/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran geometry tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

# EQUATIONS

ADD_EXECUTABLE(test_Equation_DD_2D          test_Equation_DD_2D.cc)
TARGET_LINK_LIBRARIES(test_Equation_DD_2D   transport geometry materials utilities)

ADD_EXECUTABLE(test_Equation_SC_MOC         test_Equation_SC_MOC.cc)
TARGET_LINK_LIBRARIES(test_Equation_SC_MOC  transport geometry materials utilities)

# BOUNDARY

ADD_EXECUTABLE(test_BoundaryMOC             test_BoundaryMOC.cc)
TARGET_LINK_LIBRARIES(test_BoundaryMOC      transport angle geometry materials utilities)

# STATE 

ADD_EXECUTABLE(test_State                   test_State.cc)
TARGET_LINK_LIBRARIES(test_State            transport geometry materials utilities)

# SWEEPERS

ADD_EXECUTABLE(test_Sweeper2D               test_Sweeper2D.cc)
TARGET_LINK_LIBRARIES(test_Sweeper2D        transport geometry materials utilities)

ADD_EXECUTABLE(test_Sweeper3D               test_Sweeper3D.cc)
TARGET_LINK_LIBRARIES(test_Sweeper3D        transport ioutils)

# ACCELERATION

ADD_EXECUTABLE(test_CoarseMesh              test_CoarseMesh.cc)
TARGET_LINK_LIBRARIES(test_CoarseMesh       transport)

ADD_EXECUTABLE(test_CurrentTally            test_CurrentTally.cc)
TARGET_LINK_LIBRARIES(test_CurrentTally     transport)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_Equation_DD_2D_basic test_Equation_DD_2D  0)
ADD_TEST(test_Equation_SC_MOC      test_Equation_SC_MOC 0)
ADD_TEST(test_State_basic          test_State           0)
ADD_TEST(test_Sweeper2D_basic      test_Sweeper2D       0)
ADD_TEST(test_Sweeper3D_basic      test_Sweeper3D       0)
ADD_TEST(test_BoundaryMOC          test_BoundaryMOC     0)
ADD_TEST(test_CoarseMesh           test_CoarseMesh      0)
ADD_TEST(test_CurrentTally_1D      test_CurrentTally    0)
ADD_TEST(test_CurrentTally_2D      test_CurrentTally    1)
ADD_TEST(test_CurrentTally_3D      test_CurrentTally    2)
