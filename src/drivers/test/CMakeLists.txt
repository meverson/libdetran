#------------------------------------------------------------------------------#
# libdetran/src/drivers/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran drivers tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_StupidParser test_StupidParser.cc)
TARGET_LINK_LIBRARIES(test_StupidParser 
                      detran)

INSTALL(FILES ${CMAKE_SOURCE_DIR}/drivers/test/test.inp
        DESTINATION ${CMAKE_BINARY_DIR}/drivers/test/)
        
INSTALL(FILES ${CMAKE_SOURCE_DIR}/drivers/test/test.h5
        DESTINATION ${CMAKE_BINARY_DIR}/drivers/test/)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

# StupidParser
ADD_TEST(test_StupidParser       test_StupidParser 0)
if (DETRAN_ENABLE_HDF5)
ADD_TEST(test_StupidParser_hdf5  test_StupidParser 1)
endif (DETRAN_ENABLE_HDF5)