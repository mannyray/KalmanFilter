cmake_minimum_required(VERSION 3.0)

#Locate GTest
find_package(GTest REQUIRED)
#find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${GTEST_INCLUDE_DIRS})

#Link runTests with what we want to test and the GTest and pthread library
add_executable(runTests tests.cpp)
target_link_libraries(runTests ${GTEST_LIBRARIES} pthread) #Eigen3::Eigen)


add_executable(runTests2 tests2.cpp)
target_link_libraries(runTests2 ${GTEST_LIBRARIES} pthread) #Eigen3::Eigen)


add_executable(runTests3 tests3.cpp)
target_link_libraries(runTests3 ${GTEST_LIBRARIES} pthread) #Eigen3::Eigen)


add_executable(runTests4 tests4.cpp)
target_link_libraries(runTests4 ${GTEST_LIBRARIES} pthread) #Eigen3::Eigen)
