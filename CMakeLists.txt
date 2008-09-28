cmake_minimum_required(VERSION 2.4)
FIND_PACKAGE(Boost REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS} ../parsegenbank)

ADD_EXECUTABLE(tfnetbuilder TFNetBuilder.cpp)
ADD_EXECUTABLE(tfnetperturber TFNetPerturber.cpp)
TARGET_LINK_LIBRARIES(tfnetbuilder boost_system boost_program_options boost_filesystem GenBankParser boost_regex)
TARGET_LINK_LIBRARIES(tfnetperturber boost_system boost_program_options boost_filesystem GenBankParser boost_regex)
