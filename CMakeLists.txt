FIND_PACKAGE(Boost REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS} ../parsegenbank)

ADD_EXECUTABLE(tfnetbuilder TFNetBuilder.cpp)
TARGET_LINK_LIBRARIES(tfnetbuilder boost_system boost_program_options boost_filesystem GenBankParser boost_regex)
