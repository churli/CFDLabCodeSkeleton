cmake_minimum_required(VERSION 3.9)
project(sim)

set(CMAKE_C_STANDARD 99)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Werror")

set(SOURCE_FILES main.c boundary_val.c helper.c init.c sor.c uvp.c visual.c logger.c)
add_executable(sim ${SOURCE_FILES})
target_link_libraries(sim m)

# The below is to always get an updated copy of cavity100.dat inside the cmake-build-debug folder where the binary is.
add_custom_target(copy_cavity100 COMMAND cp cavity100.dat ${sim_BINARY_DIR}/ WORKING_DIRECTORY ${sim_SOURCE_DIR})
add_dependencies(sim copy_cavity100)
#add_custom_command(TARGET sim POST_BUILD COMMAND cp cavity100.dat ${sim_BINARY_DIR}/ WORKING_DIRECTORY ${sim_SOURCE_DIR})
#add_custom_command(OUTPUT execute_always COMMAND cp cavity100.dat ${sim_BINARY_DIR}/
#        WORKING_DIRECTORY ${sim_SOURCE_DIR} DEPENDS ${sim_SOURCE_DIR}/cavity100.dat)
