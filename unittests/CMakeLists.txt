
# packers.cpp matrix1.cpp 
add_executable(UTmatrix23 main.cpp packers.cpp matrix1.cpp matrix_algebra.cpp) #vector.cpp matrix.cpp)
set_property(TARGET UTmatrix23 PROPERTY CXX_STANDARD 23)
target_compile_options(UTmatrix23 PRIVATE -Wall 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)
target_compile_definitions(UTmatrix23 PRIVATE 
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)

target_link_options(UTmatrix23 PRIVATE )
target_include_directories(UTmatrix23 PRIVATE ../include )
target_link_libraries(UTmatrix23 gtest)