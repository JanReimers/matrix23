
# 
add_executable(UTmatrix23 main.cpp initvm.cpp packers.cpp matrix1.cpp matrix_algebra.cpp blas.cpp benchmarks.cpp vector.cpp ../src/blas.cpp ../src/ran250.cpp) 
#add_executable(UTmatrix23 main.cpp vector.cpp  ../src/blas.cpp ../src/ran250.cpp) 
set_property(TARGET UTmatrix23 PROPERTY CXX_STANDARD 23)
target_compile_options(UTmatrix23 PRIVATE -Wall 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2 -fconcepts-diagnostics-depth=2 > #O3 -ftree-vectorize -mavx -march=native
)
target_compile_definitions(UTmatrix23 PRIVATE 
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)

target_link_options(UTmatrix23 PRIVATE )
target_include_directories(UTmatrix23 PRIVATE ../include )
target_link_libraries(UTmatrix23 blas gtest)